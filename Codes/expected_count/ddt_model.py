"""
ddt_model.py -- collision ratio rho_j of Section 2.4 computed from the DDT of the S-box.

For a pair whose second-round input difference Delta X_2 lies in Col(0), rho_j(Delta X_2)
is 2^{4w} times the probability that the ciphertexts of the pair (5-round AES without the
final MixColumns) have a zero difference on the inverse diagonal InvDia(j), under the
assumption of independent and uniformly distributed round keys:

    pi_j(D)  = sum_{Delta X_3} sum_{Delta X_4 : Delta X_4[Dia(j)] = 0}
               Pr[D -> Delta X_3] Pr[Delta X_3 -> Delta X_4],
    rho_j(D) = pi_j(D) / 2^{-4w}.

Delta X_3 is determined by the second-round S-box output differences Delta Y_2 of the m
nonzero bytes of D.  For a fixed Delta X_3 the probability that Delta X_4 is zero on Dia(j)
is a product over the four columns of Z_3, each evaluated exactly with the Walsh-Hadamard
transform of the DDT rows.  rho_exact sums over all possible values of Delta Y_2;
rho_sample draws Delta Y_2 at random with the probabilities of the DDT.

Usage: import ddt_model as dm; dm.init(w)   # w = 4 (small-AES) or w = 8 (AES)
"""
import numpy as np, itertools

M = [[2, 3, 1, 1], [1, 2, 3, 1], [1, 1, 2, 3], [3, 1, 1, 2]]


def init(w_):
    """Build the field, S-box, DDT and Walsh tables for cell size w_ (4 or 8)."""
    global w, q, POLY, MUL, INV, S, P, had, H, MI
    w = w_
    q = 1 << w
    POLY = {4: 0x13, 8: 0x11B}[w]
    MUL = np.array([[gm(a, b) for b in range(q)] for a in range(q)], dtype=np.int64)
    INV = np.array([0] + [int(np.where(MUL[a] == 1)[0][0]) for a in range(1, q)])
    if w == 8:
        Sl = [0] * 256; p_ = q_ = 1
        while True:
            p_ = p_ ^ ((p_ << 1) & 0xFF) ^ (0x1B if p_ & 0x80 else 0)
            q_ ^= q_ << 1; q_ ^= q_ << 2; q_ ^= q_ << 4; q_ &= 0xFF
            if q_ & 0x80: q_ ^= 0x09
            x = q_ ^ ((q_ << 1) | (q_ >> 7)) & 0xFF ^ ((q_ << 2) | (q_ >> 6)) & 0xFF ^ ((q_ << 3) | (q_ >> 5)) & 0xFF ^ ((q_ << 4) | (q_ >> 4)) & 0xFF
            Sl[p_] = (x ^ 0x63) & 0xFF
            if p_ == 1: break
        Sl[0] = 0x63
    else:
        Sl = [0x6, 0xB, 0x5, 0x4, 0x2, 0xE, 0x7, 0xA, 0x9, 0xD, 0xF, 0xC, 0x3, 0x1, 0x0, 0x8]
    S = np.array(Sl)
    P = np.zeros((q, q))
    for d in range(q):
        P[d] = np.bincount(S ^ S[np.arange(q) ^ d], minlength=q) / q   # DDT row / 2^w
    # Walsh-Hadamard: H[a][x][u] = sum_g P[x][g] (-1)^{<u, a g>}
    had = np.array([[(-1) ** bin(u & g).count('1') for u in range(q)] for g in range(q)], dtype=np.float64)
    H = {a: P @ had[MUL[a]] for a in (1, 2, 3)}
    MI = mat_inv()
    assert mcinv(mc([1, 2, 3, 4])) == [1, 2, 3, 4]


def gm(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & q: a ^= POLY
    return r


def mc(v):
    return [int(MUL[M[r][0]][v[0]] ^ MUL[M[r][1]][v[1]] ^ MUL[M[r][2]][v[2]] ^ MUL[M[r][3]][v[3]]) for r in range(4)]


def mat_inv():
    A = [[M[r][c] for c in range(4)] + [1 if r == c else 0 for c in range(4)] for r in range(4)]
    for i in range(4):
        pr = next(k for k in range(i, 4) if A[k][i]); A[i], A[pr] = A[pr], A[i]
        iv = INV[A[i][i]]; A[i] = [int(MUL[iv][x]) for x in A[i]]
        for k in range(4):
            if k != i and A[k][i]:
                f = A[k][i]; A[k] = [x ^ int(MUL[f][y]) for x, y in zip(A[k], A[i])]
    return [row[4:] for row in A]


def mcinv(v):
    return [int(MUL[MI[r][0]][v[0]] ^ MUL[MI[r][1]][v[1]] ^ MUL[MI[r][2]][v[2]] ^ MUL[MI[r][3]][v[3]]) for r in range(4)]


def geometry(D, j):
    rows = [r for r in range(4) if D[r]]
    cols = []   # per column c of Z_3: list of (r, x-coefficient, a-coefficient)
    for c in range(4):
        k = (c - j) % 4
        cols.append([(r, M[(-r - c) % 4][r], M[k][(-r - c) % 4]) for r in rows])
    return rows, cols


def rho_exact(D, j):
    """rho_j(D) summed over all possible values of Delta Y_2."""
    rows, cols = geometry(D, j)
    m = len(rows)
    if m == 1: return 0.0
    supp = [np.nonzero(P[D[r]])[0] for r in rows]
    grids = np.array(list(itertools.product(*supp)))             # (N, m) second-round S-box outputs
    wgt = np.ones(len(grids))
    for i, r in enumerate(rows): wgt *= P[D[r]][grids[:, i]]
    prod = np.ones(len(grids))
    for c in range(4):
        acc = np.ones((len(grids), q))
        for i, (r, xc, a) in enumerate(cols[c]):
            acc *= H[a][MUL[grids[:, i], xc]]
        prod *= acc.sum(axis=1) / q
    return float((wgt * prod).sum() * q ** 4)


def rho_sample(D, j, n, rng):
    """rho_j(D) estimated from n values of Delta Y_2 drawn with the DDT probabilities; returns (mean, standard error)."""
    rows, cols = geometry(D, j)
    if len(rows) == 1: return 0.0, 0.0
    dl = np.stack([rng.choice(q, size=n, p=P[D[r]]) for r in rows], axis=1)
    prod = np.ones(n)
    for c in range(4):
        acc = np.ones((n, q))
        for i, (r, xc, a) in enumerate(cols[c]):
            acc *= H[a][MUL[dl[:, i], xc]]
        prod *= acc.sum(axis=1) / q
    v = prod * q ** 4
    return float(v.mean()), float(v.std() / np.sqrt(n))


def wt(v):
    return sum(1 for x in v if x)


def mixture(dab, dac):
    return all(a == 0 or c == 0 or a == c for a, c in zip(dab, dac))


def mixture_bases(nz):
    """Delta_{a,b} X_2 of the candidates whose four X_2 states form a mixture quartet (Lemma 7):
    dict base difference (tuple) -> number of partner pairs of one base pair with that difference
    that give a mixture quartet.  Byte positions as in Section 3.3."""
    mix = {}
    if nz == 2:   # base pair active on bytes 0,5,10,15; partner exchanges bytes 0 and 10
        for T in itertools.combinations(range(4), 2):
            free = [k for k in range(4) if k not in T]
            for a, b in itertools.product(range(1, q), repeat=2):
                D = [0] * 4; D[free[0]] = a; D[free[1]] = b
                dl = mcinv(D)
                if 0 in dl: continue
                if mixture(D, mc([dl[0], 0, dl[2], 0])): mix[tuple(D)] = mix.get(tuple(D), 0) + 1
    elif nz == 1:  # base pair active on bytes 0,10; partner changes bytes 5 and 15 by nonzero amounts
        for x, y in itertools.product(range(1, q), repeat=2):
            D = mc([x, 0, y, 0])
            if wt(D) != 3: continue
            z = D.index(0); n = 0
            for bits in range(8):
                nzk = [k for k in range(4) if k != z]
                base = [0] * 4
                for b, k in enumerate(nzk):
                    if (bits >> b) & 1: base[k] = D[k]
                for zv in range(q):
                    dac = list(base); dac[z] = zv
                    v = mcinv(dac)
                    if v[0] == 0 and v[2] == 0 and v[1] and v[3]: n += 1
            if n: mix[tuple(D)] = n
    else:          # base pair active on bytes 5,15; partner exchanges byte 15 and changes bytes 0,10
        for x, y in itertools.product(range(1, q), repeat=2):
            D = mc([0, x, 0, y])
            if wt(D) != 4: continue
            n = 0
            for bits in range(16):
                dac = [D[k] if (bits >> k) & 1 else 0 for k in range(4)]
                v = mcinv(dac)
                if v[1] == 0 and v[3] == y and (v[0] or v[2]): n += 1
            if n: mix[tuple(D)] = n
    return mix


def rho_uniform(m, w_):
    """(2^w theta_m)^4 of Section 2.4 (S-box output differences uniform over the nonzero values)."""
    qq = 1 << w_
    th = 0.0
    for _ in range(1, m):
        th = (1 - th) / (qq - 1)
    return (qq * th) ** 4
