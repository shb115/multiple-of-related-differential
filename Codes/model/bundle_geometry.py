"""
bundle_geometry.py -- exact geometry of the bundles of Theorem 10 (AES, w = 8) and the
probability that a reduced structure contains at least one member of a bundle.

For each of the four pairings
    a  : n_z=1, base active on plaintext bytes {0,10}, partner shifts {5,15}   (n_z=1 on {0,10}, Yan et al.)
    c  : n_z=1, base active on {5,15}, partner shifts {0,10}                   (n_z=1 on {5,15})
    b  : n_z=0, base active on {5,15}, partner swaps byte 15, moves {0,10}
    n2 : n_z=2, base active on {0,5,10,15}, partner swaps bytes {0,10}
the script
  1. enumerates the (Delta_ab Z1, delta Z1) values for which the four X2 states form a
     mixture quartet (Lemma 7; counts must be 1020/1020/765/1020),
  2. builds bundles: the F2-affine family of X2 states of the mixture counterparts of a base
     pair, its partition into member quartets by the translation of Lemma 9, and maps every
     state to Z1 = MC^{-1}(X2) (plaintext bytes 0,5,10,15 are per-byte bijective images of the
     Z1 rows 0,1,2,3, so a random plaintext subset is a random Z1 subset),
  3. checks the disjointness of the members' coordinate values, and
  4. computes the bundle-hit probability h = Pr[at least one member lies in the structure]
     for the product grids S0 x S1 x S2 ({0,10} x {5} x {15}) for the pairings a, c, b, and
     for the A0 x A1 grids ({0,10} x {5,15}, s = 12417, 12418) for the pairings a, c, and
     compares h with the independent-thinning value 1 - exp(-k q); for n2 it records the
     number of members present on the A0 x A1 grid s = 11664 (4 bundles); for the ACP
     setting it records the number of distinct (byte0,byte10) and (byte5,byte15) values per
     bundle, from which cpoisson_model.py derives h for one ACP base structure
     (1024/2^16 for n_z = 1, 8/2^16 for n_z = 0).
Usage: python3 bundle_geometry.py [bundles=40] [draws=4000]
Writes data/bundle_geometry.json.
"""
import numpy as np, itertools, json, math, sys, os
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'expected_count'))   # ddt_model.py
sys.path.insert(0, _here)
import ddt_model as dm

HERE = os.path.dirname(os.path.abspath(__file__))
dm.init(8)
MUL = dm.MUL.astype(np.int64)
M = dm.M
MI = dm.MI


def mc_arr(v):            # v: (n,4) int array -> MC(v)
    return np.stack([MUL[M[r][0]][v[:, 0]] ^ MUL[M[r][1]][v[:, 1]] ^ MUL[M[r][2]][v[:, 2]] ^ MUL[M[r][3]][v[:, 3]] for r in range(4)], axis=1)


def mci_arr(v):
    return np.stack([MUL[MI[r][0]][v[:, 0]] ^ MUL[MI[r][1]][v[:, 1]] ^ MUL[MI[r][2]][v[:, 2]] ^ MUL[MI[r][3]][v[:, 3]] for r in range(4)], axis=1)


def mixture_ok(D, A):     # per byte: D=0 or A=0 or D=A   (D = Delta_ab X2, A = Delta_ac X2)
    return np.all((D == 0) | (A == 0) | (D == A), axis=-1)


# ---------------------------------------------------------------- 1. solutions of Lemma 7
def solutions(pairing):
    sols = []
    xs = np.array(list(itertools.product(range(1, 256), repeat=2)))
    if pairing in ('a', 'c', 'b'):
        rows_b = (0, 2) if pairing == 'a' else (1, 3)
        Z = np.zeros((len(xs), 4), dtype=np.int64); Z[:, rows_b[0]] = xs[:, 0]; Z[:, rows_b[1]] = xs[:, 1]
        D = mc_arr(Z)
        w = (D != 0).sum(axis=1)

        def accept(dz, v):
            if pairing == 'a':
                return (v[:, 0] == 0) & (v[:, 2] == 0) & (v[:, 1] != 0) & (v[:, 3] != 0)
            if pairing == 'c':
                return (v[:, 1] == 0) & (v[:, 3] == 0) & (v[:, 0] != 0) & (v[:, 2] != 0)
            return (v[:, 1] == 0) & (v[:, 3] == dz[:, 3]) & ((v[:, 0] != 0) | (v[:, 2] != 0))
        # every partner shift Delta_ac X2 of a mixture quartet lies in the family space V of D:
        # Delta_ac X2[k] in {0, D[k]} at the active bytes of D, arbitrary at its inactive bytes
        for wt in (4, 3):
            sel = w == wt
            Zs, Ds = Z[sel], D[sel]
            if wt == 4:
                for bits in range(16):
                    A = Ds * np.array([(bits >> r) & 1 for r in range(4)])[None, :]
                    v = mci_arr(A); ok = accept(Zs, v)
                    sols += [(tuple(map(int, a)), tuple(map(int, b))) for a, b in zip(Zs[ok], v[ok])]
            else:
                zpos = np.argmax(Ds == 0, axis=1)
                for bits in range(16):
                    for zv in range(256):
                        msk = np.array([[(bits >> r) & 1 for r in range(4)]] * len(Ds))
                        A = Ds * msk
                        A[np.arange(len(Ds)), zpos] = zv
                        # avoid double counting: the bit of the inactive byte is irrelevant
                        if np.any(msk[np.arange(len(Ds)), zpos] == 1):
                            keep = msk[np.arange(len(Ds)), zpos] == 0
                        else:
                            keep = np.ones(len(Ds), bool)
                        v = mci_arr(A); ok = accept(Zs, v) & keep
                        sols += [(tuple(map(int, a)), tuple(map(int, b))) for a, b in zip(Zs[ok], v[ok])]
    elif pairing == 'n2':
        # Delta_ab X2 of weight 2; Delta_ab Z1 = MC^{-1}(D) must be nonzero on all four rows
        for T in itertools.combinations(range(4), 2):
            Dl = np.zeros((len(xs), 4), dtype=np.int64); Dl[:, T[0]] = xs[:, 0]; Dl[:, T[1]] = xs[:, 1]
            Zl = mci_arr(Dl)
            keep = np.all(Zl != 0, axis=1)
            Dl, Zl = Dl[keep], Zl[keep]
            P = np.zeros_like(Zl); P[:, 0] = Zl[:, 0]; P[:, 2] = Zl[:, 2]      # swap of bytes 0 and 10
            A = mc_arr(P)
            ok = mixture_ok(Dl, A)
            for dz, p in zip(Zl[ok], P[ok]):
                sols.append((tuple(int(t) for t in dz), tuple(int(t) for t in p)))
    return sols


# ---------------------------------------------------------------- 2. bundle construction
def pack(v):
    return v[..., 0] | (v[..., 1] << 8) | (v[..., 2] << 16) | (v[..., 3] << 24)


def unpack(x):
    return np.stack([(x >> (8 * r)) & 0xFF for r in range(4)], axis=-1)


def build_bundle(dz, pz, rng):
    """Return (members, info): members is an int array (k, 4, 4) of the Z1 columns of the four
    states of each member, ordered (base_a, base_b, partner_c, partner_d)."""
    dz = np.array(dz, dtype=np.int64); pz = np.array(pz, dtype=np.int64)
    DX = mc_arr(dz[None, :])[0]; AX = mc_arr(pz[None, :])[0]
    basis = []
    for r in range(4):
        if DX[r]:
            basis.append(int(DX[r]) << (8 * r))
        else:
            basis += [(1 << b) << (8 * r) for b in range(8)]
    span = np.zeros(1, dtype=np.int64)
    for bv in basis:
        span = np.concatenate([span, span ^ bv])
    dX = int(pack(DX)); aX = int(pack(AX))
    xa = int(pack(rng.integers(0, 256, size=4)))
    st = xa ^ span                                          # all X2 states of the family
    stset = set(st.tolist())
    assert (xa ^ aX) in stset, "partner shift not in the family"
    orb = np.stack([st, st ^ dX, st ^ aX, st ^ dX ^ aX], axis=1)
    rep = orb.min(axis=1)
    _, first = np.unique(rep, return_index=True)
    mem = orb[first]                                        # (k, 4) packed X2 states
    Zm = mci_arr(unpack(mem.reshape(-1))).reshape(len(mem), 4, 4)
    return Zm, dict(dim=len(basis), k=len(mem))


# ---------------------------------------------------------------- 3./4. grid inclusion
def comp_values(Zm, rows):
    """For each member: the distinct values of the projection on Z1 rows `rows` (packed),
    as an array (k, 2); asserts that every member has exactly two distinct values."""
    pr = np.zeros(Zm.shape[:2], dtype=np.int64)
    for i, r in enumerate(rows):
        pr |= Zm[:, :, r] << (8 * i)
    srt = np.sort(pr, axis=1)
    nd = 1 + (np.diff(srt, axis=1) != 0).sum(axis=1)
    assert np.all(nd == 2), "member with %s distinct values on rows %s" % (set(nd.tolist()), rows)
    lo = srt[:, 0]; hi = srt[:, 3]
    return np.stack([lo, hi], axis=1)


def p_no_pair(Mpairs, s, n):
    """Pr[a uniformly random s-subset of an n-set contains none of Mpairs disjoint pairs]."""
    tot = 0.0
    t = 1.0
    for i in range(Mpairs + 1):
        if i > 0:
            t *= (s - 2 * i + 2) * (s - 2 * i + 1) / ((n - 2 * i + 2) * (n - 2 * i + 1))
        tot += (-1) ** i * math.comb(Mpairs, i) * t
        if t * math.comb(Mpairs, i) < 1e-18 and i > 4:
            break
    return tot


def frac2(s, n):
    return s * (s - 1) / (n * (n - 1))


def h_product(Zm, s0, s1, s2, trials, rng):
    """Bundle-hit probability on S0 x S1 x S2 (rows {0,2} x {1} x {3}); Rao-Blackwellised over S0.
    Requires the members' {0,2}-value pairs to be pairwise disjoint (checked)."""
    U = comp_values(Zm, (0, 2)); B5 = comp_values(Zm, (1,)); B15 = comp_values(Zm, (3,))
    assert len(np.unique(U)) == U.size, "u-values shared"
    acc = 0.0
    for _ in range(trials):
        in5 = np.zeros(256, bool); in5[rng.permutation(256)[:s1]] = True
        in15 = np.zeros(256, bool); in15[rng.permutation(256)[:s2]] = True
        sel = in5[B5[:, 0]] & in5[B5[:, 1]] & in15[B15[:, 0]] & in15[B15[:, 1]]
        acc += p_no_pair(int(sel.sum()), s0, 65536)
    return 1.0 - acc / trials


def h_A0A1(Zm, s0, s1, trials, rng):
    """Bundle-hit probability on A0 x A1 (rows {0,2} x {1,3}); RB over A0; A1 membership of the
    bundle's w-values drawn exactly (hypergeometric)."""
    U = comp_values(Zm, (0, 2)); W = comp_values(Zm, (1, 3))
    assert len(np.unique(U)) == U.size, "u-values shared"
    wv, winv = np.unique(W, return_inverse=True); winv = winv.reshape(W.shape)
    nw = len(wv)
    acc = 0.0
    for _ in range(trials):
        K = rng.hypergeometric(s1, 65536 - s1, nw)
        inw = np.zeros(nw, bool); inw[rng.permutation(nw)[:K]] = True
        sel = inw[winv[:, 0]] & inw[winv[:, 1]]
        acc += p_no_pair(int(sel.sum()), s0, 65536)
    return 1.0 - acc / trials


def present_stats_A0A1(Zm, s0, s1, trials, rng):
    """n_z=2 on A0 x A1: plain MC of the number of members present (u-values may be shared)."""
    U = comp_values(Zm, (0, 2)); W = comp_values(Zm, (1, 3))
    uv, uinv = np.unique(U, return_inverse=True); uinv = uinv.reshape(U.shape)
    wv, winv = np.unique(W, return_inverse=True); winv = winv.reshape(W.shape)
    cnt = []
    for _ in range(trials):
        Ku = rng.hypergeometric(s0, 65536 - s0, len(uv)); inu = np.zeros(len(uv), bool); inu[rng.permutation(len(uv))[:Ku]] = True
        Kw = rng.hypergeometric(s1, 65536 - s1, len(wv)); inw = np.zeros(len(wv), bool); inw[rng.permutation(len(wv))[:Kw]] = True
        cnt.append(int((inu[uinv[:, 0]] & inu[uinv[:, 1]] & inw[winv[:, 0]] & inw[winv[:, 1]]).sum()))
    cnt = np.array(cnt)
    return dict(mean=float(cnt.mean()), var=float(cnt.var()), p0=float((cnt == 0).mean()),
                distinct_u=int(len(uv)), distinct_w=int(len(wv)))


def main():
    rng = np.random.default_rng(20261002)
    out = {}
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 40          # bundles per pairing
    tr = int(sys.argv[2]) if len(sys.argv) > 2 else 4000        # grid draws per bundle
    grids = {'2p27.02': (11664, 108, 108), '2p27.20': (12515, 111, 111),
             '2p27.6': None, '2p27.8': None, '2p28.0': None, '2p28.3': None}
    for x in (27.6, 27.8, 28.0, 28.3):
        s1 = int(round(256 * 2 ** ((x - 32) / 4))); s0 = int(round(2 ** x / s1 ** 2))
        grids['2p%.1f' % x] = (s0, s1, s1)
    yan = {'s12417': 12417, 's12418': 12418}
    for pairing in ('a', 'c', 'b', 'n2'):
        sols = solutions(pairing)
        info = dict(n_solutions=len(sols))
        # check weights of Delta_ab X2 over the solutions
        DXs = mc_arr(np.array([s[0] for s in sols]))
        info['wt_DeltaX2'] = sorted(set(int((d != 0).sum()) for d in DXs))
        info['distinct_base_DeltaZ1'] = len(set(s[0] for s in sols))
        idx = rng.choice(len(sols), size=min(nb, len(sols)), replace=False)
        disj = dict(u_disjoint=0, w_disjoint=0, b5_shared=0, b15_shared=0)
        hp = {g: [] for g in grids}; hy = {g: [] for g in yan}; ks = set(); acpu = []; acpw = []
        st2 = []
        for t, i in enumerate(idx):
            if pairing == 'n2' and t >= 4:
                break
            Zm, bi = build_bundle(sols[i][0], sols[i][1], rng)
            ks.add(bi['k'])
            if pairing == 'n2':
                st2.append(present_stats_A0A1(Zm, 11664, 11664, 300, rng))
                continue
            U = comp_values(Zm, (0, 2))
            if pairing in ('a', 'c'):
                W = comp_values(Zm, (1, 3))
            else:   # n_z=0: a member has four distinct (byte5, byte15) values
                W = np.sort(Zm[:, :, 1] | (Zm[:, :, 3] << 8), axis=1)
                assert np.all(np.diff(W, axis=1) != 0)
            B5 = comp_values(Zm, (1,)); B15 = comp_values(Zm, (3,))
            disj['u_disjoint'] += int(len(np.unique(U)) == U.size)
            disj['w_disjoint'] += int(len(np.unique(W)) == W.size)
            disj['b5_shared'] += int(len(np.unique(B5)) < B5.size)
            disj['b15_shared'] += int(len(np.unique(B15)) < B15.size)
            acpu.append(len(np.unique(U))); acpw.append(len(np.unique(W)))
            for g, sz in grids.items():
                hp[g].append(h_product(Zm, *sz, trials=tr, rng=rng))
            if pairing in ('a', 'c'):
                for g, s in yan.items():
                    hy[g].append(h_A0A1(Zm, s, s, trials=tr, rng=rng))
        info['k'] = sorted(ks)
        if pairing == 'n2':
            info['A0xA1_11664_members_present'] = st2
        else:
            info['disjointness_over_bundles'] = disj
            info['n_bundles'] = len(acpu)
            info['distinct_u_values_per_bundle'] = sorted(set(acpu))
            info['distinct_w_values_per_bundle'] = sorted(set(acpw))
            k = info['k'][0]
            res = {}
            for g, sz in grids.items():
                s0, s1, s2 = sz
                q = frac2(s0, 65536) * frac2(s1, 256) * frac2(s2, 256)
                v = np.array(hp[g])
                res[g] = dict(s=sz, N=s0 * s1 * s2, log2N=math.log2(s0 * s1 * s2), q=q, log2q=math.log2(q),
                              kq=k * q, h_indep_exp=1 - math.exp(-k * q), h_indep_binom=1 - (1 - q) ** k,
                              h_geom=float(v.mean()), h_geom_se=float(v.std(ddof=1) / math.sqrt(len(v))),
                              h_geom_min=float(v.min()), h_geom_max=float(v.max()))
            info['product_grid'] = res
            if pairing in ('a', 'c'):
                ry = {}
                for g, s in yan.items():
                    q = frac2(s, 65536) ** 2
                    v = np.array(hy[g])
                    ry[g] = dict(s=s, N=s * s, log2N=math.log2(s * s), q=q, log2q=math.log2(q), kq=k * q,
                                 h_indep_exp=1 - math.exp(-k * q), h_indep_binom=1 - (1 - q) ** k,
                                 h_geom=float(v.mean()), h_geom_se=float(v.std(ddof=1) / math.sqrt(len(v))))
                info['A0xA1_grid'] = ry
        out[pairing] = info
        print(pairing, json.dumps(info, indent=1)[:3000], flush=True)
    with open(os.path.join(HERE, 'data', 'bundle_geometry.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
