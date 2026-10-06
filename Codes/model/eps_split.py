"""
eps_split.py -- split of epsilon(n_z) of Theorem 11 by the weight m of Delta_ab X2.

Case-2 valid quartets (Theorem 10) come in groups of g_m = 2^{m-1} (2^w)^{4-m}
(8, 2^{w+2}, 2^{2w+1} for m = 4, 3, 2), so the Case-2 part of the count contributes
    sum_m g_m eps_m    to Var(N_q)  and   sum_m eps_m / g_m   to the expected number of groups.
With delta the base-pair difference at Z1 (each candidate has two pairs with the same delta),
    eps(n_z) = 2^{-4w-2} sum_delta (n_p - n_mix(delta)) sum_j rho_j(MC delta)^2,
n_p = partners per base pair (1, (2^w-1)^2, 2^{2w}-1 for n_z = 2, 1, 0), n_mix = 1 for the
Lemma-7 bases and 0 otherwise.
  w = 4: rho_j exact from the DDT (ddt_model.rho_exact); checked against table5_w4.json.
  w = 8: rho_j by the uniform model of Section 2.4 for each weight (DDT values differ by <1e-4
         for m = 3, 4; m = 2 uses the DDT mean 1.0166 of Table 3); checked against table5_w8.json.
Writes data/eps_split.json.
"""
import numpy as np, itertools, json, os, sys, math
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'expected_count'))   # ddt_model.py
sys.path.insert(0, _here)
import ddt_model as dm

HERE = os.path.dirname(os.path.abspath(__file__))


def mc_np(MUL, Mx, v):
    return np.stack([MUL[Mx[r][0]][v[:, 0]] ^ MUL[Mx[r][1]][v[:, 1]] ^ MUL[Mx[r][2]][v[:, 2]] ^ MUL[Mx[r][3]][v[:, 3]] for r in range(4)], axis=1)


def mixture_count(MUL, MI, q, nz, deltas):
    """n_mix(delta) for each base delta (Z1 column), by the same search as bundle_geometry."""
    D = mc_np(MUL, dm.M, deltas)
    nmix = np.zeros(len(deltas), dtype=np.int64)
    for i, (dz, dx) in enumerate(zip(deltas, D)):
        act = [r for r in range(4) if dx[r]]
        ina = [r for r in range(4) if not dx[r]]
        cands = []
        for bits in itertools.product((0, 1), repeat=len(act)):
            for fr in itertools.product(range(q), repeat=len(ina)):
                a = [0] * 4
                for b, r in zip(bits, act):
                    if b: a[r] = int(dx[r])
                for f, r in zip(fr, ina):
                    a[r] = f
                cands.append(a)
        A = np.array(cands, dtype=np.int64)
        v = mc_np(MUL, MI, A)
        if nz == 1:
            ok = (v[:, 0] == 0) & (v[:, 2] == 0) & (v[:, 1] != 0) & (v[:, 3] != 0)
        elif nz == 0:
            ok = (v[:, 1] == 0) & (v[:, 3] == dz[3]) & ((v[:, 0] != 0) | (v[:, 2] != 0))
        else:
            ok = (v[:, 0] == dz[0]) & (v[:, 2] == dz[2]) & (v[:, 1] == 0) & (v[:, 3] == 0)
        nmix[i] = int(ok.sum())
    return nmix


def run_w4():
    dm.init(4)
    q = 16; MUL = dm.MUL.astype(np.int64); MI = dm.MI
    res = {}
    for nz in (1, 0, 2):
        if nz == 2:
            deltas = np.array(list(itertools.product(range(1, q), repeat=4)), dtype=np.int64)
            n_p = 1
        else:
            xs = np.array(list(itertools.product(range(1, q), repeat=2)), dtype=np.int64)
            deltas = np.zeros((len(xs), 4), dtype=np.int64)
            rows = (0, 2) if nz == 1 else (1, 3)
            deltas[:, rows[0]] = xs[:, 0]; deltas[:, rows[1]] = xs[:, 1]
            n_p = (q - 1) ** 2 if nz == 1 else q * q - 1
        D = mc_np(MUL, dm.M, deltas)
        nmix = mixture_count(MUL, MI, q, nz, deltas) if nz != 2 else None
        if nz == 2:
            # n_z=2 partner: swap of rows 0,2 of delta; n_mix computed vectorised
            P = np.zeros_like(deltas); P[:, 0] = deltas[:, 0]; P[:, 2] = deltas[:, 2]
            A = mc_np(MUL, dm.M, P)
            nmix = np.all((D == 0) | (A == 0) | (D == A), axis=1).astype(np.int64)
        cache = {}
        epsm = {}; firstm = {}
        for dz, dx, nm in zip(deltas, D, nmix):
            key = tuple(int(t) for t in dx)
            if key not in cache:
                cache[key] = [dm.rho_exact(list(key), j) for j in range(4)]
            r = cache[key]
            m = int((dx != 0).sum())
            epsm[m] = epsm.get(m, 0.0) + (n_p - nm) * sum(x * x for x in r)
            firstm[m] = firstm.get(m, 0.0) + nm * sum(r)
        epsm = {m: v * 2.0 ** (-4 * 4 - 2) for m, v in epsm.items()}
        firstm = {m: v * 2.0 ** -2 for m, v in firstm.items()}   # 2^{4w-2} n_mix candidates x 2^{-4w} rho_j
        res['nz%d' % nz] = dict(eps_by_m=epsm, eps_total=sum(epsm.values()), n_mix_total=int(nmix.sum()),
                               first_term_check=sum(firstm.values()))
        print('w4 nz', nz, res['nz%d' % nz], flush=True)
    return res


def count_weights_nz2_w8(MUL, MI):
    """number of 4-active delta (Z1) with wt(MC delta) = m, m=1..4, by enumerating D of weight <=3."""
    cnt = {}
    xs = np.arange(1, 256, dtype=np.int64)
    for m in (1, 2, 3):
        c = 0
        for T in itertools.combinations(range(4), m):
            if m == 3:
                g = np.array(list(itertools.product(range(1, 256), repeat=2)), dtype=np.int64)
                for a in range(1, 256):
                    Dv = np.zeros((len(g), 4), dtype=np.int64)
                    Dv[:, T[0]] = a; Dv[:, T[1]] = g[:, 0]; Dv[:, T[2]] = g[:, 1]
                    c += int(np.all(mc_np(MUL, MI, Dv) != 0, axis=1).sum())
            else:
                g = np.array(list(itertools.product(range(1, 256), repeat=m)), dtype=np.int64)
                Dv = np.zeros((len(g), 4), dtype=np.int64)
                for i, t in enumerate(T):
                    Dv[:, t] = g[:, i]
                c += int(np.all(mc_np(MUL, MI, Dv) != 0, axis=1).sum())
        cnt[m] = c
    cnt[4] = 255 ** 4 - sum(cnt.values())
    return cnt


def run_w8():
    dm.init(8)
    MUL = dm.MUL.astype(np.int64); MI = dm.MI
    rho = {2: 1.0166, 3: dm.rho_uniform(3, 8), 4: dm.rho_uniform(4, 8)}
    res = {}
    # n_z = 1, 0: delta 2-active (rows {0,2} or {1,3}); MC delta has weight 3 for 4*255 values
    for nz, nbase_mix in ((1, 1020), (0, 765)):
        n_p = 255 ** 2 if nz == 1 else 65535
        xs = np.array(list(itertools.product(range(1, 256), repeat=2)), dtype=np.int64)
        deltas = np.zeros((len(xs), 4), dtype=np.int64)
        rows = (0, 2) if nz == 1 else (1, 3)
        deltas[:, rows[0]] = xs[:, 0]; deltas[:, rows[1]] = xs[:, 1]
        wts = (mc_np(MUL, dm.M, deltas) != 0).sum(axis=1)
        n3 = int((wts == 3).sum()); n4 = int((wts == 4).sum())
        # Lemma-7 bases: weight 3 for n_z=1 (all 1020 weight-3 deltas), weight 4 for n_z=0
        mix3 = nbase_mix if nz == 1 else 0; mix4 = 0 if nz == 1 else nbase_mix
        e3 = 2.0 ** -34 * (n3 * n_p - mix3) * 4 * rho[3] ** 2
        e4 = 2.0 ** -34 * (n4 * n_p - mix4) * 4 * rho[4] ** 2
        res['nz%d' % nz] = dict(n_delta_by_wt={3: n3, 4: n4}, eps_by_m={3: e3, 4: e4}, eps_total=e3 + e4)
        print('w8 nz', nz, res['nz%d' % nz], flush=True)
    cnt = count_weights_nz2_w8(MUL, MI)
    e = {m: 2.0 ** -34 * (cnt[m] - (1020 if m == 2 else 0)) * 4 * rho[m] ** 2 for m in (2, 3, 4)}
    res['nz2'] = dict(n_delta_by_wt=cnt, eps_by_m=e, eps_total=sum(e.values()))
    print('w8 nz 2', res['nz2'], flush=True)
    return res


if __name__ == '__main__':
    out = {'w4': run_w4(), 'w8': run_w8()}
    for w in ('w4', 'w8'):
        ref = json.load(open(os.path.join(HERE, 'data', 'table5_%s.json' % w)))
        for nz in ('nz2', 'nz1', 'nz0'):
            out[w][nz]['eps_json'] = ref[nz]['eps']
            g = {m: 2 ** (m - 1) * (2 ** int(w[1:])) ** (4 - m) for m in (2, 3, 4)}
            em = {int(k): v for k, v in out[w][nz]['eps_by_m'].items() if int(k) >= 2}
            out[w][nz]['group_size'] = g
            out[w][nz]['var_case2'] = sum(g[m] * em[m] for m in em)
            out[w][nz]['var_case2_8eps'] = 8 * ref[nz]['eps']
            out[w][nz]['groups_case2'] = sum(em[m] / g[m] for m in em)
            out[w][nz]['groups_case2_eps_over_8'] = ref[nz]['eps'] / 8
    json.dump(out, open(os.path.join(HERE, 'data', 'eps_split.json'), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
