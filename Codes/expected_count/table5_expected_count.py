"""
table5_expected_count.py -- values of Theorem 11 (Table 5): rho-bar(n_z), epsilon(n_z), E[N_q].

    E[N_q] = 4(2^w-1) rho-bar(n_z) + epsilon(n_z)   (n_z = 2, 1)
    E[N_q] = 3(2^w-1) rho-bar(0)   + epsilon(0)     (n_z = 0)

rho-bar(n_z): mean of rho_j over j and over the candidates whose four X_2 states form a mixture
quartet (all such candidates; rho_j exact for w = 4 and for n_z = 2 at w = 8, otherwise from
1,000 values of Delta Y_2 drawn with the DDT probabilities).
epsilon(n_z): sum of 2^{-8w} rho_j^2 over j and over the candidates whose four X_2 states do not
form a mixture quartet (all candidates for w = 4; 1,000 random base pairs with all partner pairs
of each base pair for w = 8, rho_j^2 estimated unbiasedly as the product of two independent
500-sample estimates of rho_j).

Usage:  python3 table5_expected_count.py 4      (about 4 minutes)
        python3 table5_expected_count.py 8      (about 13 minutes; n_z = 2 dominates)
Writes table5_w<w>.json.  Seeds 2026 (rho-bar) and 7 (epsilon) reproduce the paper.
"""
import sys, time, json, itertools
import numpy as np
import ddt_model as dm

w = int(sys.argv[1])
dm.init(w)
q = dm.q
C = {2: 4, 1: 4, 0: 3}
out = {'w': w}

# ---- rho-bar(n_z) and the first term ----
rng = np.random.default_rng(2026)
for nz in (2, 1, 0):
    t0 = time.time()
    mix = dm.mixture_bases(nz)
    assert sum(mix.values()) == C[nz] * (q - 1), (nz, sum(mix.values()))   # Lemma 7
    vals = []
    for D, n in mix.items():
        for j in range(4):
            if w == 4 or nz == 2:
                r = dm.rho_exact(list(D), j); e = 0.0
            else:
                r, e = dm.rho_sample(list(D), j, 1000, rng)
            vals.append((r, e, n))
    den = sum(n for r, e, n in vals)
    rhobar = sum(r * n for r, e, n in vals) / den
    err = (sum((e * n) ** 2 for r, e, n in vals) ** 0.5) / den
    out[f'nz{nz}'] = dict(rhobar=rhobar, rhobar_se=err, first_term=C[nz] * (q - 1) * rhobar)
    print(f"w={w} n_z={nz}: rho-bar = {rhobar:.6f} (se {err:.1e}), first term = {C[nz]*(q-1)*rhobar:.4f}  ({time.time()-t0:.0f}s)", flush=True)

# ---- epsilon(n_z) ----
rng = np.random.default_rng(7)
def rho2_mean_j(D, exact, n_in=1000):
    """(1/4) sum_j rho_j(D)^2: exact at w = 4, product of two independent estimates otherwise."""
    if dm.wt(D) == 1: return 0.0
    tot = 0.0
    for j in range(4):
        if exact:
            r = dm.rho_exact(D, j); tot += r * r
        else:
            a, _ = dm.rho_sample(D, j, n_in // 2, rng); b, _ = dm.rho_sample(D, j, n_in // 2, rng); tot += a * b
    return tot / 4

for nz in (2, 1, 0):
    t0 = time.time()
    mix = dm.mixture_bases(nz)
    npart = {2: 1, 1: (q - 1) ** 2, 0: (q - 1) * (q + 1)}[nz]    # partner pairs per base pair
    if w == 4:
        if nz == 2:
            bases = [dm.mc(list(d)) for d in itertools.product(range(1, q), repeat=4)]
        elif nz == 1:
            bases = [dm.mc([x, 0, y, 0]) for x, y in itertools.product(range(1, q), repeat=2)]
        else:
            bases = [dm.mc([0, x, 0, y]) for x, y in itertools.product(range(1, q), repeat=2)]
        s = 0.0
        for D in bases:
            s += (npart - mix.get(tuple(D), 0)) * rho2_mean_j(D, True)
        eps = s / q ** 4; err = 0.0
    else:
        k = 4 if nz == 2 else 2
        vals = []
        for _ in range(1000):
            d = list(rng.integers(1, q, size=k))
            if nz == 2: D = dm.mc(d)
            elif nz == 1: D = dm.mc([d[0], 0, d[1], 0])
            else: D = dm.mc([0, d[0], 0, d[1]])
            vals.append((npart - mix.get(tuple(D), 0)) * rho2_mean_j(D, False))
        vals = np.array(vals)
        eps = vals.mean() * (q - 1) ** k / q ** 4
        err = vals.std() / np.sqrt(len(vals)) * (q - 1) ** k / q ** 4
    first = out[f'nz{nz}']['first_term']
    out[f'nz{nz}'].update(eps=float(eps), eps_se=float(err), E=float(first + eps))
    print(f"w={w} n_z={nz}: epsilon = {eps:.5f} (se {err:.1e}), E[N_q] = {first + eps:.2f}  ({time.time()-t0:.0f}s)", flush=True)

json.dump(out, open(f'table5_w{w}.json', 'w'), indent=1)
