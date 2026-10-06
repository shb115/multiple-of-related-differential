"""
table3_rho_m.py -- Table 3: mean of rho_j(Delta X_2) over j and over the differences Delta X_2
in Col(0) with m nonzero bytes, for m = 2, 3, 4, with the uniform value (2^w theta_m)^4.

w = 4: mean over all differences, rho_j over all possible values of Delta Y_2
       (m = 4 takes about 6 minutes).
w = 8: mean over 1,000 random differences; rho_j over all possible values of Delta Y_2 for m = 2
       (about 13 minutes) and from 1,000 values of Delta Y_2 drawn with the DDT probabilities
       for m = 3, 4.  Seed 100 + m reproduces the paper.

Usage:  python3 table3_rho_m.py 4      or      python3 table3_rho_m.py 8
Writes table3_w<w>.json.
"""
import sys, time, json, itertools
import numpy as np
import ddt_model as dm

w = int(sys.argv[1])
dm.init(w)
q = dm.q
out = {'w': w}
for m in (2, 3, 4):
    t0 = time.time()
    if w == 4:
        tot = 0.0; n = 0
        for rows in itertools.combinations(range(4), m):
            for nzv in itertools.product(range(1, q), repeat=m):
                D = [0] * 4
                for r, v in zip(rows, nzv): D[r] = v
                for j in range(4):
                    tot += dm.rho_exact(D, j); n += 1
        mean, se = tot / n, 0.0
    else:
        rng = np.random.default_rng(100 + m)
        vals, errs = [], []
        for _ in range(1000):
            rows = sorted(rng.choice(4, size=m, replace=False))
            D = [0] * 4
            for r in rows: D[r] = int(rng.integers(1, q))
            for j in range(4):
                if m == 2:
                    vals.append(dm.rho_exact(D, j)); errs.append(0.0)
                else:
                    r_, e_ = dm.rho_sample(D, j, 1000, rng); vals.append(r_); errs.append(e_)
        vals = np.array(vals); errs = np.array(errs)
        mean = float(vals.mean())
        se = float(np.sqrt(vals.var() / len(vals) + (errs ** 2).sum() / len(vals) ** 2))
    uni = dm.rho_uniform(m, w)
    out[f'm{m}'] = dict(ddt=mean, ddt_se=se, uniform=uni)
    print(f"w={w} m={m}: uniform {uni:.4f}  DDT {mean:.6f} (se {se:.1e})  ({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open(f'table3_w{w}.json', 'w'), indent=1)
