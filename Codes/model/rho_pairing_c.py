"""
rho_pairing_c.py -- rho-bar of the pairing n_z=1 on {5,15} (pairing c: base active on {5,15},
partner shifts {0,10}) against the pairing n_z=1 on {0,10} (pairing a, Yan et al.: base {0,10},
partner {5,15}), AES w=8.
Both are averaged over the 1020 base differences Delta_ab X2 (weight 3) of the Lemma-7
solutions, over j = 0..3, with rho_j estimated from n random Delta Y2 drawn from the DDT
(ddt_model.rho_sample, the estimator used for Table 5).  Writes data/rho_pairing_c.json.
Usage: python3 rho_pairing_c.py [n]   (n = 1000 for the committed JSON)
"""
import numpy as np, json, os, sys
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'expected_count'))   # ddt_model.py
sys.path.insert(0, _here)                                         # bundle_geometry.py
import bundle_geometry as bg
import ddt_model as dm

HERE = os.path.dirname(os.path.abspath(__file__))
n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
rng = np.random.default_rng(7)
out = {}
for p in ('a', 'c'):
    sols = bg.solutions(p)
    Ds = bg.mc_arr(np.array([s[0] for s in sols]))
    vals = []; ses = []
    for D in Ds:
        for j in range(4):
            m, se = dm.rho_sample([int(x) for x in D], j, n, rng)
            vals.append(m); ses.append(se)
    vals = np.array(vals); ses = np.array(ses)
    rb = float(vals.mean()); rb_se = float(np.sqrt((ses ** 2).sum()) / len(ses))
    out[p] = dict(n_bases=len(Ds), n_samples=n, rhobar=rb, rhobar_se=rb_se,
                  first_term=4 * 255 * rb, first_term_se=4 * 255 * rb_se)
    print(p, out[p], flush=True)
out['note'] = ('Results/expected_count/table5_w8.json (Table 5) gives rhobar(1)=0.9999414188 (SE 4.1e-6), first term 1019.9402 '
               'for pairing a; uniform-model rho for m=3 (Results/expected_count/table3_w8.json, Table 3) is 0.9999384866 = 1-2^-13.989')
json.dump(out, open(os.path.join(HERE, 'data', 'rho_pairing_c.json'), 'w'), indent=1)
