"""
acp_budget_sim.py -- the budget-limited ACP procedure of acp_all.c (n_z = 1 and n_z = 0),
simulated at the level of base collisions, and the resulting compound-Poisson P_find.

acp_all.c (nz = 1 or 0): while data < D: encrypt a fresh base structure (fresh random 16-byte
base, hence a fresh coset of D_0) of min(2^16, D - data) plaintexts; stop if data >= D; for each
inverse diagonal d and each colliding base pair whose two base-active bytes both differ, query
partner pairs one after another (2 queries each) until the partner set (65025 for n_z = 1,
65535 for n_z = 0, dead offsets included) is exhausted or data >= D.

Per base structure and d, the number of colliding 2-active base pairs is Poisson with mean
(2^16 * 65025 / 2) * 2^-32 = 0.49611 (rho = 1 to within 1e-4 for m = 3, 4).

In distinct cosets the members of a bundle are found one at a time (bundle_geometry.json: the
1024 (byte5,byte15) values of an n_z=1 bundle are distinct, as are the 8 (byte0,byte10) values
of an n_z=0 bundle), so a colliding base pair with n_c queried partners yields a valid quartet
with probability r * n_c, r = E[N_q] / (number of (colliding base pair, partner) checks of a full
structure) = E[N_q] / 255^4 (n_z = 1) and E[N_q] / (255^2 * 65535) (n_z = 0).
P_find = 1 - E[ prod_c (1 - r n_c) ],  mean = r E[N_p].
Writes data/acp_budget_sim.json.
Usage: python3 acp_budget_sim.py [trials]   (default 40000)
The committed data/acp_budget_sim.json was produced with 40000 trials (the default) and rng seed 1.
"""
import numpy as np, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
T5 = json.load(open(os.path.join(HERE, 'data', 'table5_w8.json')))
LAM_D = (65536 * 65025 / 2) * 2.0 ** -32


def simulate(nz, Dexp, trials, rng):
    D = int(2.0 ** Dexp + 0.5)
    npart = 65025 if nz == 1 else 65535
    E = T5['nz%d' % nz]['E']
    r = E / (255 ** 4) if nz == 1 else E / (255 ** 2 * 65535)
    logq = np.zeros(trials); Np = np.zeros(trials); nstruct = np.zeros(trials); ncoll = np.zeros(trials)
    for t in range(trials):
        data = 0; lq = 0.0; n_p = 0; ns = 0; nc = 0
        while data < D:
            lim = min(D - data, 65536); data += lim; ns += 1
            if data >= D:
                break
            c = rng.poisson(4 * LAM_D)
            for _ in range(c):
                avail = (D - data) // 2 if data < D else 0
                if avail <= 0:
                    break
                n_c = min(npart, avail + (1 if (D - data) % 2 else 0))   # loop checks data<D before each pair
                n_c = min(n_c, npart)
                data += 2 * n_c; n_p += n_c; nc += 1
                lq += math.log1p(-r * n_c)
                if data >= D:
                    break
        logq[t] = lq; Np[t] = n_p; nstruct[t] = ns; ncoll[t] = nc
    P = 1 - np.exp(logq)
    return dict(nz=nz, Dexp=Dexp, D=D, trials=trials, r=r, log2r=math.log2(r),
                P_find=float(P.mean()), P_find_mc_se=float(P.std() / math.sqrt(trials)),
                mean_valid=float(r * Np.mean()), E_Np=float(Np.mean()), log2_E_Np=math.log2(Np.mean()),
                E_structures=float(nstruct.mean()), E_collisions_used=float(ncoll.mean()),
                base_plaintexts=float(np.minimum(nstruct * 65536, D).mean()))


def fixed_structures(nz, n_struct):
    """Yan et al.'s ACP: a fixed number of base structures in distinct cosets, all partners queried."""
    npart = 65025 if nz == 1 else 65535
    E = T5['nz%d' % nz]['E']
    r = E / (255 ** 4) if nz == 1 else E / (255 ** 2 * 65535)
    p = r * npart
    mu_c = n_struct * 4 * LAM_D
    lam = mu_c * p                        # E[(1-p)^C] = exp(-mu_c p) for C ~ Poisson(mu_c)
    return dict(nz=nz, n_struct=n_struct, collisions=mu_c, p_per_collision=p,
                P_find=1 - math.exp(-lam), mean_valid=mu_c * p,
                data=n_struct * 65536 + 2 * mu_c * npart, log2_data=math.log2(n_struct * 65536 + 2 * mu_c * npart))


def one_coset(n_struct):
    """All n_struct base structures of n_z=1 in ONE coset of D_0 (distinct (byte5,byte15) values):
    a bundle is hit when one of its 1024 distinct (byte5,byte15) values is among the n_struct values."""
    T1 = T5['nz1']['first_term']; mu_b = T1 / 512
    h = 1 - math.exp(sum(math.log((65536 - 1024 - i) / (65536 - i)) for i in range(n_struct)))
    eps = T5['nz1']['eps']
    lam = mu_b * h + eps * n_struct * 2.0 ** -15
    return dict(n_struct=n_struct, h=h, P_find=1 - math.exp(-lam))


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    out = {'budget': [], 'fixed': [], 'one_coset': []}
    for Dexp in (22.0, 23.32):
        for nz in (1, 0):
            res = simulate(nz, Dexp, trials, rng)
            out['budget'].append(res); print(res, flush=True)
    for nz in (1, 0):
        for ns in (32, 2 ** 5.75, 53, 54):
            res = fixed_structures(nz, ns); out['fixed'].append(res); print(res)
    for ns in (32, 54):
        res = one_coset(ns); out['one_coset'].append(res); print(res)
    json.dump(out, open(os.path.join(HERE, 'data', 'acp_budget_sim.json'), 'w'), indent=1)
