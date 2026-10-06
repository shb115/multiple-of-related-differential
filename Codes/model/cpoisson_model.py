"""
cpoisson_model.py -- compound-Poisson model of the success probability of the 5-round AES
related-differential distinguishers, and every number the paper needs from it.

Model.  In a full diagonal structure the Case-1 valid quartets of a pairing come in bundles of
k = 2^{2-n_z} (2^w)^{n_z}; the number of bundles is Poisson with mean mu_b = T1 / k, T1 the first
term of Theorem 11.  The Case-2 valid quartets come in groups of g_m = 2^{m-1} (2^w)^{4-m}
(m = weight of Delta_ab X2), the number of groups of size g_m being Poisson with mean eps_m / g_m,
sum_m eps_m = eps(n_z).  A structure with fewer plaintexts keeps each member of a bundle (group)
independently with probability q, the probability that the four plaintexts of a candidate lie in
the structure.  Hence
    P_find = 1 - exp( - mu_b (1 - e^{-k q}) - sum_m (eps_m / g_m)(1 - e^{-g_m q}) )
           ~ 1 - exp( - mu_b (1 - e^{-k q}) - eps q ),
    Var(N_q) (full structure) = k T1 + sum_m g_m eps_m  ~  k T1 + 8 eps.
Inputs: data/table5_w8.json, data/table5_w4.json (verbatim copies of Results/expected_count/table5_w{8,4}.json,
Table 5), data/table3_w8.json (verbatim copy of Results/expected_count/table3_w8.json, Table 3),
data/bundle_geometry.json, data/eps_split.json, data/acp_budget_sim.json, data/rho_pairing_c.json
(written by the other scripts of this directory), measured data in data/5r_nz{0,1,2}.txt,
data/raw_nz{0,1,2}_5r.txt, data/cp_nz01_2p27{,2}.txt, data/cp_nz2_2p27{,2}.txt and data/acp_all_2p2332.txt
(verbatim copies of the files of the same name under Results/), and the Section 4.3 chunk files under
Results/5r_AES_distinguishers/combined (read through Codes/tools/aggregate_combined.py).
Outputs: data/cpoisson_results.json and data/cpoisson_results.txt.
"""
import json, math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
import aggregate_combined
DATA = os.path.join(HERE, 'data')
J = lambda f: json.load(open(os.path.join(DATA, f)))
T = {8: J('table5_w8.json'), 4: J('table5_w4.json')}
BG = J('bundle_geometry.json')
ES = J('eps_split.json')
ACP = J('acp_budget_sim.json')
try:
    RC = J('rho_pairing_c.json')
except FileNotFoundError:
    RC = None
K = {8: {2: 65536, 1: 512, 0: 4}, 4: {2: 256, 1: 32, 0: 4}}
OUT = {}
LINES = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LINES.append(s)


def frac2(s, n):
    return s * (s - 1) / (n * (n - 1))


def q_prod(s0, s1, s2):
    return frac2(s0, 65536) * frac2(s1, 256) * frac2(s2, 256)


def q_A(s0, s1):
    return frac2(s0, 65536) * frac2(s1, 65536)


def eps_m(w, nz):
    e = ES['w%d' % w]['nz%d' % nz]['eps_by_m']
    return {int(m): v for m, v in e.items() if int(m) >= 2}


def gsize(w, m):
    return 2 ** (m - 1) * (2 ** w) ** (4 - m)


def T1(nz, w=8, pairing=None):
    if pairing == 'c' and RC is not None:
        return RC['c']['first_term']
    return T[w]['nz%d' % nz]['first_term']


def EPS(nz, w=8):
    return T[w]['nz%d' % nz]['eps']


def lam(nz, q, w=8, h=None, pairing=None):
    """Poisson rate of 'at least one valid quartet' contributions in a structure with keep-prob q."""
    k = K[w][nz]
    mu_b = T1(nz, w, pairing) / k
    if h is None:
        h = -math.expm1(-k * q)
    c2 = sum(e / gsize(w, m) * -math.expm1(-gsize(w, m) * q) for m, e in eps_m(w, nz).items())
    return mu_b * h + c2


def Pf(*lams):
    return -math.expm1(-sum(lams))


def binse(p, n):
    return math.sqrt(max(p * (1 - p), 1e-12) / n)


def z(meas, pred, n):
    return (meas - pred) / binse(pred, n)


def hgeom(pairing, grid_key, kind='product_grid'):
    try:
        return BG[pairing][kind][grid_key]['h_geom']
    except KeyError:
        return None


# ------------------------------------------------------------------ measured data
def parse_cp():
    m = {}
    for f in ('cp_nz01_2p27.txt', 'cp_nz01_2p272.txt'):
        s = open(os.path.join(DATA, f)).read()
        g = re.search(r's0=(\d+) s1=(\d+) s2=(\d+)', s)
        sz = tuple(int(x) for x in g.groups())
        for nz in (0, 1):
            r = re.search(r'n_z=%d: P_find=(\d+)/(\d+)=[\d.]+ mean=([\d.]+)' % nz, s)
            m[(nz, sz)] = dict(hits=int(r.group(1)), n=int(r.group(2)), mean=float(r.group(3)))
    for f in ('cp_nz2_2p27.txt', 'cp_nz2_2p272.txt'):
        s = open(os.path.join(DATA, f)).read()
        g = re.search(r's0=(\d+) s1=(\d+)', s)
        sz = tuple(int(x) for x in g.groups())
        r = re.search(r'P_find=(\d+)/(\d+)=[\d.]+ mean=([\d.]+)', s)
        m[(2, sz)] = dict(hits=int(r.group(1)), n=int(r.group(2)), mean=float(r.group(3)))
    return m


def parse_acp():
    s = open(os.path.join(DATA, 'acp_all_2p2332.txt')).read()
    m = {}
    for r in re.finditer(r'ACP nz=(\d) D=2\^([\d.]+) .*?P_find=(\d+)/(\d+)=[\d.]+ mean=([\d.]+)', s):
        m[(int(r.group(1)), float(r.group(2)))] = dict(hits=int(r.group(3)), n=int(r.group(4)), mean=float(r.group(5)))
    return m


def freq(f):
    d = {}
    for r in re.finditer(r'Count\s+(\d+)\s*:\s*(\d+)\s*times', open(os.path.join(DATA, f)).read()):
        d[int(r.group(1))] = d.get(int(r.group(1)), 0) + int(r.group(2))
    return d


def raw(f):
    v = []
    for line in open(os.path.join(DATA, f)):
        t = re.findall(r'-?\d+', line)
        if t:
            v.append(int(t[-1]))
    return np.array(v, dtype=np.float64)


def boot_sd(x, B=4000, seed=5):
    rng = np.random.default_rng(seed)
    n = len(x)
    s = [np.std(x[rng.integers(0, n, n)], ddof=1) for _ in range(B)]
    return float(np.std(s, ddof=1))


# ================================================================== (i)-(ii) q
say('=' * 100)
say('(i)-(iv)  inclusion probability q')
say('=' * 100)
QQ = {}
for sz in ((11664, 108, 108), (12515, 111, 111)):
    q = q_prod(*sz); N = sz[0] * sz[1] * sz[2]
    QQ['product_%d_%d_%d' % sz] = dict(N=N, log2N=math.log2(N), q=q, log2q=math.log2(q), ideal_log2q=2 * (math.log2(N) - 32))
    say('product grid S0xS1xS2 =', sz, 'N=2^%.4f' % math.log2(N), ' q = s0(s0-1)/(2^16(2^16-1)) * s1(s1-1)/(256*255) * s2(s2-1)/(256*255) = %.6e = 2^%.4f  [(N/2^32)^2 = 2^%.4f]' % (q, math.log2(q), 2 * (math.log2(N) - 32)),
        ' -- same q for pairings a, b, c (each candidate uses 2 {0,10} values, 2 byte-5 values, 2 byte-15 values)')
for s in (11664, 12417, 12418):
    q = q_A(s, s)
    QQ['A0xA1_%d' % s] = dict(N=s * s, log2N=math.log2(s * s), q=q, log2q=math.log2(q))
    say('A0xA1 grid s=%d N=2^%.4f  q = [s(s-1)/(2^16(2^16-1))]^2 = %.6e = 2^%.4f' % (s, math.log2(s * s), q, math.log2(q)))
say('ACP base structure (one fresh coset, all 2^16 values of the base-active bytes, every partner queried):')
say('  per member q = 2/2^16 = 2^-15 (a member is found iff one of its two pairs lies in the slice);')
say('  n_z=1: the 512 members of a bundle use 1024 distinct (byte5,byte15) values (checked on %d bundles), so the members are mutually exclusive and h = 1024/2^16 = 2^-6 exactly (independent thinning would give 1-e^{-1/64} = %.6f)' % (BG['a']['n_bundles'], -math.expm1(-1 / 64)))
say('  n_z=0: the 4 members use 8 distinct (byte0,byte10) values, h = 8/2^16 = 2^-13 exactly')
say('  per base structure: L = mu_b h + eps 2^-15;   2^5 structures in distinct cosets: Lambda = 32 L')
OUT['q'] = QQ

# ================================================================== model statement constants
say('')
say('constants (AES, Table 5, data/table5_w8.json):')
CONST = {}
for nz in (2, 1, 0):
    k = K[8][nz]
    CONST[nz] = dict(k=k, T1=T1(nz), eps=EPS(nz), E=T[8]['nz%d' % nz]['E'], mu_b=T1(nz) / k,
                     eps_by_m=eps_m(8, nz))
    say('  n_z=%d k=%d T1=%.4f eps=%.4f mu_b=T1/k=%.6f  eps_m=%s' % (nz, k, T1(nz), EPS(nz), T1(nz) / k,
                                                                    {m: round(v, 6) for m, v in eps_m(8, nz).items()}))
if RC is not None:
    say('  pairing c (n_z=1 on {5,15}): rhobar = %.6f +- %.6f -> T1 = %.3f (pairing a recomputed: %.6f +- %.6f, T1 %.3f; Table 5: 1019.940)' % (
        RC['c']['rhobar'], RC['c']['rhobar_se'], RC['c']['first_term'], RC['a']['rhobar'], RC['a']['rhobar_se'], RC['a']['first_term']))
T3 = J('table3_w8.json')
S41 = {}
for nz in (2, 1, 0):
    lam_nz = (3 if nz == 0 else 4) / 255 ** 3
    rb = T[8]['nz%d' % nz]['rhobar']
    S41[nz] = dict(lam=lam_nz, rhobar=rb, log2_factor=math.log2(2 ** 32 * lam_nz * rb))
    say('  n_z=%d: lambda = %d/(2^8-1)^3 = 2^%.2f; a candidate is valid with probability 2^-30 lambda rhobar = 2^-62 x 2^%.2f (Section 4.1)' % (
        nz, 3 if nz == 0 else 4, math.log2(lam_nz), S41[nz]['log2_factor']))
for nz in (2, 1, 0):
    wt = {int(m): c for m, c in ES['w8']['nz%d' % nz]['n_delta_by_wt'].items()}
    w2 = {m: c for m, c in wt.items() if m >= 2}
    rho_bp = sum(c * T3['m%d' % m]['ddt'] for m, c in w2.items()) / sum(w2.values())
    S41[nz]['mean_rho_base_pairs'] = rho_bp
    say('  n_z=%d: mean rho over the base pairs = %.7f (weights %s of the base-pair differences at X2 from data/eps_split.json, rho_m of Table 3 from data/table3_w8.json%s) (Section 4.1)' % (
        nz, rho_bp, w2, '; the %d base pairs with m = 1, a fraction %.1e, have no Table 3 entry and are left out' % (wt[1], wt[1] / sum(wt.values())) if 1 in wt else ''))
OUT['section41'] = {str(k): v for k, v in S41.items()}
OUT['constants'] = {str(k): v for k, v in CONST.items()}

# ================================================================== (2) Table 12
say('')
say('=' * 100)
say('(2) Table 12 (Section 4.2): CP, one structure, 300 trials')
say('=' * 100)
MCP = parse_cp()
rows = []
hdr = '%-5s %-18s %-9s %-8s %-8s | %-7s %-7s %-7s | %-7s %-6s | %-6s %-6s' % (
    'n_z', 'grid', 'log2N', 'mean_pr', 'mean_ms', 'P_cpois', 'P_geom', 'Poisson', 'P_meas', 'SE', 'z_cp', 'z_Poi')
say(hdr)
for (nz, sz), md in sorted(MCP.items(), key=lambda t: (-t[0][0], t[0][1])):
    if nz == 2:
        q = q_A(*sz); key = None
    else:
        q = q_prod(*sz); key = '2p27.02' if sz[0] == 11664 else '2p27.20'
    E = T[8]['nz%d' % nz]['E']
    pc = Pf(lam(nz, q))
    hg = hgeom('a' if nz == 1 else 'b', key) if nz != 2 else None
    pg = Pf(lam(nz, q, h=hg)) if hg is not None else pc
    pp = -math.expm1(-E * q)
    pm = md['hits'] / md['n']
    r = dict(nz=nz, grid=sz, N=int(math.prod(sz)), log2N=math.log2(math.prod(sz)), q=q, log2q=math.log2(q),
             kq=K[8][nz] * q, mean_pred=E * q, mean_meas=md['mean'], P_cpois=pc, P_cpois_geom=pg, P_poisson=pp,
             P_meas=pm, hits=md['hits'], trials=md['n'], SE_meas=binse(pm, md['n']), z_cpois=z(pm, pc, md['n']),
             z_cpois_geom=z(pm, pg, md['n']), z_poisson=z(pm, pp, md['n']),
             cond_mean_members=(K[8][nz] * q) / -math.expm1(-K[8][nz] * q))
    rows.append(r)
    say('%-5d %-18s %-9.4f %-8.4f %-8.4f | %-7.4f %-7.4f %-7.4f | %-7.4f %-6.4f | %+6.2f %+6.2f' % (
        nz, str(sz), r['log2N'], r['mean_pred'], r['mean_meas'], pc, pg, pp, pm, r['SE_meas'], r['z_cpois'], r['z_poisson']))
OUT['table12'] = rows
say('  (P_cpois: compound-Poisson model of Section 4.2 with h = 1-e^{-kq}; P_geom: h from the exact bundle geometry on the product grid '
    '(n_z=2 on the A0 x A1 grid has no geometry entry, so P_geom = P_cpois there);')
say('   z against the model/Poisson value with the binomial SE of the prediction; E[members present | bundle hit] = kq/(1-e^{-kq}) = ' +
    ', '.join('%.2f' % r['cond_mean_members'] for r in rows) + ')')

# Yan CP
say('')
say("Yan et al.'s CP (A0 x A1, s = 2^13.6, 1000 trials, 614 successes):")
YAN = []
for s in (12417, 12418):
    q = q_A(s, s); E = T[8]['nz1']['E']
    pc = Pf(lam(1, q)); hg = hgeom('a', 's%d' % s, 'A0xA1_grid'); pg = Pf(lam(1, q, h=hg)); pp = -math.expm1(-E * q)
    r = dict(s=s, N=s * s, log2N=math.log2(s * s), q=q, log2q=math.log2(q), kq=512 * q, mean_pred=E * q,
             P_cpois=pc, P_cpois_geom=pg, P_poisson=pp, P_meas=0.614, trials=1000, SE_meas=binse(0.614, 1000),
             z_cpois=z(0.614, pc, 1000), z_poisson=z(0.614, pp, 1000))
    YAN.append(r)
    say('  s=%d N=2^%.4f q=2^%.4f kq=%.4f mean=%.4f  P_cpois=%.4f P_geom=%.4f Poisson=%.4f | meas 0.614 (SE %.4f) z_cp=%+.2f z_Poisson=%+.2f' % (
        s, r['log2N'], r['log2q'], r['kq'], r['mean_pred'], pc, pg, pp, r['SE_meas'], r['z_cpois'], r['z_poisson']))
# N at which n_z=1 reaches 63.2% on Yan's grid and on the product grid
def solve(fun, target, lo=24.0, hi=32.0):
    if fun(hi) < target:
        return None
    for _ in range(80):
        mid = (lo + hi) / 2
        if fun(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi
yanN = solve(lambda x: Pf(lam(1, q_A(2 ** (x / 2), 2 ** (x / 2)))), 1 - math.exp(-1))
say('  Yan grid: n_z=1 reaches 1-1/e = 63.2%% at N = 2^%.3f (compound Poisson); Poisson would give 2^%.3f' % (
    yanN, solve(lambda x: -math.expm1(-T[8]['nz1']['E'] * q_A(2 ** (x / 2), 2 ** (x / 2))), 1 - math.exp(-1))))
OUT['yan_cp'] = YAN
OUT['yan_cp_N_63'] = yanN

# ================================================================== (2b) Table 14 (Section 5.2) / Yan ACP
say('')
say('=' * 100)
say('(2b) Table 14 (Section 5.2): ACP at 2^23.32 (and the unreported 2^22 block), 300 trials')
say('=' * 100)
MACP = parse_acp()
T14 = []
for b in ACP['budget']:
    nz, Dx = b['nz'], b['Dexp']
    md = MACP[(nz, Dx)]; pm = md['hits'] / md['n']
    pp = -math.expm1(-b['mean_valid'])
    r = dict(nz=nz, Dexp=Dx, P_cpois=b['P_find'], mean_pred=b['mean_valid'], P_poisson_at_mean=pp, E_structures=b['E_structures'],
             E_partner_pairs=b['E_Np'], P_meas=pm, mean_meas=md['mean'], SE_meas=binse(pm, md['n']),
             z_cpois=z(pm, b['P_find'], md['n']))
    T14.append(r)


def acp_nz2(D):
    reserve = int(D * D / 1073741824.0)
    M = D - reserve if D > reserve else D // 2
    ps = -math.expm1(M * math.log1p(-2.0 ** -32))
    q = 2 * ps * ps - ps ** 4
    k = 65536
    h = -math.expm1(k * math.log1p(-q))
    L = lam(2, q, h=h)
    coll = M * M / 2 * 2.0 ** -30 * (255 / 256) ** 4
    return dict(M=M, log2M=math.log2(M), p_sample=ps, q=q, log2q=math.log2(q), kq=k * q, h=h, P_cpois=-math.expm1(-L),
                mean_pred=T[8]['nz2']['E'] * q, swap_queries=2 * coll, reserve=reserve)


for Dx in (22.0, 23.32):
    D = int(2.0 ** Dx + 0.5)
    a2 = acp_nz2(D); md = MACP[(2, Dx)]; pm = md['hits'] / md['n']
    r = dict(nz=2, Dexp=Dx, P_cpois=a2['P_cpois'], mean_pred=a2['mean_pred'], P_poisson_at_mean=-math.expm1(-a2['mean_pred']),
             P_meas=pm, mean_meas=md['mean'], SE_meas=binse(pm, md['n']), z_cpois=z(pm, a2['P_cpois'], md['n']), detail=a2)
    T14.append(r)
say('%-4s %-6s | %-8s %-8s %-8s | %-8s %-8s %-6s %-6s' % ('n_z', 'log2D', 'P_cpois', 'Poi@mean', 'mean_pr', 'P_meas', 'mean_ms', 'SE', 'z'))
for r in sorted(T14, key=lambda t: (t['Dexp'], -t['nz'])):
    say('%-4d %-6.2f | %-8.4f %-8.4f %-8.4f | %-8.4f %-8.4f %-6.4f %+6.2f' % (r['nz'], r['Dexp'], r['P_cpois'], r['P_poisson_at_mean'],
                                                                          r['mean_pred'], r['P_meas'], r['mean_meas'], r['SE_meas'], r['z_cpois']))
say('  (n_z=1,0: budget-limited procedure of acp_all.c simulated; n_z=2: M = D - D^2/2^30 i.i.d. 4-active plaintexts of one coset,')
say('   q = 2 p^2 - p^4 with p = 1-(1-2^-32)^M; at 2^23.32: kq = %.3f, h = %.3f)' % (T14[-1]['detail']['kq'], T14[-1]['detail']['h']))
d2 = T14[-1]['detail']
say('   n_z=2 at 2^23.32: M = 2^%.3f sampled plaintexts, D^2/2^30 = 2^%.2f queries reserved for adaptive swaps, expected swaps 2^%.2f (Section 5.2)' % (
    d2['log2M'], math.log2(d2['reserve']), math.log2(d2['swap_queries'])))
BUD = {}
for nz, ncoll in ((0, 255 / 3),):
    nst = ncoll / 2
    qry = ncoll * 255 ** 2 * 2
    BUD[nz] = dict(collisions=ncoll, structures=nst, log2_structures=math.log2(nst), log2_base=math.log2(nst * 2 ** 16),
                   log2_queries=math.log2(qry), log2_total=math.log2(nst * 2 ** 16 + qry))
    say('  Section 5.1 budget, n_z=%d: %.1f base collisions, about 2 per base structure, need %.1f = 2^%.2f base structures (2^%.2f plaintexts); %.1f x (2^8-1)^2 x 2 = 2^%.2f partner queries; total 2^%.2f' % (
        nz, ncoll, nst, BUD[nz]['log2_structures'], BUD[nz]['log2_base'], ncoll, BUD[nz]['log2_queries'], BUD[nz]['log2_total']))
OUT['acp_budget_51'] = {str(k): v for k, v in BUD.items()}
OUT['table14_acp'] = T14

say('')
say("Yan et al.'s ACP (fixed number of base structures in distinct cosets, all (2^8-1)^2 partners queried):")
YA = []
for f in ACP['fixed']:
    if f['nz'] != 1:
        continue
    meas = {32: (0.633, 2000)}.get(f['n_struct'])
    if abs(f['n_struct'] - 2 ** 5.75) < 1e-6 or f['n_struct'] in (53, 54):
        meas = (0.825, 2000)
    r = dict(n_struct=f['n_struct'], log2_data=f['log2_data'], P_cpois=f['P_find'], mean=f['mean_valid'])
    if meas:
        r.update(P_meas=meas[0], trials=meas[1], SE_meas=binse(meas[0], meas[1]), z=z(meas[0], f['P_find'], meas[1]))
    YA.append(r)
    say('  %6.2f structures (data 2^%.3f): P_cpois = %.4f  mean %.4f' % (f['n_struct'], f['log2_data'], f['P_find'], f['mean_valid']) +
        ('  | measured %.3f (SE %.4f, %d trials) z = %+.2f' % (meas[0], binse(meas[0], meas[1]), meas[1], r['z']) if meas else ''))
for f in ACP['fixed']:
    if f['nz'] == 0 and f['n_struct'] == 32:
        say('  n_z=0 with 32 structures: P_cpois = %.4f (mean %.4f)' % (f['P_find'], f['mean_valid']))
for o in ACP['one_coset']:
    say('  n_z=1 with all %d base structures in ONE coset of D_0: h = %.4f per bundle, P_find = %.4f (distinct cosets: see above)' % (o['n_struct'], o['h'], o['P_find']))
OUT['yan_acp'] = YA
OUT['acp_one_coset'] = ACP['one_coset']

# ================================================================== (3) variance
say('')
say('=' * 100)
say('(3) Var(N_q) = k T1 + sum_m g_m eps_m  (~ k T1 + 8 eps) over a full structure')
say('=' * 100)
VAR = []
small = {2: freq('5r_nz2.txt'), 1: freq('5r_nz1.txt'), 0: freq('5r_nz0.txt')}
for w in (4, 8):
    for nz in (2, 1, 0):
        k = K[w][nz]; t1 = T1(nz, w); e = EPS(nz, w); em = eps_m(w, nz)
        v_case2 = sum(gsize(w, m) * x for m, x in em.items())
        var_split = k * t1 + v_case2; var_8 = k * t1 + 8 * e
        if w == 4:
            fr = small[nz]
            x = np.repeat(np.array(list(fr.keys()), dtype=float), list(fr.values()))
        else:
            x = raw('raw_nz%d_5r.txt' % nz)
        sd = float(np.std(x, ddof=1)); sdse = boot_sd(x)
        mu_b = t1 / k
        psat = Pf(mu_b, sum(v / gsize(w, m) for m, v in em.items()))
        r = dict(w=w, nz=nz, k=k, T1=t1, eps=e, kT1=k * t1, var_case2_split=v_case2, var_case2_8eps=8 * e,
                 SD_model=math.sqrt(var_split), SD_model_8eps=math.sqrt(var_8), SD_meas=sd, SD_meas_bootSE=sdse,
                 z_SD=(sd - math.sqrt(var_split)) / sdse, n=len(x), mean_meas=float(x.mean()),
                 mean_model=T[w]['nz%d' % nz]['E'], P_ge1_model=psat, P_ge1_bundles_only=-math.expm1(-mu_b),
                 P_ge1_meas=float((x > 0).mean()), P_ge1_meas_SE=binse(float((x > 0).mean()), len(x)))
        VAR.append(r)
        say('w=%d n_z=%d n=%4d  kT1=%10.1f  case2: split %.1f / 8eps %.1f   SD model %.1f (8eps: %.1f)  SD meas %.1f +- %.1f  (z %+.2f) | Pr[N>=1] model %.4f (bundles only %.4f) meas %.4f +- %.4f' % (
            w, nz, len(x), k * t1, v_case2, 8 * e, r['SD_model'], r['SD_model_8eps'], sd, sdse, r['z_SD'], psat, r['P_ge1_bundles_only'],
            r['P_ge1_meas'], r['P_ge1_meas_SE']))
OUT['variance'] = VAR
# dispersion of the bundle count (AES n_z=1: the 8B parts, i.e. the residues mod 512 printed below, are at most 24 < 512,
# so floor(N/512) is the bundle count)
x1 = raw('raw_nz1_5r.txt'); A1 = np.floor(x1 / 512)
disp = dict(n=len(A1), mean=float(A1.mean()), var=float(A1.var(ddof=1)), index=float(A1.var(ddof=1) / A1.mean()),
            mu_b_model=CONST[1]['mu_b'], residues=sorted(set((x1 % 512).astype(int).tolist())))
# small-AES n_z=2: floor(N/256) (one Case-2 group of 512 would be counted as 2 bundles; eps_2 is small)
fr = small[2]; xs2 = np.repeat(np.array(list(fr.keys()), dtype=float), list(fr.values())); A2 = np.floor(xs2 / 256)
disp2 = dict(n=len(A2), mean=float(A2.mean()), var=float(A2.var(ddof=1)), index=float(A2.var(ddof=1) / A2.mean()),
             mu_b_model=T1(2, 4) / 256)
# chi-square-type check of the dispersion index: SE ~ sqrt(2/(n-1))
say('bundle-count dispersion: AES n_z=1 floor(N/512): mean %.3f (mu_b %.3f), var %.3f, index %.3f (Poisson 1, SE ~ %.3f); residues mod 512: %s' % (
    disp['mean'], disp['mu_b_model'], disp['var'], disp['index'], math.sqrt(2 / (disp['n'] - 1)), disp['residues']))
say('                         small-AES n_z=2 floor(N/256): mean %.3f (mu_b %.3f), var %.3f, index %.3f (SE ~ %.3f)' % (
    disp2['mean'], disp2['mu_b_model'], disp2['var'], disp2['index'], math.sqrt(2 / (disp2['n'] - 1))))
# bundle totals over the 700 AES structures (n_z=2: the residues mod 65536 are at most 16, far below 65536,
# so floor(N/65536) is the bundle count)
x2 = raw('raw_nz2_5r.txt'); A2a = np.floor(x2 / 65536)
disp['total'] = int(A1.sum()); disp['expected_total'] = len(A1) * CONST[1]['mu_b']
disp_nz2 = dict(n=len(A2a), total=int(A2a.sum()), expected_total=len(A2a) * CONST[2]['mu_b'])
say('bundles in the %d AES structures (Section 3.4.2): n_z=1 %d (expected %d mu_b = %.1f), n_z=2 %d (expected %.1f)' % (
    len(A1), disp['total'], len(A1), disp['expected_total'], disp_nz2['total'], disp_nz2['expected_total']))
OUT['dispersion'] = dict(aes_nz1=disp, aes_nz2=disp_nz2, small_nz2=disp2)

# ================================================================== (4) product-grid configurations (Table 13, Section 4.3)
say('')
say('=' * 100)
say('(4) Table 13 (Section 4.3) model columns: product-grid configurations of the pairings a (n_z=1 on {0,10}, Yan et al.), '
    'b (n_z=0), c (n_z=1 on {5,15}) and their unions')
say('=' * 100)
# grids of the Table 13 runs, in ascending N; the keys 2p27.6 ... 2p28.3 index data/bundle_geometry.json,
# which has no entry for the three smallest grids (their geom column is printed as '-')
GR = {}
for x in (26.05, 26.33, 26.71):
    s1 = int(round(256 * 2 ** ((x - 32) / 4))); s0 = int(round(2 ** x / s1 ** 2))
    GR['2p%.2f' % x] = (s0, s1, s1)
GR['2p27.02'] = (11664, 108, 108)
GR['2p27.20'] = (12515, 111, 111)
for x in (27.6, 27.8, 28.0, 28.3):
    s1 = int(round(256 * 2 ** ((x - 32) / 4))); s0 = int(round(2 ** x / s1 ** 2))
    GR['2p%.1f' % x] = (s0, s1, s1)


def ncand(s0, s1, s2, pairing):
    w_pairs = s1 * (s1 - 1) * s2 * (s2 - 1) / 2
    a_pairs = s0 * (s0 - 1) / 2
    if pairing in ('a', 'c'):
        a_pairs *= 255 * 255 / 65535
    return a_pairs * w_pairs


CONF = []
say('grid rule: s1 = s2 = round(256 (N/2^32)^{1/4}), s0 = round(N / s1^2) for N = 2^26.05, 2^26.33, 2^26.71, 2^27.6, 2^27.8, 2^28.0, 2^28.3;'
    ' 2^27.02 and 2^27.20 are the grids of Table 12')
say('%-8s %-18s %-8s | %-14s %-14s %-14s | %-7s %-7s %-7s %-7s | %-9s' % ('N', 'grid', 'log2q', 'P_a (geom)', 'P_b (geom)', 'P_c (geom)',
                                                                     'a|b', 'a|c', 'b|c', 'a|b|c', 'typeI abc'))
for gk, sz in GR.items():
    q = q_prod(*sz)
    La = lam(1, q); Lb = lam(0, q); Lc = lam(1, q, pairing='c')
    ha = hgeom('a', gk); hb = hgeom('b', gk); hc = hgeom('c', gk)
    geom = None not in (ha, hb, hc)
    if geom:
        Lag = lam(1, q, h=ha); Lbg = lam(0, q, h=hb); Lcg = lam(1, q, h=hc, pairing='c')
    fa = 4 * 2.0 ** -64 * ncand(*sz, 'a'); fb = 4 * 2.0 ** -64 * ncand(*sz, 'b'); fc = fa
    r = dict(grid_key=gk, grid=sz, N=int(math.prod(sz)), log2N=math.log2(math.prod(sz)), q=q, log2q=math.log2(q),
             mean_a=T[8]['nz1']['E'] * q, mean_b=T[8]['nz0']['E'] * q,
             P_a=Pf(La), P_b=Pf(Lb), P_c=Pf(Lc), P_ab=Pf(La, Lb), P_ac=Pf(La, Lc), P_bc=Pf(Lb, Lc), P_abc=Pf(La, Lb, Lc),
             P_a_geom=Pf(Lag) if geom else None, P_b_geom=Pf(Lbg) if geom else None, P_c_geom=Pf(Lcg) if geom else None,
             P_ab_geom=Pf(Lag, Lbg) if geom else None, P_ac_geom=Pf(Lag, Lcg) if geom else None,
             P_bc_geom=Pf(Lbg, Lcg) if geom else None, P_abc_geom=Pf(Lag, Lbg, Lcg) if geom else None,
             P_a_poisson=-math.expm1(-T[8]['nz1']['E'] * q), P_b_poisson=-math.expm1(-T[8]['nz0']['E'] * q),
             typeI_a=fa, typeI_b=fb, typeI_c=fc, typeI_abc=-math.expm1(-(fa + fb + fc)),
             log2_typeI_a=math.log2(fa), log2_typeI_b=math.log2(fb), memory_GB=math.prod(sz) * 16 / 1e9,
             log2_time=math.log2(1.8 * math.prod(sz)))
    CONF.append(r)
    gs = lambda v: '%.4f' % v if v is not None else '  -   '
    say('2^%-6.2f %-18s %-8.3f | %.4f (%s) %.4f (%s) %.4f (%s) | %.4f %.4f %.4f %.4f | 2^%.3f' % (
        r['log2N'], str(sz), r['log2q'], r['P_a'], gs(r['P_a_geom']), r['P_b'], gs(r['P_b_geom']), r['P_c'], gs(r['P_c_geom']),
        r['P_ab'], r['P_ac'], r['P_bc'], r['P_abc'], math.log2(r['typeI_abc'])))
    say('         means: a %.4f  b %.4f ; Poisson P: a %.4f  b %.4f ; type-I per pairing: a 2^%.3f  b 2^%.3f  c 2^%.3f ; ciphertext memory %.2f GB' % (
        r['mean_a'], r['mean_b'], r['P_a_poisson'], r['P_b_poisson'], r['log2_typeI_a'], r['log2_typeI_b'], math.log2(r['typeI_c']), r['memory_GB']))
OUT['configs'] = CONF


def q_cont(x):
    """q on the proportional grid with real-valued sizes: s0 = 2^16 (N/2^32)^{1/2}, s1 = s2 = 2^8 (N/2^32)^{1/4}."""
    s0 = 65536 * 2 ** ((x - 32) / 2); s1 = 256 * 2 ** ((x - 32) / 4)
    return frac2(s0, 65536) * frac2(s1, 256) ** 2


TGT = {}
funcs = {
    'a': lambda x: Pf(lam(1, q_cont(x))),
    'b': lambda x: Pf(lam(0, q_cont(x))),
    'a|c': lambda x: Pf(lam(1, q_cont(x)), lam(1, q_cont(x), pairing='c')),
    'a|b': lambda x: Pf(lam(1, q_cont(x)), lam(0, q_cont(x))),
    'a|b|c': lambda x: Pf(lam(1, q_cont(x)), lam(0, q_cont(x)), lam(1, q_cont(x), pairing='c')),
    'a Poisson': lambda x: -math.expm1(-T[8]['nz1']['E'] * q_cont(x)),
}
say('')
say('data N (proportional grid, one structure) for a target P_find:')
for name, f in funcs.items():
    TGT[name] = {}
    row = []
    for tg in (0.5, 1 - math.exp(-1), 0.8, 0.85, 0.9, 0.95):
        x = solve(f, tg)
        TGT[name]['%.4f' % tg] = x
        row.append('%.3f: %s' % (tg, ('2^%.3f' % x) if x else 'unreachable (max %.4f)' % f(32.0)))
    say('  %-10s ' % name + ' | '.join(row) + '   [P at 2^32 = %.4f]' % f(32.0))
OUT['target_N'] = TGT

# ================================================================== (4b) Table 13: model vs measured
say('')
say('=' * 100)
say('(4b) Table 13 (Section 4.3): model vs measured, 300 trials per size')
say('=' * 100)
AGG = aggregate_combined.aggregate()
say('measured: Results/5r_AES_distinguishers/combined, read by Codes/tools/aggregate_combined.py (hits / 300);')
say('model: section (4); z = (meas - model) / sqrt(model (1 - model) / 300), the binomial SE of the prediction;')
say('false positive = 1 - exp(-4 (Q_a + Q_b + Q_c) / 2^64) for the union; time = N_s + 4 N_s / 5 = 1.8 N_s encryptions.')
say('columns in the order of Table 13: a = n_z=1 on {0,10} (Yan et al.), c = n_z=1 on {5,15}, b = n_z=0, union a|b|c')
say('%-8s %-18s | %-21s | %-21s | %-21s | %-21s | %-9s %-8s' % ('N_s', 'grid', 'a meas (model) z', 'c meas (model) z',
                                                                'b meas (model) z', 'union meas (model) z', 'false pos', 'time'))
T13 = []
for r in CONF:
    g = [v for v in AGG.values() if v['rounds'] == 5 and v['grid'] == tuple(r['grid'])]
    assert len(g) == 1, 'no measured run for %s' % (r['grid'],)
    g = g[0]; n = g['n']
    row = dict(grid=r['grid'], log2N=r['log2N'], group=g['group'], trials=n, log2_false_positive=math.log2(r['typeI_abc']),
               log2_time=r['log2_time'])
    cells = []
    for p, key in (('a', 'P_a'), ('c', 'P_c'), ('b', 'P_b'), ('a|b|c', 'P_abc')):
        pm = g['hits'][p] / n; pr = r[key]; zz = z(pm, pr, n)
        row[p] = dict(hits=g['hits'][p], P_meas=pm, P_model=pr, z=zz)
        cells.append('%.4f (%.4f) %+5.2f' % (pm, pr, zz))
    ca, cb, cc = (g['counts'][p] for p in 'abc')
    both = lambda x, y: sum(1 for i in range(n) if x[i] > 0 and y[i] > 0)
    row['both'] = {'a&c': both(ca, cc), 'a&b': both(ca, cb), 'b&c': both(cb, cc)}
    row['both_indep'] = {'a&c': g['hits']['a'] * g['hits']['c'] / n, 'a&b': g['hits']['a'] * g['hits']['b'] / n,
                         'b&c': g['hits']['b'] * g['hits']['c'] / n}
    T13.append(row)
    say('2^%-6.2f %-18s | %s | %s | %s | %s | 2^%-7.3f 2^%.3f' % ((r['log2N'], str(tuple(r['grid']))) + tuple(cells) +
                                                                (row['log2_false_positive'], r['log2_time'])))
    say('         both succeed: a&c %d, a&b %d, b&c %d ; under independence of the measured marginals %.1f, %.1f, %.1f' % (
        row['both']['a&c'], row['both']['a&b'], row['both']['b&c'], row['both_indep']['a&c'], row['both_indep']['a&b'],
        row['both_indep']['b&c']))
big = [t for t in T13 if t['log2N'] > 27]; sml = [t for t in T13 if t['log2N'] < 27]
zmax = lambda rows, ps: max(abs(t[p]['z']) for t in rows for p in ps)
ZM = dict(nz1_large=zmax(big, 'ac'), nz0_large=zmax(big, 'b'), union_large=zmax(big, ['a|b|c']),
          small_pairings=zmax(sml, 'acb'), small_all=zmax(sml, ['a', 'c', 'b', 'a|b|c']))
say('largest |z|: 2^27.02 .. 2^28.30: n_z=1 (a, c) %.2f, n_z=0 (b) %.2f, union %.2f ; 2^26.05 .. 2^26.71: a, c, b %.2f, with the union %.2f' % (
    ZM['nz1_large'], ZM['nz0_large'], ZM['union_large'], ZM['small_pairings'], ZM['small_all']))
r10 = [v for v in AGG.values() if v['rounds'] == 10]
for g in r10:
    cf = [r for r in CONF if tuple(r['grid']) == g['grid']][0]
    say('10 rounds (%s, %s, N = 2^%.2f): union hits %d/%d ; expected false positives %d x 2^%.3f = %.2f' % (
        g['group'], g['grid'], g['log2N'], g['hits']['a|b|c'], g['n'], g['n'], math.log2(cf['typeI_abc']), g['n'] * cf['typeI_abc']))
    ZM['r10_' + g['group']] = dict(union_hits=g['hits']['a|b|c'], trials=g['n'], expected_false_positives=g['n'] * cf['typeI_abc'])
OUT['table13'] = dict(rows=T13, largest_abs_z=ZM)

# ================================================================== (5) saturation
say('')
say('=' * 100)
say('(5) saturation of n_z=1 in one full diagonal structure')
say('=' * 100)
v8 = [r for r in VAR if r['w'] == 8 and r['nz'] == 1][0]; v4 = [r for r in VAR if r['w'] == 4 and r['nz'] == 1][0]
say('  AES:       1 - exp(-mu_b - sum eps_m/g_m) = 1 - exp(-%.4f - %.4f) = %.4f ; bundles only 1 - e^{-%.4f} = %.4f ; measured %.4f +- %.4f (700 structures)' % (
    CONST[1]['mu_b'], sum(v / gsize(8, m) for m, v in eps_m(8, 1).items()), v8['P_ge1_model'], CONST[1]['mu_b'], v8['P_ge1_bundles_only'],
    v8['P_ge1_meas'], v8['P_ge1_meas_SE']))
say('  small-AES: %.4f (bundles only %.4f) ; measured %.4f +- %.4f (1000 structures)' % (v4['P_ge1_model'], v4['P_ge1_bundles_only'], v4['P_ge1_meas'], v4['P_ge1_meas_SE']))
OUT['saturation'] = dict(aes=dict(model=v8['P_ge1_model'], bundles_only=v8['P_ge1_bundles_only'], meas=v8['P_ge1_meas'], se=v8['P_ge1_meas_SE']),
                         small=dict(model=v4['P_ge1_model'], bundles_only=v4['P_ge1_bundles_only'], meas=v4['P_ge1_meas'], se=v4['P_ge1_meas_SE']))

# ================================================================== (6) ACP n_z=2 data
say('')
say('=' * 100)
say('(6) ACP n_z=2 data requirement')
say('=' * 100)
mu2 = CONST[2]['mu_b']


def acp2_coset(M):
    ps = M / 2 ** 32
    q = 2 * ps * ps - ps ** 4
    h = -math.expm1(-65536 * q)
    L = lam(2, q, h=h)
    coll = M * M / 2 * 2.0 ** -30 * (255 / 256) ** 4
    return L, M + 2 * coll, q, coll


L1, D1, q1, c1 = acp2_coset(2 ** 26.5)
say('  one structure of 2^26.5 random 4-active plaintexts of one coset: q = 2^%.3f, kq = %.1f, mean = %.4f, P_find = %.4f ; '
    'data incl. the %.3g adaptive swap queries = 2^%.3f' % (math.log2(q1), 65536 * q1, T[8]['nz2']['E'] * q1, -math.expm1(-L1), 2 * c1, math.log2(D1)))
say('  bundles alone cap the single-coset success at 1 - e^{-mu_b} = %.4f; a FULL coset (2^32) also catches the Case-2 groups of 8: P = %.4f (Table 10: 0.104)' % (-math.expm1(-mu2), Pf(lam(2, 1.0))))
for xx in (27.0, 28.0, 29.0, 30.0, 31.0):
    Lx, Dx, qx, cx = acp2_coset(2 ** xx)
    say('    single coset, M = 2^%.1f: P_find = %.4f (bundle part %.4f), data 2^%.3f' % (xx, -math.expm1(-Lx), mu2 * -math.expm1(-65536 * qx), math.log2(Dx)))
ACP6 = dict(single_2p26p5=dict(q=q1, kq=65536 * q1, mean=T[8]['nz2']['E'] * q1, P_find=-math.expm1(-L1), data=D1, log2_data=math.log2(D1)),
            single_full=Pf(lam(2, 1.0)))
for tg in (0.5, 1 - math.exp(-1), 0.8, 0.9):
    need = -math.log1p(-tg)
    best = None
    for x in np.arange(18.0, 30.0, 0.005):
        L, Dc, q, c = acp2_coset(2 ** x)
        n = need / L
        tot = n * Dc
        if best is None or tot < best[0]:
            best = (tot, x, n, L)
    ACP6['acp_multicoset_%.4f' % tg] = dict(total=best[0], log2_total=math.log2(best[0]), log2_M=best[1], cosets=best[2])
    say('  ACP, many cosets, target %.3f: total data 2^%.3f with %.1f cosets of M = 2^%.2f sampled plaintexts (+ swaps)' % (tg, math.log2(best[0]), best[2], best[1]))
    # CP analogue: A0 x A1 product grid of N per coset (partners inside the grid), no adaptive queries
    bestc = None
    for x in np.arange(18.0, 32.0, 0.005):
        s = 2 ** (x / 2)
        q = q_A(s, s)
        L = lam(2, q)
        n = need / L
        if bestc is None or n * 2 ** x < bestc[0]:
            bestc = (n * 2 ** x, x, n)
    ACP6['cp_multicoset_%.4f' % tg] = dict(total=bestc[0], log2_total=math.log2(bestc[0]), log2_N=bestc[1], structures=bestc[2])
    say('  CP  (A0xA1 grid per coset), target %.3f: total data 2^%.3f with %.1f structures of 2^%.2f' % (tg, math.log2(bestc[0]), bestc[2], bestc[1]))
xb = ACP6['acp_multicoset_%.4f' % (1 - math.exp(-1))]['log2_M']
Lb, Db, qb, cb = acp2_coset(2 ** xb)
hb = -math.expm1(-65536 * qb)
ACP6['acp_optimum'] = dict(log2_M=xb, q=qb, kq=65536 * qb, h=hb, L=Lb, inv_L=1 / Lb)
say('  at the optimum M = 2^%.3f (the same for every target): q = %.3e, 2^16 q = %.3f, h = 1 - e^{-2^16 q} = %.3f, per-structure L = mu_b h + sum_m (eps_m/g_m)(1 - e^{-g_m q}) = %.5f = 1/%.2f' % (
    xb, qb, 65536 * qb, hb, Lb, 1 / Lb))
q27 = q_A(11664, 11664)
n63 = 1 / lam(2, q27)
say('  CP with repeated 2^27.02 structures (distinct cosets): P per structure %.4f, %.1f structures = 2^%.2f for 63.2%%' % (
    Pf(lam(2, q27)), n63, math.log2(n63 * 11664 * 11664)))
ACP6['cp_repeat_2p27'] = dict(P_per_structure=Pf(lam(2, q27)), structures_63=n63, log2_total=math.log2(n63 * 11664 * 11664))
OUT['acp_nz2'] = ACP6

json.dump(OUT, open(os.path.join(DATA, 'cpoisson_results.json'), 'w'), indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))
open(os.path.join(DATA, 'cpoisson_results.txt'), 'w').write('\n'.join(LINES) + '\n')
