/* Exhaustive search of byte-related differentials over MC = circ(2,3,1,1) in GF(2^w)^4.
 * Verifies Lemma 7 / Table 2 (DCC paper) for w = 4..8.
 *
 * usage: lemma7_search w poly mode [raw]
 *   mode 0 = "exactly one of u,v,u+v zero at every byte" on input AND output (the form of
 *            Table 2; stricter than Definition 1)
 *   mode 1 = "at least one zero" (paper Definition 1, admits all-zero bytes) on input AND output
 *   raw   = 1: additionally enumerate WITHOUT scalar normalization (cross-check, small w only)
 *
 * Input parameterization per byte i: type t_i
 *   0: u_i=0, v_i=eta      1: v_i=0, u_i=eta      2: u_i=v_i=eta      3 (mode 1 only): u_i=v_i=0
 * eta_i != 0. Normalization: eta of the first byte with type != 3 is 1 (every projective
 * class {tau*(u,v)} has exactly one representative, since tau acts freely on eta).
 * Patterns where u=0, v=0 or u=v are excluded.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int W, Q, POLY, MODE;
static unsigned char mul[256][256], inv_[256];
static const int MC[4][4] = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};

static int gmul_slow(int a, int b) {
    int r = 0;
    while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a & Q) a ^= POLY; }
    return r;
}
static void mc(const unsigned char *x, unsigned char *y) {
    for (int j = 0; j < 4; j++) {
        int s = 0;
        for (int k = 0; k < 4; k++) s ^= mul[MC[j][k]][x[k]];
        y[j] = (unsigned char)s;
    }
}
/* per-byte related condition */
static int byte_ok(int a, int b) {
    int c = a ^ b;
    int nz = (a != 0) + (b != 0) + (c != 0);
    if (MODE == 0) return nz == 2;          /* exactly one zero (two zero => all zero) */
    return nz <= 2;                          /* at least one zero */
}
static int out_ok(const unsigned char *u, const unsigned char *v) {
    unsigned char mu[4], mv[4];
    mc(u, mu); mc(v, mv);
    for (int j = 0; j < 4; j++) if (!byte_ok(mu[j], mv[j])) return 0;
    return 1;
}
static int zmask(const unsigned char *x) { int m = 0; for (int i = 0; i < 4; i++) if (!x[i]) m |= 1 << i; return m; }
static int wt(const unsigned char *x) { int c = 0; for (int i = 0; i < 4; i++) c += x[i] != 0; return c; }

/* ---------------- solution storage ---------------- */
typedef struct { unsigned char u[4], v[4]; } Sol;
static Sol *sols; static int nsols, capsols;
static void add_sol(const unsigned char *u, const unsigned char *v) {
    if (nsols == capsols) { capsols = capsols ? 2 * capsols : 1024; sols = realloc(sols, capsols * sizeof(Sol)); }
    memcpy(sols[nsols].u, u, 4); memcpy(sols[nsols].v, v, 4); nsols++;
}
/* normalize: scale so that eta at first nonzero position of (u or v) is 1 */
static void normalize(const unsigned char *u, const unsigned char *v, unsigned char *nu, unsigned char *nv) {
    int e = 0;
    for (int i = 0; i < 4 && !e; i++) { if (u[i]) e = u[i]; else if (v[i]) e = v[i]; }
    int t = inv_[e];
    for (int i = 0; i < 4; i++) { nu[i] = mul[t][u[i]]; nv[i] = mul[t][v[i]]; }
}
static int sol_eq(const unsigned char *a, const unsigned char *b, const unsigned char *c, const unsigned char *d) {
    return !memcmp(a, c, 4) && !memcmp(b, d, 4);
}

/* ---------------- reference set from Table 2 ---------------- */
static const int ROWS[4][2][4] = {
    {{0,1,4,7},{5,1,0,7}},
    {{2,1,1,3},{0,1,0,3}},
    {{7,0,7,7},{7,7,7,0}},
    {{2,3,2,3},{0,3,2,0}}};
static const int ROWS_Y[4][2][4] = {
    {{0x0,0x9,0x0,0xB},{0xE,0x0,0xD,0x0}},
    {{0x5,0x0,0x4,0x0},{0x0,0x1,0x4,0x7}},
    {{0xE,0x9,0x0,0x0},{0x0,0x0,0xE,0x9}},
    {{0x0,0x1,0x0,0x1},{0x7,0x0,0x7,0x1}}};
typedef struct { Sol s; int row, rot, ord; } Ref;
static Ref refs[4 * 4 * 6]; static int nrefs;

static void build_refs(void) {
    nrefs = 0;
    for (int r = 0; r < 4; r++) for (int k = 0; k < 4; k++) {
        unsigned char A[4], B[4], C[4];
        for (int i = 0; i < 4; i++) { /* rot_k: byte at i moves to (i+k) mod 4 ; rot_1 [a,b,c,d]->[d,a,b,c] */
            A[(i + k) % 4] = (unsigned char)ROWS[r][0][i];
            B[(i + k) % 4] = (unsigned char)ROWS[r][1][i];
        }
        for (int i = 0; i < 4; i++) C[i] = A[i] ^ B[i];
        const unsigned char *P[6][2] = {{A,B},{B,A},{A,C},{C,A},{B,C},{C,B}};
        for (int o = 0; o < 6; o++) {
            unsigned char nu[4], nv[4];
            normalize(P[o][0], P[o][1], nu, nv);
            int dup = 0;
            for (int t = 0; t < nrefs; t++) if (sol_eq(refs[t].s.u, refs[t].s.v, nu, nv)) { dup = 1; break; }
            if (!dup) { memcpy(refs[nrefs].s.u, nu, 4); memcpy(refs[nrefs].s.v, nv, 4); refs[nrefs].row = r + 1; refs[nrefs].rot = k; refs[nrefs].ord = o; nrefs++; }
        }
    }
}
static int find_ref(const unsigned char *u, const unsigned char *v) {
    for (int t = 0; t < nrefs; t++) if (sol_eq(refs[t].s.u, refs[t].s.v, u, v)) return t;
    return -1;
}

/* class key: sorted triple of zero masks */
static int class_key(int a, int b, int c) {
    int t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    return (a << 8) | (b << 4) | c;
}
static void print_mask(FILE *f, int m) {
    fprintf(f, "{"); int first = 1;
    for (int i = 0; i < 4; i++) if (m >> i & 1) { fprintf(f, first ? "%d" : ",%d", i); first = 0; }
    fprintf(f, "}");
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s w poly mode [raw]\n", argv[0]); return 1; }
    W = atoi(argv[1]); POLY = (int)strtol(argv[2], NULL, 0); MODE = atoi(argv[3]);
    int RAW = argc > 4 ? atoi(argv[4]) : 0;
    Q = 1 << W;
    clock_t t0 = clock();

    /* field tables + sanity */
    for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) mul[a][b] = (unsigned char)gmul_slow(a, b);
    for (int a = 1; a < Q; a++) {
        int seen[256] = {0}, cnt = 0;
        for (int b = 0; b < Q; b++) { if (!seen[mul[a][b]]) cnt++; seen[mul[a][b]] = 1; if (mul[a][b] == 1) inv_[a] = (unsigned char)b; }
        if (cnt != Q) { printf("FIELD CHECK FAILED: poly 0x%X not irreducible (row %d)\n", POLY, a); return 2; }
    }
    printf("=== w=%d poly=0x%X mode=%d (%s) ===\n", W, POLY, MODE, MODE == 0 ? "exactly one zero per byte" : "at least one zero per byte");
    printf("field check: OK (every nonzero row of the mult table is a permutation)\n");

    /* Table 2 MC outputs */
    int tab_ok = 1;
    for (int r = 0; r < 4; r++) for (int s = 0; s < 2; s++) {
        unsigned char x[4], y[4];
        for (int i = 0; i < 4; i++) x[i] = (unsigned char)ROWS[r][s][i];
        mc(x, y);
        for (int i = 0; i < 4; i++) if (y[i] != ROWS_Y[r][s][i]) tab_ok = 0;
    }
    printf("Table 2 dy=MC(dx), dy'=MC(dx') reproduced in this field: %s\n", tab_ok ? "YES" : "NO");

    build_refs();
    printf("reference set (Table 2 rows x rot_k x 6 orderings, normalized, deduplicated): %d ordered pairs\n", nrefs);

    /* ---------------- normalized search ---------------- */
    int NT = MODE == 0 ? 3 : 4;
    int npat = NT * NT * NT * NT;
    long long ncand = 0;
    for (int p = 0; p < npat; p++) {
        int t[4], pp = p;
        for (int i = 0; i < 4; i++) { t[i] = pp % NT; pp /= NT; }
        /* exclusion: u=0 (all t in {0,3}), v=0 (all t in {1,3}), u=v (all t in {2,3}) */
        int allu0 = 1, allv0 = 1, alleq = 1, first = -1;
        for (int i = 0; i < 4; i++) {
            if (!(t[i] == 0 || t[i] == 3)) allu0 = 0;
            if (!(t[i] == 1 || t[i] == 3)) allv0 = 0;
            if (!(t[i] == 2 || t[i] == 3)) alleq = 0;
            if (t[i] != 3 && first < 0) first = i;
        }
        if (allu0 || allv0 || alleq) continue;
        /* active bytes other than 'first' get free eta */
        int act[4], na = 0;
        for (int i = 0; i < 4; i++) if (t[i] != 3 && i != first) act[na++] = i;
        long long tot = 1; for (int i = 0; i < na; i++) tot *= (Q - 1);
        int eta[4] = {0, 0, 0, 0}; eta[first] = 1;
        for (long long c = 0; c < tot; c++) {
            long long cc = c;
            for (int i = 0; i < na; i++) { eta[act[i]] = 1 + (int)(cc % (Q - 1)); cc /= (Q - 1); }
            unsigned char u[4], v[4];
            for (int i = 0; i < 4; i++) {
                switch (t[i]) {
                case 0: u[i] = 0; v[i] = (unsigned char)eta[i]; break;
                case 1: u[i] = (unsigned char)eta[i]; v[i] = 0; break;
                case 2: u[i] = v[i] = (unsigned char)eta[i]; break;
                default: u[i] = v[i] = 0;
                }
            }
            ncand++;
            if (out_ok(u, v)) add_sol(u, v);
        }
    }
    double tsearch = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("normalized candidates tested: %lld ; normalized solutions (ordered (u,v), eta_first=1): %d\n", ncand, nsols);
    printf("=> ordered (u,v) solutions in total = %d x (2^w-1) = %lld\n", nsols, (long long)nsols * (Q - 1));

    /* scalar orbit check: tau*(u,v) is a solution for every tau and normalizes back */
    int orbit_ok = 1;
    for (int s = 0; s < nsols; s++) for (int tau = 1; tau < Q; tau++) {
        unsigned char u[4], v[4], nu[4], nv[4];
        for (int i = 0; i < 4; i++) { u[i] = mul[tau][sols[s].u[i]]; v[i] = mul[tau][sols[s].v[i]]; }
        if (!out_ok(u, v)) orbit_ok = 0;
        normalize(u, v, nu, nv);
        if (!sol_eq(nu, nv, sols[s].u, sols[s].v)) orbit_ok = 0;
    }
    printf("scalar-orbit check (every tau*(u,v) is a solution; orbit size 2^w-1): %s\n", orbit_ok ? "OK" : "FAILED");

    /* (i) membership both directions */
    int notin = 0, refmiss = 0;
    for (int s = 0; s < nsols; s++) if (find_ref(sols[s].u, sols[s].v) < 0) {
        notin++;
        printf("  SOLUTION NOT IN TABLE-2 SPAN: u=[%X,%X,%X,%X] v=[%X,%X,%X,%X]\n",
               sols[s].u[0], sols[s].u[1], sols[s].u[2], sols[s].u[3], sols[s].v[0], sols[s].v[1], sols[s].v[2], sols[s].v[3]);
    }
    for (int t = 0; t < nrefs; t++) {
        int f = 0;
        for (int s = 0; s < nsols; s++) if (sol_eq(sols[s].u, sols[s].v, refs[t].s.u, refs[t].s.v)) { f = 1; break; }
        if (!f) refmiss++;
    }
    printf("(i) solutions outside {tau*rot_k(row)} : %d ; reference pairs not found by search: %d ; set equality: %s\n",
           notin, refmiss, (notin == 0 && refmiss == 0 && nsols == nrefs) ? "YES" : "NO");

    /* (b),(c) classes */
    int keys[4096], kcnt[4096], nkeys = 0;
    for (int s = 0; s < nsols; s++) {
        unsigned char sm[4];
        for (int i = 0; i < 4; i++) sm[i] = sols[s].u[i] ^ sols[s].v[i];
        int k = class_key(zmask(sols[s].u), zmask(sols[s].v), zmask(sm));
        int f = -1;
        for (int j = 0; j < nkeys; j++) if (keys[j] == k) f = j;
        if (f < 0) { keys[nkeys] = k; kcnt[nkeys] = 0; f = nkeys++; }
        kcnt[f]++;
    }
    printf("number of classes Z(D): %d\n", nkeys);
    for (int j = 0; j < nkeys; j++) {
        int k = keys[j], a = k >> 8 & 15, b = k >> 4 & 15, c = k & 15;
        printf("class Z = {"); print_mask(stdout, a); printf(", "); print_mask(stdout, b); printf(", "); print_mask(stdout, c); printf("}\n");
        printf("   normalized ordered pairs: %d ; ordered (u,v) pairs total: %lld ; /6 (unordered triples = Lemma 7 count): %lld = %g x (2^w-1)\n",
               kcnt[j], (long long)kcnt[j] * (Q - 1), (long long)kcnt[j] * (Q - 1) / 6, kcnt[j] / 6.0);
        /* distinct zero sets? */
        printf("   three zero sets distinct: %s\n", (a != b && b != c && a != c) ? "yes" : "NO");
        /* weights per zero set, consistency */
        int wout[16]; for (int m = 0; m < 16; m++) wout[m] = -1;
        int wconsistent = 1;
        int rotlist[64], nrot = 0;
        for (int s = 0; s < nsols; s++) {
            unsigned char sm[4], mu[4], mv[4], ms[4];
            for (int i = 0; i < 4; i++) sm[i] = sols[s].u[i] ^ sols[s].v[i];
            if (class_key(zmask(sols[s].u), zmask(sols[s].v), zmask(sm)) != k) continue;
            mc(sols[s].u, mu); mc(sols[s].v, mv); for (int i = 0; i < 4; i++) ms[i] = mu[i] ^ mv[i];
            int zm[3] = {zmask(sols[s].u), zmask(sols[s].v), zmask(sm)}, wo[3] = {wt(mu), wt(mv), wt(ms)};
            for (int q = 0; q < 3; q++) { if (wout[zm[q]] < 0) wout[zm[q]] = wo[q]; else if (wout[zm[q]] != wo[q]) wconsistent = 0; }
            int r = find_ref(sols[s].u, sols[s].v);
            if (r >= 0) {
                int lab = refs[r].row * 10 + refs[r].rot, dup = 0;
                for (int q = 0; q < nrot; q++) if (rotlist[q] == lab) dup = 1;
                if (!dup) rotlist[nrot++] = lab;
            }
        }
        int ms3[3] = {a, b, c};
        printf("   (input zero set: wt_in -> wt_out of MC):");
        for (int q = 0; q < 3; q++) { printf("  "); print_mask(stdout, ms3[q]); printf(": %d -> %d", 4 - __builtin_popcount(ms3[q]), wout[ms3[q]]); }
        printf("   [constant within class: %s]\n", wconsistent ? "yes" : "NO");
        printf("   table-2 rotations in class (row D_r, rot_k; a rotation may appear under several labels if it equals another rotation up to reordering):");
        for (int q = 0; q < nrot; q++) printf(" rot%d(D%d)", rotlist[q] % 10, rotlist[q] / 10);
        printf("\n");
    }

    /* (a) per ordered pattern (zu,zv): Lemma 7 quantity for a fixed assignment of zero sets */
    printf("(a) per ordered zero-set assignment (Z(u),Z(v)) [Z(u+v) determined]: count of ordered (u,v) pairs\n");
    int pc[16][16]; memset(pc, 0, sizeof pc);
    for (int s = 0; s < nsols; s++) pc[zmask(sols[s].u)][zmask(sols[s].v)]++;
    int npc = 0;
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) if (pc[a][b]) {
        npc++;
        printf("   Z(u)="); print_mask(stdout, a); printf(" Z(v)="); print_mask(stdout, b);
        printf(" : %lld = %d x (2^w-1)\n", (long long)pc[a][b] * (Q - 1), pc[a][b]);
    }
    printf("   number of ordered zero-set assignments with solutions: %d\n", npc);

    /* raw cross-check without normalization */
    if (RAW) {
        long long rawcnt = 0;
        for (int p = 0; p < npat; p++) {
            int t[4], pp = p;
            for (int i = 0; i < 4; i++) { t[i] = pp % NT; pp /= NT; }
            int allu0 = 1, allv0 = 1, alleq = 1;
            for (int i = 0; i < 4; i++) {
                if (!(t[i] == 0 || t[i] == 3)) allu0 = 0;
                if (!(t[i] == 1 || t[i] == 3)) allv0 = 0;
                if (!(t[i] == 2 || t[i] == 3)) alleq = 0;
            }
            if (allu0 || allv0 || alleq) continue;
            int act[4], na = 0;
            for (int i = 0; i < 4; i++) if (t[i] != 3) act[na++] = i;
            long long tot = 1; for (int i = 0; i < na; i++) tot *= (Q - 1);
            int eta[4] = {0};
            for (long long c = 0; c < tot; c++) {
                long long cc = c;
                for (int i = 0; i < na; i++) { eta[act[i]] = 1 + (int)(cc % (Q - 1)); cc /= (Q - 1); }
                unsigned char u[4], v[4];
                for (int i = 0; i < 4; i++) {
                    switch (t[i]) {
                    case 0: u[i] = 0; v[i] = (unsigned char)eta[i]; break;
                    case 1: u[i] = (unsigned char)eta[i]; v[i] = 0; break;
                    case 2: u[i] = v[i] = (unsigned char)eta[i]; break;
                    default: u[i] = v[i] = 0;
                    }
                }
                if (out_ok(u, v)) rawcnt++;
            }
        }
        printf("RAW (no normalization) ordered (u,v) solutions: %lld ; matches normalized x (2^w-1): %s\n",
               rawcnt, rawcnt == (long long)nsols * (Q - 1) ? "YES" : "NO");
    }
    printf("search time: %.2f s ; total time: %.2f s\n\n", tsearch, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
