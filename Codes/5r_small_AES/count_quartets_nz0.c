/*
 * count_quartets_nz0.c
 * Small-AES valid-quartet count of the pairing with base pair (P^0, P^3) in a
 * Family-B quartet of Table 4 (n_z = 0), base pair active on nibbles {5,15}.
 * Counts valid quartets per trial over a structure of 2^16 plaintexts and
 * thereby verifies the n_z = 0 instance of the multiple-of property (Theorem 10 of the paper),
 *      N = 2^{2-n_z}(2^w)^{n_z} A + 8B = 4A + 8B   (w = 4, n_z = 0),
 * i.e. the bundle size is 2^{2-0}(2^4)^0 = 4.
 * Paper: Table 6 and Table 9 (Section 3.4.1) and the small-AES column of
 * Table 11 (Section 3.4.2), from Results/5r_small_AES/{5r,10r}_nz0.txt.
 *
 * Structure (3D, full diagonal of 2^16 plaintexts):
 *   i in [0,256) : value of (nibble 0, nibble 10)   = n0=i>>4, n10=i&0xF
 *   j in [0,16)  : value of nibble 5
 *   k in [0,16)  : value of nibble 15
 *
 * n_z = 0 quartet (base pair 2-active on {5,15}):
 *   base    : (i, jP, kP), (i, jQ, kQ)     same i, differ in both j and k
 *   partner : (i', jP, kQ), (i', jQ, kP)   other i', byte-15 cross-swapped
 * Valid when both pairs collide (zero difference) on the SAME inverse
 * diagonal.  Each quartet is enumerated twice (base / partner i-slice roles),
 * so the raw pair count is halved, matching the convention of count_quartets_nz2.c.
 *
 * Usage:
 *   ./count_quartets_nz0 <rounds> [test_count] [seed]
 * Defaults: test_count = 1000, seed = 42 (env SEED also honoured).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "small_aes.h"
#include "../common/rmd_random.h"

#define S0 256
#define S1 16
#define S2 16
#define SLICE       (S1 * S2)        /* 256 plaintexts per i-slice */
#define STRUCT_SIZE (S0 * SLICE)     /* 2^16 */

typedef struct { uint8_t ct[8]; } Element;

static const int INV_DIAG_INDICES[4][4] = {
    { 0,  7, 10, 13},
    { 4,  1, 14, 11},
    { 8,  5,  2, 15},
    {12,  9,  6,  3}
};

static inline uint8_t get_nibble(const uint8_t *d, int idx) {
    return (idx & 1) ? (d[idx / 2] & 0xF) : ((d[idx / 2] >> 4) & 0xF);
}
static inline void set_nibble(uint8_t *d, int idx, uint8_t v) {
    int b = idx / 2;
    if (idx & 1) d[b] = (d[b] & 0xF0) | (v & 0xF);
    else         d[b] = (d[b] & 0x0F) | ((v & 0xF) << 4);
}
static inline uint16_t inv_diag(const uint8_t *ct, int d) {
    return (uint16_t)((get_nibble(ct, INV_DIAG_INDICES[d][0]) << 12) |
                      (get_nibble(ct, INV_DIAG_INDICES[d][1]) <<  8) |
                      (get_nibble(ct, INV_DIAG_INDICES[d][2]) <<  4) |
                       get_nibble(ct, INV_DIAG_INDICES[d][3]));
}

typedef struct { uint16_t v; uint8_t j, k; } SV;
static int cmp_sv(const void *a, const void *b) {
    int va = ((const SV *)a)->v, vb = ((const SV *)b)->v;
    return (va > vb) - (va < vb);
}
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static uint64_t run_single_trial(Element *arr, SV *sv,
                                 const uint8_t *mk, const uint8_t *base_pt,
                                 int num_rounds) {
    /* Phase 1: encrypt the full 2^16 structure. */
    for (int i = 0; i < S0; i++)
        for (int j = 0; j < S1; j++)
            for (int k = 0; k < S2; k++) {
                uint8_t pt[8];
                memcpy(pt, base_pt, 8);
                set_nibble(pt, 0,  (i >> 4) & 0xF);
                set_nibble(pt, 10,  i       & 0xF);
                set_nibble(pt, 5,   j       & 0xF);
                set_nibble(pt, 15,  k       & 0xF);
                small_aes_enc_ttable(arr[(long)i * SLICE + j * S2 + k].ct,
                                     pt, mk, num_rounds);
            }

    /* Phase 2: per inverse diagonal, find {5,15}-active colliding base pairs
     * inside each i-slice, then test the byte-15-swapped partner in i'. */
    uint64_t count = 0;
    for (int d = 0; d < 4; d++) {
        for (int iP = 0; iP < S0; iP++) {
            Element *sl = &arr[(long)iP * SLICE];
            for (int j = 0; j < S1; j++)
                for (int k = 0; k < S2; k++) {
                    int idx = j * S2 + k;
                    sv[idx].v = inv_diag(sl[idx].ct, d);
                    sv[idx].j = (uint8_t)j;
                    sv[idx].k = (uint8_t)k;
                }
            qsort(sv, SLICE, sizeof(SV), cmp_sv);
            for (int a = 0; a < SLICE; ) {
                int b = a + 1;
                while (b < SLICE && sv[b].v == sv[a].v) b++;
                for (int x = a; x < b; x++)
                    for (int y = x + 1; y < b; y++) {
                        int jP = sv[x].j, kP = sv[x].k, jQ = sv[y].j, kQ = sv[y].k;
                        if (jP == jQ || kP == kQ) continue;   /* need {5,15} 2-active */
                        for (int ip = 0; ip < S0; ip++) {
                            if (ip == iP) continue;
                            Element *sl2 = &arr[(long)ip * SLICE];
                            if (inv_diag(sl2[jP * S2 + kQ].ct, d) ==
                                inv_diag(sl2[jQ * S2 + kP].ct, d))
                                count++;
                        }
                    }
                a = b;
            }
        }
    }
    return count / 2;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "Usage: %s <rounds> [test_count] [seed]\n", argv[0]);
        return 1;
    }
    int num_rounds = atoi(argv[1]);
    int test_count = (argc >= 3) ? atoi(argv[2]) : 1000;
    unsigned int base_seed = 42u;
    const char *env_seed = getenv("SEED");
    if (env_seed) base_seed = (unsigned int)strtoul(env_seed, NULL, 10);
    if (argc >= 4) base_seed = (unsigned int)strtoul(argv[3], NULL, 10);
    if (num_rounds < 1 || num_rounds > 10 || test_count <= 0) { fprintf(stderr, "bad args\n"); return 1; }

    printf("Config: n_z=0, Tests=%d, StructureSize=%d, Rounds=%d, Seed=%u\n\n",
           test_count, STRUCT_SIZE, num_rounds, base_seed);
    fflush(stdout);

    Element  *arr     = (Element  *)malloc(sizeof(Element)  * STRUCT_SIZE);
    SV       *sv      = (SV       *)malloc(sizeof(SV)       * SLICE);
    uint64_t *results = (uint64_t *)malloc(sizeof(uint64_t) * test_count);
    if (!arr || !sv || !results) { perror("malloc"); return 1; }

    for (int t = 0; t < test_count; t++) {
        unsigned int trial_seed = rmd_derive_small_trial_seed(base_seed, t);
        uint8_t mk[8], base_pt[8];
        for (int x = 0; x < 8; x++) mk[x]      = (uint8_t)rmd_rand_byte(&trial_seed);
        for (int x = 0; x < 8; x++) base_pt[x] = (uint8_t)rmd_rand_byte(&trial_seed);
        results[t] = run_single_trial(arr, sv, mk, base_pt, num_rounds);
        if ((t + 1) % 10 == 0) {
            fprintf(stderr, "\r[Progress] %d/%d (last: %llu)   ",
                    t + 1, test_count, (unsigned long long)results[t]);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\n");

    qsort(results, test_count, sizeof(uint64_t), cmp_u64);
    printf("--- %d-Round Frequency Results (n_z=0) ---\n", num_rounds);
    uint64_t cur = results[0]; int freq = 0;
    for (int i = 0; i < test_count; i++) {
        if (results[i] == cur) freq++;
        else { printf("Count %-8llu : %4d times\n", (unsigned long long)cur, freq);
               cur = results[i]; freq = 1; }
    }
    printf("Count %-8llu : %4d times\n", (unsigned long long)cur, freq);
    printf("\n(Total: %d trials)\n", test_count);

    free(arr); free(sv); free(results);
    return 0;
}
