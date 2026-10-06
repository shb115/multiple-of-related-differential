/*
 * count_quartets_nz1.c
 * Small-AES valid-quartet count of the pairing n_z = 1 on {0,10}: base pair
 * (P^0, P^2) in a Family-A quartet of Table 4, 2-active on nibbles {0,10}
 * (the pairing of Yan et al.).  Counts valid quartets per
 * trial over a structure of 2^16 plaintexts and verifies the n_z = 1 instance
 * of the multiple-of property (Theorem 10 of the paper),
 *      N = 2^{2-n_z}(2^w)^{n_z} A + 8B = 2 * 2^4 A + 8B = 32A + 8B   (w=4, n_z=1),
 * i.e. the bundle size is 2^{2-1}(2^4)^1 = 32.
 * Paper: Table 6 and Table 8 (Section 3.4.1) and the small-AES column of
 * Table 11 (Section 3.4.2), from Results/5r_small_AES/{5r,10r}_nz1.txt.
 *
 * Index encoding within the 2^16 structure (same as count_quartets_nz2.c):
 *   i  ->  nibble0=(i>>12), nibble5=(i>>8), nibble10=(i>>4), nibble15=i, all &0xF.
 *   mask 0xF0F0 selects (nibble0, nibble10); 0x0F0F selects (nibble5, nibble15).
 *
 * n_z=1 quartet: base pair (i0,i1) 2-active on {0,10} (equal on {5,15});
 * partner pair (i0^g, i1^g) for a 2-active offset g on {5,15}.  Valid when both
 * pairs collide on the SAME inverse diagonal.  Each quartet is enumerated twice
 * (base / partner roles), so the raw count is halved.
 *
 * Usage:
 *   ./count_quartets_nz1 <rounds> [test_count] [seed]
 * Defaults: test_count = 1000, seed = 42 (env SEED also honoured).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "small_aes.h"
#include "../common/rmd_random.h"

#define STRUCT_SIZE 65536            /* 2^16 */
#define BUCKET_COUNT 65536           /* 2^16: ct[d] is 16 bits */
#define SENTINEL 0xFFFFFFFFu

typedef struct { uint8_t ct[8]; } Element;

static inline uint8_t get_nibble(const uint8_t *d, int idx) {
    return (idx & 1) ? (d[idx / 2] & 0xF) : ((d[idx / 2] >> 4) & 0xF);
}
static inline void set_nibble(uint8_t *d, int idx, uint8_t v) {
    int b = idx / 2;
    if (idx & 1) d[b] = (d[b] & 0xF0) | (v & 0xF);
    else         d[b] = (d[b] & 0x0F) | ((v & 0xF) << 4);
}
static void build_pt_from_index(uint8_t *pt, const uint8_t *base_pt, uint32_t i) {
    memcpy(pt, base_pt, 8);
    set_nibble(pt, 0,  (i >> 12) & 0xF);
    set_nibble(pt, 5,  (i >>  8) & 0xF);
    set_nibble(pt, 10, (i >>  4) & 0xF);
    set_nibble(pt, 15,  i        & 0xF);
}

static const int INV_DIAG_INDICES[4][4] = {
    { 0,  7, 10, 13},
    { 4,  1, 14, 11},
    { 8,  5,  2, 15},
    {12,  9,  6,  3}
};
static inline uint16_t get_inv_diag_val(const uint8_t *ct, int d) {
    return (uint16_t)((get_nibble(ct, INV_DIAG_INDICES[d][0]) << 12) |
                      (get_nibble(ct, INV_DIAG_INDICES[d][1]) <<  8) |
                      (get_nibble(ct, INV_DIAG_INDICES[d][2]) <<  4) |
                       get_nibble(ct, INV_DIAG_INDICES[d][3]));
}
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/* Offsets that are 2-active on {5,15}: g = (g5<<8)|g15, g5,g15 in 1..15. */
static uint16_t g_offsets[225];
static int      g_count;
static void build_offsets(void) {
    g_count = 0;
    for (int g5 = 1; g5 < 16; g5++)
        for (int g15 = 1; g15 < 16; g15++)
            g_offsets[g_count++] = (uint16_t)((g5 << 8) | g15);
}

static uint64_t run_single_trial(Element *arr, uint32_t *head, uint32_t *next,
                                 const uint8_t *mk, const uint8_t *base_pt,
                                 int num_rounds) {
    for (uint32_t i = 0; i < STRUCT_SIZE; i++) {
        uint8_t pt[8];
        build_pt_from_index(pt, base_pt, i);
        small_aes_enc_ttable(arr[i].ct, pt, mk, num_rounds);
    }

    uint64_t count = 0;
    for (int d = 0; d < 4; d++) {
        memset(head, 0xFF, sizeof(uint32_t) * BUCKET_COUNT);
        for (uint32_t i = 0; i < STRUCT_SIZE; i++) {
            uint16_t v = get_inv_diag_val(arr[i].ct, d);
            next[i] = head[v];
            head[v] = i;
        }
        for (uint32_t v = 0; v < BUCKET_COUNT; v++) {
            for (uint32_t i0 = head[v]; i0 != SENTINEL; i0 = next[i0]) {
                for (uint32_t i1 = next[i0]; i1 != SENTINEL; i1 = next[i1]) {
                    uint16_t diff = (uint16_t)(i0 ^ i1);
                    /* base pair must be 2-active on {0,10} and 0 on {5,15}. */
                    if ((diff & 0x0F0Fu) != 0) continue;          /* {5,15} must match */
                    if ((diff & 0xF000u) == 0 || (diff & 0x00F0u) == 0) continue; /* {0,10} both active */
                    /* partner = base translated by a 2-active offset g on {5,15}. */
                    for (int gi = 0; gi < g_count; gi++) {
                        uint32_t g = g_offsets[gi];
                        uint32_t p2 = i0 ^ g;
                        uint32_t p3 = i1 ^ g;
                        if (get_inv_diag_val(arr[p2].ct, d) ==
                            get_inv_diag_val(arr[p3].ct, d))
                            count++;
                    }
                }
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

    build_offsets();
    printf("Config: n_z=1, Tests=%d, StructureSize=%d, Rounds=%d, Seed=%u\n\n",
           test_count, STRUCT_SIZE, num_rounds, base_seed);
    fflush(stdout);

    Element  *arr     = (Element  *)malloc(sizeof(Element)  * STRUCT_SIZE);
    uint32_t *head    = (uint32_t *)malloc(sizeof(uint32_t) * BUCKET_COUNT);
    uint32_t *next    = (uint32_t *)malloc(sizeof(uint32_t) * STRUCT_SIZE);
    uint64_t *results = (uint64_t *)malloc(sizeof(uint64_t) * test_count);
    if (!arr || !head || !next || !results) { perror("malloc"); return 1; }

    for (int t = 0; t < test_count; t++) {
        unsigned int trial_seed = rmd_derive_small_trial_seed(base_seed, t);
        uint8_t mk[8], base_pt[8];
        for (int x = 0; x < 8; x++) mk[x]      = (uint8_t)rmd_rand_byte(&trial_seed);
        for (int x = 0; x < 8; x++) base_pt[x] = (uint8_t)rmd_rand_byte(&trial_seed);
        results[t] = run_single_trial(arr, head, next, mk, base_pt, num_rounds);
        if ((t + 1) % 10 == 0) {
            fprintf(stderr, "\r[Progress] %d/%d (last: %llu)   ",
                    t + 1, test_count, (unsigned long long)results[t]);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\n");

    qsort(results, test_count, sizeof(uint64_t), cmp_u64);
    printf("--- %d-Round Frequency Results (n_z=1) ---\n", num_rounds);
    uint64_t cur = results[0]; int freq = 0;
    for (int i = 0; i < test_count; i++) {
        if (results[i] == cur) freq++;
        else { printf("Count %-8llu : %4d times\n", (unsigned long long)cur, freq);
               cur = results[i]; freq = 1; }
    }
    printf("Count %-8llu : %4d times\n", (unsigned long long)cur, freq);
    printf("\n(Total: %d trials)\n", test_count);

    free(arr); free(head); free(next); free(results);
    return 0;
}
