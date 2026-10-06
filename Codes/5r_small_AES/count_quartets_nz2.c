/*
 * count_quartets_nz2.c
 * Small-AES valid-quartet count of the pairing with base pair (P^0, P^1) in a
 * Family-A quartet of Table 4 (n_z = 2; base pair 4-active on all nibbles
 * {0,5,10,15}).  Counts valid quartets per trial
 * over a structure of 2^16 plaintexts and verifies the n_z = 2 instance
 * of the multiple-of property (Theorem 10 of the paper),
 *      N = 2^{2-n_z}(2^w)^{n_z} A + 8B = 256 A + 8B   (w = 4, n_z = 2),
 * i.e. the bundle size is 2^{2-2}(2^4)^2 = 256.
 * Paper: Table 6 and Table 7 (Section 3.4.1) and the small-AES column of
 * Table 11 (Section 3.4.2), from Results/5r_small_AES/{5r,10r}_nz2.txt.
 *
 * Implementation notes:
 *   - Pair search uses direct-addressed linked-list bucketing on the
 *     16-bit inverse-diagonal value (O(N) instead of qsort's O(N log N)).
 *   - struct_arr.pt is dropped: the active nibbles are encoded directly
 *     in the array index (i = (n0<<12)|(n5<<8)|(n10<<4)|n15), and the
 *     plaintext is reconstructed on-the-fly only during encryption.
 *   - The swap pair (P^2, P^3) is not re-encrypted: their indices are
 *     computed by bit-masking from (i_0, i_1) and their ciphertexts are
 *     looked up in struct_arr.ct directly.
 *
 * Per trial:
 *   1. Pick a random base plaintext and master key.
 *   2. Encrypt all 2^16 plaintexts; only ct is stored (8 bytes per slot).
 *   3. For each inverse diagonal d, build head[1<<16] / next[1<<16]
 *      linked lists keyed by the 16-bit ct[d] value.
 *   4. For each colliding pair (i_0, i_1) in a bucket: skip pairs with
 *      any zero active-nibble difference (so the base pair is 4-active),
 *      compute the swap pair indices via bit-mask, and check whether
 *      their ciphertexts collide on the same inverse diagonal.
 *
 * Usage:
 *   ./count_quartets_nz2 <rounds> [test_count] [seed]
 *
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

typedef struct {
    uint8_t ct[8];
} Element;

static inline uint8_t get_nibble(const uint8_t *data, int idx) {
    if (idx % 2 == 0) return (data[idx / 2] >> 4) & 0xF;
    return data[idx / 2] & 0xF;
}

static inline void set_nibble(uint8_t *data, int idx, uint8_t val) {
    int byte_idx = idx / 2;
    if (idx % 2 == 0)
        data[byte_idx] = (data[byte_idx] & 0x0F) | ((val & 0xF) << 4);
    else
        data[byte_idx] = (data[byte_idx] & 0xF0) |  (val & 0xF);
}

/* Index encoding within a structure: i in [0, 2^16) maps to nibbles
 *   nibble 0  = (i >> 12) & 0xF
 *   nibble 5  = (i >>  8) & 0xF
 *   nibble 10 = (i >>  4) & 0xF
 *   nibble 15 =  i        & 0xF
 * Bit mask 0xF0F0 selects (nibble 0, nibble 10), 0x0F0F selects
 * (nibble 5, nibble 15).  The byte-swap of nibbles 0 and 10 between two
 * indices i0, i1 maps to:
 *   p2_idx = (i1 & 0xF0F0) | (i0 & 0x0F0F)
 *   p3_idx = (i0 & 0xF0F0) | (i1 & 0x0F0F)
 */
static inline void build_pt_from_index(uint8_t *pt, const uint8_t *base_pt, uint32_t i) {
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
    uint16_t v = 0;
    v |= (uint16_t)get_nibble(ct, INV_DIAG_INDICES[d][0]) << 12;
    v |= (uint16_t)get_nibble(ct, INV_DIAG_INDICES[d][1]) <<  8;
    v |= (uint16_t)get_nibble(ct, INV_DIAG_INDICES[d][2]) <<  4;
    v |= (uint16_t)get_nibble(ct, INV_DIAG_INDICES[d][3]);
    return v;
}

static int compare_uint64(const void *a, const void *b) {
    uint64_t ua = *(const uint64_t *)a;
    uint64_t ub = *(const uint64_t *)b;
    if (ua < ub) return -1;
    if (ua > ub) return  1;
    return 0;
}

static uint64_t run_single_trial(Element *struct_arr,
                                 uint32_t *head,
                                 uint32_t *next,
                                 const uint8_t *masterKey,
                                 const uint8_t *base_pt,
                                 int num_rounds) {
    uint64_t mixture_count = 0;

    /* Phase 1: encrypt all 2^16 plaintexts; store only ct. */
    for (uint32_t i = 0; i < STRUCT_SIZE; i++) {
        uint8_t pt[8];
        build_pt_from_index(pt, base_pt, i);
        small_aes_enc_ttable(struct_arr[i].ct, pt, masterKey, num_rounds);
    }

    /* Phase 2: per inverse diagonal, bucket by ct[d] and enumerate pairs. */
    for (int d = 0; d < 4; d++) {
        memset(head, 0xFF, sizeof(uint32_t) * BUCKET_COUNT);
        for (uint32_t i = 0; i < STRUCT_SIZE; i++) {
            uint16_t v = get_inv_diag_val(struct_arr[i].ct, d);
            next[i] = head[v];
            head[v] = i;
        }

        for (uint32_t v = 0; v < BUCKET_COUNT; v++) {
            for (uint32_t i0 = head[v]; i0 != SENTINEL; i0 = next[i0]) {
                for (uint32_t i1 = next[i0]; i1 != SENTINEL; i1 = next[i1]) {
                    /* Skip if any active-nibble difference is zero. */
                    uint16_t diff = (uint16_t)(i0 ^ i1);
                    if ((diff & 0xF000) == 0 ||
                        (diff & 0x0F00) == 0 ||
                        (diff & 0x00F0) == 0 ||
                        (diff & 0x000F) == 0) continue;

                    /* Swap nibbles 0 and 10 between (P^0, P^1) -> (P^2, P^3). */
                    uint32_t p2 = (i1 & 0xF0F0u) | (i0 & 0x0F0Fu);
                    uint32_t p3 = (i0 & 0xF0F0u) | (i1 & 0x0F0Fu);

                    if (get_inv_diag_val(struct_arr[p2].ct, d) ==
                        get_inv_diag_val(struct_arr[p3].ct, d))
                        mixture_count++;
                }
            }
        }
    }
    /* Each unique quartet is enumerated twice (once with (P^0,P^1) as base
     * pair, once with the swap pair (P^2,P^3) as base); halve to obtain
     * the unique-quartet count, matching Theorem 10 of the paper. */
    return mixture_count / 2;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "Usage: %s <rounds> [test_count] [seed]\n", argv[0]);
        fprintf(stderr, "  Default: test_count=1000, seed=42 (env SEED also honoured).\n");
        return 1;
    }
    int num_rounds = atoi(argv[1]);
    int test_count = (argc >= 3) ? atoi(argv[2]) : 1000;
    unsigned int base_seed = 42u;
    const char *env_seed = getenv("SEED");
    if (env_seed) base_seed = (unsigned int)strtoul(env_seed, NULL, 10);
    if (argc >= 4)  base_seed = (unsigned int)strtoul(argv[3], NULL, 10);
    if (num_rounds < 1 || num_rounds > 10) {
        fprintf(stderr, "Invalid rounds: %s\n", argv[1]);
        return 1;
    }
    if (test_count <= 0) {
        fprintf(stderr, "Invalid test_count: %s\n", argv[2]);
        return 1;
    }

    printf("Config: Tests=%d, StructureSize=%d, Rounds=%d, Seed=%u\n\n",
           test_count, STRUCT_SIZE, num_rounds, base_seed);
    fflush(stdout);

    Element  *struct_arr = (Element  *)malloc(sizeof(Element)  * STRUCT_SIZE);
    uint32_t *head       = (uint32_t *)malloc(sizeof(uint32_t) * BUCKET_COUNT);
    uint32_t *next       = (uint32_t *)malloc(sizeof(uint32_t) * STRUCT_SIZE);
    uint64_t *results    = (uint64_t *)malloc(sizeof(uint64_t) * test_count);
    if (!struct_arr || !head || !next || !results) { perror("malloc"); return 1; }

    for (int t = 0; t < test_count; t++) {
        unsigned int trial_seed = rmd_derive_small_trial_seed(base_seed, t);
        uint8_t mk[8], base_pt[8];
        for (int k = 0; k < 8; k++) mk[k]      = (uint8_t)rmd_rand_byte(&trial_seed);
        for (int k = 0; k < 8; k++) base_pt[k] = (uint8_t)rmd_rand_byte(&trial_seed);
        results[t] = run_single_trial(struct_arr, head, next, mk, base_pt, num_rounds);
        if ((t + 1) % 10 == 0) {
            fprintf(stderr, "\r[Progress] %d/%d done (last count: %llu)   ",
                    t + 1, test_count, (unsigned long long)results[t]);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\n");

    qsort(results, test_count, sizeof(uint64_t), compare_uint64);

    printf("--- %d-Round Frequency Results ---\n", num_rounds);
    uint64_t cur_val = results[0];
    int cur_freq = 0;
    for (int i = 0; i < test_count; i++) {
        if (results[i] == cur_val) {
            cur_freq++;
        } else {
            printf("Count %-8llu : %4d times\n",
                   (unsigned long long)cur_val, cur_freq);
            cur_val = results[i];
            cur_freq = 1;
        }
    }
    printf("Count %-8llu : %4d times\n",
           (unsigned long long)cur_val, cur_freq);
    printf("\n(Total: %d trials)\n", test_count);

    free(struct_arr);
    free(head);
    free(next);
    free(results);
    return 0;
}
