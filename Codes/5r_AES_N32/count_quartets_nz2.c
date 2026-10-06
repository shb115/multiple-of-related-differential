/*
 * count_quartets_nz2.c
 * AES-128 valid-quartet count of the pairing with base pair (P^0, P^1) in a
 * Family-A quartet of Table 4 (n_z = 2; base pair 4-active on all bytes
 * {0,5,10,15}).  Counts valid quartets per trial
 * over a full structure of 2^STRUCT_BITS plaintexts (default 2^32;
 * AES-NI + OpenMP).  Used for the n_z = 2 instance of Theorem 10 (multiple-of property)
 * verification on the full AES (Section 3.4.1 of the paper), with
 *      N = 2^{2-n_z}(2^w)^{n_z} A + 8B = 65536 A + 8B   (w = 8, n_z = 2),
 * i.e. the bundle size is 2^{2-2}(2^8)^2 = 65536.
 * Paper: Table 10 (Section 3.4.1) and the AES column of Table 11 (Section 3.4.2),
 * from Results/5r_AES_N32/raw_nz2_{5,10}r.txt; Codes/tools/summarize_counts.py prints
 * the mean, the detection rate and the sample standard deviation quoted in Section 3.4.1.
 *
 * Implementation notes:
 *   - Pair search uses direct-addressed linked-list bucketing on the
 *     32-bit inverse-diagonal value (O(N) instead of qsort's O(N log N)).
 *   - struct_arr.pt is dropped: the 4 active bytes are encoded directly
 *     in the array index (i = byte0 | byte10<<8 | byte5<<16 |
 *     byte15<<24), and the plaintext is built on-the-fly only during
 *     encryption.
 *   - The swap pair (P^2, P^3) is not re-encrypted: their indices are
 *     computed by bit-mask from (i_0, i_1) and their ciphertexts are
 *     looked up in struct_arr.ct directly.
 *
 * For a quick pipeline smoke-test on a small machine, build with e.g.
 * -DSTRUCT_BITS=16 -DBUCKET_BITS=16 (this does NOT reproduce the exact
 * multiple-of bundle, which requires the full 2^32 diagonal space).
 *
 * Usage:
 *   ./count_quartets_nz2 <rounds> [test_count] [raw_file] [seed]
 *
 * <rounds> is 5 (5-round AES) or 10 (random-permutation baseline).
 * Defaults: test_count = 100, raw_file = raw_counts_<rounds>r.txt,
 * seed = 42 (env SEED also honoured).
 * The committed results (Results/5r_AES_N32/raw_nz2_5r.txt, raw_nz2_10r.txt) use
 * 700 trials at 5 rounds and 200 trials at 10 rounds, seed 42.
 *
 * Memory (per OpenMP thread, with default STRUCT_BITS=32):
 *   struct_arr.ct : 2^32 * 16 B = 64 GB
 *   head          : 2^32 *  8 B = 32 GB
 *   next          : 2^32 *  8 B = 32 GB
 *   total         : ~128 GB.  Use OMP_NUM_THREADS=1 unless the host has
 *   well over 128 GB RAM.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <wmmintrin.h>
#include <omp.h>
#include "../common/rmd_aesni.h"
#include "../common/rmd_random.h"

#ifndef STRUCT_BITS
#define STRUCT_BITS 32
#endif
#ifndef BUCKET_BITS
#define BUCKET_BITS 32
#endif
#define STRUCT_SIZE  (1ULL << STRUCT_BITS)
#define BUCKET_COUNT (1ULL << BUCKET_BITS)
#define BUCKET_MASK  (BUCKET_COUNT - 1ULL)
#define SENTINEL     (~(uint64_t)0)
#define BYTE_BITS    (STRUCT_BITS / 4)
#define BYTE_MASK    ((1u << BYTE_BITS) - 1u)

/* Index encoding:
 *   bits [0,         BYTE_BITS)   = byte0
 *   bits [BYTE_BITS, 2*BYTE_BITS) = byte10
 *   bits [2*BB,      3*BB)        = byte5
 *   bits [3*BB,      4*BB)        = byte15
 * The byte-swap of (byte0, byte10) between two indices i0, i1 yields:
 *   p2_idx = (i1 & MASK_LOW)  | (i0 & MASK_HIGH)
 *   p3_idx = (i0 & MASK_LOW)  | (i1 & MASK_HIGH)
 * MASK_LOW selects byte0 and byte10 bits; MASK_HIGH selects byte5 and byte15.
 */
#define MASK_LOW  ((1ULL << (2*BYTE_BITS)) - 1ULL)
#define MASK_HIGH (MASK_LOW << (2*BYTE_BITS))
#define BYTE0_MASK   ((uint64_t)BYTE_MASK)
#define BYTE10_MASK  ((uint64_t)BYTE_MASK <<   BYTE_BITS)
#define BYTE5_MASK   ((uint64_t)BYTE_MASK << (2*BYTE_BITS))
#define BYTE15_MASK  ((uint64_t)BYTE_MASK << (3*BYTE_BITS))

typedef struct {
    uint8_t ct[16];
} Element;

static inline void build_pt_from_index(uint8_t *pt, const uint8_t *base_pt, uint64_t idx) {
    memcpy(pt, base_pt, 16);
    pt[0]  = (uint8_t)( idx                  & BYTE_MASK);
    pt[10] = (uint8_t)((idx >>  BYTE_BITS   ) & BYTE_MASK);
    pt[5]  = (uint8_t)((idx >> (BYTE_BITS*2)) & BYTE_MASK);
    pt[15] = (uint8_t)((idx >> (BYTE_BITS*3)) & BYTE_MASK);
}

static uint64_t run_single_trial(Element *struct_arr,
                                 uint64_t *head,
                                 uint64_t *next,
                                 const __m128i *round_keys, int rounds,
                                 const uint8_t *base_pt) {
    /* Phase 1: encrypt all 2^STRUCT_BITS plaintexts; store only ct. */
    #pragma omp parallel for schedule(static)
    for (uint64_t idx = 0; idx < STRUCT_SIZE; idx++) {
        uint8_t pt[16];
        build_pt_from_index(pt, base_pt, idx);
        rmd_aes128_encrypt(pt, struct_arr[idx].ct, round_keys, rounds);
    }

    /* Phase 2: per inverse diagonal, bucket by ct[d] and enumerate pairs. */
    uint64_t count = 0;
    for (int d = 0; d < 4; d++) {
        /* Reset head buckets to SENTINEL (-1). */
        #pragma omp parallel for schedule(static)
        for (uint64_t v = 0; v < BUCKET_COUNT; v++) head[v] = SENTINEL;

        /* Sequential bucket build: bucket index = ct[d] & BUCKET_MASK so
         * that the same code works at full BUCKET_BITS=32 (no hash
         * collisions, the equality check below is always true) or with a
         * smaller BUCKET_BITS (e.g. for testing on small machines). */
        for (uint64_t i = 0; i < STRUCT_SIZE; i++) {
            uint32_t v = rmd_aes_inv_diag_val(struct_arr[i].ct, d);
            uint64_t b = (uint64_t)v & BUCKET_MASK;
            next[i] = head[b];
            head[b] = i;
        }

        for (uint64_t b = 0; b < BUCKET_COUNT; b++) {
            for (uint64_t i0 = head[b]; i0 != SENTINEL; i0 = next[i0]) {
                uint32_t v0 = rmd_aes_inv_diag_val(struct_arr[i0].ct, d);
                for (uint64_t i1 = next[i0]; i1 != SENTINEL; i1 = next[i1]) {
                    /* Real ct[d] equality (always true at BUCKET_BITS=32). */
                    if (rmd_aes_inv_diag_val(struct_arr[i1].ct, d) != v0) continue;

                    /* Skip if any active-byte difference is zero. */
                    uint64_t diff = i0 ^ i1;
                    if ((diff & BYTE0_MASK)  == 0 ||
                        (diff & BYTE10_MASK) == 0 ||
                        (diff & BYTE5_MASK)  == 0 ||
                        (diff & BYTE15_MASK) == 0) continue;

                    uint64_t p2 = (i1 & MASK_LOW) | (i0 & MASK_HIGH);
                    uint64_t p3 = (i0 & MASK_LOW) | (i1 & MASK_HIGH);

                    if (rmd_aes_inv_diag_val(struct_arr[p2].ct, d) ==
                        rmd_aes_inv_diag_val(struct_arr[p3].ct, d))
                        count++;
                }
            }
        }
    }
    /* Each unique quartet is enumerated twice (once with (P^0,P^1) as base
     * pair, once with the swap pair (P^2,P^3) as base); halve to obtain
     * the unique-quartet count, matching Theorem 10 of the paper. */
    return count / 2;
}

int main(int argc, char **argv) {
    setbuf(stdout, NULL);

    if (argc < 2 || argc > 5) {
        fprintf(stderr, "Usage: %s <rounds> [test_count] [raw_file] [seed]\n", argv[0]);
        fprintf(stderr, "  Default: test_count=100, raw_file=raw_counts_<rounds>r.txt, seed=42 (env SEED also honoured).\n");
        fprintf(stderr, "  Committed results: 700 trials (5 rounds), 200 trials (10 rounds).\n");
        return 1;
    }
    int rounds = atoi(argv[1]);
    if (rounds != 5 && rounds != 10) {
        fprintf(stderr, "Rounds must be 5 or 10; got %s\n", argv[1]);
        return 1;
    }
    int test_count = (argc >= 3) ? atoi(argv[2]) : 100;
    if (test_count <= 0) {
        fprintf(stderr, "Invalid test_count: %s\n", argv[2]);
        return 1;
    }

    char default_raw[64];
    snprintf(default_raw, sizeof(default_raw), "raw_counts_%dr.txt", rounds);
    const char *raw_file = (argc >= 4) ? argv[3] : default_raw;

    unsigned int base_seed = 42u;
    const char *env_seed = getenv("SEED");
    if (env_seed)  base_seed = (unsigned int)strtoul(env_seed, NULL, 10);
    if (argc >= 5) base_seed = (unsigned int)strtoul(argv[4], NULL, 10);

    int prior_trials = 0;
    {
        FILE *fr = fopen(raw_file, "r");
        if (fr) {
            int v;
            while (fscanf(fr, "%d", &v) == 1) prior_trials++;
            fclose(fr);
        }
    }

    fprintf(stderr, "[*] %d-Round AES Experiment, n_z=2 (2^%d structure, bundle 65536)\n",
            rounds, STRUCT_BITS);
    fprintf(stderr, "[*] Seed=%u, Raw counts -> %s (prior: %d, adding: %d)\n",
            base_seed, raw_file, prior_trials, test_count);

    /* Collect counts in trial-index order so the raw file is deterministic
     * regardless of thread completion order. */
    uint64_t *trial_counts = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)test_count);
    if (!trial_counts) { perror("malloc trial_counts"); return 1; }

    #pragma omp parallel
    {
        Element  *struct_arr = (Element  *)malloc(sizeof(Element)  * STRUCT_SIZE);
        uint64_t *head       = (uint64_t *)malloc(sizeof(uint64_t) * BUCKET_COUNT);
        uint64_t *next       = (uint64_t *)malloc(sizeof(uint64_t) * STRUCT_SIZE);
        if (!struct_arr || !head || !next) {
            #pragma omp critical
            fprintf(stderr, "OOM in thread %d\n", omp_get_thread_num());
            exit(1);
        }

        #pragma omp for schedule(static)
        for (int t = 0; t < test_count; t++) {
            unsigned int seed = rmd_derive_trial_seed(base_seed, prior_trials + t, rounds);
            uint8_t mk[16], base_pt[16];
            for (int k = 0; k < 16; k++) mk[k]      = (uint8_t)rmd_rand_byte(&seed);
            for (int k = 0; k < 16; k++) base_pt[k] = (uint8_t)rmd_rand_byte(&seed);
            __m128i round_keys[11];
            rmd_aes128_generate_round_keys(mk, round_keys);

            uint64_t c = run_single_trial(struct_arr, head, next,
                                          round_keys, rounds, base_pt);
            trial_counts[t] = c;

            #pragma omp critical
            fprintf(stderr, "Trial %3d/%d: count %llu\n",
                    prior_trials + t + 1,
                    prior_trials + test_count,
                    (unsigned long long)c);
        }
        free(struct_arr);
        free(head);
        free(next);
    }

    /* Append in deterministic trial-index order. */
    {
        FILE *fa = fopen(raw_file, "a");
        if (!fa) { perror(raw_file); free(trial_counts); return 1; }
        for (int t = 0; t < test_count; t++)
            fprintf(fa, "%llu\n", (unsigned long long)trial_counts[t]);
        fclose(fa);
    }
    free(trial_counts);

    /* Re-read raw_file (this run + prior) for cumulative summary. */
    int   total_trials = 0;
    int   max_count    = 0;
    int   cap          = 256;
    int  *all_counts   = (int *)malloc(sizeof(int) * cap);
    if (!all_counts) { perror("malloc all_counts"); return 1; }
    {
        FILE *fr = fopen(raw_file, "r");
        if (!fr) { perror(raw_file); free(all_counts); return 1; }
        int v;
        while (fscanf(fr, "%d", &v) == 1) {
            if (total_trials >= cap) {
                cap *= 2;
                int *tmp = (int *)realloc(all_counts, sizeof(int) * cap);
                if (!tmp) { perror("realloc"); free(all_counts); fclose(fr); return 1; }
                all_counts = tmp;
            }
            all_counts[total_trials++] = v;
            if (v > max_count) max_count = v;
        }
        fclose(fr);
    }

    int *freq = (int *)calloc((size_t)max_count + 1, sizeof(int));
    if (!freq) { perror("calloc"); free(all_counts); return 1; }
    for (int t = 0; t < total_trials; t++) freq[all_counts[t]]++;

    printf("Config: Tests=%d, StructureSize=2^%d, Rounds=%d\n\n",
           total_trials, STRUCT_BITS, rounds);
    printf("--- %d-Round Frequency Results ---\n", rounds);
    for (int i = 0; i <= max_count; i++) {
        if (freq[i] > 0)
            printf("Count %3d : %4d times\n", i, freq[i]);
    }
    printf("\n(Total: %d trials)\n", total_trials);

    free(freq);
    free(all_counts);
    return 0;
}
