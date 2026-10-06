/*
 * count_quartets_nz0.c
 * AES-128 valid-quartet count of the pairing with base pair (P^0, P^3) in a
 * Family-B quartet of Table 4 (n_z = 0).  Standard-AES
 * (w = 8) counterpart of 5r_small_AES/count_quartets_nz0.c, built on the
 * 2^32 array + linked-list bucketing engine of
 * 5r_AES_N32/count_quartets_nz2.c.
 *
 * n_z = 0 is the Family-B pairing (P^0, P^3) whose base pair is 2-active on
 * bytes {5,15} and equal on {0,10}.  Counts valid quartets per trial over a
 * full structure of 2^STRUCT_BITS plaintexts (default 2^32; AES-NI +
 * OpenMP), verifying the n_z = 0 instance of Theorem 10 (multiple-of property) on the full AES:
 *      N = 2^{2-n_z}(2^w)^{n_z} A + 8B = 4 A + 8B     (w = 8, n_z = 0),
 *   i.e. the bundle size is 2^{2-0}(2^8)^0 = 4.
 * Paper: Table 10 (Section 3.4.1) and the AES column of Table 11 (Section 3.4.2),
 * from Results/5r_AES_N32/raw_nz0_{5,10}r.txt; Codes/tools/summarize_counts.py prints
 * the mean, the detection rate and the sample standard deviation quoted in Section 3.4.1.
 *
 * Index encoding (identical to count_quartets_nz2.c):
 *   idx = byte0 | byte10<<BB | byte5<<2BB | byte15<<3BB,   BB = STRUCT_BITS/4
 *   MASK_LOW  = low  2*BB bits = {byte0, byte10}   (the {0,10} half)
 *   MASK_HIGH = high 2*BB bits = {byte5, byte15}   (the {5,15} half)
 *
 * n_z = 0 quartet (mirrors 5r_small_AES/count_quartets_nz0.c):
 *   base pair    (i0, i1) : equal on {0,10}, 2-active on {5,15}
 *                           (i0 ^ i1 lies entirely in MASK_HIGH, with both
 *                            byte5 and byte15 active).
 *   partner pair (i0^g, i1^g) with g = sigma | byte15diff, where
 *                  byte15diff = (i0 ^ i1) & BYTE15_MASK   (fixed by the base
 *                               pair -- byte15 is cross-swapped, byte5 kept),
 *                  sigma      ranges over all NONZERO {0,10} offsets
 *                               (1 .. MASK_LOW; i.e. every other i-slice).
 * Valid when BOTH pairs collide (zero difference) on the SAME inverse
 * diagonal.  Only the two "horizontal" pairings (i0,i1) and (i0^g,i1^g) are
 * valid n_z=0 base pairs of the 4 plaintexts, so each unique quartet is
 * enumerated twice; the raw pair count is halved, matching
 * count_quartets_nz2.c and the paper's unique-quartet convention.
 *
 * Number of partner offsets per base pair: MASK_LOW = 2^(2*BB) - 1 = 65535
 * at BB = 8.  Expected work per trial ~ 2-3x count_quartets_nz2.c.
 *
 * For a quick pipeline smoke-test on a small machine, build with e.g.
 * -DSTRUCT_BITS=16 -DBUCKET_BITS=16 (this does NOT reproduce the exact
 * multiple-of bundle, which requires the full 2^32 diagonal space).
 *
 * Usage:
 *   ./count_quartets_nz0 <rounds> [test_count] [raw_file] [seed]
 * <rounds> is 5 (5-round AES) or 10 (random-permutation baseline).
 * Defaults: test_count = 100, raw_file = raw_counts_nz0_<rounds>r.txt,
 * seed = 42 (env SEED also honoured).  The committed results use 700 trials
 * at 5 rounds and 100 trials at 10 rounds.
 *
 * Memory (per OpenMP thread, default STRUCT_BITS=32):
 *   struct_arr.ct : 2^32 * 16 B = 64 GB
 *   dv (inv-diag) : 2^32 *  4 B = 16 GB   (compact per-diagonal memoization)
 *   head          : 2^32 *  8 B = 32 GB
 *   next          : 2^32 *  8 B = 32 GB
 *   total         : ~144 GB per thread.  On a 2 TB host set
 *   OMP_NUM_THREADS to ~13 to run ~13 trials concurrently.
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

/* Field layout inside the 32-bit index (see header):
 *   byte0  = bits [0,        BYTE_BITS)
 *   byte10 = bits [BYTE_BITS, 2*BYTE_BITS)
 *   byte5  = bits [2*BB,      3*BB)
 *   byte15 = bits [3*BB,      4*BB)
 */
#define MASK_LOW   ((1ULL << (2*BYTE_BITS)) - 1ULL)   /* {byte0, byte10} */
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
                                 uint32_t *dv,
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

    /* Phase 2: per inverse diagonal, bucket by the precomputed inv-diag value. */
    uint64_t count = 0;
    for (int d = 0; d < 4; d++) {
        /* Memoize this diagonal's 32-bit inverse-diagonal value into the
         * compact dv[] array (4 B/elt = 16 GB).  All subsequent random
         * accesses then touch dv instead of the 64 GB ct array -- 4x smaller
         * working set, far fewer TLB misses.  dv[i] == inv_diag(ct[i],d)
         * exactly, so results are identical to the ct-based version. */
        #pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < STRUCT_SIZE; i++)
            dv[i] = rmd_aes_inv_diag_val(struct_arr[i].ct, d);

        #pragma omp parallel for schedule(static)
        for (uint64_t v = 0; v < BUCKET_COUNT; v++) head[v] = SENTINEL;

        for (uint64_t i = 0; i < STRUCT_SIZE; i++) {
            uint64_t b = (uint64_t)dv[i] & BUCKET_MASK;
            next[i] = head[b];
            head[b] = i;
        }

        for (uint64_t b = 0; b < BUCKET_COUNT; b++) {
            for (uint64_t i0 = head[b]; i0 != SENTINEL; i0 = next[i0]) {
                uint32_t v0 = dv[i0];
                for (uint64_t i1 = next[i0]; i1 != SENTINEL; i1 = next[i1]) {
                    /* Real inv-diag equality (always true at BUCKET_BITS=32). */
                    if (dv[i1] != v0) continue;

                    /* n_z = 0 base pair: equal on {0,10}, 2-active on {5,15}. */
                    uint64_t diff = i0 ^ i1;
                    if ((diff & MASK_LOW) != 0)     continue;   /* {0,10} must match */
                    if ((diff & BYTE5_MASK)  == 0 ||
                        (diff & BYTE15_MASK) == 0)  continue;   /* {5,15} both active */

                    /* Partner: byte5 kept, byte15 cross-swapped (fixed offset),
                     * {0,10} translated by any nonzero sigma (every other slice). */
                    uint64_t byte15diff = diff & BYTE15_MASK;
                    for (uint64_t sigma = 1; sigma <= MASK_LOW; sigma++) {
                        uint64_t g = sigma | byte15diff;
                        if (dv[i0 ^ g] == dv[i1 ^ g]) count++;
                    }
                }
            }
        }
    }
    /* Each unique quartet is enumerated twice; halve to match Theorem 10 (multiple-of property). */
    return count / 2;
}

int main(int argc, char **argv) {
    setbuf(stdout, NULL);

    if (argc < 2 || argc > 5) {
        fprintf(stderr, "Usage: %s <rounds> [test_count] [raw_file] [seed]\n", argv[0]);
        fprintf(stderr, "  Default: test_count=100, raw_file=raw_counts_nz0_<rounds>r.txt, seed=42 (env SEED also honoured).\n");
        fprintf(stderr, "  Committed results: 700 trials (5 rounds), 100 trials (10 rounds).\n");
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
    snprintf(default_raw, sizeof(default_raw), "raw_counts_nz0_%dr.txt", rounds);
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

    fprintf(stderr, "[*] %d-Round AES Experiment, n_z=0 (2^%d structure, bundle 4)\n",
            rounds, STRUCT_BITS);
    fprintf(stderr, "[*] Seed=%u, Raw counts -> %s (prior: %d, adding: %d)\n",
            base_seed, raw_file, prior_trials, test_count);

    uint64_t *trial_counts = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)test_count);
    if (!trial_counts) { perror("malloc trial_counts"); return 1; }

    int n_threads = omp_get_max_threads();
    double wall_start = omp_get_wtime();

    #pragma omp parallel
    {
        Element  *struct_arr = (Element  *)malloc(sizeof(Element)  * STRUCT_SIZE);
        uint32_t *dv         = (uint32_t *)malloc(sizeof(uint32_t) * STRUCT_SIZE);
        uint64_t *head       = (uint64_t *)malloc(sizeof(uint64_t) * BUCKET_COUNT);
        uint64_t *next       = (uint64_t *)malloc(sizeof(uint64_t) * STRUCT_SIZE);
        if (!struct_arr || !dv || !head || !next) {
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

            double tt0 = omp_get_wtime();
            uint64_t c = run_single_trial(struct_arr, dv, head, next,
                                          round_keys, rounds, base_pt);
            double tt1 = omp_get_wtime();
            trial_counts[t] = c;

            #pragma omp critical
            fprintf(stderr, "Trial %3d/%d: count %llu  (%.1f s/trial)\n",
                    prior_trials + t + 1,
                    prior_trials + test_count,
                    (unsigned long long)c, tt1 - tt0);
        }
        free(struct_arr);
        free(dv);
        free(head);
        free(next);
    }

    double wall_end = omp_get_wtime();
    double wall  = wall_end - wall_start;
    int    waves = (test_count + n_threads - 1) / n_threads;   /* parallel waves run */
    double per_wave = wall / (waves > 0 ? waves : 1);          /* ~ one trial's duration */
    fprintf(stderr,
            "[*] Done: %d trials on %d threads in %.1f s (%.2f min); ~%.1f s per trial (per wave).\n"
            "[*] At this rate: 100 trials ~%.2f h, 200 trials ~%.2f h (per pairing/round cell).\n",
            test_count, n_threads, wall, wall / 60.0, per_wave,
            ((100 + n_threads - 1) / n_threads) * per_wave / 3600.0,
            ((200 + n_threads - 1) / n_threads) * per_wave / 3600.0);

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

    printf("Config: n_z=0, Tests=%d, StructureSize=2^%d, Rounds=%d, bundle=4\n\n",
           total_trials, STRUCT_BITS, rounds);
    printf("--- %d-Round Frequency Results (n_z=0) ---\n", rounds);
    for (int i = 0; i <= max_count; i++) {
        if (freq[i] > 0)
            printf("Count %6d : %4d times\n", i, freq[i]);
    }
    printf("\n(Total: %d trials)\n", total_trials);

    free(freq);
    free(all_counts);
    return 0;
}
