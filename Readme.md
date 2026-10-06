# A Multiple-of Property and Bundle Structure for Related-Differential Distinguishers on 5-Round AES

Reference implementation and experiment outputs for the paper
*"A Multiple-of Property and Bundle Structure for Related-Differential
Distinguishers on 5-Round AES"* (submitted to Designs, Codes and Cryptography).
A preliminary version is available as IACR Cryptology ePrint Archive
Report 2026/1286 (https://eprint.iacr.org/2026/1286).

All reported counts are unique-quartet counts.  Each quartet is enumerated
twice as ordered pairs during the collision search, so the C binaries divide
the raw pair count by two before printing summaries.  Pair search uses
direct-addressed linked-list bucketing on the inverse-diagonal value (O(N)),
and the partner pair is recovered from the index without re-encryption.

## Layout

```text
Codes/
  common/                 Shared AES-NI and deterministic-seed helpers.
  expected_count/         DDT computation of Sections 2.4 and 3.3 (Python):
                          ddt_model.py, table3_rho_m.py, table5_expected_count.py.
  5r_small_AES/           Small-AES (4-bit cells, small-AES key schedule), structure 2^16.
                          count_quartets_nz2 (n_z=2), _nz1 (n_z=1), _nz0 (n_z=0).
  5r_AES_N32/             AES full diagonal structure, N = 2^32.
                          count_quartets_nz2 (n_z=2), _nz1 (n_z=1), _nz0 (n_z=0).
  5r_AES_distinguishers/  5-round AES CP and ACP distinguishers:
                          cp_nz01 (CP, n_z=0 and n_z=1), cp_nz2 (CP, n_z=2),
                          acp_all (ACP, one pairing per run),
                          cp_nz012, cp_nz012f (CP, the three pairings on one structure).
  model/                  Compound-Poisson model of Sections 3.4, 4.1 to 4.3 and 5 (Python):
                          bundle_geometry.py, eps_split.py, rho_pairing_c.py,
                          acp_budget_sim.py, cpoisson_model.py; inputs and outputs in data/.
  table2_search/          Exhaustive check that Table 2 is complete (Section 2.2), with the
                          trail counts of Lemma 7 and the weights of Table 4:
                          lemma7_search.c, run_all.sh.
  tools/                  summarize_counts.py (summary of a count file),
                          aggregate_combined.py (summary of the Section 4.3 runs).
Results/                  Experiment outputs reported in the paper.
  expected_count/         Tables 3 and 5 (JSON).
  5r_small_AES/           Frequency tables for n_z = 2, 1, 0 (5 and 10 rounds).
  5r_AES_N32/             Raw per-structure counts for n_z = 2, 1, 0 (5 and 10 rounds).
  5r_AES_distinguishers/  CP and ACP distinguisher runs of Tables 12 and 14.
    combined/             Per-chunk runs of cp_nz012 / cp_nz012f for Table 13, their
                          summary AGGREGATE.txt, and the cross-checks verify_*.txt.
  table2_search/          Output of the Table 2 completeness search.
```

## Paper tables and the files that produce them

Paths in the Code column are relative to `Codes/`.

| Paper | Content | Code | Result |
|---|---|---|---|
| Table 1 (Section 1) | rows "Related diff. (three pairings)": model 63% at 2^26.33 CP with time 2^27.18 E and 90% at 2^27.02 CP with time 2^27.87 E, measured 65% and 91% (footnote 2), also stated in the abstract (with the false-positive probability 2^-9.8), Section 1 and Section 6; 61.8% of the model for Yan et al. (footnote 1); the measurements 63.3% and 82.5% of Yan et al. (footnote 3) against the model | `model/cpoisson_model.py`; `tools/aggregate_combined.py` | `Codes/model/data/cpoisson_results.txt`: section (4b), rows 2^26.33 and 2^27.02 (union model, union measured, false positive, time), block "Yan et al.'s CP" and block "Yan et al.'s ACP" of section (2b); `Results/5r_AES_distinguishers/combined/AGGREGATE.txt` (line `P_find` of the union of a, b and c) |
| Section 2.2 | exhaustive check that Table 2 is complete for 4 <= w <= 8 | `table2_search/lemma7_search.c` (driver `run_all.sh`) | `Results/table2_search/results_mode*.txt` |
| Lemma 7 (Section 3.1) | trail counts 4(2^w-1), 3(2^w-1) and 2^w-1, checked exhaustively for 4 <= w <= 8 | `table2_search/lemma7_search.c` (driver `run_all.sh`) | `Results/table2_search/results_mode{0,1}.txt`, lines "(unordered triples = Lemma 7 count)" (w = 8: 1020, 765, 255) |
| Table 4 (Section 3.2) | weights at P/Z_1 and at X_2 of the three base pairs of Families A and B | `table2_search/lemma7_search.c` | the same files, lines "(input zero set: wt_in -> wt_out of MC)": Family A {}: 4 -> 2, {0,2} and {1,3}: 2 -> 3; Family B {0} and {2}: 3 -> 2, {1,3}: 2 -> 4 |
| Table 3 (Section 2.4) | mean of rho_j over all differences with m nonzero bytes | `expected_count/table3_rho_m.py` | `Results/expected_count/table3_w{4,8}.json` |
| Table 5 (Section 3.3) | rho-bar(n_z), epsilon(n_z), E[N_q] of Theorem 11 | `expected_count/table5_expected_count.py` | `Results/expected_count/table5_w{4,8}.json` |
| Tables 6-9 (Section 3.4.1) | small-AES bundles, detection rates, frequency tables (Theorem 10) | `5r_small_AES/count_quartets_nz*` | `Results/5r_small_AES/*.txt` |
| Table 10 (Section 3.4.1) | AES bundles and detection rates (Theorem 10) | `5r_AES_N32/count_quartets_nz*` | `Results/5r_AES_N32/raw_nz*_*.txt` |
| Section 3.4.1 text | expected bundles per structure 75.71/256, 59.02/32, 45.05/4; standard deviations 14, 44, 139 (small-AES) and 723, 55 (AES) of the model k T_1 + 8 eps against the measured 15, 47, 143 and 752, 57 (the model values are the `8eps:` values of section (3), not its main `SD model` column, which gives 55.5 for AES n_z = 0; the measured values are the sample standard deviations with ddof = 1, printed as 15.46, 46.83, 143.06, 752.49 and 56.77 by `tools/summarize_counts.py` and to one decimal in section (3)) | `expected_count/table5_expected_count.py`; `model/cpoisson_model.py`; `tools/summarize_counts.py` | `Results/expected_count/table5_w4.json` (`first_term`); `Codes/model/data/cpoisson_results.txt`, section (3); line `std deviation` of `summarize_counts.py` |
| Table 11 (Section 3.4.2) | small-AES and AES means (Theorem 11) and the random-permutation expected row | the count programs above; `tools/summarize_counts.py` | the result files above; the line `random-perm E[N_q]` of `summarize_counts.py`; `Results/expected_count/table5_w{4,8}.json`, field `E`, for the row "Theorem 11" |
| Section 3.4.2 text | bundles in the 700 AES structures: 1,461 against 700 T_1/512 = 1,394 for n_z = 1, and 6 against 11 for n_z = 2; the deviations 3.6%, 3.0%, 0.2%, 4.8% and 0.08 follow in one line from the Table 11 means and the E values | `tools/summarize_counts.py`; `model/cpoisson_model.py` | line `bundles` of `summarize_counts.py`; `Codes/model/data/cpoisson_results.txt`, section (3), line "bundles in the 700 AES structures"; field `E` of `Results/expected_count/table5_w{4,8}.json` |
| Section 4.1 text | a candidate of 5-round AES is valid with probability 2^-62 times 2^10.04, 2^10.02, 2^9.60 for n_z = 2, 1, 0; the mean of rho over all base pairs is 1 to within 10^-5 | `model/cpoisson_model.py` | `Codes/model/data/cpoisson_results.txt`, block `constants`, lines "(Section 4.1)" |
| Table 12 (Section 4.2) | CP success probabilities and model predictions | `5r_AES_distinguishers/cp_nz01`, `cp_nz2`; `model/cpoisson_model.py` | `Results/5r_AES_distinguishers/cp_*.txt`; `Codes/model/data/cpoisson_results.txt`, section (2) |
| Section 4.2 text | q = 2^-9.98 and 2^-9.61 of the two product grids; all six predictions of Table 12 within 1.4 standard errors; Poisson 0.64 and 0.73 for n_z = 1 more than three standard errors above the measured values; 61.8% of the model on the grid of Yan et al., against 73% for a Poisson count | `model/cpoisson_model.py` | `Codes/model/data/cpoisson_results.txt`, block "(i)-(iv)", section (2) (columns P_cpois, Poisson, z_cp, z_Poi) and block "Yan et al.'s CP" |
| Table 13 (Section 4.3) | three pairings on one structure, their union, larger structures, false positives, 10-round run | `5r_AES_distinguishers/cp_nz012` (`cp_nz012f`); `tools/aggregate_combined.py`; `model/cpoisson_model.py` | `Results/5r_AES_distinguishers/combined/` (summary in `AGGREGATE.txt`); `Codes/model/data/cpoisson_results.txt`, section (4) for the model columns and false positives and section (4b) for model against measured, z-scores and time |
| Section 4.3 text | T_1 = 1019.94 of the pairing n_z = 1 on {5, 15}; time 2^27.18 and 2^27.87; 1.1 false positives expected in the 10-round run; bound 0.88 of a pairing n_z = 1 in one structure; predictions within 1.3, 2.3 and 1.7 standard errors (largest absolute z-scores 1.23, 2.26 and 1.65 in section (4b)); about 2^27.8 plaintexts for 90% with n_z = 0 (2^27.821 on the proportional grid; the 2^27.80 run grid of Table 13 gives the model value 0.893), and no 90% for a pairing n_z = 1 in one structure (the model saturates at 0.879) | `model/rho_pairing_c.py`; `model/cpoisson_model.py` | `Codes/model/data/rho_pairing_c.json`; `Codes/model/data/cpoisson_results.txt`, sections (4b) and (5); section (4), block "data N (proportional grid, one structure) for a target P_find" |
| Table 14 (Section 5.2) | ACP success probabilities and model predictions; about 2^16.6 adaptive swaps for n_z = 2 at 2^23.32 (Section 5.2 text) | `5r_AES_distinguishers/acp_all`; `model/acp_budget_sim.py` (n_z = 1, 0) and `acp_nz2()` in `model/cpoisson_model.py` (n_z = 2) | `Results/5r_AES_distinguishers/acp_all_2p2332.txt`; model values in `Codes/model/data/cpoisson_results.txt`, section (2b) |
| Section 5.1 text | n_z = 0 with 2^5.4 base structures, 2^23.4 adaptive partner queries and 2^23.7 in total; n_z = 2 in the ACP setting: 1.7% at 2^26.5, 89 base structures of about 2^23.65 plaintexts, 2^30.15 in total, with 2^16 q of about 1.2, 0.71 and about 1/89 per base structure | `model/cpoisson_model.py` | `Codes/model/data/cpoisson_results.txt`, sections (2b) and (6) |

## Building

Each C experiment is self-contained under its own `Codes/<exp>/` directory and
builds with `make`.  Each Makefile sets its own flags: `-O3 -maes -mpclmul
-fopenmp -Wall` for N32, `-O3 -march=native` for the distinguishers
(`cp_nz012f` uses AES-NI through `common/rmd_aesni.h`), and `-O3 -Wall` for
small-AES.  The Table 2 search has no Makefile; `run_all.sh` builds it in a
temporary directory.  Native Windows builds are intentionally unsupported
because MinGW-style builds create `.exe` artifacts; build under Linux or WSL.

```sh
(cd Codes/5r_small_AES          && make)
(cd Codes/5r_AES_N32            && make)
(cd Codes/5r_AES_distinguishers && make)
gcc -O3 -march=native Codes/table2_search/lemma7_search.c -o Codes/table2_search/lemma7_search
```

The Python scripts need Python 3.8 or later with NumPy (`summarize_counts.py`
and `aggregate_combined.py` need only the standard library).

## Running

The count programs of `5r_small_AES` and `5r_AES_N32` take the seed as the
last, optional argument and default to `42` (the environment variable `SEED`
is also honoured).  Per-trial seeds are derived from
`(base_seed, trial_index, rounds)` for N32, so results do not depend on the
OpenMP completion order, and from `(base_seed, trial_index)` for small-AES.
The distinguishers `cp_nz01`, `cp_nz2`, `acp_all`, `cp_nz012` and
`cp_nz012f` are single-threaded and seed trial `t` with
`seed ^ 0x9E3779B97F4A7C15*(t+1)`, so each trial depends only on the seed
and `t`.  `cp_nz012` and `cp_nz012f` also take a first trial index
`t_start`, so their runs can be split into chunks of trials without changing
any trial; `cp_nz01`, `cp_nz2` and `acp_all` always start at `t = 0`.  The
seed is optional (default `42`) for `cp_nz01`, `cp_nz2` and `acp_all` and
required for `cp_nz012` and `cp_nz012f`.  `lemma7_search` is deterministic
and takes no seed.
`<rounds>` is `5` for 5-round AES (the distinguisher target) or `10` for the
random-permutation baseline.  Progress is printed to `stderr`; summaries are
printed to `stdout`.  The 5-round and 10-round ciphers omit MixColumns in the
last round, so a zero difference on an inverse diagonal of the ciphertext is
tested directly.

### Table 2 completeness (Section 2.2)

```sh
bash Codes/table2_search/run_all.sh <out_dir>                # about 1 minute
./Codes/table2_search/lemma7_search <w> <poly> <mode> [raw]  # one field
```

`run_all.sh` builds `lemma7_search` in a temporary directory and runs it for
(w, poly) = (4, 0x13), (5, 0x25), (6, 0x43), (7, 0x83) and (8, 0x11B).  Mode 0
requires exactly one of `dx`, `dx'`, `dx + dx'` to be zero at every byte of the
input and the output, the form of Table 2; mode 1 requires at least one zero,
as in Definition 1.  `raw = 1` adds an enumeration without scalar
normalization as a cross-check (used for w <= 6 in mode 0).  The script writes
`results_mode0.txt` and `results_mode1.txt` to `<out_dir>`; they agree with
`Results/table2_search/results_mode{0,1}.txt` except for the timing lines
(`diff <(grep -vE 'time|wall' a) <(grep -vE 'time|wall' b)`).  The line `set equality: YES`
means that the solutions found are exactly the scalar multiples of byte
rotations of the rows of Table 2, in any of the six orders of `dx`, `dx'` and
`dx + dx'`, so Table 2 is complete for that w.  For each class of solutions
the files also give the number of unordered triples, which is the trail count
of Lemma 7 (`4 x (2^w-1)`, `3 x (2^w-1)` or `1 x (2^w-1)`), and the weights
`wt_in -> wt_out` of the three differences at the input and the output of
MixColumns, which are the weights at P/Z_1 and X_2 of Table 4.

### Expected number of valid quartets (Sections 2.4 and 3.3)

`ddt_model.py` computes the collision ratio rho_j of Section 2.4 from the DDT
of the S-box, assuming independent and uniformly distributed round keys.
The sum over the second-round S-box output differences Delta Y_2 runs over all
possible values, or over 1,000 values drawn with the DDT probabilities where
the paper says so.

```sh
cd Codes/expected_count
python3 table5_expected_count.py 4    # about 4 minutes
python3 table5_expected_count.py 8    # about 13 minutes
python3 table3_rho_m.py 4             # about 6 minutes
python3 table3_rho_m.py 8             # about 13 minutes
```

Each script writes `table5_w<w>.json` or `table3_w<w>.json` into the current
directory.  The fixed seeds in the scripts reproduce the committed files and
the values printed in the paper; compare, for example, with
`cmp table5_w4.json ../../Results/expected_count/table5_w4.json`.

### Small-AES (Section 3.4)

Small-AES is SR(n,4,4,4) of Cid, Murphy and Robshaw: the AES round structure
with 4-bit cells, the MixColumns matrix `[2,1,1,3]` over GF(2^4) and the
small-AES key schedule.  The whitening key is the 64-bit master key and the
round keys follow the key schedule.  (An earlier version of this artifact
added the master key in every round; the committed results were regenerated
with the key schedule.)  The programs count the valid quartets of one pairing
over a full diagonal structure of 2^16 plaintexts and verify the multiple-of
form `N = 2^{2-n_z}(2^w)^{n_z} A + 8B` at `w = 4` (bundle 256, 32, 4 for
n_z = 2, 1, 0).  The pairing n_z = 1 of `count_quartets_nz1` is the pairing
n_z = 1 on {0, 10} of Yan et al. (base pairs differ at cells 0 and 10,
partner pairs change cells 5 and 15).

```sh
./Codes/5r_small_AES/count_quartets_nz2 <rounds> [test_count] [seed]
./Codes/5r_small_AES/count_quartets_nz1 <rounds> [test_count] [seed]
./Codes/5r_small_AES/count_quartets_nz0 <rounds> [test_count] [seed]
```

The committed results under `Results/5r_small_AES/` use 1,000 structures at
seed 42 for each pairing and each number of rounds.

### Full AES, N = 2^32 (Section 3.4)

The programs use the AES-128 key schedule with a fresh random key for every
structure, count the valid quartets of one pairing over the full diagonal
structure of 2^32 plaintexts, append one count per structure to `<raw_file>`,
and print the cumulative summary.  The pairing n_z = 1 of
`count_quartets_nz1` is the pairing n_z = 1 on {0, 10} of Yan et al. (base
pairs differ at bytes 0 and 10, partner pairs change bytes 5 and 15).

```sh
./Codes/5r_AES_N32/count_quartets_nz2 <rounds> [test_count] [raw_file] [seed]
./Codes/5r_AES_N32/count_quartets_nz1 <rounds> [test_count] [raw_file] [seed]
./Codes/5r_AES_N32/count_quartets_nz0 <rounds> [test_count] [raw_file] [seed]
```

Runs are resumable: the trial index starts at the number of counts already in
`<raw_file>`, so rerunning with the same file and seed extends one seed
sequence.  The committed files under `Results/5r_AES_N32/` use seed 42:

| File | Pairing | Rounds | Structures |
|---|---|---|---|
| `raw_nz2_5r.txt`  | n_z = 2 | 5  | 700 |
| `raw_nz1_5r.txt`  | n_z = 1 | 5  | 700 |
| `raw_nz0_5r.txt`  | n_z = 0 | 5  | 700 |
| `raw_nz2_10r.txt` | n_z = 2 | 10 | 200 |
| `raw_nz1_10r.txt` | n_z = 1 | 10 | 100 |
| `raw_nz0_10r.txt` | n_z = 0 | 10 | 100 |

`test_count` (default 100) is the number of structures added to
`<raw_file>`; without `raw_file` the counts go to `raw_counts_<rounds>r.txt`
(n_z = 2) or `raw_counts_nz<k>_<rounds>r.txt` (n_z = 1, 0).  To reproduce the
committed files, run in an empty directory, because the programs append to
`<raw_file>` and resume from its length, then `cmp` each file against
`Results/5r_AES_N32/`:

```sh
N32=<path to>/Codes/5r_AES_N32
$N32/count_quartets_nz2 5  700 raw_nz2_5r.txt  42
$N32/count_quartets_nz1 5  700 raw_nz1_5r.txt  42
$N32/count_quartets_nz0 5  700 raw_nz0_5r.txt  42
$N32/count_quartets_nz2 10 200 raw_nz2_10r.txt 42
$N32/count_quartets_nz1 10 100 raw_nz1_10r.txt 42
$N32/count_quartets_nz0 10 100 raw_nz0_10r.txt 42
```

Memory at the default `STRUCT_BITS=32` and `BUCKET_BITS=32` is about 128 GB
(n_z = 2) and 144 GB (n_z = 1, 0) per OpenMP thread; set `OMP_NUM_THREADS`
to the installed memory divided by the per-thread size.  For small-machine
testing, rebuild with reduced compile-time parameters (this does not
reproduce the bundle, which needs the full 2^32 structure), e.g.

```sh
cd Codes/5r_AES_N32
gcc -O3 -maes -mpclmul -fopenmp -Wall \
    -DSTRUCT_BITS=24 -DBUCKET_BITS=24 \
    count_quartets_nz2.c -o count_quartets_nz2 -lm
```

### 5-Round AES distinguishers (Sections 4 and 5)

`Codes/5r_AES_distinguishers/` holds the chosen-plaintext (CP) and
adaptively chosen-plaintext (ACP) distinguishers that compare the three
pairings `n_z = 0, 1, 2` on AES (`w = 8`) with the AES-128 key schedule.

| Program | Setting | Pairings | Paper |
|---|---|---|---|
| `cp_nz01` | CP, one fixed structure | n_z=0, n_z=1 on {0,10} | Section 4.2 |
| `cp_nz2`  | CP, one fixed structure | n_z=2 | Section 4.2 |
| `cp_nz012` (`cp_nz012f`: AES-NI build; same `T` lines and `# P_find` footer, only the program name in the header line differs) | CP, one fixed structure | n_z=1 on {0,10} (Yan et al.), n_z=1 on {5,15}, n_z=0, and their union | Section 4.3 |
| `acp_all` | ACP (adaptive queries) | n_z=0, 1, 2 (one per run) | Section 5.2 |

```sh
./cp_nz01   s0 s1 s2 rounds [trials] [seed]     # structure {0,10}=s0 x {5}=s1 x {15}=s2; trials default 200
./cp_nz2    s0 s1 rounds [trials] [seed]        # structure {0,10}=s0 x {5,15}=s1; trials default 100
./cp_nz012  s0 s1 s2 rounds trials seed [t_start] [verify]   # trials t_start .. t_start+trials-1
./cp_nz012f s0 s1 s2 rounds trials seed [t_start]            # same T lines and footer as cp_nz012 (header names cp_nz012f), faster
./acp_all   nz Dexp rounds [trials] [seed]      # ACP, total data complexity 2^Dexp; trials default 200
```

`cp_nz01`, `cp_nz2` and `acp_all` print the per-structure success
probability `P_find` (fraction of trials with at least one valid quartet) and
the mean valid-quartet count.  `cp_nz01` and `acp_all` also print diagnostic
counters (base collisions, partner-pair checks, and the per-quartet rate
`RV/checks`); `acp_all` adds a histogram of the counts.  `cp_nz012` and
`cp_nz012f` print one `T` line per trial, described in the next subsection.

| File | Experiment |
|---|---|
| `cp_nz01_2p27.txt`   | CP, n_z=0 and n_z=1 on {0,10}, one structure of about 2^27.02 plaintexts, 300 trials |
| `cp_nz01_2p272.txt`  | CP, n_z=0 and n_z=1 on {0,10}, one structure of about 2^27.20 plaintexts, 300 trials |
| `cp_nz2_2p27.txt`    | CP, n_z=2, one structure of about 2^27.02 plaintexts, 300 trials |
| `cp_nz2_2p272.txt`   | CP, n_z=2, one structure of about 2^27.20 plaintexts, 300 trials |
| `acp_all_2p2332.txt` | ACP, n_z=0/1/2 at data complexities 2^22 and 2^23.32, 300 trials, seed 2024 |

The CP files of Table 12 were produced by

```sh
./cp_nz01 11664 108 108 5 300 42 > cp_nz01_2p27.txt
./cp_nz01 12515 111 111 5 300 42 > cp_nz01_2p272.txt
./cp_nz2  11664 11664 5 300      > cp_nz2_2p27.txt
./cp_nz2  12418 12418 5 300      > cp_nz2_2p272.txt
```

The seed 42 of `cp_nz01` is confirmed by the seed-42 runs of `cp_nz012`,
which give the same `P_find` and means for both structures; the `cp_nz2`
files do not record the seed, and the commands above use the default 42.
The ACP results use seed 2024, not the default 42:

```sh
for D in 22 23.32; do for nz in 1 0 2; do ./acp_all $nz $D 5 300 2024; done; done
```

The committed files keep the summary lines of standard output.  The final
`ALL_DONE` line, and in the ACP file the run banner
(`=== ACP CLEAN (2-active) T=300 seed 2024 ===`, `-- D=... --`), were written
by the shell loop that ran the programs, not by the programs.  The `[diag]`
and `[hist]` lines of `acp_all` are omitted, as are the `[diag]` lines of
`cp_nz01_2p272.txt`, so compare a rerun with a committed file on the summary
lines only (`ACP nz=...`, `AES ...`, `  n_z=...: P_find=...`,
`  P_find=...`).  On the machine used for the model timings below, one
single-threaded ACP run of 300 trials takes under a minute for n_z = 1 at
2^22 and about 6 minutes for n_z = 2 at 2^23.32.

#### Three pairings on one structure (Section 4.3)

`cp_nz012` extends `cp_nz01` by the pairing n_z = 1 on {5, 15} (base pairs
that differ at bytes 5 and 15, partner pairs that change bytes 0 and 10) and
reports, per trial, the counts of the three pairings.  In the outputs,
pairing `a` is n_z = 1 on {0, 10} (Yan et al.), `b` is n_z = 0 and `c` is
n_z = 1 on {5, 15}; Table 13 lists them in the order a, c, b.  A `T` line
reads `T t ca cb cc rawA rawB rawC baseA baseB P2both candC`: the three
counts, the raw pair counts (twice the counts), the base collisions of the
pairings a and b (b and c share their base pairs), the number of pairs of
{0,10} values that differ at both bytes, and the number of (base collision,
partner) checks of pairing c.  `cp_nz012` draws keys and plaintexts exactly as
`cp_nz01`, so with seed 42 its pairings a and b give the results of
`cp_nz01_2p27.txt` and `cp_nz01_2p272.txt`.

The runs of Table 13 use seed 42, trials 0-299, 5 rounds, and the
structures below; `r10_s2720` is the 10-round run at 2^27.20.  The chunk
files are named `<group>_t<first trial>.txt`, where the group `s<100 log2 N>`
encodes the structure size.

| Group | Structure s0 x s1 x s2 | N_s | Program |
|---|---|---|---|
| `s2605` | 8390 x 91 x 91 | 2^26.05 | `cp_nz012f` |
| `s2633` | 9153 x 96 x 96 | 2^26.33 | `cp_nz012f` |
| `s2671` | 10551 x 102 x 102 | 2^26.71 | `cp_nz012f` |
| `s2702` | 11664 x 108 x 108 | 2^27.02 | `cp_nz012` |
| `s2720` | 12515 x 111 x 111 | 2^27.20 | `cp_nz012` |
| `s2760` | 14366 x 119 x 119 | 2^27.60 | `cp_nz012f` |
| `s2780` | 15198 x 124 x 124 | 2^27.80 | `cp_nz012f` |
| `s2800` | 16384 x 128 x 128 | 2^28.00 | `cp_nz012f` |
| `s2830` | 18133 x 135 x 135 | 2^28.30 | `cp_nz012f` |
| `r10_s2720` | 12515 x 111 x 111, 10 rounds | 2^27.20 | `cp_nz012f` |

A chunk is one run such as

```sh
./cp_nz012f 9153 96 96 5 30 42 0 > s2633_t0.txt
```

(most chunks request 30 trials, those of `s2702` and `s2720` 60 trials, and
`r10_s2720` uses `rounds = 10`).  Interrupted chunks were resumed from the
next trial index, so a chunk can hold fewer trials than its header requests
and then has no `# P_find` footer; every group covers the trials 0-299
exactly once.

```sh
python3 Codes/tools/aggregate_combined.py > Results/5r_AES_distinguishers/combined/AGGREGATE.txt
```

checks the chunks (headers, trial coverage, raw = 2 x count, footers) and
writes the measured values of Table 13: `P_find` with standard errors for
a, b, c and their unions, the count distributions, and the number of trials
in which two pairings both succeed.  `Codes/model/cpoisson_model.py` reads
the same chunks through this script for section (4b) of its output.

`combined/verify_r3.txt` and `combined/verify_mask.txt` check the fast counts
of `cp_nz012` against a brute force written from the plaintext-byte
definitions (8th argument `verify = 1`, small grids only); every `V` line
(one per trial) ends with `OK`, and every run ends with `VERIFY_ALL_OK` in its
`# P_find` footer (one run in `verify_r3.txt`, three in `verify_mask.txt`).
They are not used in any table.

```sh
./cp_nz012 30 10 10 3 3 7 0 1 > verify_r3.txt      # 3-round AES, normal build
gcc -O3 -march=native -DINVD_MASK=0xFF000000u cp_nz012.c aes.c -o cp_nz012_m8  -lm
gcc -O3 -march=native -DINVD_MASK=0xFFFF0000u cp_nz012.c aes.c -o cp_nz012_m16 -lm
{ ./cp_nz012_m8  30 10 10 5 3 11 0 1
  ./cp_nz012_m8  30 10 10 4 2 12 0 1
  ./cp_nz012_m16 40 12 12 5 2 13 0 1; } > verify_mask.txt
```

The test builds compare only 8 or 16 bits of the inverse-diagonal value, so
that collisions occur on these small 4- and 5-round grids.

### Success-probability model (Sections 3.4, 4.1 to 4.3 and 5)

`Codes/model/` computes the compound-Poisson model of Section 4.2 and every
model value of the paper.  Run from `Codes/model`, in this order:

| Step | Command | Time | Writes |
|---|---|---|---|
| 1 | `python3 bundle_geometry.py` (defaults: 40 bundles, 4000 draws) | about 2 minutes | `data/bundle_geometry.json` (bundle geometry, h_geom) |
| 2 | `python3 eps_split.py` | more than 5 minutes | `data/eps_split.json` (epsilon(n_z) split by the weight m) |
| 3 | `python3 rho_pairing_c.py` (default n = 1000) | about 80 seconds | `data/rho_pairing_c.json` (T_1 of the pairing n_z = 1 on {5, 15}) |
| 4 | `python3 acp_budget_sim.py` (default 40000 trials) | about 5 seconds | `data/acp_budget_sim.json` (budget-limited ACP of `acp_all`) |
| 5 | `python3 cpoisson_model.py` | a few seconds | `data/cpoisson_results.txt` and `.json` |

The times were measured on a 20-core WSL machine.  `cpoisson_model.py` reads
the four JSON files above; `data/table5_w{4,8}.json` and `data/table3_w8.json`;
the measured files `data/5r_nz{0,1,2}.txt`, `data/raw_nz{0,1,2}_5r.txt`,
`data/cp_nz01_2p27.txt`, `data/cp_nz01_2p272.txt`, `data/cp_nz2_2p27.txt`,
`data/cp_nz2_2p272.txt` and `data/acp_all_2p2332.txt`; and the Section 4.3
chunks under `Results/5r_AES_distinguishers/combined/`.  (`eps_split.py` and
`acp_budget_sim.py` also read `data/table5_w*.json`.)  Every
`data/table*.json` file and every measured `data/*.txt` file is a verbatim
copy of the file of the same name under `Results/`.  The copies
`data/table3_w4.json`, `data/10r_nz{0,1,2}.txt` and `data/raw_nz{0,1,2}_10r.txt`
are kept for completeness, and no script reads them; `eps_split.py` uses the
m = 2 DDT mean 1.0166 of Table 3 (`table3_w8.json`) as a constant.  The files
`data/*.log` are the standard output of steps 1, 2 and 4.  A rerun can differ
from the committed JSON files in the last digit of a float.  The sections of
`cpoisson_results.txt` give:

| Section of `cpoisson_results.txt` | Paper |
|---|---|
| `(i)-(iv)` | inclusion probability q of the product grids of Table 12 (Section 4.2: 2^-9.98 at 2^27.02 and 2^-9.61 at 2^27.20), q of the A_0 x A_1 grids of Yan et al., and the per-member q and the bundle-hit probability h of one ACP base structure for n_z = 1 and 0 (used by the ACP model of section (2b)) |
| `constants` | T_1, epsilon and mu_b of Table 5 and Theorem 11; T_1 = 1019.94 of the pairing n_z = 1 on {5, 15} (Section 4.3); the Section 4.1 factors over 2^-62 and the mean rho over the base pairs |
| (2) Table 12 | model column of Table 12 (P_cpois), the Poisson column and the z-scores (Section 4.2: all six predictions within 1.33 standard errors; Poisson 0.64 and 0.73 for n_z = 1 at z = -3.85 and -3.56) |
| block "Yan et al.'s CP" | 61.8% on the grid of Yan et al., Poisson 73% (Section 4.2) |
| (2b) Table 14 | model column of Table 14 and the ACP of Yan et al. (Section 5.2); 2^16.6 adaptive swaps for n_z = 2 at 2^23.32 (Section 5.2); n_z = 0 budget of 2^5.4 base structures, 2^23.4 partner queries, 2^23.7 in total (Section 5.1) |
| (3) | standard deviations of the count over a full structure (Section 3.4.1); bundle totals of Section 3.4.2 |
| (4) | model columns and false positives of Table 13 (Section 4.3); block "data N (proportional grid, one structure) for a target P_find": the size at which the model reaches a target success probability on a proportional grid (Section 4.3: 2^27.82 for 90% with n_z = 0, a pairing n_z = 1 saturates at 0.879; the union reaches 63% at 2^26.33 and exactly 90% at 2^27.00, against 0.904 on the 2^27.02 run grid) |
| (4b) | Table 13 model against measured with z-scores, largest deviations, time 1.8 N_s, false positives of the 10-round run (Section 4.3) |
| (5) | saturation bound 0.88 of a pairing n_z = 1 in one structure, measured 0.89 (Section 4.3) |
| (6) | n_z = 2 in the ACP setting: 1.7% at 2^26.5, 89 base structures, 2^30.15 (Section 5.1); at the optimum M = 2^23.655: 2^16 q = 1.24, h = 0.71 and 1/88.8 per base structure (Section 5.1, which rounds M to 2^23.7) |

### Summaries

```sh
python3 Codes/tools/summarize_counts.py Results/5r_AES_N32/raw_nz1_5r.txt --nz 1 --w 8
python3 Codes/tools/summarize_counts.py Results/5r_small_AES/5r_nz0.txt   --nz 0 --w 4
```

prints the number of structures, the mean with its standard error, the
sample standard deviation, the detection rate `Pr[N_q >= 1]`, the number of
structures with at least one bundle, the sum of floor(N_q / bundle) (the
number of bundles when the 8B part of every count is below the bundle size,
as for the AES files with n_z = 1 and 2; not for n_z = 0), the number of
counts that are not multiples of 8 or do not have the form `bundle * A + 8B`,
and the expected count of a random permutation over a full diagonal
structure (the "random, expected" row of Table 11), to be compared with the
mean of a 10-round file.

## Tested platform

The committed results under `Results/` were measured on the following
machine.  Other Linux x86-64 hosts with AES-NI should build and run
identically.

| Component | Spec |
|---|---|
| CPU | 2× Intel Xeon Gold 6230R @ 2.10 GHz (52 cores / 104 threads, AES-NI + PCLMUL + AVX-512) |
| RAM | 1.9 TiB DDR4 (2 NUMA nodes) |
| OS | Ubuntu 22.04 LTS, kernel 6.8 |
| Compiler | gcc 11.4.0 |
| Storage | 1.8 TiB NVMe |

For N32 multi-thread runs on a 2-socket machine, `numactl --interleave=all`
spreads memory across both sockets and is recommended (install via
`sudo apt install numactl`).

## License

MIT (see `LICENSE`).
