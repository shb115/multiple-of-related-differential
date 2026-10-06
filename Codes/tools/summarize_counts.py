#!/usr/bin/env python3
"""
summarize_counts.py -- summary of valid-quartet counts N_q, one value per structure.

Input: either a raw per-trial file (one integer count per line, as written by
Codes/5r_AES_N32/count_quartets_nz*) or a frequency file ("Count <N> : <k> times" lines,
as printed by Codes/5r_small_AES/count_quartets_nz*).

Output: number of structures, mean with its standard error, sample standard deviation
(ddof = 1), detection rate Pr[N_q >= 1], the number of structures with N_q >= bundle, the
sum of floor(N_q / bundle) (the number of bundles when the 8B part of every count is below
the bundle size, as for the AES files with n_z = 1 and 2, whose counts are at most 24 and 16
above a multiple of the bundle; not for n_z = 0, where the bundle is 4), the number of
counts that are not multiples of 8, the number of counts that do not have the form
bundle * A + 8B, and the expected count E[N_q] of a random permutation over a full diagonal
structure (Section 3.4; the "random, expected" row of Table 11), to be compared with the mean
of a 10-round file.

Usage:
    python3 summarize_counts.py <file> --nz <0|1|2> --w <4|8>
"""
import argparse, math, re
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument('file')
ap.add_argument('--nz', type=int, required=True, choices=(0, 1, 2))
ap.add_argument('--w', type=int, required=True, choices=(4, 8))
args = ap.parse_args()

text = open(args.file).read()
freq = Counter()
pairs = re.findall(r'Count\s+(\d+)\s*:\s*(\d+)\s+times', text)
if pairs:
    for a, b in pairs:
        freq[int(a)] += int(b)
else:
    for line in text.split():
        if line.strip().lstrip('-').isdigit():
            freq[int(line)] += 1

n = sum(freq.values())
mean = sum(k * v for k, v in freq.items()) / n
var = sum(v * (k - mean) ** 2 for k, v in freq.items()) / (n - 1) if n > 1 else 0.0
bundle = 2 ** (2 - args.nz) * (2 ** args.w) ** args.nz
q = 2 ** args.w
e_rand = (q - 1) ** 4 / q ** 4 if args.nz in (1, 2) else (q - 1) ** 3 * (q + 1) / q ** 4

def has_form(k):
    """k = bundle * A + 8 B with A, B >= 0."""
    return any((k - bundle * a) % 8 == 0 for a in range(k // bundle + 1))

print(f'file            : {args.file}')
print(f'n_z, w, bundle  : {args.nz}, {args.w}, {bundle}')
print(f'structures      : {n}')
print(f'mean N_q        : {mean:.3f}  (standard error {math.sqrt(var / n):.3f})')
print(f'std deviation   : {math.sqrt(var):.2f}')
print(f'Pr[N_q >= 1]    : {1 - freq.get(0, 0) / n:.3f}')
print(f'N_q >= bundle   : {sum(v for k, v in freq.items() if k >= bundle)} structures')
print(f'bundles         : {sum(v * (k // bundle) for k, v in freq.items())} (sum of floor(N_q / bundle))')
print(f'not mult. of 8  : {sum(v for k, v in freq.items() if k % 8)} structures')
print(f'not bundle*A+8B : {sum(v for k, v in freq.items() if not has_form(k))} structures')
print('frequency       : ' + ', '.join(f'{k}:{v}' for k, v in sorted(freq.items())))
print(f'random-perm E[N_q]: {e_rand:.4f}  (Section 3.4; compare with the mean of a 10-round file)')
