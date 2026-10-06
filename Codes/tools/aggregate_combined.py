#!/usr/bin/env python3
"""
aggregate_combined.py -- measured success probabilities of the three pairings on one product
structure (Section 4.3, Table 13) from the per-chunk outputs of cp_nz012 / cp_nz012f.

Input: the chunk files <group>_t<first trial>.txt in Results/5r_AES_distinguishers/combined
(verify_*.txt are skipped).  Only the '# cp_nz012[f] s0=... t_start=...' header, the 'T' lines
and the optional '# P_find' footer of a chunk are read.  Groups: s2605 ... s2830 (5 rounds,
N = 2^26.05 ... 2^28.30) and r10_s2720 (10 rounds, N = 2^27.20).

Pairings: a = n_z=1 on {0,10} (Yan et al.), b = n_z=0, c = n_z=1 on {5,15}.
Table 13 lists them in the order a, c, b.
T line: T t ca cb cc rawA rawB rawC baseA baseB P2both candC  (a hit is a count > 0).

Checks (an AssertionError stops the script):
  - the headers of a group agree apart from trials= and t_start=, and t_start equals the
    trial index in the file name;
  - every trial index t = 0 .. 299 occurs exactly once in a group (a repeated T line of an
    overlapping chunk would have to be identical);
  - rawA = 2 ca, rawB = 2 cb, rawC = 2 cc (each quartet is found twice);
  - a '# P_find' footer, where present, agrees with the T lines of its chunk.
A chunk that was interrupted holds fewer T lines than its header requests and has no footer;
the run was resumed with t_start = next trial index (the per-trial seed is
seed ^ 0x9E3779B97F4A7C15*(t+1), so the split into chunks does not change any trial).

Output (stdout): per group the header of its first chunk, the chunk files, the number of
candidates per trial with the random-permutation expectation 4Q/2^64, P_find with binomial
standard error for a, b, c and their unions, count mean / sample sd (n-1) / se / distribution,
and the number of trials in which two pairings both succeed.

Usage:
    python3 Codes/tools/aggregate_combined.py [combined_dir] > Results/5r_AES_distinguishers/combined/AGGREGATE.txt
"""
import glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DIR = os.path.join(HERE, '..', '..', 'Results', '5r_AES_distinguishers', 'combined')
NTRIALS = 300

HDR = re.compile(r'^# (cp_nz012f?) s0=(\d+) s1=(\d+) s2=(\d+) N=(\d+)\(2\^([\d.]+)\) rounds=(\d+) '
                 r'trials=(\d+) seed=(\d+) t_start=(\d+)\s*$')
FOOT = re.compile(r'^# P_find a=(\d+) b=(\d+) c=(\d+) a\|b=(\d+) a\|c=(\d+) b\|c=(\d+) a\|b\|c=(\d+) of (\d+)')
FNAME = re.compile(r'^(.+)_t(\d+)\.txt$')
UNIONS = ('a', 'b', 'c', 'a|b', 'a|c', 'b|c', 'a|b|c')


def hits_of(rows):
    """rows: iterable of (ca, cb, cc) -> dict of hit counts for a, b, c and the unions."""
    h = dict.fromkeys(UNIONS, 0)
    for ca, cb, cc in rows:
        A, B, C = ca > 0, cb > 0, cc > 0
        h['a'] += A; h['b'] += B; h['c'] += C
        h['a|b'] += A or B; h['a|c'] += A or C; h['b|c'] += B or C; h['a|b|c'] += A or B or C
    return h


def read_chunk(path, t_first):
    hdr = None; trials = {}; foot = None
    for line in open(path):
        line = line.rstrip('\r\n')
        if line.startswith('# cp_nz012'):
            m = HDR.match(line)
            assert m and hdr is None, 'bad or repeated header in %s' % path
            hdr = (line, m.groups())
        elif line.startswith('# P_find'):
            m = FOOT.match(line)
            assert m, 'bad footer in %s' % path
            foot = tuple(int(x) for x in m.groups())
        elif line.startswith('T '):
            v = [int(x) for x in line.split()[1:]]
            assert len(v) == 11, 'bad T line in %s: %s' % (path, line)
            t = v[0]
            assert t not in trials, 'trial %d twice in %s' % (t, path)
            assert v[4] == 2 * v[1] and v[5] == 2 * v[2] and v[6] == 2 * v[3], 'raw != 2 x count in %s: %s' % (path, line)
            trials[t] = tuple(v[1:])
    assert hdr is not None, 'no header in %s' % path
    assert int(hdr[1][9]) == t_first, 't_start in %s differs from the file name' % path
    if foot is not None:
        h = hits_of((x[0], x[1], x[2]) for x in trials.values())
        assert foot[:7] == tuple(h[u] for u in UNIONS) and foot[7] == len(trials), 'footer mismatch in %s' % path
    return hdr, trials, foot


def aggregate(combined_dir=None):
    """Return {group: record}; record holds the header, grid, rounds, trials and the summary."""
    d = combined_dir or DEFAULT_DIR
    files = {}
    for path in sorted(glob.glob(os.path.join(d, '*_t*.txt'))):
        m = FNAME.match(os.path.basename(path))
        if m:
            files.setdefault(m.group(1), []).append((int(m.group(2)), path))
    out = {}
    for grp, lst in files.items():
        lst.sort()
        key = None; first_hdr = None; trials = {}; interrupted = []
        for t_first, path in lst:
            hdr, tr, foot = read_chunk(path, t_first)
            g = hdr[1]
            k = (g[0],) + g[1:7] + (g[8],)          # program, s0, s1, s2, N, log2N, rounds, seed
            if key is None:
                key = k; first_hdr = hdr[0]
            assert k == key, 'header of %s differs from the other chunks of %s' % (path, grp)
            if len(tr) < int(g[7]):
                assert foot is None, 'short chunk with footer: %s' % path
                interrupted.append('%s (%d of %d)' % (os.path.basename(path), len(tr), int(g[7])))
            for t, v in tr.items():
                assert t not in trials or trials[t] == v, 'trial %d of %s differs between chunks' % (t, grp)
                trials[t] = v
        assert sorted(trials) == list(range(NTRIALS)), '%s: trial indices are not exactly 0..%d' % (grp, NTRIALS - 1)
        s0, s1, s2, N = (int(x) for x in key[1:5])
        rows = [trials[t] for t in range(NTRIALS)]
        out[grp] = dict(group=grp, header=first_hdr, program=key[0], grid=(s0, s1, s2), N=N, log2N=math.log2(N),
                        rounds=int(key[6]), seed=int(key[7]), n=NTRIALS, chunks=[os.path.basename(p) for _, p in lst],
                        interrupted=interrupted, counts={p: [r[i] for r in rows] for i, p in enumerate('abc')},
                        hits=hits_of((r[0], r[1], r[2]) for r in rows))
    return out


def ncand(s0, s1, s2):
    """Candidates per trial: Q_b = C(s0,2) W and Q_a = Q_c = Q_b 255^2/65535, W = s1(s1-1)s2(s2-1)/2."""
    W = s1 * (s1 - 1) * s2 * (s2 - 1) / 2
    qb = s0 * (s0 - 1) / 2 * W
    return qb * 255 * 255 / 65535, qb


def binse(p, n):
    return math.sqrt(p * (1 - p) / n)


def report(agg):
    say = print
    say('# Section 4.3 (Table 13): three pairings on one product structure, measured in 300 trials per size.')
    say('# Generated by Codes/tools/aggregate_combined.py from the chunk files of this directory.')
    say('# a = pairing n_z=1 on {0,10} (Yan et al.), b = pairing n_z=0, c = pairing n_z=1 on {5,15};')
    say('# Table 13 lists them in the order a, c, b.  T line: T t ca cb cc rawA rawB rawC baseA baseB P2both candC.')
    for grp, r in sorted(agg.items(), key=lambda kv: (kv[1]['rounds'], kv[1]['N'])):
        n = r['n']
        say('== %s: %s' % (grp, r['header']))
        say('   trials=%d (t=0..%d)' % (n, n - 1))
        say('   chunks: %d files; interrupted and resumed: %s' % (len(r['chunks']), ', '.join(r['interrupted']) or 'none'))
        qa, qb = ncand(*r['grid'])
        ea = 4 * qa * 2.0 ** -64; eb = 4 * qb * 2.0 ** -64
        say('   candidates/trial: Q_a=Q_c=%.4e (2^%.3f), Q_b=%.4e (2^%.3f); random-perm E[valid]=4Q/2^64: a,c 2^%.3f, b 2^%.3f' % (
            qa, math.log2(qa), qb, math.log2(qb), math.log2(ea), math.log2(eb)))
        for u in UNIONS:
            k = r['hits'][u]; p = k / n
            say('   P_find %-7s= %d/%d = %.4f +- %.4f' % (u, k, n, p, binse(p, n)))
        for p in 'abc':
            x = r['counts'][p]
            mean = sum(x) / n
            sd = math.sqrt(sum((v - mean) ** 2 for v in x) / (n - 1))
            dist = {}
            for v in x:
                dist[v] = dist.get(v, 0) + 1
            say('   counts %s: mean=%.4f sd=%.3f se=%.4f  dist=%s' % (p, mean, sd, sd / math.sqrt(n), dict(sorted(dist.items()))))
        ca, cb, cc = (r['counts'][p] for p in 'abc')
        both = lambda x, y: sum(1 for i in range(n) if x[i] > 0 and y[i] > 0)
        say('   both a&c=%d, a&b=%d, b&c=%d' % (both(ca, cc), both(ca, cb), both(cb, cc)))


if __name__ == '__main__':
    report(aggregate(sys.argv[1] if len(sys.argv) > 1 else None))
