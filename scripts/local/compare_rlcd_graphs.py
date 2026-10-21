# -*- coding: utf-8 -*-
"""
compare_rlcd_graphs.py — compare causal graphs output by RLCD (.dot) with a reference graph and check the physical prior

Usage (any Python 3, no extra packages):
    python3 scripts/local/compare_rlcd_graphs.py <reference.dot> <graph1.dot> [<graph2.dot> ...]

Notes:
  - Latent labels (names starting with L) are arbitrary (L1/L2 may be swapped); all relabelings are tried
    and the one closest to the reference graph is used
  - Prior check: age must have no observed parent (paper Sec. 3.2)
  - Undirected edges (edges in the "Undirected" subgraph) count as both directions and are listed separately
  - The reference may be given as git:<repo dir>:<relative path> to read the version at git HEAD
"""
import re
import sys
import itertools
import subprocess


def parse_dot(text):
    directed, undirected, nodes = set(), set(), []
    section = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith('subgraph Undirected'):
            section = 'u'
        elif s.startswith('subgraph Directed'):
            section = 'd'
        m = re.match(r'^(\w+)\s*->\s*(\w+)', s)
        if m:
            a, b = m.group(1), m.group(2)
            if section == 'u':
                undirected.add(frozenset((a, b)))
            else:
                directed.add((a, b))
        for n in re.findall(r'(\w+)\s*\[color', s):
            nodes.append(n)
    return {'nodes': nodes, 'directed': directed, 'undirected': undirected}


def latents(g):
    return sorted({n for n in g['nodes'] if re.match(r'^L\d+$', n)})


def relabel(g, mapping):
    f = lambda x: mapping.get(x, x)
    return {'nodes': [f(n) for n in g['nodes']],
            'directed': {(f(a), f(b)) for a, b in g['directed']},
            'undirected': {frozenset(f(x) for x in e) for e in g['undirected']}}


def best_match(ref, g):
    """Try all latent relabelings and return the relabeled graph closest to the reference"""
    lr, lg = latents(ref), latents(g)
    best, best_d = g, None
    if len(lr) == len(lg) and lg:
        for perm in itertools.permutations(lr):
            h = relabel(g, dict(zip(lg, perm)))
            d = len(ref['directed'] ^ h['directed']) + len(ref['undirected'] ^ h['undirected'])
            if best_d is None or d < best_d:
                best, best_d = h, d
    return best


def fmt(edges):
    return ', '.join(f'{a}→{b}' for a, b in sorted(edges)) or 'none'


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ref_path, paths = sys.argv[1], sys.argv[2:]
    if ref_path.startswith('git:'):          # git:<repo dir>:<relative path>  read the version at git HEAD
        _, repo, rel = ref_path.split(':', 2)
        ref_text = subprocess.check_output(['git', '-C', repo, 'show', f'HEAD:{rel}'], text=True)
    else:
        ref_text = open(ref_path).read()
    ref = parse_dot(ref_text)
    print(f"Reference graph: {ref_path}")
    print(f"  latents {latents(ref)}, {len(ref['directed'])} directed edges, {len(ref['undirected'])} undirected edges")
    for L in latents(ref):
        ch = sorted(b for a, b in ref['directed'] if a == L)
        pa = sorted(a for a, b in ref['directed'] if b == L)
        print(f"    {L}: parents {pa or 'none'}, children {ch}")

    for p in paths:
        g = best_match(ref, parse_dot(open(p).read()))
        print('\n' + '=' * 90)
        print(f"Graph: {p}")
        print(f"  latents {latents(g)} (labels aligned to the reference), {len(g['directed'])} directed edges, {len(g['undirected'])} undirected edges")
        for L in latents(g):
            ch = sorted(b for a, b in g['directed'] if a == L)
            pa = sorted(a for a, b in g['directed'] if b == L)
            print(f"    {L}: parents {pa or 'none'}, children {ch}")
        missing = ref['directed'] - g['directed']
        extra = g['directed'] - ref['directed']
        # list reversed edges separately
        reversed_ = {(a, b) for a, b in extra if (b, a) in missing}
        extra_only = extra - reversed_
        missing_only = missing - {(b, a) for a, b in reversed_}
        print(f"  directed edges shared with the reference: {len(ref['directed'] & g['directed'])} / {len(ref['directed'])}")
        print(f"  missing edges: {fmt(missing_only)}")
        print(f"  extra edges: {fmt(extra_only)}")
        print(f"  reversed edges: {fmt(reversed_)}")
        if g['undirected'] or ref['undirected']:
            print(f"  undirected edges: this graph {[tuple(sorted(e)) for e in g['undirected']]}, reference {[tuple(sorted(e)) for e in ref['undirected']]}")
        age_par = sorted(a for a, b in g['directed'] if b == 'X_age' and not re.match(r'^L\d+$', a))
        print(f"  prior check (age has no observed parent): {'pass' if not age_par else 'FAIL, observed parents of age: ' + str(age_par)}")
        same = not missing and not extra and g['undirected'] == ref['undirected']
        print(f"  verdict: {'identical to the reference' if same else 'different from the reference'}")


if __name__ == '__main__':
    main()
