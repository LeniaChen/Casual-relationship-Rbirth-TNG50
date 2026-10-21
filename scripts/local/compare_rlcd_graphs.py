# -*- coding: utf-8 -*-
"""
compare_rlcd_graphs.py —— 比较 RLCD 输出的因果图 (.dot) 与参考图, 并检查物理先验

用法 (任何 Python 3 都可以, 不需要额外的包):
    python3 scripts/local/compare_rlcd_graphs.py <参考图.dot> <图1.dot> [<图2.dot> ...]

说明:
  - 隐变量 (名字以 L 开头) 的编号是任意的 (L1/L2 可能互换), 比较前会尝试所有编号对应方式,
    取与参考图差别最小的那一种
  - 先验检查: age 在观测变量里不能有父节点 (论文 3.2 节)
  - 无向边 (Undirected 子图里的边) 按"两个方向都算"处理, 单独列出
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
    """尝试所有隐变量编号对应, 返回差别最小的重新编号后的图"""
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
    return ', '.join(f'{a}→{b}' for a, b in sorted(edges)) or '无'


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ref_path, paths = sys.argv[1], sys.argv[2:]
    if ref_path.startswith('git:'):          # git:<仓库目录>:<相对路径>  读 git HEAD 中的版本
        _, repo, rel = ref_path.split(':', 2)
        ref_text = subprocess.check_output(['git', '-C', repo, 'show', f'HEAD:{rel}'], text=True)
    else:
        ref_text = open(ref_path).read()
    ref = parse_dot(ref_text)
    print(f"参考图: {ref_path}")
    print(f"  隐变量 {latents(ref)}, 有向边 {len(ref['directed'])} 条, 无向边 {len(ref['undirected'])} 条")
    for L in latents(ref):
        ch = sorted(b for a, b in ref['directed'] if a == L)
        pa = sorted(a for a, b in ref['directed'] if b == L)
        print(f"    {L}: 父节点 {pa or '无'}, 子节点 {ch}")

    for p in paths:
        g = best_match(ref, parse_dot(open(p).read()))
        print('\n' + '=' * 90)
        print(f"图: {p}")
        print(f"  隐变量 {latents(g)} (已按参考图对齐编号), 有向边 {len(g['directed'])} 条, 无向边 {len(g['undirected'])} 条")
        for L in latents(g):
            ch = sorted(b for a, b in g['directed'] if a == L)
            pa = sorted(a for a, b in g['directed'] if b == L)
            print(f"    {L}: 父节点 {pa or '无'}, 子节点 {ch}")
        missing = ref['directed'] - g['directed']
        extra = g['directed'] - ref['directed']
        # 方向相反的边单独标出
        reversed_ = {(a, b) for a, b in extra if (b, a) in missing}
        extra_only = extra - reversed_
        missing_only = missing - {(b, a) for a, b in reversed_}
        print(f"  与参考图相同的有向边: {len(ref['directed'] & g['directed'])} / {len(ref['directed'])}")
        print(f"  缺少的边: {fmt(missing_only)}")
        print(f"  多出的边: {fmt(extra_only)}")
        print(f"  方向相反的边: {fmt(reversed_)}")
        if g['undirected'] or ref['undirected']:
            print(f"  无向边: 本图 {[tuple(sorted(e)) for e in g['undirected']]}, 参考图 {[tuple(sorted(e)) for e in ref['undirected']]}")
        age_par = sorted(a for a, b in g['directed'] if b == 'X_age' and not re.match(r'^L\d+$', a))
        print(f"  先验检查 (age 没有观测父节点): {'通过' if not age_par else '不通过, age 的观测父节点: ' + str(age_par)}")
        same = not missing and not extra and g['undirected'] == ref['undirected']
        print(f"  结论: {'与参考图完全相同' if same else '与参考图不同'}")


if __name__ == '__main__':
    main()
