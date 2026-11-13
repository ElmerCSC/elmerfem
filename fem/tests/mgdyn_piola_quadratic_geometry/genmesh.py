#!/usr/bin/env python3
"""Mesh for mgdyn_piola_quadratic_geometry: three disconnected unit boxes made of
quadratic 820 bricks (x in [0,1]), 715 wedges (x in [2,3]), 718 wedges (x in [4,5]) and
613 pyramids (x in [6,7]) and 827 bricks (x in [8,9]), each 2x2x2 cells (wedges: each cell split into two prisms,
pyramids: six per cell with the apex at the cell center). Boundary 1 = all outer faces.
Usage: genmesh.py <meshdir>
"""
import os, sys
from collections import defaultdict

N = 2           # cells per direction
L = 2 * N       # quadratic lattice points per direction minus one

HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
HEX_MID = {(0, 1): 8, (1, 2): 9, (2, 3): 10, (3, 0): 11, (0, 4): 12, (1, 5): 13, (2, 6): 14,
           (3, 7): 15, (4, 5): 16, (5, 6): 17, (6, 7): 18, (7, 4): 19}
HEX_FC = [24, 25, 20, 21, 22, 23]   # 827: face-center node (0-based) per HEX_FACES entry
PRI_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
PRI_MID = {(0, 1): 6, (1, 2): 7, (2, 0): 8, (3, 4): 9, (4, 5): 10, (5, 3): 11, (0, 3): 12, (1, 4): 13, (2, 5): 14}
PRI_QC = {(0, 1, 4, 3): 15, (1, 2, 5, 4): 16, (2, 0, 3, 5): 17}


def mid(tab, a, b):
    return tab.get((a, b), tab.get((b, a)))


nodes, elems, bnd = [], [], []


def box(x0, typ):
    ids, lat = {}, {}
    for k in range(L + 1):
        for j in range(L + 1):
            for i in range(L + 1):
                nodes.append((x0 + i / L, j / L, k / L))
                ids[i, j, k] = len(nodes)
                lat[len(nodes)] = (i, j, k)
    first = len(elems)
    for k in range(N):
        for j in range(N):
            for i in range(N):
                I, J, K = 2 * i, 2 * j, 2 * k
                P = lambda a, b, c: ids[I + a, J + b, K + c]
                if typ in (820, 827):
                    q = [(0, 0), (2, 0), (2, 2), (0, 2)]
                    c = [P(a, b, 0) for a, b in q] + [P(a, b, 2) for a, b in q]
                    e = [P(1, 0, 0), P(2, 1, 0), P(1, 2, 0), P(0, 1, 0),
                         P(0, 0, 1), P(2, 0, 1), P(2, 2, 1), P(0, 2, 1),
                         P(1, 0, 2), P(2, 1, 2), P(1, 2, 2), P(0, 1, 2)]
                    fc = [P(1, 0, 1), P(2, 1, 1), P(1, 2, 1), P(0, 1, 1), P(1, 1, 0), P(1, 1, 2), P(1, 1, 1)]
                    elems.append((typ, c + e + (fc if typ == 827 else [])))
                else:
                    for t in (((0, 0), (2, 0), (2, 2)), ((0, 0), (2, 2), (0, 2))):
                        m = lambda p, r: ((p[0] + r[0]) // 2, (p[1] + r[1]) // 2)
                        c = [P(*p, 0) for p in t] + [P(*p, 2) for p in t]
                        e = [P(*m(t[0], t[1]), 0), P(*m(t[1], t[2]), 0), P(*m(t[2], t[0]), 0),
                             P(*m(t[0], t[1]), 2), P(*m(t[1], t[2]), 2), P(*m(t[2], t[0]), 2),
                             P(*t[0], 1), P(*t[1], 1), P(*t[2], 1)]
                        f = [P(*m(t[0], t[1]), 1), P(*m(t[1], t[2]), 1), P(*m(t[2], t[0]), 1)]
                        elems.append((typ, c + e + (f if typ == 718 else [])))
    for ei in range(first, len(elems)):
        typ_e, nn = elems[ei]
        faces = HEX_FACES if typ_e in (820, 827) else PRI_FACES
        tab = HEX_MID if typ_e in (820, 827) else PRI_MID
        for fi, f in enumerate(faces):
            cs = [lat[nn[q]] for q in f]
            if not any(all(p[ax] == v for p in cs) for ax in range(3) for v in (0, L)):
                continue
            fn = [nn[q] for q in f] + [nn[mid(tab, f[q], f[(q + 1) % len(f)])] for q in range(len(f))]
            if typ_e == 718 and len(f) == 4:
                fn.append(nn[PRI_QC[f]])
            if typ_e == 827:
                fn.append(nn[HEX_FC[fi]])
            bnd.append((ei + 1, {6: 306, 8: 408, 9: 409}[len(fn)], fn))


def pyramid_box(x0):
    key = {}

    def node(x, y, z):
        k = (round(x * 64), round(y * 64), round(z * 64))
        if k not in key:
            nodes.append((x0 + x, y, z))
            key[k] = len(nodes)
        return key[k]

    def midpt(a, b):
        return tuple((a[q] + b[q]) / 2 for q in range(3))

    h = 1.0 / N
    for k in range(N):
        for j in range(N):
            for i in range(N):
                o = (i * h, j * h, k * h)
                cr = [(o[0] + a * h, o[1] + b * h, o[2] + c * h) for a, b, c in
                      ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))]
                ctr = (o[0] + h / 2, o[1] + h / 2, o[2] + h / 2)
                for f in HEX_FACES:
                    b = [cr[q] for q in f]
                    # orient the base so that its right-hand normal points to the apex
                    e1 = [b[1][q] - b[0][q] for q in range(3)]
                    e2 = [b[3][q] - b[0][q] for q in range(3)]
                    nrm = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
                    if sum(nrm[q] * (ctr[q] - b[0][q]) for q in range(3)) < 0:
                        b = [b[0], b[3], b[2], b[1]]
                    pts = b + [ctr]
                    pts += [midpt(b[0], b[1]), midpt(b[1], b[2]), midpt(b[2], b[3]), midpt(b[3], b[0])]
                    pts += [midpt(b[q], ctr) for q in range(4)]
                    nn = [node(*p) for p in pts]
                    elems.append((613, nn))
                    if any(all(abs(p[ax] - v) < 1e-12 for p in b) for ax in range(3) for v in (0.0, 1.0)):
                        bnd.append((len(elems), 408, nn[0:4] + nn[5:9]))


box(0.0, 820)
box(2.0, 715)
box(4.0, 718)
pyramid_box(6.0)
box(8.0, 827)

d = sys.argv[1]
os.makedirs(d, exist_ok=True)
with open(f'{d}/mesh.nodes', 'w') as f:
    for i, (x, y, z) in enumerate(nodes, 1):
        f.write(f'{i} -1 {x:.16g} {y:.16g} {z:.16g}\n')
with open(f'{d}/mesh.elements', 'w') as f:
    for i, (t, nn) in enumerate(elems, 1):
        f.write(f'{i} 1 {t} ' + ' '.join(map(str, nn)) + '\n')
with open(f'{d}/mesh.boundary', 'w') as f:
    for i, (p, t, nn) in enumerate(bnd, 1):
        f.write(f'{i} 1 {p} 0 {t} ' + ' '.join(map(str, nn)) + '\n')
types = defaultdict(int)
for t, _ in elems:
    types[t] += 1
for _, t, _ in bnd:
    types[t] += 1
with open(f'{d}/mesh.header', 'w') as f:
    f.write(f'{len(nodes)} {len(elems)} {len(bnd)}\n{len(types)}\n')
    for t in sorted(types):
        f.write(f'{t} {types[t]}\n')
