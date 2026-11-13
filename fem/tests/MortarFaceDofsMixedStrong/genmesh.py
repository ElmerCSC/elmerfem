#!/usr/bin/env python3
"""Generate Elmer meshes for the mortar regression tests.

Usage: genmesh.py mixed|circle <meshdir>

mixed : box, lower body 3x3x2 hexes, upper body 4x4x2 hexes with prisms in the
        right half (x > 0.5), with disconnected nodes at z=0.5.
        BC1 = slave (upper bottom), BC2 = master (lower top), BC3 = outer.
        Environment options:
          SLAVE_FIRST=404|303  which slave face type is listed first (default 404)
          LOWER_TRI=1          also put prisms in the right half of the lower body
          TRIX=<x>             prisms where the cell center x > TRIX (default 0.5)
        MortarFaceDofsMixed:       defaults
        MortarFaceDofsMixedStrong: LOWER_TRI=1 SLAVE_FIRST=303
circle: full-circle quadratic (820) rotor r in [0.5,1] (24 sectors), stator
        r in [1,1.5] (32 sectors, offset 4 deg), z in [0,0.5], disconnected at
        r=1. BC1 = slave (rotor outer), BC2 = master (stator inner), BC3 = outer.
        LINEAR=1 writes 808 elements instead (corner nodes only).
"""
import math, os, sys
from collections import defaultdict


class Mesh:
    def __init__(self):
        self.nodes = []          # (x,y,z)
        self.elems = []          # (body, type, [nodes 1-based])

    def node(self, x, y, z):
        self.nodes.append((x, y, z))
        return len(self.nodes)


HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
PRI_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
# 820: corner face -> midside nodes (0-based local) for the 408 face
MID820 = {(0, 1): 8, (1, 2): 9, (2, 3): 10, (3, 0): 11, (0, 4): 12, (1, 5): 13, (2, 6): 14,
          (3, 7): 15, (4, 5): 16, (5, 6): 17, (6, 7): 18, (7, 4): 19}


def mid820(a, b):
    return MID820.get((a, b), MID820.get((b, a)))


def write(mesh, bnd, d):
    os.makedirs(d, exist_ok=True)
    with open(f'{d}/mesh.nodes', 'w') as f:
        for i, (x, y, z) in enumerate(mesh.nodes, 1):
            f.write(f'{i} -1 {x:.16g} {y:.16g} {z:.16g}\n')
    with open(f'{d}/mesh.elements', 'w') as f:
        for i, (b, t, n) in enumerate(mesh.elems, 1):
            f.write(f'{i} {b} {t} ' + ' '.join(map(str, n)) + '\n')
    with open(f'{d}/mesh.boundary', 'w') as f:
        for i, (bc, p, t, n) in enumerate(bnd, 1):
            f.write(f'{i} {bc} {p} 0 {t} ' + ' '.join(map(str, n)) + '\n')
    types = defaultdict(int)
    for _, t, _ in mesh.elems:
        types[t] += 1
    for _, _, t, _ in bnd:
        types[t] += 1
    with open(f'{d}/mesh.header', 'w') as f:
        f.write(f'{len(mesh.nodes)} {len(mesh.elems)} {len(bnd)}\n{len(types)}\n')
        for t in sorted(types):
            f.write(f'{t} {types[t]}\n')


def boundary_faces(mesh, classify):
    """Return list of (bc, parent, type, nodes) for faces where classify gives a bc."""
    out = []
    for ei, (b, t, n) in enumerate(mesh.elems, 1):
        if t == 808:
            faces = [[n[i] for i in f] for f in HEX_FACES]
        elif t == 706:
            faces = [[n[i] for i in f] for f in PRI_FACES]
        elif t == 820:
            faces = []
            for f in HEX_FACES:
                c = [n[i] for i in f]
                m = [n[mid820(f[k], f[(k + 1) % 4])] for k in range(4)]
                faces.append(c + m)
        for fn in faces:
            ncorner = 3 if len(fn) == 3 else 4
            bc = classify(b, [mesh.nodes[i - 1] for i in fn[:ncorner]])
            if bc:
                ft = {3: 303, 4: 404, 8: 408}[len(fn)]
                out.append((bc, ei, ft, fn))
    return out


def gen_mixed(d):
    m = Mesh()
    tol = 1e-12

    def body(bodyid, nx, z0, z1, nz, tri_right):
        ids = {}
        for k in range(nz + 1):
            z = z0 + (z1 - z0) * k / nz
            for j in range(nx + 1):
                for i in range(nx + 1):
                    ids[i, j, k] = m.node(i / nx, j / nx, z)
        for k in range(nz):
            for j in range(nx):
                for i in range(nx):
                    q = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                    lo = [ids[a, b, k] for a, b in q]
                    hi = [ids[a, b, k + 1] for a, b in q]
                    if tri_right and (i + 0.5) / nx > float(os.environ.get("TRIX", "0.5")):
                        for t in ((0, 1, 2), (0, 2, 3)):
                            m.elems.append((bodyid, 706, [lo[t[0]], lo[t[1]], lo[t[2]],
                                                          hi[t[0]], hi[t[1]], hi[t[2]]]))
                    else:
                        m.elems.append((bodyid, 808, lo + hi))

    # LOWER_TRI=1: also the master (lower) body gets prisms in its right half
    body(1, 3, 0.0, 0.5, 2, bool(os.environ.get('LOWER_TRI')))
    body(2, 4, 0.5, 1.0, 2, True)

    def classify(b, xs):
        if all(abs(p[2] - 0.5) < tol for p in xs):
            return 1 if b == 2 else 2
        for c in range(3):
            for v in (0.0, 1.0):
                if all(abs(p[c] - v) < tol for p in xs):
                    return 3
        return 0

    bnd = boundary_faces(m, classify)
    # Slave faces: put quads first, then triangles. With second-order Piola
    # quads have 4 and triangles 2 face DOFs, so the old bdofs*(ind-1) row
    # numbering makes triangle rows collide with quad rows.
    first = int(os.environ.get('SLAVE_FIRST', '404'))  # 303: triangles first
    bnd.sort(key=lambda f: (f[0], 0 if (f[0] == 1 and f[2] == first) else 1))
    write(m, bnd, d)


def gen_circle(d):
    m = Mesh()
    tol = 1e-9

    def ring(bodyid, r0, r1, nr, nth, th0, z0, z1):
        # quadratic nodes on a (2nr+1) x (2nth) x 3 lattice, periodic in theta
        ids = {}
        for k in range(3):
            z = z0 + (z1 - z0) * k / 2
            for j in range(2 * nth):
                th = th0 + math.pi * j / nth
                for i in range(2 * nr + 1):
                    r = r0 + (r1 - r0) * i / (2 * nr)
                    ids[i, j, k] = m.node(r * math.cos(th), r * math.sin(th), z)

        def L(i, j, k):
            return ids[i, j % (2 * nth), k]

        for j in range(nth):
            for i in range(nr):
                I, J = 2 * i, 2 * j
                c = [L(I, J, 0), L(I + 2, J, 0), L(I + 2, J + 2, 0), L(I, J + 2, 0),
                     L(I, J, 2), L(I + 2, J, 2), L(I + 2, J + 2, 2), L(I, J + 2, 2)]
                e = [L(I + 1, J, 0), L(I + 2, J + 1, 0), L(I + 1, J + 2, 0), L(I, J + 1, 0),
                     L(I, J, 1), L(I + 2, J, 1), L(I + 2, J + 2, 1), L(I, J + 2, 1),
                     L(I + 1, J, 2), L(I + 2, J + 1, 2), L(I + 1, J + 2, 2), L(I, J + 1, 2)]
                if os.environ.get('LINEAR'):
                    m.elems.append((bodyid, 808, c))
                else:
                    m.elems.append((bodyid, 820, c + e))

    ring(1, 0.5, 1.0, 2, 24, 0.0, 0.0, 0.5)
    ring(2, 1.0, 1.5, 2, 32, math.radians(4.0), 0.0, 0.5)

    def classify(b, xs):
        rs = [math.hypot(p[0], p[1]) for p in xs]
        if all(abs(r - 1.0) < tol for r in rs):
            return 1 if b == 1 else 2
        if all(abs(r - 0.5) < tol for r in rs) or all(abs(r - 1.5) < tol for r in rs):
            return 3
        if all(abs(p[2]) < tol for p in xs) or all(abs(p[2] - 0.5) < tol for p in xs):
            return 3
        return 0

    write(m, boundary_faces(m, classify), d)


if __name__ == '__main__':
    which, d = sys.argv[1], sys.argv[2]
    {'mixed': gen_mixed, 'circle': gen_circle}[which](d)
