#!/usr/bin/env python3
"""Map an MPI process count ``np`` to an OpenFOAM ``simple`` decomposition.

``factor3(np, dims)`` returns ``(nx, ny, nz)`` whose product is exactly ``np``,
balanced as evenly as possible.  For 1D/2D cases the out-of-plane direction(s)
are pinned to 1 so the mesh is only split in-plane (a z-split of a one-cell-thick
2D case would put empty processors on the empty front/back patches).

Examples (dims=3): 1->(1,1,1) 2->(2,1,1) 4->(2,2,1) 8->(2,2,2) 16->(4,2,2)
Examples (dims=2): 4->(2,2,1) 8->(4,2,1) 16->(4,4,1)
"""
import sys


def _prime_factors(n):
    factors, d = [], 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def factor3(np_, dims=3):
    np_ = int(np_)
    if np_ < 1:
        raise ValueError("np must be >= 1")
    nbuckets = 3 if dims >= 3 else (2 if dims == 2 else 1)
    buckets = [1, 1, 1]
    for f in sorted(_prime_factors(np_), reverse=True):
        i = min(range(nbuckets), key=lambda b: buckets[b])
        buckets[i] *= f
    return tuple(buckets)


def _main(argv):
    targets = [int(x) for x in argv] if argv else [1, 2, 4, 8, 16]
    for dims in (3, 2, 1):
        print(f"dims={dims}")
        for n in targets:
            t = factor3(n, dims)
            assert t[0] * t[1] * t[2] == n, (n, dims, t)
            print(f"  np={n:3d} -> {t}")


if __name__ == "__main__":
    _main(sys.argv[1:])
