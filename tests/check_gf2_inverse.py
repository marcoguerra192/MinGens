"""
Check the exact inverse over Z2 used by CySPT, and show where the previous
floating point formula fails.

Usage, from the MinGens folder, after make:

    python tests/check_gf2_inverse.py
"""
import sys

import numpy as np

sys.path.insert(0, 'src')
import CySPT  # noqa: E402


def random_invertible(n: int, rng: np.random.Generator) -> np.ndarray:
    """Random 0/1 matrix of size n x n, invertible over Z2."""
    while True:
        a = rng.integers(0, 2, (n, n))
        try:
            CySPT.gf2_inverse(a)
            return a
        except np.linalg.LinAlgError:
            pass


def main() -> None:
    rng = np.random.default_rng(0)
    failures = 0
    for n in (1, 2, 5, 20, 60, 120, 250):
        a = random_invertible(n, rng)
        identity = np.eye(n, dtype=np.int64)
        new = CySPT.gf2_inverse(a)
        new_ok = np.array_equal(a @ new % 2, identity) and np.array_equal(new @ a % 2, identity)
        old = np.mod(np.round(np.linalg.inv(a) * np.linalg.det(a)), 2).astype(np.int64)
        old_ok = np.array_equal(a @ old % 2, identity)
        print(f'n = {n:3d}: exact inverse {"ok" if new_ok else "WRONG"}, '
              f'previous formula {"ok" if old_ok else "wrong"}')
        failures += not new_ok
    if failures:
        raise SystemExit(1)
    print('ALL OK')


if __name__ == '__main__':
    main()
