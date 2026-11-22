# Changes

## 2026-11 — cycle vectors in uint8 (`fast-v4`)

`closeCycle` stores each cycle as a uint8 vector over all edges instead of
int64. The fundamental cycles and the cache of closed cycles keep one such
vector per cycle, so on large complexes (thousands of cycles over tens of
thousands of edges) this took gigabytes; it now takes 8 times less. The
cycles of the minimal basis are still returned as int, so outputs, and
their types, are unchanged.

## 2026-11 — exact Z2 inverse, less memory, faster (`fast-v3`)

One correctness fix; otherwise same algorithm, same outputs.

### Fix: `inv_z2` is now exact

`inv_z2` computed the inverse over Z2 as `mod(round(inv(A) * det(A)), 2)` in
floating point. That is exact only while `det(A)` and the adjugate stay
below 2^53; on random 0/1 matrices it already fails at 60 x 60. The matrix
inverted in `UpdateMehlhorn` has about half the dimension of H1, so with
H1 of dimension ~100 or more at one scale the support vectors, and hence
the generators, could be wrong. It is replaced by Gauss-Jordan elimination
over Z2 (`gf2_inverse`), exact at any size (`tests/check_gf2_inverse.py`).
The unused `sympy` import is removed, so sympy is no longer needed.

### Speed and memory (`src/CySPT.pyx`)

- `getAnnotation` reduces rows stored as bit vectors (Python ints) and no
  longer builds identity matrices of size (#triangles + #cycles)^2. Same
  reduction order, same annotation. It accepts the dense `d2` as before,
  or the rows of `d2` as bit vectors.
- `UpdateMehlhorn` recomputes only the rows of B = L C that changed
  (rows lam..mu-1 of L; C does not change there), instead of all of B.
- `mulz2` multiplies in float64 (BLAS); the entries are small integers,
  exact in float64.
- Candidate cycle lengths come from the shortest-path distances,
  d(r,u) + d(r,v) + l(uv) - 2 d(r, lca(u,v)), once per graph; a cycle
  vector is built only for the candidates within 1e-9 (relative) of the
  shortest, and among those the choice is made with the original lengths
  `WEdges . C` and the original rules, so rounding cannot change it.
- `computeLabels` fills all shortest-path trees at once, level by level.
- `cycleAnnotation` is vectorised; with draws, identical tied cycles are
  compared once.

Checked with `tests/check_equivalence.py` against `fast-v2`: identical
outputs on `data/` (theta 0.05-0.4), on random clouds of 150 and 300
points and a noisy ring of 200 (theta 0.04-0.1), and on a random cloud of
1500 points at theta 0.012 and 0.016 (H1 of dimension up to 116; 78 s
before, 16 s now). At theta 0.02 (H1 of dimension 187) the previous
version ran out of memory on an 8 GB machine; this one takes 28 s.

## 2026-11 — faster minimal generators and death points

Same algorithm, same outputs. Verified identical (generators, edge lists,
death points) against the previous version on 22 real point clouds of
400-450 points at theta = 6, and on 3 of them at theta = 4, 6, 8.
On those clouds the time per cloud went from ~25 s to ~1.5-2 s.
The previous version is kept in the `old_draft` branch.

### `src/CySPT.pyx` (`CYSPT_IMPLEMENTATION = 'fast-v2'`)

- Edges are looked up in a dictionary instead of `list.index`.
- `closeCycle` builds each cycle once from the indices of its edges.
- The Z2 label of every sentinel edge with respect to a support vector is
  computed with one matrix product, instead of one dot product per use.
- `computeLabels` uses a `deque` for its breadth-first traversal.
- The orthogonality test in `findShortestNonOrtho` is evaluated for all
  (root, non-tree edge) candidates at once, keeping their order.
- `fixNTE` uses `nonTreeEdges`, which builds each component's edge set once
  and reproduces the candidate order of `shortestPathTree` exactly. The
  candidate order decides which cycle is chosen among equally short ones.

After editing `CySPT.pyx`, run `make`: Python loads the compiled module in
preference to the `.pyx`.

### `src/ComputeMinGensDP.py`

- New argument `n_jobs` (default 1): number of csv files processed in
  parallel. Outputs do not depend on it.
- Death points by incremental reduction over Z2: reduced boundary columns
  are kept by pivot, each cycle keeps its residual, and a residual is
  updated only when a new pivot equals its lowest edge. Columns are bit
  vectors (Python ints). Same death values as before.
- A file with no cycles now writes an empty death-point file and processing
  continues with the next file (previously the whole run stopped).
- Files are processed in sorted order; figures are closed after saving.

### Units

`Thresh` (theta) is a radius. Gudhi's Alpha filtration values are squared
radii; the complex is cut at `Thresh**2`, and death points are returned as
radii (or squared radii with `square=True`).

### `tests/check_equivalence.py`

Runs two MinGens folders on the same csv files at several theta and
requires identical outputs. See the docstring for usage with a worktree of
`old_draft`.
