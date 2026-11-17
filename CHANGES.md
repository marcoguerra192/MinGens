# Changes

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
