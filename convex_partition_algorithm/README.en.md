# Convex Partition Algorithm

[한국어](README.md) | English

A stack-based algorithm that partitions a simple polygon into convex polygons **using diagonals only** (no new points are added).
It also includes a dynamic programming solver for the optimal partition, so the results can be compared.

## Files

| File | Description |
|---|---|
| `convex_partition.py` | Stack-based convex partition algorithm, examples, random tests, plotting |
| `optimal_partition.py` | Dynamic programming for the minimum number of pieces (for comparison) |
| `requirements.txt` | Required packages (`matplotlib`, for plotting) |

## Setup and Usage

```bash
python3 -m venv myenv
source myenv/bin/activate
pip install -r requirements.txt
```

```bash
python3 convex_partition.py                  # run built-in examples (L-shape, U-shape, Comb, Star, Notch)
python3 convex_partition.py polygon.txt      # read a polygon from a file ("x y" per line, counterclockwise)
python3 convex_partition.py --random 50      # test on 100 random polygons with 50 vertices
python3 convex_partition.py --compare 50     # compare piece counts with the optimal solution
```

Options: `--quiet` (turn off step-by-step logs), `--save out.png` (save the figure), `--count N`, `--seed S`.
In environments that cannot open a plot window (e.g. WSL), the figure is saved to `convex_partition.png` automatically.

The step-by-step logs printed to the terminal are in Korean.

## Algorithm

Input: a list of vertices in counterclockwise (CCW) order

**Step 1–4: Build the stack**
1. Walk the polygon counterclockwise, find the first reflex vertex, set it as `vt`, and push it onto the stack.
2. Continue counterclockwise from `vt`. For each vertex `w`:
   - If the segment `vt`–`w` is not a diagonal of the polygon, set `vt` to the vertex before `w` and push it.
   - If `w` is a reflex vertex, set `vt = w` and push it.
   - After a full loop, go to the next step.

**Step 5–8: Process stack / stack2**
- Pop from the stack into stack2, then compare the stack top `a` with the stack2 top `b`:
  - `b` is already convex → discard `b`
  - Drawing diagonal `a`–`b` makes `b` convex → draw it and discard `b`
  - `b` stays reflex even after drawing `a`–`b` → draw it and keep `b` in stack2
  - Then move `a` to stack2. Repeat until the stack is empty.

**Step 9–13: Fix the remaining reflex vertices**
- Pop vertices from stack2 one by one. If a vertex is still reflex:
  1. Shoot the **bisector** of its largest angle. Starting where it hits the boundary, check vertices
     alternately to the right and left, and draw the first diagonal that makes the vertex convex.
  2. If no single diagonal works, keep adding the diagonal that reduces the largest angle the most,
     until the vertex becomes convex.

(The previous approach, "scan counterclockwise", is kept as a comment inside `_phase3`.)

### Implementation details

- **Visibility**: For each vertex, the set of visible vertices is computed with a rotational sweep (once, when first needed), so the diagonal test takes O(1).
- **Face tracking**: Each new diagonal splits a face, and the faces are recorded. "Does it cross an existing diagonal?" becomes "Are both endpoints in the same face?".
- `ConvexPartition(points, fast=False)` uses the previous approach (checking every edge for intersection each time).

## Results

Compared with the optimal solution on random polygons (seed=0, 100 polygons each):

| Vertices | Avg. pieces | Optimal | Ratio | Worst case |
|---|---|---|---|---|
| 20 | 9.80 | 8.00 | 1.23× | 1.67× |
| 50 | 28.17 | 22.82 | 1.23× | 1.41× |

- In every test, all pieces are convex.
- It produces about 23% more pieces than optimal on average (similar to Hertel–Mehlhorn).

## Time Complexity

| | Time |
|---|---|
| This algorithm (with visibility) | about O(n² log n) worst case |
| Previous approach (`fast=False`) | O(n³) worst case |
| Hertel–Mehlhorn | O(n log n) |
| Optimal (Keil–Snoeyink) | O(n + r² min(r², n)) |

n: number of vertices, r: number of reflex vertices

## Optimal Solution (`optimal_partition.py`)

A Keil-style dynamic programming algorithm. For each sub-polygon cut off by a chord (i, j), it computes
the minimum number of pieces, using only diagonals that touch a reflex vertex as candidates.
Its results were verified against a brute-force search over all combinations on 400 polygons with 4–10 vertices.
