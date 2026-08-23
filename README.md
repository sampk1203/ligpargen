# LAMMPS force field pipeline (EMC → OPLS-AA/ligpargen → LAMMPS data file)

Turns a combined multi-molecule PSF+PDB (polymer chains + ions, e.g. from EMC)
into a single full-system LAMMPS data file with OPLS-AA parameters, by:
finding unique molecule species, running ligpargen per species (splitting
large ones into ligpargen/BOSS-sized fragments first), stitching fragment
parameters back together, then reassembling the full system with real
coordinates.

This pipeline ships as part of the `ligpargen` package (`ligpargen/pipeline/`)
and is exposed as its own console script, `molpargen`, installed
alongside the normal single-molecule `ligpargen` command by the same
`pip install -e .`. No separate scripts to keep track of, no `--scripts-dir`,
no `sys.path` tricks — `molpargen` is on `PATH` the same way
`ligpargen` is.

## Quick start — one command

Two entry points, pick whichever matches your input:

```bash
# PSF + PDB (e.g. from EMC)
molpargen run input.psf input.pdb species_dir/ out.lmp --box-from-pdb

# single mol2 file (any packing tool that can emit mol2 -- packmol, moltemplate,
# or a hand-built system; works for any packed multi-molecule system, not just EMC)
molpargen run-mol2 input.mol2 species_dir/ out.lmp --box-from FILE
```

Both run the whole pipeline end to end and block until `out.lmp` exists (or
raise on the first failing stage — neither silently continues past a broken
step). Nothing else below is required for normal use; the rest of this file
documents the individual stages for debugging, re-running a single step, or
understanding what `run`/`run-mol2` are doing under the hood.

**Requires:** `ligpargen` on `PATH` (installed by this same package).

**mol2 input**: mol2 already carries coordinates, elements (via SYBYL atom
type), connectivity, and partial charges all in one file — `run-mol2` reads
those directly, no PSF/PDB synthesis and no `obabel` round-trip. A packed
multi-molecule system as mol2 is just several `@<TRIPOS>MOLECULE` blocks
concatenated back to back in one file (what packing tools normally emit);
`run-mol2` handles that natively. mol2 files generally don't carry box
dimensions, so pass `--box-from a-reference.data` (from your packing
tool's own output, if it wrote one) — `--box-from-pdb` only works if
whatever wrote your mol2 also produced a `CRYST1`-bearing PDB alongside it.

**If any species has a `PLACEHOLDER` (-9999) charge/LJ entry in `out.lmp`**
(bare ions — Li⁺, Na⁺, ... — that ligpargen can't parametrize):
```bash
grep PLACEHOLDER out.lmp
```
and fill those in by hand from `$BOSSdir/oplsaa.par` before using the file.

---

## Pipeline stages

```
input.psf + input.pdb  ─┐                    input.mol2  ─┐
                         ├─ extract ──────────────────────┤
                         ▼                                 ▼
              species_dir/manifest.json, species_N.pdb, (+ system_from_mol2.pdb for mol2)
        │
        ▼
  (per species, if > --max-atoms)
split    →  fragments_N/fragment_N_K.pdb, fragment_manifest_N.json
        │
        ▼
  ligpargen (per species_N.pdb or fragment_N_K.pdb)  →  *.lammps.lmp
        │
        ▼
  (per split species)
stitch   →  species_dir/species_N.lammps.lmp
        │
        ▼
combine  →  out.lmp   (full system, real coords, real box)
```

`molpargen run` (PSF+PDB) and `molpargen run-mol2` (mol2) only differ in
stage 0 — `extract_unique.run()` vs `extract_unique_mol2.run()`. Both hand
off to the exact same `run_pipeline.run_from_manifest()` for everything
past that: split/ligpargen/stitch/combine don't know or care which format
the species came from. `run-mol2` calls the `ligpargen` binary itself where
BOSS actually needs to run, same as `run`. Everything past this point is
for running a stage manually, via the matching subcommand.

---

## `molpargen run` / `molpargen run-mol2` — orchestrators (normal entry point)

```bash
molpargen run input.psf input.pdb species_dir/ out.lmp [options]
molpargen run-mol2 input.mol2 species_dir/ out.lmp [options]
```

| Option | Default | Meaning |
|---|---|---|
| `--max-atoms N` | 220 | BOSS/ligpargen atom cap; species above this get split |
| `--max-orbitals N` | 480 | BOSS orbital-count budget passed to the splitter (hard BOSS limit: 550) |
| `--cgen CM1A\|CM1A-LBCC` | CM1A | charge model passed to ligpargen |
| `--dry-run` | off | print stages without running them; stops after stage 0 (extraction) since `manifest.json` doesn't exist yet to plan the rest against |
| `--force-split` | off | pass `--force` to the splitter (skip unsaturated-backbone-atom warnings) |
| `--box-from FILE` | — | copy box bounds from this LAMMPS data file (passed to combine) |
| `--box-from-pdb` | off | read box from `input.pdb`'s `CRYST1` record instead (passed to combine) |

Species with `n_atoms == 1` (bare ions) are skipped by ligpargen entirely —
`combine` synthesizes a placeholder entry for them (see above).

---

## Running stages individually (debugging / re-running one step)

### `molpargen extract`
Finds unique molecule species (by Weisfeiler-Leman graph isomorphism on
element-labeled connectivity, not just atom/bond counts) and writes one
representative single-molecule PDB per species.

```bash
molpargen extract input.psf input.pdb species_dir/
```
Writes `species_dir/manifest.json`, `species_dir/species_N.pdb`.

### `molpargen extract-mol2`
Same as `extract`, from a single mol2 file instead. Also writes
`species_dir/system_from_mol2.pdb` (full system, derived directly from the
mol2's own coordinates) — this is what `combine`/`run-mol2` use in place of
the original PDB for real system coordinates.

```bash
molpargen extract-mol2 input.mol2 species_dir/
```
Writes `species_dir/manifest.json`, `species_dir/species_N.pdb`,
`species_dir/system_from_mol2.pdb`. The manifest's schema is identical to
`extract`'s (same `species`/`instances` keys) — every stage past this one
reads either manifest the same way.

### `molpargen split`
For one species over the atom/orbital cap: splits along the backbone
(graph-topology BFS tree-diameter, side branches never cut) into
overlapping fragments small enough for ligpargen/BOSS, capping cut bonds
with H.

```bash
molpargen split manifest.json species_dir/ SPECIES_ID outdir/ \
    [--max-atoms N] [--max-orbitals N] [--force]
```
Writes `outdir/fragment_SPECIES_ID_K.pdb`, `outdir/fragment_manifest_SPECIES_ID.json`.
Prints "no split needed" and writes nothing if the species fits in one
fragment after all — feed `species_N.pdb` to ligpargen directly in that case.

### ligpargen (per species or fragment PDB)
```bash
ligpargen -i species_N.pdb -n species_N -c CHARGE -o 0 -cgen CM1A
```
Run from inside `species_dir/` (or `fragments_N/` for fragments) so the
`.lammps.lmp` output lands next to its input. Fragment charge is always `0`
regardless of the whole species' net charge — the stitch step's
redistribution corrects any sum drift afterward.

### `molpargen stitch`
Merges per-fragment ligpargen outputs back into one `species_N.lammps.lmp`,
in the exact format `combine` expects. Drops any bond/angle/dihedral/improper
touching a cap atom; keeps every other term exactly once even when it spans
two adjacent fragments' overlap region (deduplicated by canonical atom-index
form, not gated on single-fragment ownership — a term that straddles a cut
point is real in whichever fragment(s) can see it, so ownership alone would
silently drop it). Redistributes residual charge equally across all atoms to
hit the manifest's `expected_net_charge`.

```bash
molpargen stitch manifest.json fragment_manifest_N.json \
    fragment_dir/ species_dir/ [--charge-tol F]
```
`--charge-tol F` (default 0.5): warn, don't fail, if residual charge after
redistribution exceeds this.

### `molpargen combine`
Reassembles the full system: real per-instance coordinates from the
original PDB, ligpargen/stitched parameters per species, correct
per-species type offsets, box from a reference data file or the PDB's
`CRYST1` record, per-atom periodic wrapping, and an overlap sanity scan.

```bash
molpargen combine manifest.json original.pdb species_dir/ out.lmp \
    [box_from.data] [--box-from-pdb]
```
Without either box option, the box falls back to atom-coordinate bbox + 10 Å
padding — **not physically meaningful** if the source PDB has chains
unwrapped across the periodic boundary (common for EMC-packed systems); use
one of the box flags whenever you have a real box available.

Bare-ion (1-atom) species get a synthesized entry: real mass, but
`charge`/`epsilon`/`sigma` all set to the sentinel `-9999` and flagged
`# PLACEHOLDER` in the output — fill from `oplsaa.par` before running.

---

## `molpargen check-connectivity` — verify a combined system against a reference

Compares atom counts, bond/angle/dihedral/improper sets, per-atom bond
degree, and element-pair bond histogram between two full-system LAMMPS data
files (e.g. `out.lmp` vs a known-good reference).

```bash
molpargen check-connectivity reference.lmp out.lmp
molpargen check-connectivity reference.lmp out.lmp --map coord --tol 0.05
molpargen check-connectivity reference.lmp out.lmp --pdb system.pdb --verbose
```

| Option | Default | Meaning |
|---|---|---|
| `--map id\|coord\|sequential` | id | how to pair reference atoms to output atoms: same IDs, nearest-coordinate match within `--tol`, or 1↔1 by sorted-ID file order |
| `--tol F` | 0.05 | coordinate-matching tolerance in Å, only used with `--map coord` |
| `--pdb FILE` | — | source PDB for element lookup (element-pair histogram); falls back to mass-based element guessing without it |
| `--verbose` | off | print identical-set confirmations and full per-key histogram, not just mismatches |

Exit code 0 = topologies consistent, 1 = differences found (details printed
either way). Run this after every `combine` (or after patching `stitch`/`split`)
against a trusted reference before using `out.lmp` in production — it's the
only step that catches silently-dropped bonds/angles/dihedrals at fragment
boundaries or elsewhere.

---

## Using the pipeline as a library

Every stage's core logic lives in `ligpargen.pipeline.<stage>.run(...)`
(`extract_unique`, `extract_unique_mol2`, `split_species`, `stitch_fragments`,
`combine_lammps`, `check_connectivity`), independent of the CLI parsing in
each module's `cli(argv)`. `ligpargen.pipeline.run_pipeline.run(...)` /
`run_pipeline_mol2.run(...)` call all of them directly — no subprocess, no
separate script files — so you can import and call any stage (or the whole
pipeline) from your own Python code:

```python
from ligpargen.pipeline import run_pipeline, run_pipeline_mol2

run_pipeline.run("input.psf", "input.pdb", "species_dir", "out.lmp",
                  box_from_pdb=True)

run_pipeline_mol2.run("input.mol2", "species_dir", "out.lmp",
                       box_from="reference.data")
```

Both orchestrators share the same format-independent tail via
`run_pipeline.run_from_manifest(...)`, so a custom format-3 input source
only ever needs its own `extract_*.run()` writing the same manifest schema
— everything past stage 0 is already generic.

---

## Common failure points

- **`ligpargen ran but output not found`** — check ligpargen's own stdout
  (printed above the error) for a BOSS failure; usually an orbital/integral
  cap issue the splitter's pre-check estimator missed, or a genuinely bad
  input geometry.
- **`PLACEHOLDER` values reach a production run** — your LAMMPS input script
  should have a hard stop checking for this after `read_data`; don't remove
  it. Fill from `oplsaa.par` first.
- **Bonds/angles/dihedrals missing at fragment boundaries** —
  `check-connectivity` will show this as missing terms clustered at a
  handful of atom-index pairs with degree mismatches. Fixed as of the
  current `stitch_fragments.py` (terms are no longer gated on
  single-fragment ownership); if you're on an older copy, update it.
- **Non-orthogonal box** — `combine --box-from-pdb` only supports orthogonal
  cells (α=β=γ=90°); use `--box-from` with a LAMMPS data file for triclinic
  cells instead (tilt factors aren't written by this pipeline either way —
  check your reference file supports it).
