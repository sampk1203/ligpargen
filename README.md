# **LigParGen v2.1**

**Author:** &nbsp;&nbsp;Israel Cabeza de Vaca Lopez</br>
**Email:**  &nbsp;&nbsp;&nbsp;israel.cabezadevaca@yale.edu // israel.cabezadevaca@icm.uu.se </br>
**Place:** &nbsp;&nbsp;&nbsp; William L. Jorgensen Lab at Yale University // Jens Carlsson Lab at Uppsala university</br>
**Date:** &nbsp;&nbsp;  2020-2021

**Description:** An automatic OPLS-AA parameter generator for small organic molecules using CM1A, 1.14CM1A and CM1A-LBCC charge models. LigParGen accepts any Open Babel molecular format including SMILES, PDB, MOL, MOL2, among others. Final OPLSAA parameter outputs will be written in topology/coordinate input files for BOSS, Q, Tinker, PQR, openMM, CHARMM/NAMD, Gromacs, LAMMPS, Desmond, and xplor softwares.

**Note:** This new version has been written from scratch but it is based on the **Leela Dodda** initial ligpargen python code.


**Note 2:** LigParGen is a wrapper around BOSS that reads and transforms OPLS atom types and force field parameters for use with different software packages. Therefore, issues related to molecule parameterization (incorrect atom types, incorrect torsional parameters,...) are due to BOSS and should be reported to Bill Jorgensen. On the other hand, issues related to output generation (incorrect formats, FEP problems,...) should be reported here.. 

**Note 3:** **Do NOT report** issues related to the LigParGen server version here. As mentioned earlier, this is a new implementation, and it has not yet been deployed on the server. If you encounter issues with the server version, please try rerunning your task with this improved version to see if the problem persists or if it has already been fixed.

New LigParGen features:

- The order and the name of the atoms will remain the same in the output files.
- This new version of the ligpargen includes a robust version of the alchemical transformation method to generate single and dual topologies for four different molecular mechanics softwares (BOSS, CHARMM/NAMD, Gromacs, and Tinker).
- Sanity checks to detect incorrect inputs have been implemented (incompatible charge model, net molecule charge, input format, ...).
- A log file to check the inputs, outputs, and intermediate processes has been created. This allows tracking the input information (charge model used, input molecule,...) and also provides useful warnings in case of error.
- Automatic net charge detection in the input molecule. If the user specifies a different charge than one automatically estimated, the log file will include a warnning.
- Atom XYZ positions in the molecule input remain unchanged in the output files.
- Hydrogen atoms will be added automatically just for SMILES inputs. Molecules in any other input format require to have all hydrogens to provide more flexibility of the protonation states. In this way, the user can avoid parameterization problems with tautomers or stereoisomers.
- Different molecule input format, net charges, charge models, and optimization steps can be used for each molecule in alchemical transformations.
- Bugs fixed (Q wrong torsion parameters, alchemical molecule overlap failure,...)
- Additional default input parameters have been included such as residue name, charge model,...

## **INSTALATION**

LigParGen requires the free BOSS software to generate the OPLSAA parameters.

1 - Download and install BOSS software from the official William L. Jorgensen lab website: http://zarbi.chem.yale.edu/software.html

BOSS is compiled for linux using 32 bits libraries so it can not run in windows using WSL. Alternatively, you can use a virtual machine in windows such as virtualBox to install a linux distro (ubuntu, centos, ...) and run LigParGen.

1.1 - Set the BOSSdir enviromental variable:

  - bashrc
  
            export BOSSdir=PATH_TO_BOSS_DIRECTORY
  - cshrc
        
            setenv BOSSdir PATH_TO_BOSS_DIRECTORY

    **TIP:** add this command line in your ~/.bashrc or ~/.cshrc file.

2 - Download and install conda (anaconda or minicoda):  

    wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh

or

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

3 - Create and activate an enviroment for python3.7 (Everything was tested with 3.7):

    conda create --name py37 python=3.7

    conda activate py37  
**Optional TIP:** add this line to your .bashrc/.cshrc,....

4- Install rdkit and openbabel in py37 enviroment

    conda install -c rdkit rdkit
    conda install -c conda-forge openbabel

5 - Download and install LigParGen

    git clone https://github.com/Isra3l/ligpargen.git

    pip install -e ligpargen

**Optional:** Check your installation by runing the tests included in the ligpargen folder.

    cd ligpargen;python -m unittest test_ligpargen/test_ligpargen.py

**TIP:** Do not forget to activate your py37 enviroment before using LigParGen.

    conda activate py37

## **USAGE**

Input arguments list:

- ' -s ': SMILES input  (Ex. -s 'CCCC' )
- ' -i ': PDB, MOL and MOL2 files + any input format supported by Open Babel (Ex. -i benzene.pdb )
- ' -n ': Molecule file name (Ex. -n benzene - it will produce benzene.gmx.gro, benzene.charmm.rtf,...)
- ' -p ': Folder for the output files (Ex. -p bnz )
- ' -r ': Residue name for the output files (Ex. -r BNZ)
- ' -c ': Molecule net charge (Ex. -c +1)
- ' -o ': Number of optimizations (Ex. -o 3)
- ' -cgen ': Charge model to be used - CM1A or CM1A-LBCC (Ex. -cgen CM1A)

For alchemical transformations:

- ' -sb ': SMILES input for molecule B (Ex. -sb 'CCCCO' )
- ' -ib ': PDB, MOL and MOL2 files + any input format supported by Open Babel (Ex. -ib phenol.pdb )
- ' -cb ': Molecule net charge (Ex. -cb +1)
- ' -ob ': Number of optimizations (Ex. -ob 3)
- ' -cgenb ': Charge model to be used - CM1A or CM1A-LBCC (Ex. -cgenb CM1A)

**Notes:**

- CM1A is automatically scaled by 1.14 in neutral molecules.
- CM1A-LBCC is just for neutral molecules and it is also scaled by 1.14.
  
## **Examples**

Molecule template generation:

    ligpargen -i phenol.pdb -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A
    ligpargen -s 'c1ccc(cc1)O' -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A
    ligpargen -s 'c1ccc(cc1)O'
    ligpargen -s 'c1ccc(cc1)O' -n phenol -cgen CM1A-LBCC

Alchemical transformations:

    ligpargen -i phenol.pdb -ib benzene.pdb -n phenolToBenzene -p phenol2bnz -r A2B -c 0 -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC
    ligpargen -i phenol.pdb -ib benzene.pdb -n phenolToBenzene
    ligpargen -s 'c1ccc(cc1)O' -sb 'c1ccccc1' -n phenol_benzene
    ligpargen -s 'c1ccc(cc1)O' -sb 'c1ccccc1' -n phenol_benzene -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC
    ligpargen -i phenol.pdb -sb 'c1ccccc1' -n phenol_benzene -o 0 -cgen CM1A -cb 0 -ob 1 -cgenb CM1A-LBCC

Default values:

- -n : molecule
- -p : Current working directory
- -r : MOL
- -c,-cb : Automatically determined from input molecule
- -o,-ob : 0
- -cgen,-cgenb: CM1A-LBCC for neutral molecules and CM1A for charged molecules.

For help use the -h flag:

    ligpargen -h


## **REPORTING ISSUES**

If you encounter an issue related to parameters (weak improper torsions,...), check the Zmat file provided by LigParGen. This file is directly generated by BOSS and contains all the relevant OPLS parameter information.

Example of atom type information:

    Final Non-Bonded Parameters for QM (AM1 CM1Ax1.14) Atoms:                      
    800  6 CM  -0.250555  3.550000  0.076000                                                                     
    801  6 CM  -0.250555  3.550000  0.076000                                       
    802  1 HC   0.125278  2.500000  0.030000                                       
    803  1 HC   0.125278  2.500000  0.030000 

Example of force-field interaction types:

                    Variable Dihedrals follow     (3I4,F12.6)                   
    6 162 162    2.000000
    7 222 222    2.000000
    8 162 162    2.000000

You can lookup the interaction types in the original boss file (BOSSDir/oplsaa.par and BOSSDir/oplsaa.sb).

Dihedral 162 is defined here:

    162   0.0       5.0       0.0       0.0        Z -CA-X -Y      improper torsion 9/08 was 2.2  

This line indicates the atom types ("CA") used for selecting this interaction type.


Additionally, use the LigParGen debug flag to retain the intermediate BOSS-generated files. In the file named out, you will find the OPLS parameters and atom types directly assigned by BOSS.

For more details, I strongly recommend reading the BOSS manual (PDF) located in the BOSS folder.

If you want to report an issue, please include all the necessary files, command lines, and information required to reproduce it quickly. Do not provide just an image of your molecule; include the SMILES, PDB files, and any other essential data needed to run it efficiently. 

**Note:** This section was suggested by Andrew Jewett. Thank you for your feedback!  

## **REFERENCES**

Please do not forget to cite the following references:

1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
        Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
        Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336

2. 1.14*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations

## What this fork does
# LAMMPS force field pipeline (EMC → OPLS-AA/ligpargen → LAMMPS data file)

Turns a combined multi-molecule PSF+PDB (polymer chains + ions, e.g. from EMC)
into a single full-system LAMMPS data file with OPLS-AA parameters, by:
finding unique molecule species, running ligpargen per species (splitting
large ones into ligpargen/BOSS-sized fragments first), stitching fragment
parameters back together, then reassembling the full system with real
coordinates. A final coarsening step (on by default) collapses LigParGen's
per-atom unique type numbering into a compact, physically meaningful type
set suitable for GNN force-field training and CG-MD workflows.

This pipeline ships as part of the `ligpargen` package (`ligpargen/pipeline/`)
and is exposed as its own console script, `molpargen`, installed
alongside the normal single-molecule `ligpargen` command by the same
`pip install -e .`. The coarsening script (`00_coarsen_datafile.py`) lives at
the repo root and is loaded dynamically by `ligpargen/tools/coarsen.py` — keep
it there. No separate scripts to track, no `--scripts-dir`, no `sys.path`
tricks — `molpargen` is on `PATH` the same way `ligpargen` is.

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
        │
        ▼
  (default, unless --no-coarsen)
coarsen  →  out_coarse.data + out_topology.json
             (mol × element-flavor types; verified clean before write)
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
| `--no-coarsen` | off | skip the post-combine coarsening step (coarsening is **on by default**) |

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

## `00_coarsen_datafile.py` — collapse per-atom types to mol × element-flavor

After `combine` produces `out.lmp`, the file still carries LigParGen's original
per-atom unique type numbering (every atom gets its own type). This step
collapses those into a much smaller set of physically meaningful types, making
the file suitable for GNN force-field training and coarse-grained MD workflows.

**Runs automatically** as the final stage of `molpargen run` / `molpargen
run-mol2` (pass `--no-coarsen` to skip). Can also be run standalone:

```bash
python 00_coarsen_datafile.py out.lmp [output_prefix]
```

Outputs:
- `<prefix>_coarse.data` — coarsened LAMMPS data file
- `<prefix>_topology.json` — per-molecule atom/bond/angle/dihedral lists + element-sequence dihedral pattern map

### Type-collapsing scheme

| Level | Key |
|---|---|
| Atoms | `(mol_id, element)` |
| Bonds | `(mol_id, canonical elem-pair)` |
| Angles | `(mol_id, canonical elem-triple)` |
| Dihedrals | `(mol_id, canonical flavor-quad, coeff)` — see below |

All `*Coeffs` sections are remapped to the new type indices. No values are
averaged or silently discarded.

### Dihedral keying — why it's richer than atom/bond/angle

A bare element-quad key (`C-C-C-C`) collapses genuinely different OPLS torsion
types that share the same four elements but different local bonded environments
(e.g. a backbone C-C-C-C vs one adjacent to a CF₂ group). The fix is a
two-level key:

1. **1-hop neighbor fingerprint ("flavor")**: each atom is represented as
   `(element, sorted-tuple-of-neighbor-elements)` rather than element alone.
   A backbone carbon bonded to `[C, C, H, H]` and a CF₂-adjacent carbon bonded
   to `[C, C, F, F]` get distinct flavors even though both are `C`.

2. **Coeff tiebreak**: if two dihedral entries still share the same
   `(mol, flavor-quad)` (rare geometric ambiguity), their Fourier coefficient
   tuple is appended to the key. This guarantees each new type index uniquely
   identifies the FF parameters — the same invariant that atoms (mol × element),
   bonds (mol × elem-pair), and angles (mol × elem-triple) already satisfy.

### Self-check before writing

`verify_coarsening()` runs automatically before the output file is written. For
every level it checks:

- **(a) no drop** — every old type that had a Coeffs entry maps to a new type
  that also has one.
- **(b) no silent merge** — every old type folded into the same new type carried
  identical coeff values.

On failure it prints itemised diffs and raises `RuntimeError` — the coarsened
file is never written with bad data. With the flavor + coeff key, the dihedral
check always passes by construction; the check remains as a guard for future
changes to the collapsing scheme.

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
- **`coarsening self-check failed`** — `verify_coarsening()` found that two
  or more old types with *different* FF parameters were being collapsed into
  the same new type. With the current dihedral flavor + coeff key this should
  not occur; if it does on a new system, check whether the `Dihedral Coeffs`
  section in `out.lmp` contains `None`/missing entries for some types (a BOSS
  edge case for certain improper torsions). Pass `--no-coarsen` to bypass
  while investigating, but don't leave it bypassed in production.
- **Atoms section parse crash (`invalid literal for int: '#'`)** — the LAMMPS
  data file has inline comments on atom lines (e.g. `... 0 0 0 # mol1`). This
  is now handled: comment text after `#` is stripped before column parsing.
  If you see this on a freshly combined file, ensure you're running the current
  `00_coarsen_datafile.py`, not an older copy.
