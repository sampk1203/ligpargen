"""
Dispatcher for the multi-molecule PSF+PDB -> LAMMPS pipeline built on top of
ligpargen (extract unique species, split oversized ones, run ligpargen per
piece, stitch, combine). Installed as the `molpargen` console script.

Usage:
    molpargen run <psf> <pdb> <species_dir> <out.lmp> [options]
    molpargen run-mol2 <mol2> <species_dir> <out.lmp> [options]
    molpargen extract <psf> <pdb> <outdir>
    molpargen extract-mol2 <mol2> <outdir>
    molpargen split <manifest> <species_dir> <species_id> <outdir> [options]
    molpargen stitch <manifest> <frag_manifest> <fragment_dir> <species_dir> [options]
    molpargen combine <manifest> <original.pdb> <species_dir> <out.lmp> [options]
    molpargen check-connectivity <reference.lmp> <out.lmp> [options]

`run` (PSF+PDB) and `run-mol2` (single mol2 file) are the two full-pipeline
entry points -- pick whichever matches your input; both call ligpargen
itself as needed and produce the same out.lmp. The other subcommands are
the individual stages, exposed for debugging/re-running a single step.
"""
import sys

from ..pipeline import (
    run_pipeline,
    run_pipeline_mol2,
    extract_unique,
    extract_unique_mol2,
    split_species,
    stitch_fragments,
    combine_lammps,
    check_connectivity,
)

SUBCOMMANDS = {
    'run': run_pipeline.cli,
    'run-mol2': run_pipeline_mol2.cli,
    'extract': extract_unique.cli,
    'extract-mol2': extract_unique_mol2.cli,
    'split': split_species.cli,
    'stitch': stitch_fragments.cli,
    'combine': combine_lammps.cli,
    'check-connectivity': check_connectivity.cli,
}


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ('-h', '--help'):
        print(__doc__)
        return 0 if args and args[0] in ('-h', '--help') else 1
    sub, rest = args[0], args[1:]
    if sub not in SUBCOMMANDS:
        print(f"Unknown subcommand: {sub!r}")
        print(f"Available: {', '.join(SUBCOMMANDS)}")
        return 1
    return SUBCOMMANDS[sub](rest)


if __name__ == "__main__":
    sys.exit(main())
