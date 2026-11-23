"""
coarsen.py
==========
Importable wrapper around 00_coarsen_datafile.coarsen().

Exposes:
    coarsen(datafile, prefix=None) -> (coarse_data_path, topology_json_path)

The coarsen function collapses a LigParGen-combined LAMMPS data file from
per-atom unique types into (mol_id, element, coeff) grouped types for atoms,
(mol_id, elem-pair, coeff) for bonds, (mol_id, elem-triple, coeff) for
angles, and (mol_id, element-flavor-quad, coeff) for dihedrals — where
"flavor" is a 1-hop neighbor-element fingerprint that distinguishes OPLS
torsion types sharing the same 4-element sequence but different local
environments (e.g. backbone C-C-C-C vs branch C-C-C-C adjacent to CF₂).
Every level appends its FF coeff tuple as a tiebreak so the coarsened type
index always uniquely identifies the FF parameters, never merging two
distinct OPLS types that happen to share the same element/elem-sequence.
All *Coeffs sections are remapped to the new type indices so force field
values are preserved without averaging or silent data loss.
"""

import os
import importlib.util
import sys


def _load_coarsen_module():
    from . import coarsen_datafile
    return coarsen_datafile


_mod = None


def coarsen(datafile, prefix=None):
    """
    Collapse a LigParGen LAMMPS data file to coarsened types,
    preserving all force field coefficients without loss or averaging.

    Args:
        datafile: path to LAMMPS data file produced by combine_lammps
        prefix:   output prefix (default: datafile with .data/.lmp stripped)

    Returns:
        (coarse_data_path, topology_json_path)
    """
    global _mod
    if _mod is None:
        _mod = _load_coarsen_module()
    if prefix is None:
        for ext in (".data", ".lmp"):
            if datafile.endswith(ext):
                prefix = datafile[:-len(ext)]
                break
        else:
            prefix = datafile
    return _mod.coarsen(datafile, prefix)
