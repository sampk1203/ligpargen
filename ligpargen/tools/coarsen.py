"""
coarsen.py
==========
Importable wrapper around 00_coarsen_datafile.coarsen().

Exposes:
    coarsen(datafile, prefix=None) -> (coarse_data_path, topology_json_path)

The coarsen function collapses a LigParGen-combined LAMMPS data file from
per-atom unique types into (mol_id, element) grouped types, remapping all
*Coeffs sections to the new type indices so force field values are preserved.
"""

import os
import importlib.util
import sys


def _load_coarsen_module():
    """Locate and load 00_coarsen_datafile.py relative to this package root."""
    # Package tree: <repo>/ligpargen/tools/coarsen.py
    # 00_coarsen_datafile.py lives at <repo> root (tracked in git). Walk up from here.
    here = os.path.dirname(os.path.abspath(__file__))
    # here = <repo>/ligpargen/tools
    candidate = os.path.join(here, os.pardir, os.pardir,
                             "00_coarsen_datafile.py")
    candidate = os.path.normpath(candidate)
    if not os.path.isfile(candidate):
        raise FileNotFoundError(
            f"00_coarsen_datafile.py not found at expected path: {candidate}\n"
            f"Place it at the repo root, alongside setup.py.")
    spec = importlib.util.spec_from_file_location("_coarsen_datafile", candidate)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_mod = None


def coarsen(datafile, prefix=None):
    """
    Collapse a LigParGen LAMMPS data file to (mol × element) types,
    preserving all force field coefficients.

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
