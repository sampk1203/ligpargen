#!/usr/bin/env python3
"""
Compare connectivity of two full-system LAMMPS data files
(e.g. a known-good reference vs out.lmp from 2_combine_lammps.py).

Checks (after optional atom mapping):
  - atom counts
  - bond sets (undirected pairs)
  - angle sets (canonical A-B-C)
  - dihedral sets (canonical A-B-C-D)
  - improper sets
  - per-atom bond degree
  - element-pair bond histogram (if masses allow element guess, or --pdb)

Atom mapping modes:
  --map id          assume same atom IDs / same order (default)
  --map coord       match atoms by coordinates within --tol (Angstrom)
  --map sequential  1->1, 2->2, ... by file order (ignore stored IDs)

Usage:
  python compare_lammps_topology.py reference.lmp out.lmp
  python compare_lammps_topology.py reference.lmp out.lmp --map coord --tol 0.05
  python compare_lammps_topology.py reference.lmp out.lmp --pdb system.pdb --verbose
"""
import sys
import os
import math
from collections import Counter, defaultdict


def parse_lammps_data(path):
    """Minimal full-system LAMMPS data reader.
    Returns dict with:
      atoms: {id: (mol, type, charge, x, y, z)}
      bonds: list of (type, a, b)
      angles: list of (type, a, b, c)
      dihedrals: list of (type, a, b, c, d)
      impropers: list of (type, a, b, c, d)
      masses: {type: mass}
      n_atoms, n_bonds, ... from header if present
    """
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]

    def section(name):
        start = None
        for i, l in enumerate(lines):
            if l.strip() == name:
                start = i + 1
                break
        if start is None:
            return []
        i = start
        while i < len(lines) and lines[i].strip() == '':
            i += 1
        out = []
        while i < len(lines):
            row = lines[i].strip()
            if row == '' or not row.split()[0].lstrip('-').isdigit():
                break
            out.append(lines[i].split())
            i += 1
        return out

    # header counts (best-effort)
    header = {}
    for l in lines:
        parts = l.split()
        if len(parts) >= 2 and parts[0].lstrip('-').isdigit():
            key = ' '.join(parts[1:]).lower()
            if key in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                       'atom types', 'bond types', 'angle types',
                       'dihedral types', 'improper types'):
                header[key] = int(parts[0])

    masses = {}
    for r in section('Masses'):
        masses[int(r[0])] = float(r[1])

    atoms = {}
    for r in section('Atoms'):
        # id mol type charge x y z  [optional image flags]
        aid = int(r[0])
        mol = int(r[1])
        atype = int(r[2])
        charge = float(r[3])
        x, y, z = float(r[4]), float(r[5]), float(r[6])
        atoms[aid] = (mol, atype, charge, x, y, z)

    bonds = [(int(r[1]), int(r[2]), int(r[3])) for r in section('Bonds')]
    angles = [(int(r[1]), int(r[2]), int(r[3]), int(r[4])) for r in section('Angles')]
    dihedrals = [(int(r[1]), int(r[2]), int(r[3]), int(r[4]), int(r[5]))
                 for r in section('Dihedrals')]
    impropers = [(int(r[1]), int(r[2]), int(r[3]), int(r[4]), int(r[5]))
                 for r in section('Impropers')]

    return dict(path=path, header=header, masses=masses, atoms=atoms,
                bonds=bonds, angles=angles, dihedrals=dihedrals,
                impropers=impropers)


def undirected(a, b):
    return frozenset((a, b))


def canon_angle(a, b, c):
    ends = tuple(sorted((a, c)))
    return (ends[0], b, ends[1])


def canon_torsion(a, b, c, d):
    t1 = (a, b, c, d)
    t2 = (d, c, b, a)
    return t1 if t1 <= t2 else t2


def mass_to_element(mass, tol=0.3):
    """Rough element guess from atomic mass (OPLS-style)."""
    table = [
        (1.008, 'H'), (12.011, 'C'), (14.007, 'N'), (15.999, 'O'),
        (18.998, 'F'), (22.990, 'NA'), (24.305, 'MG'), (28.086, 'SI'),
        (30.974, 'P'), (32.06, 'S'), (35.45, 'CL'), (39.098, 'K'),
        (79.904, 'BR'), (126.90, 'I'), (6.94, 'LI'),
    ]
    best, best_d = '?', 1e9
    for m, el in table:
        d = abs(mass - m)
        if d < best_d:
            best, best_d = el, d
    return best if best_d <= tol else '?'


def elements_from_data(data, pdb_elems=None):
    """Map atom id -> element symbol."""
    out = {}
    for aid, (mol, atype, ch, x, y, z) in data['atoms'].items():
        if pdb_elems and aid in pdb_elems:
            out[aid] = pdb_elems[aid]
        elif atype in data['masses']:
            out[aid] = mass_to_element(data['masses'][atype])
        else:
            out[aid] = '?'
    return out


def read_pdb_elements(pdb_path):
    elems = {}
    idx = 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                idx += 1
                elems[idx] = l.split()[-1].upper()
    return elems


def build_coord_map(ref_atoms, act_atoms, tol):
    """Map ref atom id -> act atom id by nearest coordinate within tol.
    Greedy 1-1 matching. Returns (map, unmatched_ref, unmatched_act)."""
    # spatial hash
    inv = 1.0 / max(tol, 1e-6)
    buckets = defaultdict(list)
    for aid, (_, _, _, x, y, z) in act_atoms.items():
        key = (int(math.floor(x * inv)), int(math.floor(y * inv)),
               int(math.floor(z * inv)))
        buckets[key].append((aid, x, y, z))

    mapping = {}
    used_act = set()
    unmatched_ref = []

    for rid, (_, _, _, rx, ry, rz) in sorted(ref_atoms.items()):
        rkey = (int(math.floor(rx * inv)), int(math.floor(ry * inv)),
                int(math.floor(rz * inv)))
        best, best_d2 = None, tol * tol
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for aid, ax, ay, az in buckets.get(
                            (rkey[0] + dx, rkey[1] + dy, rkey[2] + dz), []):
                        if aid in used_act:
                            continue
                        d2 = (rx - ax) ** 2 + (ry - ay) ** 2 + (rz - az) ** 2
                        if d2 < best_d2:
                            best, best_d2 = aid, d2
        if best is None:
            unmatched_ref.append(rid)
        else:
            mapping[rid] = best
            used_act.add(best)

    unmatched_act = [a for a in act_atoms if a not in used_act]
    return mapping, unmatched_ref, unmatched_act


def remap_pair(a, b, mp):
    return undirected(mp[a], mp[b])


def remap_angle(a, b, c, mp):
    return canon_angle(mp[a], mp[b], mp[c])


def remap_torsion(a, b, c, d, mp):
    return canon_torsion(mp[a], mp[b], mp[c], mp[d])


def compare_sets(label, expected, actual, verbose, limit=25):
    missing = expected - actual
    extra = actual - expected
    ok = not missing and not extra
    tag = "OK" if ok else "MISMATCH"
    print(f"  {label}: ref={len(expected)}  out={len(actual)}  [{tag}]")
    if missing:
        sample = sorted(missing, key=lambda x: tuple(sorted(x)) if isinstance(x, frozenset)
                        else x)[:limit]
        print(f"    missing in out ({len(missing)}): {sample}"
              + (" ..." if len(missing) > limit else ""))
    if extra:
        sample = sorted(extra, key=lambda x: tuple(sorted(x)) if isinstance(x, frozenset)
                        else x)[:limit]
        print(f"    extra in out ({len(extra)}): {sample}"
              + (" ..." if len(extra) > limit else ""))
    if verbose and ok:
        print(f"    identical sets")
    return ok


def bond_element_hist(bond_set, elems):
    hist = Counter()
    for pair in bond_set:
        a, b = tuple(pair)
        ea, eb = elems.get(a, '?'), elems.get(b, '?')
        key = tuple(sorted((ea, eb)))
        hist[key] += 1
    return hist


def run(ref_path, act_path, map_mode='id', tol=0.05, pdb_path=None, verbose=False):
    print(f"reference: {ref_path}")
    print(f"actual:    {act_path}")
    print(f"map mode:  {map_mode}" + (f"  tol={tol}" if map_mode == 'coord' else ""))

    ref = parse_lammps_data(ref_path)
    act = parse_lammps_data(act_path)

    print(f"\n[counts]")
    print(f"  atoms:     ref={len(ref['atoms']):6d}  out={len(act['atoms']):6d}")
    print(f"  bonds:     ref={len(ref['bonds']):6d}  out={len(act['bonds']):6d}")
    print(f"  angles:    ref={len(ref['angles']):6d}  out={len(act['angles']):6d}")
    print(f"  dihedrals: ref={len(ref['dihedrals']):6d}  out={len(act['dihedrals']):6d}")
    print(f"  impropers: ref={len(ref['impropers']):6d}  out={len(act['impropers']):6d}")

    # --- atom mapping ref_id -> act_id ---
    if map_mode == 'id':
        common = set(ref['atoms']) & set(act['atoms'])
        only_ref = set(ref['atoms']) - set(act['atoms'])
        only_act = set(act['atoms']) - set(ref['atoms'])
        mp = {i: i for i in common}
        if only_ref or only_act:
            print(f"\n[map id] ID sets differ: only_in_ref={len(only_ref)} "
                  f"only_in_out={len(only_act)}")
            if verbose:
                print(f"  only_in_ref sample: {sorted(only_ref)[:20]}")
                print(f"  only_in_out sample: {sorted(only_act)[:20]}")
        else:
            print(f"\n[map id] all atom IDs match ({len(mp)} atoms)")
    elif map_mode == 'sequential':
        rids = sorted(ref['atoms'])
        aids = sorted(act['atoms'])
        if len(rids) != len(aids):
            print(f"\n[map sequential] FAIL: different atom counts "
                  f"{len(rids)} vs {len(aids)}")
            return 1
        mp = {r: a for r, a in zip(rids, aids)}
        print(f"\n[map sequential] paired {len(mp)} atoms by sorted ID order")
    elif map_mode == 'coord':
        mp, un_ref, un_act = build_coord_map(ref['atoms'], act['atoms'], tol)
        print(f"\n[map coord] matched={len(mp)}  unmatched_ref={len(un_ref)} "
              f"unmatched_out={len(un_act)}  (tol={tol} A)")
        if un_ref or un_act:
            if verbose:
                print(f"  unmatched_ref sample: {un_ref[:15]}")
                print(f"  unmatched_out sample: {un_act[:15]}")
            print(f"  WARN: incomplete mapping — topology compare uses matched atoms only")
    else:
        raise ValueError(f"unknown map mode: {map_mode}")

    if not mp:
        print("No atom mapping — cannot compare topology")
        return 1

    # Drop terms that touch unmapped atoms
    def atoms_mapped(ids):
        return all(i in mp for i in ids)

    # --- topology sets in *actual* ID space ---
    ref_bonds = {remap_pair(a, b, mp)
                 for _, a, b in ref['bonds'] if atoms_mapped((a, b))}
    act_bonds = {undirected(a, b)
                 for _, a, b in act['bonds']
                 if a in set(mp.values()) and b in set(mp.values())}

    ref_angles = {remap_angle(a, b, c, mp)
                  for _, a, b, c in ref['angles'] if atoms_mapped((a, b, c))}
    act_angles = {canon_angle(a, b, c)
                  for _, a, b, c in act['angles']
                  if a in set(mp.values()) and b in set(mp.values())
                  and c in set(mp.values())}

    ref_dih = {remap_torsion(a, b, c, d, mp)
               for _, a, b, c, d in ref['dihedrals'] if atoms_mapped((a, b, c, d))}
    act_dih = {canon_torsion(a, b, c, d)
               for _, a, b, c, d in act['dihedrals']
               if all(x in set(mp.values()) for x in (a, b, c, d))}

    ref_imp = {remap_torsion(a, b, c, d, mp)
               for _, a, b, c, d in ref['impropers'] if atoms_mapped((a, b, c, d))}
    act_imp = {canon_torsion(a, b, c, d)
               for _, a, b, c, d in act['impropers']
               if all(x in set(mp.values()) for x in (a, b, c, d))}

    print(f"\n[connectivity] (mapped into out atom-ID space)")
    ok_b = compare_sets("bonds", ref_bonds, act_bonds, verbose)
    ok_a = compare_sets("angles", ref_angles, act_angles, verbose)
    ok_d = compare_sets("dihedrals", ref_dih, act_dih, verbose)
    ok_i = True
    if ref_imp or act_imp:
        ok_i = compare_sets("impropers", ref_imp, act_imp, verbose)
    else:
        print(f"  impropers: both empty [OK]")

    # degree
    deg_ref = Counter()
    for p in ref_bonds:
        for x in p:
            deg_ref[x] += 1
    deg_act = Counter()
    for p in act_bonds:
        for x in p:
            deg_act[x] += 1
    mapped_act_ids = set(mp.values())
    deg_bad = []
    for aid in sorted(mapped_act_ids):
        if deg_ref[aid] != deg_act[aid]:
            deg_bad.append((aid, deg_ref[aid], deg_act[aid]))
    print(f"\n[degree]")
    if deg_bad:
        print(f"  mismatches: {len(deg_bad)}")
        for aid, er, ea in deg_bad[:20]:
            print(f"    out atom {aid}: ref_degree={er}  out_degree={ea}")
        if len(deg_bad) > 20:
            print(f"    ... and {len(deg_bad) - 20} more")
    else:
        print(f"  OK: per-atom bond degrees match for all mapped atoms")

    # element-pair bond histogram
    pdb_elems = read_pdb_elements(pdb_path) if pdb_path else None
    # elems keyed by *out* atom ids for act; for ref we remap
    act_elems = elements_from_data(act, pdb_elems)
    # build ref elems in out-id space
    ref_elems_native = elements_from_data(ref, pdb_elems)
    ref_elems_out = {mp[r]: ref_elems_native[r] for r in mp if r in ref_elems_native}

    print(f"\n[bond element pairs]")
    h_ref = bond_element_hist(ref_bonds, ref_elems_out)
    h_act = bond_element_hist(act_bonds, act_elems)
    all_keys = sorted(set(h_ref) | set(h_act))
    hist_ok = True
    for k in all_keys:
        nr, na = h_ref.get(k, 0), h_act.get(k, 0)
        mark = "  " if nr == na else "!!"
        if nr != na:
            hist_ok = False
        if verbose or nr != na:
            print(f"  {mark} {k[0]:>3s}-{k[1]:<3s}  ref={nr:5d}  out={na:5d}")
    if hist_ok:
        print(f"  OK: element-pair bond counts match")
    else:
        print(f"  MISMATCH in element-pair bond counts")

    print(f"\n{'=' * 55}")
    all_ok = (ok_b and ok_a and ok_d and ok_i and not deg_bad and hist_ok
              and len(ref['atoms']) == len(act['atoms']))
    if all_ok:
        print("RESULT: topologies consistent")
        return 0
    print("RESULT: differences found (see above)")
    return 1


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    verbose = '--verbose' in args
    args = [a for a in args if a != '--verbose']

    map_mode = 'id'
    if '--map' in args:
        i = args.index('--map')
        map_mode = args[i + 1]
        del args[i:i + 2]

    tol = 0.05
    if '--tol' in args:
        i = args.index('--tol')
        tol = float(args[i + 1])
        del args[i:i + 2]

    pdb_path = None
    if '--pdb' in args:
        i = args.index('--pdb')
        pdb_path = args[i + 1]
        del args[i:i + 2]

    if len(args) != 2:
        print(__doc__)
        print("Usage: python compare_lammps_topology.py reference.lmp out.lmp "
              "[--map id|coord|sequential] [--tol 0.05] [--pdb system.pdb] [--verbose]")
        return 2

    return run(args[0], args[1], map_mode=map_mode, tol=tol,
               pdb_path=pdb_path, verbose=verbose)


if __name__ == "__main__":
    sys.exit(cli())
