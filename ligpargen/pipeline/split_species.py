#!/usr/bin/env python3
"""
Step 1b (optional), between extract_unique_molecules.py and running
ligpargen. For any species over the BOSS/ligpargen atom cap, split it
into overlapping fragments small enough for ligpargen, capping cut
bonds with H. Run ligpargen on each fragment_N_M.pdb, then run
stitch_fragments.py (step 1c) to reassemble species_N.lammps.lmp,
which plugs into combine_lammps.py unmodified.

No bond-order or aromaticity input exists anywhere in this pipeline
(PSF/PDB carry connectivity only) -- so cuts are restricted to the
molecule's main-chain (backbone) path, found purely from graph
topology (BFS tree-diameter). Side branches (nitrile groups, CF3
groups, whatever) hang off backbone atoms and are never cut into,
regardless of what bond order they actually contain -- this is what
makes triple/double bonds safe without ever needing to know they're
there. See conversation for the full reasoning; short version: BOSS
perceives bond order fine from fragment geometry, but only if the
fragment's geometry hasn't already had that bond deleted by an
ill-placed cut. Backbone-only cutting prevents that, unconditionally.

Remaining blind spot: a double/triple bond ON the backbone itself
(not in a branch) is not detected. Backbone atoms are checked against
expected saturated (sp3-style) total connection count for their
element; any mismatch aborts rather than silently cutting near it.

Usage: python 1b_split_large_species.py manifest.json species_dir/ species_id outdir/ [--max-atoms N] [--force]
Writes: outdir/fragment_<species_id>_<k>.pdb, outdir/fragment_manifest_<species_id>.json
"""
import sys
import os
import json
import math

CHUNK_TARGET = 200      # atoms per fragment core budget (BOSS/ligpargen hard cap is 250)
TRUST_RADIUS = 3        # bonds -- an atom's params are only trusted >=3 bonds from any cap
CAP_BOND_LENGTH = 1.09  # Angstrom, C-H; fine as a generic cap-H placement distance

# BOSS's semi-empirical charge step (run 1 of ligpargen) has a hard cap on total
# valence orbitals in the system -- NOT on atom count. Observed directly from a
# BOSS stderr banner: "MAX. NUMBER OF ORBITALS: 550". Orbital cost per atom is
# element-dependent (minimal valence basis): H = 1 (1s), second-row elements
# (C/N/O/F) = 4 (2s + 3x2p), third-row-and-heavier default to 9 (ns + 3xnp +
# 5xnd) as a conservative guess when the element isn't in this table -- better
# to over-estimate and over-split than crash deep inside a Fortran binary.
# A flat atom-count budget (the old behavior) silently assumes ~1 orbital/atom,
# which is fine for a hydrocarbon but badly wrong for an F-heavy polymer like
# PVDF-HFP, where the true average is closer to 3-3.5 orbitals/atom -- that's
# exactly what caused a 236-atom fragment to hit 803 orbitals (550 cap) and
# fail with a cryptic "end-of-file during read, unit 26, file optzmat" deep in
# BOSS instead of an understandable error.
ORBITAL_WEIGHTS = {'H': 1, 'HE': 1,
                    'LI': 4, 'BE': 4, 'B': 4, 'C': 4, 'N': 4, 'O': 4, 'F': 4, 'NE': 4,
                    'NA': 4, 'MG': 4, 'AL': 4, 'SI': 4, 'CL': 4, 'BR': 4, 'I': 4}
DEFAULT_ORBITAL_WEIGHT = 9   # conservative fallback for unlisted (typically heavier) elements
BOSS_ORBITAL_CAP = 550       # hard limit, from BOSS's own "MAX. NUMBER OF ORBITALS" banner
ORBITAL_TARGET = 480         # default chunking target, with margin below the hard cap to
                              # absorb the trust_radius overlap extension added after chunking

# expected total connection count (heavy + H) for a saturated (sp3-only) atom.
# Backbone atoms outside this are flagged, since that's the one case (a double/
# triple bond sitting ON the backbone) topology-only cutting can't see coming.
EXPECTED_SATURATED_DEGREE = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'SI': 4}
NO_CHECK_ELEMENTS = {'F', 'CL', 'BR', 'I', 'H'}  # always degree 1, nothing to flag


def orbital_weight(element):
    return ORBITAL_WEIGHTS.get(element.upper(), DEFAULT_ORBITAL_WEIGHT)


def estimate_orbitals(atom_indices, elements):
    return sum(orbital_weight(elements[a]) for a in atom_indices)


BOSS_INTEGRAL_CAP = 656175   # hard limit, from BOSS's "MAX. NUMBER OF TWO-ELECTRON INTEGRALS" banner
# NDDO-type semi-empirical methods compute two-electron integrals over
# essentially every pair of atoms in the system (one-center + two-center
# terms), not just bonded neighbors -- so integral count scales with the
# *pairwise product* of per-atom orbital counts, roughly (sum(orbitals))^2,
# not linearly with atom/orbital count the way the 550-orbital cap does. A
# fragment can clear the orbital cap comfortably and still blow this one --
# that's exactly what happened: a 152-atom C/H/F fragment estimated at ~505
# orbitals (well under 550) hit 720145 integrals (over the 656175 cap).
# INTEGRAL_CALIBRATION converts an orbital-weighted pairwise-product sum into
# an integral-count estimate; backed out from that failure (true ratio ~5.68,
# rounded up for margin -- this is an empirical fit to one data point, not a
# first-principles NDDO integral count, so it's deliberately conservative).
INTEGRAL_CALIBRATION = 6.5
INTEGRAL_SAFETY = 0.85       # extra margin below the hard cap for the default target


def estimate_integrals_from_weights(weights):
    """weights: per-atom orbital weights (see orbital_weight/ORBITAL_WEIGHTS).
    Returns the calibrated two-electron-integral estimate for that atom set."""
    s = sum(weights)
    ssq = sum(w * w for w in weights)
    return INTEGRAL_CALIBRATION * (s * s - ssq) / 2.0




def read_species_pdb(pdb_path):
    """Returns {local_idx: (x,y,z)}, {local_idx: element}, in file order (1-indexed)."""
    coords, elements = {}, {}
    idx = 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                idx += 1
                coords[idx] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                elements[idx] = l.split()[-1].upper()
    return coords, elements


def build_adjacency(n_atoms, bonds):
    adj = {i: [] for i in range(1, n_atoms + 1)}
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def spanning_tree(adj, n_atoms):
    """Iterative DFS spanning tree from atom 1. Returns parent, depth, tree_edges
    (as frozensets). Works whether or not the graph has cycles."""
    parent = {1: None}
    depth = {1: 0}
    visited = {1}
    tree_edges = set()
    stack = [(1, iter(adj[1]))]
    while stack:
        node, it = stack[-1]
        advanced = False
        for nb in it:
            if nb not in visited:
                visited.add(nb)
                parent[nb] = node
                depth[nb] = depth[node] + 1
                tree_edges.add(frozenset((node, nb)))
                stack.append((nb, iter(adj[nb])))
                advanced = True
                break
        if not advanced:
            stack.pop()
    if len(visited) != n_atoms:
        raise ValueError(f"molecule graph is disconnected ({len(visited)}/{n_atoms} atoms reached) "
                          f"-- can't be a single species")
    return parent, depth, tree_edges


def find_ring_atoms(bonds, parent, depth, tree_edges):
    """Any non-tree ('chord') edge closes a cycle == tree-path between its
    endpoints (via their lowest common ancestor). Union of all such paths
    = every atom that sits on some ring. Empty for a tree (no cycles)."""
    ring_atoms = set()
    for a, b in bonds:
        if frozenset((a, b)) in tree_edges:
            continue
        u, v = a, b
        pu, pv = [u], [v]
        while depth[u] > depth[v]:
            u = parent[u]; pu.append(u)
        while depth[v] > depth[u]:
            v = parent[v]; pv.append(v)
        while u != v:
            u = parent[u]; pu.append(u)
            v = parent[v]; pv.append(v)
        ring_atoms.update(pu)
        ring_atoms.update(pv)
    return ring_atoms


def bfs_farthest(adj, start):
    """BFS from start. Returns (farthest atom, parent map for path reconstruction)."""
    dist = {start: 0}
    parent = {start: None}
    queue = [start]
    i = 0
    while i < len(queue):
        u = queue[i]; i += 1
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                parent[v] = u
                queue.append(v)
    far = max(dist, key=dist.get)
    return far, parent


def path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    path.reverse()
    return path


def reduced_backbone_graph(adj):
    """Drop every atom with total degree 1 (bonded to exactly one other
    atom, nothing else) -- that's H always, plus any monovalent heavy atom
    (terminal F/Cl/Br/I, nitrile N, carbonyl O, whatever) dynamically, no
    element list needed. A real chain-end atom is never dropped: even the
    last carbon of the chain has degree 4 (3 H + 1 backbone neighbor), so
    it stays a valid node -- only its H's (which are individually degree-1)
    get excluded, same as everywhere else along the chain."""
    nodes = {a for a in adj if len(adj[a]) > 1}
    return {a: [b for b in adj[a] if b in nodes] for a in nodes}


def find_backbone(adj, n_atoms, ring_atoms):
    """Tree-diameter trick, run on the reduced (monovalent-stripped) graph
    so the farthest-atom search can't wander onto a branch tip: BFS from
    any atom to find a farthest atom A, BFS from A to find a farthest atom
    B -- the A-B shortest path is the main chain for a linear/branched
    polymer with no backbone rings. Refuses if either endpoint sits on a
    ring, or if the reduced graph is disconnected (both mean auto-detection
    isn't safe here and the backbone needs to be given manually)."""
    reduced = reduced_backbone_graph(adj)
    if not reduced:
        raise ValueError("no atom has degree > 1 after stripping monovalent atoms -- "
                          "this doesn't look like a chain molecule, backbone needs manual input")
    a0 = next(iter(reduced))
    A, _ = bfs_farthest(reduced, a0)
    B, parent_from_A = bfs_farthest(reduced, A)
    backbone = path_between(parent_from_A, A, B)
    unreached = set(reduced) - set(parent_from_A) - {A}
    if unreached:
        raise ValueError(f"reduced (backbone-candidate) graph is disconnected -- "
                          f"{len(unreached)} atom(s) unreachable from the main chain search, "
                          f"e.g. {sorted(unreached)[:5]} -- needs manual backbone")
    if A in ring_atoms or B in ring_atoms:
        raise ValueError(f"backbone auto-detection landed on a ring atom (endpoint {A if A in ring_atoms else B}) "
                          f"-- this molecule needs --backbone specified manually, auto-detection isn't safe here")
    return backbone


def check_backbone_saturation(backbone, adj, elements, force):
    """Flag any backbone atom whose total connection count doesn't match a
    saturated (all-single-bond) valence for its element. Can't detect a
    backbone double/triple bond directly (no bond-order data exists), this
    is the closest available proxy: an sp2/sp center has fewer total
    neighbors than its saturated valence. Not proof of a problem, but the
    only signal available -- abort unless --force overrides it."""
    flagged = []
    for atom in backbone:
        el = elements[atom]
        expected = EXPECTED_SATURATED_DEGREE.get(el)
        if expected is None or el in NO_CHECK_ELEMENTS:
            continue
        actual = len(adj[atom])
        if actual < expected:
            flagged.append((atom, el, actual, expected))
    if flagged and not force:
        lines = "\n".join(f"  local atom {a} ({el}): {actual} connections, expected {exp} for saturated"
                           for a, el, actual, exp in flagged)
        raise ValueError(
            f"{len(flagged)} backbone atom(s) look unsaturated (fewer neighbors than a single-bonded "
            f"atom would have) -- possible double/triple bond ON the backbone, which topology-only "
            f"cutting cannot see. Won't auto-cut near these without confirmation:\n{lines}\n"
            f"Pass --force to proceed anyway (only safe if you know these are e.g. a valid terminus, "
            f"not a real multiple bond), or fix upstream structure generation.")
    return flagged


def backbone_components(adj, n_atoms, backbone):
    """Remove backbone-path edges, take connected components of what's left.
    Each component contains exactly one backbone atom (since the backbone
    is a single path through a tree, backbone atoms have no other route to
    each other once path edges are cut) -- that backbone atom 'owns' every
    other atom in its component (its branch: side groups, terminal caps,
    whatever -- entire subtrees, never split)."""
    backbone_set = set(backbone)
    backbone_edges = {frozenset((backbone[i], backbone[i + 1])) for i in range(len(backbone) - 1)}

    owner = {}
    for root in backbone:
        if root in owner:
            continue
        seen = {root}
        stack = [root]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v in seen:
                    continue
                if frozenset((u, v)) in backbone_edges:
                    continue  # don't cross into the next backbone atom
                seen.add(v)
                stack.append(v)
        for atom in seen:
            owner[atom] = root
    missing = set(range(1, n_atoms + 1)) - owner.keys()
    if missing:
        raise ValueError(f"{len(missing)} atom(s) not reachable from the backbone without crossing "
                          f"a ring or other unexpected structure: {sorted(missing)[:10]}...")
    return owner, backbone_set


def chunk_backbone(backbone, owner, chunk_target, elements, orbital_target, integral_target):
    """Walk the backbone in order, accumulating each backbone atom's full
    branch-component size (raw atom count, estimated BOSS orbital count, and
    estimated BOSS two-electron-integral count), and close a core segment
    (list of backbone positions, i.e. indices into `backbone`) once adding
    the next atom would exceed chunk_target (atoms), orbital_target
    (estimated orbitals), or integral_target (estimated two-electron
    integrals) -- whichever binds first. The integral estimate is a function
    of Sigma(w) and Sigma(w^2) over the running atom set (see
    estimate_integrals_from_weights), both of which are additive per branch
    component, so it can be tracked incrementally without re-scanning the
    whole running set at every step.
    Each segment always has >=1 atom even if that one atom's branch alone
    exceeds any target (can't split a branch)."""
    comp_size = {}
    comp_orbitals = {}
    comp_orbitals_sq = {}
    for atom, root in owner.items():
        w = orbital_weight(elements[atom])
        comp_size[root] = comp_size.get(root, 0) + 1
        comp_orbitals[root] = comp_orbitals.get(root, 0) + w
        comp_orbitals_sq[root] = comp_orbitals_sq.get(root, 0) + w * w

    segments = []  # list of (start_pos, end_pos) inclusive, positions into `backbone`
    seg_start = 0
    running = 0
    running_orb = 0
    running_orb_sq = 0
    for pos, atom in enumerate(backbone):
        size = comp_size[atom]
        orb = comp_orbitals[atom]
        orb_sq = comp_orbitals_sq[atom]
        integral_est = INTEGRAL_CALIBRATION * ((running_orb + orb) ** 2 - (running_orb_sq + orb_sq)) / 2.0
        if running > 0 and (running + size > chunk_target
                             or running_orb + orb > orbital_target
                             or integral_est > integral_target):
            segments.append((seg_start, pos - 1))
            seg_start = pos
            running = 0
            running_orb = 0
            running_orb_sq = 0
        running += size
        running_orb += orb
        running_orb_sq += orb_sq
    segments.append((seg_start, len(backbone) - 1))
    return segments


def atoms_in_backbone_range(backbone, owner_atoms_by_root, lo, hi):
    """All atoms (backbone + their branches) for backbone positions lo..hi inclusive."""
    out = set()
    for pos in range(lo, hi + 1):
        out.update(owner_atoms_by_root[backbone[pos]])
    return out


def build_owner_atoms_by_root(owner):
    d = {}
    for atom, root in owner.items():
        d.setdefault(root, set()).add(atom)
    return d


def make_fragments(backbone, segments, owner, coords, elements, trust_radius):
    """For each core segment, extend by trust_radius backbone positions on
    each side (clipped at the chain ends) to get the fragment's full
    backbone span; include those extended positions' branch atoms too
    (real, un-capped -- this is the overlap). Cap with H at the fragment's
    own outer boundary, only where that boundary isn't a true chain end.
    Returns list of dicts: species_atoms (set), core_atoms (set, this
    fragment's ownership claim), cap_local_ids (in fragment-local numbering,
    filled in later), backbone_span (lo,hi positions)."""
    owner_atoms_by_root = build_owner_atoms_by_root(owner)
    n_seg = len(segments)
    fragments = []
    for i, (core_lo, core_hi) in enumerate(segments):
        ext_lo = max(0, core_lo - trust_radius)
        ext_hi = min(len(backbone) - 1, core_hi + trust_radius)
        species_atoms = atoms_in_backbone_range(backbone, owner_atoms_by_root, ext_lo, ext_hi)
        core_atoms = atoms_in_backbone_range(backbone, owner_atoms_by_root, core_lo, core_hi)

        cap_specs = []  # (attach_atom, direction_source_atom) -- H added onto attach_atom,
                         # pointing away from where direction_source_atom used to be
        if ext_lo > 0:
            attach = backbone[ext_lo]
            removed_neighbor = backbone[ext_lo - 1]
            cap_specs.append((attach, removed_neighbor))
        if ext_hi < len(backbone) - 1:
            attach = backbone[ext_hi]
            removed_neighbor = backbone[ext_hi + 1]
            cap_specs.append((attach, removed_neighbor))

        fragments.append(dict(
            fragment_id=i + 1,
            species_atoms=species_atoms,
            core_atoms=core_atoms,
            backbone_span=(ext_lo, ext_hi),
            cap_specs=cap_specs,
        ))
    return fragments


def cap_coord(attach_coord, removed_coord, bond_length):
    dx = attach_coord[0] - removed_coord[0]
    dy = attach_coord[1] - removed_coord[1]
    dz = attach_coord[2] - removed_coord[2]
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm < 1e-6:
        raise ValueError("cap direction vector degenerate -- attach and removed-neighbor coords coincide")
    ux, uy, uz = dx / norm, dy / norm, dz / norm
    return (attach_coord[0] + ux * bond_length,
            attach_coord[1] + uy * bond_length,
            attach_coord[2] + uz * bond_length)


def write_fragment_pdb(out_path, species_atoms_ordered, coords, elements, cap_local_coords,
                       local_bonds=(), cap_attach_locals=()):
    """species_atoms_ordered: list of species-local indices in the order they'll
    be numbered in this fragment (real atoms first, caps appended after).
    cap_local_coords: list of (x,y,z) for the cap atoms, same count/order as
    however many caps this fragment has (appended after real atoms).
    local_bonds: pairs of fragment-local serials (real atoms only) to emit as
    CONECT records. cap_attach_locals: fragment-local serial of the real atom
    each cap bonds to, same order as cap_local_coords.

    CONECT records are written explicitly (rather than leaving RDKit to guess
    bonds from atomic distances) because distance-based perception is fragile
    on F-dense fragments like PVDF-HFP: nonbonded F...F/F...C contacts can sit
    close enough to the perceiver's cutoff to be misread as a bond, which
    then fails valence sanitization (F caps out at 1 bond) deep in ligpargen's
    RDKit step. Explicit CONECT removes the guesswork entirely -- BOSS still
    perceives bond ORDER from geometry as before, only connectivity is now
    given directly."""
    with open(out_path, 'w') as f:
        local_i = 0
        for species_idx in species_atoms_ordered:
            local_i += 1
            x, y, z = coords[species_idx]
            el = elements[species_idx]
            name = (el[:1] + el[1:].lower()) if len(el) > 1 else el
            # resName must be exactly 3 chars (cols 18-20) + blank (col 21) so
            # chain lands at col 22 and resSeq at cols 23-26, per PDB spec --
            # RDKit's PDB reader parses those columns strictly.
            f.write(f"ATOM  {local_i:5d}  {name:<3s} MOL A{1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
        for x, y, z in cap_local_coords:
            local_i += 1
            f.write(f"ATOM  {local_i:5d}  H   MOL A{1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           H\n")

        n_real = len(species_atoms_ordered)
        cap_local_ids = range(n_real + 1, n_real + len(cap_local_coords) + 1)
        all_bonds = set()
        for a, b in local_bonds:
            all_bonds.add((a, b))
            all_bonds.add((b, a))
        for cap_id, attach_local in zip(cap_local_ids, cap_attach_locals):
            all_bonds.add((attach_local, cap_id))
            all_bonds.add((cap_id, attach_local))

        neighbors = {}
        for a, b in all_bonds:
            neighbors.setdefault(a, []).append(b)
        for atom in sorted(neighbors):
            nbrs = sorted(neighbors[atom])
            # PDB CONECT caps at 4 partners per line; chunk if ever exceeded
            for k in range(0, len(nbrs), 4):
                chunk = nbrs[k:k + 4]
                f.write("CONECT" + f"{atom:5d}" + "".join(f"{n:5d}" for n in chunk) + "\n")

        f.write("END\n")
    return local_i


def run(manifest_path, species_dir, species_id, outdir, chunk_target, trust_radius, force,
        orbital_target=ORBITAL_TARGET):
    with open(manifest_path) as f:
        manifest = json.load(f)
    sp = manifest['species'][str(species_id)]
    n_atoms = sp['n_atoms']
    bonds = [tuple(b) for b in sp['representative_bonds']]

    pdb_path = os.path.join(species_dir, sp['pdb'])
    coords, elements = read_species_pdb(pdb_path)
    if len(coords) != n_atoms:
        raise ValueError(f"{pdb_path} has {len(coords)} atoms, manifest says {n_atoms}")

    adj = build_adjacency(n_atoms, bonds)
    parent, depth, tree_edges = spanning_tree(adj, n_atoms)
    ring_atoms = find_ring_atoms(bonds, parent, depth, tree_edges)
    backbone = find_backbone(adj, n_atoms, ring_atoms)

    flagged = check_backbone_saturation(backbone, adj, elements, force)
    if flagged:
        print(f"WARNING: proceeding past {len(flagged)} flagged backbone atom(s) due to --force")

    owner, backbone_set = backbone_components(adj, n_atoms, backbone)

    # Chunk, then verify every finished fragment (core + trust_radius overlap
    # + H caps) actually clears BOSS's real orbital AND two-electron-integral
    # caps. The overlap extension happens *after* chunk_backbone's budgeting,
    # so the targets that go into chunk_backbone are necessarily a guess at
    # how much headroom to leave for it -- if a fragment still comes in over
    # either cap, shrink orbital_target (the integral target is always
    # re-derived from it) in proportion to the worse overshoot and re-chunk,
    # rather than making the user hand-tune targets by trial and error.
    current_orbital_target = orbital_target
    MAX_ATTEMPTS = 8
    MIN_ORBITAL_TARGET = 50
    SAFETY = 0.92  # extra margin below the proportional estimate each retry

    def derive_integral_target(orb_target):
        # integral count ~ orbital_target^2 (pairwise), so derive a matching
        # integral budget from whatever orbital target is currently in play,
        # capped at a safety fraction of the real hard cap either way.
        modeled = INTEGRAL_CALIBRATION * (orb_target ** 2) / 2.0
        return min(modeled, BOSS_INTEGRAL_CAP * INTEGRAL_SAFETY)

    current_integral_target = derive_integral_target(current_orbital_target)

    # Post-check margin: the integral estimate is calibrated from a single
    # observed failure, not derived first-principles, so treat anything
    # within POST_CHECK_MARGIN of either hard cap as a violation too --
    # otherwise convergence can land a fragment at >99% of the integral cap,
    # which is too close given the model's uncertainty.
    POST_CHECK_MARGIN = 0.90

    for attempt in range(1, MAX_ATTEMPTS + 1):
        segments = chunk_backbone(backbone, owner, chunk_target, elements,
                                   current_orbital_target, current_integral_target)

        if len(segments) == 1:
            print(f"species {species_id}: {n_atoms} atoms fits in one fragment "
                  f"(target {chunk_target} atoms / {current_orbital_target} orbitals / "
                  f"{current_integral_target:.0f} integrals) -- no split needed, "
                  f"feed species_{species_id}.pdb to ligpargen directly.")
            return

        fragments = make_fragments(backbone, segments, owner, coords, elements, trust_radius)

        # For each fragment, find its worse-of-the-two overshoot ratio (1.0 = at cap).
        worst = None  # (fragment, total_orbitals, total_integrals, worst_ratio, which)
        for frag in fragments:
            n_caps = len(frag['cap_specs'])
            weights = [orbital_weight(elements[a]) for a in frag['species_atoms']]
            weights += [ORBITAL_WEIGHTS['H']] * n_caps
            total_orbitals = sum(weights)
            total_integrals = estimate_integrals_from_weights(weights)
            orb_ratio = total_orbitals / BOSS_ORBITAL_CAP
            int_ratio = total_integrals / BOSS_INTEGRAL_CAP
            which = 'integral' if int_ratio >= orb_ratio else 'orbital'
            ratio = max(orb_ratio, int_ratio)
            if ratio > POST_CHECK_MARGIN and (worst is None or ratio > worst[3]):
                worst = (frag, total_orbitals, total_integrals, ratio, which)

        if worst is None:
            if attempt > 1:
                print(f"species {species_id}: converged on orbital target "
                      f"{current_orbital_target} / integral target "
                      f"{current_integral_target:.0f} after {attempt} attempt(s)")
            break

        frag, total_orbitals, total_integrals, ratio, which = worst
        if attempt == MAX_ATTEMPTS or current_orbital_target <= MIN_ORBITAL_TARGET:
            raise ValueError(
                f"fragment {frag['fragment_id']} of species {species_id}: estimated "
                f"{total_orbitals} orbitals (cap {BOSS_ORBITAL_CAP}) / "
                f"{total_integrals:.0f} integrals (cap {BOSS_INTEGRAL_CAP}) still exceeds "
                f"BOSS's hard cap after {attempt} auto-retry attempt(s) (last orbital target "
                f"{current_orbital_target}). Try a smaller --max-atoms (currently "
                f"{chunk_target}) -- a single backbone atom's branch may itself be too large.")

        if which == 'orbital':
            shrink = (BOSS_ORBITAL_CAP / total_orbitals) * SAFETY
        else:
            # integral count ~ orbital_target^2, so correct with sqrt(ratio)
            shrink = math.sqrt(BOSS_INTEGRAL_CAP / total_integrals) * SAFETY
        new_orbital_target = max(int(current_orbital_target * shrink), MIN_ORBITAL_TARGET)
        new_integral_target = derive_integral_target(new_orbital_target)
        print(f"species {species_id}: fragment {frag['fragment_id']} estimated at "
              f"{total_orbitals} orbitals (cap {BOSS_ORBITAL_CAP}) / {total_integrals:.0f} "
              f"integrals (cap {BOSS_INTEGRAL_CAP}) -- {which} cap binding -- retrying with "
              f"orbital target {new_orbital_target} (was {current_orbital_target})")
        current_orbital_target = new_orbital_target
        current_integral_target = new_integral_target

    os.makedirs(outdir, exist_ok=True)
    frag_manifest = {"species_id": species_id, "n_atoms": n_atoms,
                      "trust_radius": trust_radius, "fragments": {}}

    for frag in fragments:
        species_atoms_ordered = sorted(frag['species_atoms'])  # stable, deterministic local numbering
        cap_specs = frag['cap_specs']
        cap_coords = []
        for attach, removed in cap_specs:
            cap_coords.append(cap_coord(coords[attach], coords[removed], CAP_BOND_LENGTH))

        # species-local index -> this fragment's local serial (real atoms only)
        species_to_local = {si: li for li, si in enumerate(species_atoms_ordered, start=1)}
        species_atoms_set = frag['species_atoms']
        local_bonds = [(species_to_local[a], species_to_local[b])
                       for a, b in bonds if a in species_atoms_set and b in species_atoms_set]
        cap_attach_locals = [species_to_local[attach] for attach, _removed in cap_specs]

        fid = frag['fragment_id']
        out_pdb = os.path.join(outdir, f"fragment_{species_id}_{fid}.pdb")
        n_written = write_fragment_pdb(out_pdb, species_atoms_ordered, coords, elements, cap_coords,
                                        local_bonds=local_bonds, cap_attach_locals=cap_attach_locals)

        n_real = len(species_atoms_ordered)
        cap_local_ids = list(range(n_real + 1, n_real + len(cap_coords) + 1))
        # fragment-local id -> species-local id, for real atoms only (caps have no species id)
        local_to_species = {li: si for li, si in enumerate(species_atoms_ordered, start=1)}

        frag_manifest["fragments"][fid] = {
            "pdb": os.path.basename(out_pdb),
            "n_atoms": n_written,
            "n_real": n_real,
            "n_caps": len(cap_coords),
            "cap_local_ids": cap_local_ids,
            "local_to_species": local_to_species,          # fragment-local -> species-local (real atoms)
            "core_species_atoms": sorted(frag['core_atoms']),  # this fragment is authoritative for these
            "backbone_span": frag['backbone_span'],
        }
        print(f"  fragment {fid}: {n_real} real + {len(cap_coords)} cap = {n_written} atoms, "
              f"owns {len(frag['core_atoms'])} atoms -- {out_pdb}")

    with open(os.path.join(outdir, f"fragment_manifest_{species_id}.json"), 'w') as f:
        json.dump(frag_manifest, f, indent=2)

    total_owned = sum(len(fr['core_atoms']) for fr in fragments)
    if total_owned != n_atoms:
        raise ValueError(f"ownership partition bug: {total_owned} atoms owned, species has {n_atoms}")

    print(f"\nspecies {species_id}: {n_atoms} atoms -> {len(fragments)} fragments, "
          f"ownership partition verified complete. Run ligpargen on each fragment_{species_id}_*.pdb, "
          f"then run molpargen stitch.")


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    force = '--force' in args
    if force:
        args.remove('--force')
    chunk_target = CHUNK_TARGET
    if '--max-atoms' in args:
        i = args.index('--max-atoms')
        chunk_target = int(args[i + 1])
        del args[i:i + 2]
    orbital_target = ORBITAL_TARGET
    if '--max-orbitals' in args:
        i = args.index('--max-orbitals')
        orbital_target = int(args[i + 1])
        del args[i:i + 2]
    if len(args) != 4:
        print("Usage: molpargen split manifest.json species_dir/ species_id outdir/ "
              "[--max-atoms N] [--max-orbitals N] [--force]")
        return 1
    manifest_path, species_dir, species_id, outdir = args
    run(manifest_path, species_dir, int(species_id), outdir, chunk_target, TRUST_RADIUS, force,
        orbital_target=orbital_target)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
