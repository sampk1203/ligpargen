"""
Shared molecule/species-detection graph algorithms. Used identically by the
PSF+PDB path (extract_unique.py) and the mol2 path (extract_unique_mol2.py) --
species detection only needs (natom, bonds, elements); it doesn't care
whether those came from a PSF+PDB pair or a mol2 file.
"""
from collections import Counter


def connected_components(natom, bonds):
    adj = {i: [] for i in range(1, natom + 1)}
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)

    visited = {}
    comp_id = 0
    comp_atoms = {}
    for start in range(1, natom + 1):
        if start in visited:
            continue
        comp_id += 1
        stack = [start]
        visited[start] = comp_id
        comp_atoms[comp_id] = []
        while stack:
            cur = stack.pop()
            comp_atoms[comp_id].append(cur)
            for nb in adj[cur]:
                if nb not in visited:
                    visited[nb] = comp_id
                    stack.append(nb)
    for c in comp_atoms:
        comp_atoms[c].sort()  # preserve serial order within each molecule
    return visited, comp_atoms


def component_bonds(bonds, comp_atoms_set):
    return [(a, b) for a, b in bonds if a in comp_atoms_set and b in comp_atoms_set]


def species_signature(atom_serials, elements, bonds_in_comp):
    """Weisfeiler-Leman canonical signature -- element-labeled connectivity,
    not just atom/element/bond COUNTS. Two molecules only get the same
    signature if they're actually isomorphic (same topology), not merely
    matched on aggregate composition. Composition-count matching alone
    silently merges e.g. two random-copolymer chains that share monomer
    counts but differ in sequence -- ligpargen then parametrizes only the
    representative, and every other "same species" instance gets that
    representative's bond graph stamped on positionally-wrong atoms.
    """
    adj = {s: [] for s in atom_serials}
    for a, b in bonds_in_comp:
        adj[a].append(b)
        adj[b].append(a)

    color = {s: elements.get(s, '?') for s in atom_serials}

    rounds = min(len(atom_serials), 8)  # plenty for organic-sized molecules
    for _ in range(rounds):
        new_color = {}
        for s in atom_serials:
            neighbor_colors = tuple(sorted(color[nb] for nb in adj[s]))
            new_color[s] = hash((color[s], neighbor_colors))
        color = new_color

    canonical = tuple(sorted(Counter(color.values()).items()))
    return (len(atom_serials), canonical)


def group_into_species(natom, bonds, elements):
    """Given a whole system's (natom, bonds, elements), split into connected
    components and group them into unique species by graph isomorphism.

    Returns (comp_atoms, comp_to_species, species_list) where:
      comp_atoms      : {comp_id: [sorted 1-based atom serials]}
      comp_to_species : {comp_id: species_id}
      species_list    : [(species_id, representative_comp_id, signature), ...]
    """
    visited, comp_atoms = connected_components(natom, bonds)

    sig_to_species = {}
    comp_to_species = {}
    species_list = []
    for comp_id, serials in comp_atoms.items():
        serial_set = set(serials)
        comp_b = component_bonds(bonds, serial_set)
        sig = species_signature(serials, elements, comp_b)
        if sig not in sig_to_species:
            species_id = len(species_list) + 1
            sig_to_species[sig] = species_id
            species_list.append((species_id, comp_id, sig))
        comp_to_species[comp_id] = sig_to_species[sig]

    return comp_atoms, comp_to_species, species_list


ION_ATOM_LIMIT = 30


def estimate_net_charge(species_id, n_atoms, charge_sum, verbose=True):
    """Same heuristic extract_unique.py has always used: small molecules
    with a near-integer nonzero charge sum are treated as genuine ions;
    everything else is assumed neutral (override in manifest if wrong).
    """
    nearest_int = round(charge_sum)
    deviation = abs(charge_sum - nearest_int)
    if n_atoms <= ION_ATOM_LIMIT and deviation < 0.1 and abs(nearest_int) >= 1:
        expected_net_charge = nearest_int
        if verbose:
            print(f"  species {species_id}: charge sum = {charge_sum:.4f} "
                  f"-> expected_net_charge = {expected_net_charge} (ion)")
    else:
        expected_net_charge = 0
        if verbose and abs(charge_sum) > 0.01:
            print(f"  species {species_id}: charge sum = {charge_sum:.4f} "
                  f"(polymer/neutral) -> expected_net_charge = 0. "
                  f"Override in manifest if wrong.")
    return expected_net_charge
