#!/usr/bin/env python
# Boltz Generalized Covalent Residue Modification
#
# Modifies the Boltz Chemical Component Dictionary cache with a new entry
# corresponding to an arbitrary covalent modification of a residue.
#
# original code : https://github.com/benf549/boltz-generalized-covalent-modification/

import os
import pickle
from pathlib import Path
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem

amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# THIS IS ABSOLUTELY NECESSARY TO MAKE SURE THAT ALL PROPERTIES ARE PICKLED OR ELSE BOLTZ WILL CRASH.
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

def main():
    parser = argparse.ArgumentParser(description='Add a custom covalent residue modification to Boltz CCD cache')
    parser.add_argument('smiles', type=str, help='SMILES string of the covalent modification')
    parser.add_argument('reference_residue', type=str, help='Reference residue to modify', choices=amino_acids)
    parser.add_argument('output_ccd_code', type=str, help='New CCD code for the covalent ligand (first 5 letters are written to the CIF file as the resname)')
    parser.add_argument('--force', '-f', action='store_true', help='Force overwrite of existing CCD code')
    parser.add_argument('--cache_dir', type=str, default=os.path.expanduser("~/.boltz"), help='Boltz cache directory')
    parser.add_argument('--export_sdf', action='store_true', help='Export the modified residue as an SDF file')
    args = parser.parse_args()

    ncaa_smiles_str = args.smiles
    reference_residue_to_modify = args.reference_residue
    output_ccd_code = args.output_ccd_code
    cache_dir = Path(args.cache_dir)

    # Load the Boltz Cache Chemical Component Dictionary
    ccd_path = cache_dir / 'ccd.pkl'
    with ccd_path.open("rb") as file:
        ccd = pickle.load(file)  # noqa: S301
    if ccd.get(output_ccd_code) is not None and not args.force:
        raise ValueError(f"CCD code {output_ccd_code} already exists in the CCD cache. Use --force to overwrite.")
    print(f"Loaded Boltz CCD cache from {ccd_path}")

    # Parse the SMILES string into a molecule
    smiles_mol = Chem.MolFromSmiles(ncaa_smiles_str)
    if smiles_mol is None:
        raise ValueError("Invalid SMILES string provided")

    # Load reference residue (e.g., Cysteine) from the CCD and remove hydrogens
    if reference_residue_to_modify not in ccd:
        raise ValueError(f"Reference residue {reference_residue_to_modify} not found in CCD")

    reference_mol = ccd[reference_residue_to_modify]
    reference_mol = AllChem.RemoveHs(reference_mol)

    # Search for the reference substructure in the NCAA molecule
    has_match = smiles_mol.HasSubstructMatch(reference_mol)
    if not has_match:
        raise ValueError(f"NCAA molecule does not contain {reference_residue_to_modify} substructure")

    match_indices = smiles_mol.GetSubstructMatch(reference_mol)
    substruct_to_match = {i.GetProp('name'): match_indices[idx] for idx, i in enumerate(reference_mol.GetAtoms())}

    # Construct mapping of reference atom name to atom index in the NCAA molecule
    idx_to_name = {j: i for i, j in substruct_to_match.items()}
    print(f"Found mapping between reference residue and NCAA molecule: {idx_to_name}")

    # Add the metadata boltz expects to the atoms
    for idx, atom in enumerate(smiles_mol.GetAtoms()):
        default_name = f'{atom.GetSymbol()}{str(atom.GetIdx())}'

        # If the index is a canonical reference atom, use the reference atom name
        name = idx_to_name.get(idx, default_name)

        # Set atom properties
        atom.SetProp('name', name)
        atom.SetProp('alt_name', name)
        is_leaving = False
        if name == 'OXT':
            is_leaving = True
        atom.SetBoolProp('leaving_atom', is_leaving)

    # Reorder atoms to canonical ordering (N, Ca, C, O, CB, ...)
    # Map atom name to atom index in the NCAA molecule
    curr_atom_order = {atom.GetProp('name'): idx for idx, atom in enumerate(smiles_mol.GetAtoms()) if atom.GetSymbol() != 'H'}

    # Map atom name to atom index in the reference molecule
    target_atom_order = {}
    for atom in reference_mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        target_atom_order[atom.GetProp('name')] = atom.GetIdx()

    # Add atoms not in target_atom_order that are in curr_atom_order
    remapped_atom_order = {}
    offset_idx = len(target_atom_order)
    for atom in curr_atom_order:
        if atom in target_atom_order:
            remapped_atom_order[atom] = target_atom_order[atom]
        else:
            remapped_atom_order[atom] = offset_idx
            offset_idx += 1

    # Remove hydrogens and reorder atoms according to the order in the reference residue
    trim = AllChem.RemoveHs(smiles_mol)
    remap_order = {x.GetProp('name'): (remapped_atom_order[x.GetProp('name')], x.GetIdx()) for x in trim.GetAtoms()}
    remap_idx_list = [x[1] for x in sorted(remap_order.values())]
    trim_reordered = Chem.RenumberAtoms(trim, remap_idx_list)

    # Generate a simple conformer
    trim_reordered = AllChem.AddHs(trim_reordered)
    embed_success = AllChem.EmbedMolecule(trim_reordered)
    if embed_success == -1:
        print("Warning: Could not embed molecule, trying with random coordinates...")
        AllChem.EmbedMolecule(trim_reordered, useRandomCoords=True)

    optimize_success = AllChem.UFFOptimizeMolecule(trim_reordered)
    if optimize_success == -1:
        print("Warning: UFF optimization failed, structure may not be ideal")

    # Set conformer properties
    for c in trim_reordered.GetConformers():
        c.SetProp('name', 'Ideal')

    # Verify atom ordering
    print("\nVerifying atom ordering:")
    for atom in trim_reordered.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        props = {k: v for k, v in atom.GetPropsAsDict().items() if k not in ['__computedProps']}
        print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()} - {props}")

    # Save the conformer to the CCD cache
    ccd[output_ccd_code] = trim_reordered

    with ccd_path.open("wb") as file:
        pickle.dump(ccd, file)
    if args.export_sdf:
        with Chem.SDWriter(f"{output_ccd_code}.sdf") as writer:
            writer.write(trim_reordered)

    print(f"\nSuccessfully added {output_ccd_code} to the Boltz CCD cache")

if __name__ == "__main__":
    main()
