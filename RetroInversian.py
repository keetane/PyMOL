from rdkit import Chem
from rdkit.Chem import AllChem
from pymol import cmd

# aminoacid format
laa = 'N[C@@]([H])({X})C(=O)'
daa = 'N[C@@]({X})([H])C(=O)'

# residue dict
amino_acids = {
    'G': '[H]',
    'A': 'C',
    'V': 'C(C)C',
    'L': 'CC(C)C',
    'I': 'C(C)CC',
    'M': 'CCSC',
    'F': 'CC1=CC=CC=C1',
    'W': 'CC1=CNC2=C1C=CC=C2',
    'S': 'CO',
    'C': 'CS',
    'Y': 'CC1=CC=C(C=C1)O',
    'N': 'C(C(=O)N)',
    'Q': 'C(CC(=O)N)',
    'D': 'C(C(=O)O)',
    'E': 'C(CC(=O)O)',
    'H': 'C(C1=CNC=N1)',
    'K': 'CCCCN',
    'R': 'CCCNC(=N)N',
}

# insert variable to format
laaa = {k: laa.format(X=v) for k, v in amino_acids.items()}
daaa = {k: daa.format(X=v) for k, v in amino_acids.items()}

# add prolines and threonines
laaa['P'] = 'N1[C@@]([H])(CCC1)C(=O)'
daaa['P'] = 'N1[C@@](CCC1)([H])C(=O)'
laaa['T'] = 'N[C@@]([H])([C@]([H])(O)C)C(=O)'
daaa['T'] = 'N[C@@]([C@]([H])(O)C)([H])C(=O)'

# PepType
PepType = {
    '1': 'Liner',
    '2': 'Head to Tail',
    '3': 'Cystein Cyclization',
    '4': 'Symple Disulfide'
}

def generate_peptide_smiles(fasta, mode):
    l_sec = ''.join([laaa[aa] for aa in fasta])
    d_sec = ''.join([daaa[aa] for aa in fasta[::-1]])

    # liner
    l_liner = l_sec + 'O'
    d_liner = d_sec + 'O'

    # Head to Tail
    l_ht = 'N9' + l_sec[1:-4] + '9(=O)'
    d_ht = 'N9' + d_sec[1:-4] + '9(=O)'

    # Cystein Cyclization
    l_cc = 'C9C(=O)' + l_sec + 'N[C@@]([H])(CS9)C(=O)N'
    d_cc = 'C9C(=O)' + d_sec + 'N[C@@](CS9)([H])C(=O)N'

    # Simple Disulfide
    l_ds = 'N[C@@]([H])(CS9)C(=O)' + l_sec + 'N[C@@]([H])(CS9)C(=O)O'
    d_ds = 'N([H])[C@@](CS9)C(=O)' + d_sec + 'N([H])[C@@](CS9)C(=O)O'

    if mode == '2':
        l_smiles = l_ht
        d_smiles = d_ht
    elif mode == '3':
        l_smiles = l_cc
        d_smiles = d_cc
    elif mode == '4':
        l_smiles = l_ds
        d_smiles = d_ds
    else:
        l_smiles = l_liner
        d_smiles = d_liner

    return l_smiles, d_smiles

def display_peptide_structure(fasta, mode=1):
    l_smiles, d_smiles = generate_peptide_smiles(fasta, mode)

    # generate molecules
    l_mol = Chem.MolFromSmiles(l_smiles)
    d_mol = Chem.MolFromSmiles(d_smiles)

    # Check if molecules are valid
    if not l_mol:
        raise ValueError("Invalid L-peptide SMILES")
    if not d_mol:
        raise ValueError("Invalid D-peptide SMILES")

    # 3D conformer generation
    if AllChem.EmbedMolecule(l_mol, randomSeed=42) == -1:
        raise ValueError("Embedding failed for L-peptide")
    if AllChem.EmbedMolecule(d_mol, randomSeed=42) == -1:
        # Try a different random seed
        if AllChem.EmbedMolecule(d_mol, randomSeed=43) == -1:
            raise ValueError("Embedding failed for D-peptide")
    
    AllChem.UFFOptimizeMolecule(l_mol)
    AllChem.UFFOptimizeMolecule(d_mol)

    # Convert to PDB
    l_pdb = Chem.MolToPDBBlock(l_mol)
    d_pdb = Chem.MolToPDBBlock(d_mol)

    # Load into PyMOL
    cmd.read_pdbstr(l_pdb, 'l_peptide')
    cmd.read_pdbstr(d_pdb, 'd_peptide')

    cmd.show('sticks', 'l_peptide')
    cmd.show('sticks', 'd_peptide')
    cmd.orient()

cmd.extend('retroinversian', display_peptide_structure)
