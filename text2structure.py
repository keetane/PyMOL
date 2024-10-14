from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd

def smiles(name: str, smile: str):
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    cmd.read_pdbstr(pdb_block, name)

def pc(name: str):
    sdf = pcp.get_sdf(name, 'name')
    mol = Chem.MolFromMolBlock(sdf)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    cmd.read_pdbstr(pdb_block, name)

cmd.extend('smiles', smiles)
cmd.extend('pc', pc)