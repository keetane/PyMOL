# from rdkit import Chem
# from rdkit.Chem import AllChem
# import pubchempy as pcp
# from pymol import cmd

# def smiles(name: str, smile: str):
#     mol = Chem.MolFromSmiles(smile)
#     mol = Chem.AddHs(mol)
#     AllChem.EmbedMolecule(mol)
#     AllChem.UFFOptimizeMolecule(mol)
#     pdb_block = Chem.MolToPDBBlock(mol)
#     cmd.read_pdbstr(pdb_block, name)

# def pc(name: str):
#     sdf = pcp.get_sdf(name, 'name')
#     mol = Chem.MolFromMolBlock(sdf)
#     mol = Chem.AddHs(mol)
#     AllChem.EmbedMolecule(mol)
#     AllChem.UFFOptimizeMolecule(mol)
#     pdb_block = Chem.MolToPDBBlock(mol)
#     cmd.read_pdbstr(pdb_block, name)

# cmd.extend('smiles', smiles)
# cmd.extend('pc', pc)

from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd
from openbabel import openbabel as ob

obc = ob.OBConversion()
obc.SetInAndOutFormats('smi', 'smi')

# https://iwatobipen.wordpress.com/2024/05/14/add-hydrogen-with-user-defined-ph-from-python-openbabel-cheminformatics/
def get_smi_with_pH(smi, pH=7.4):
    obmol = ob.OBMol()
    obc.ReadString(obmol, smi)
    obmol.CorrectForPH(pH)
    return obc.WriteString(obmol)

def smiles(name: str, smile: str, pH=7.4):
    smile = get_smi_with_pH(smile, pH)
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    cmd.read_pdbstr(pdb_block, name)

def pc(name: str, pH=7.4):
    try:
        sdf = pcp.get_sdf(name, 'name')
        if sdf is None:
            print(f"Error: Could not retrieve SDF for {name}")
            return
        mol = Chem.MolFromMolBlock(sdf)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        pdb_block = Chem.MolToPDBBlock(mol)
        cmd.read_pdbstr(pdb_block, name)
    except Exception as e:
        print(f"An error occurred: {e}")

cmd.extend('smiles', smiles)
cmd.extend('pc', pc)

