from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd
from openbabel import openbabel as ob

obc = ob.OBConversion()
obc.SetInAndOutFormats('smi', 'smi')

# https://iwatobipen.wordpress.com/2024/05/14/add-hydrogen-with-user-defined-ph-from-python-openbabel-cheminformatics/
def get_smi_with_pH(smi, pH:float=7.4):
    obmol = ob.OBMol()
    obc.ReadString(obmol, smi)
    obmol.CorrectForPH(float(pH))
    return obc.WriteString(obmol)

def smiles(name: str, smile:str=None, pH:float=None):
    if smile is None:
        smile=name
    if not pH is None:
        smile = get_smi_with_pH(smile, pH)
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    cmd.read_pdbstr(pdb_block, name)
    cmd.hide(f'({smile} and hydro and (elem C extend 1))')
cmd.extend('smiles', smiles)

def pc(name: str, pH:float=None):
    try:
        sdf = pcp.get_sdf(name, 'name')
        if sdf is None:
            print(f"Error: Could not retrieve SDF for {name}")
            return
        mol = Chem.MolFromMolBlock(sdf)
        smile = Chem.MolToSmiles(mol)
        if not pH is None:
            smile = get_smi_with_pH(smile, pH)
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        pdb_block = Chem.MolToPDBBlock(mol)
        cmd.read_pdbstr(pdb_block, name)
    except Exception as e:
        print(f"An error occurred: {e}")
    cmd.hide(f'({name} and hydro and (elem C extend 1))')
cmd.extend('pc', pc)