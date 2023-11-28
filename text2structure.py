from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import pubchempy as pcp

def smiles(x: str, y: str):
    mol = Chem.MolFromSmiles(y)
    Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    # Draw.MolToImage(mol)
    writer = Chem.SDWriter(x + '.sdf')
    writer.write(mol)
    writer.close()
    cmd.load(x + '.sdf')

def pc(x: str):
    sdf = pcp.get_sdf(x, 'name')
    mol = Chem.MolFromMolBlock(sdf)
    Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(x + '.sdf')
    writer.write(mol)
    writer.close()
    cmd.load(x + '.sdf')

cmd.extend('smiles', smiles)
cmd.extend('pc', pc)
