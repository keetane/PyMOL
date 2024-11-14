from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd
from openbabel import openbabel as ob

# https://iwatobipen.wordpress.com/2024/05/14/add-hydrogen-with-user-defined-ph-from-python-openbabel-cheminformatics/
obc = ob.OBConversion()
obc.SetInAndOutFormats('smi', 'smi')

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
    cmd.hide(f'({name} and hydro and (elem C extend 1))')
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


def pubchem2smi(compound_name: str):
    try:
        # PubChemから化合物情報を取得
        compound = pcp.get_compounds(compound_name, 'name')
        if not compound:
            print(f"No compounds found for {compound_name}")
            return

        # SMILESを取得
        smiles = compound[0].canonical_smiles
        if not smiles:
            print(f"SMILES not found for {compound_name}")
            return

        # ファイルに保存
        file_name = f"{compound_name.replace(' ', '_')}.smi"
        with open(file_name, 'w') as f:
            f.write(smiles)
        print(f"SMILES for {compound_name} saved in {file_name} as")
        print(f"{smiles}")

    except Exception as e:
        print(f"An error occurred: {e}")
cmd.extend('pubchem2smi', pubchem2smi)


def sele2smi(selection="sele"):
    # PDB形式で一時ファイルに保存
    tmp_pdb = "/tmp/tmp_molecule.pdb"
    cmd.save(tmp_pdb, selection)
    
    # RDKitでPDBファイルを読み込み
    molecule = Chem.MolFromPDBFile(tmp_pdb)
    
    if molecule is not None:
        # SMILES形式に変換して出力
        smiles = Chem.MolToSmiles(molecule)
        print("SMILES:", smiles)
    else:
        print("分子を読み込めませんでした。")
# コマンドをPyMOLに拡張
cmd.extend("sele2smi", sele2smi)
