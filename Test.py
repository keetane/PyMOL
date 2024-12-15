from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from pymol import cmd

# https://iwatobipen.wordpress.com/2024/05/14/add-hydrogen-with-user-defined-ph-from-python-openbabel-cheminformatics/
obc = openbabel.OBConversion()
obc.SetInAndOutFormats('smi', 'smi')

# adjusting the pH of SMILES molecules
def get_smi_with_pH(smi, pH:float=7.4):
    obmol = openbabel.OBMol()
    obc.ReadString(obmol, smi)
    obmol.CorrectForPH(float(pH))
    return obc.WriteString(obmol)

def smiles(arg_string):
    # 引数をスペースで分割
    args = arg_string.split()
    
    # 最低限の引数チェック
    if len(args) < 2:
        print("エラー: コマンドには少なくとも名前とSMILES文字列が必要です。")
        return
    
    name = args[0]
    smile = args[1]
    
    # 必要に応じてpH引数を取得
    pH = None
    if len(args) > 2:
        try:
            pH = float(args[2])
        except ValueError:
            print(f"エラー: pH値 '{args[2]}' は有効な数値ではありません。")
            return
    
    if pH is not None:
        smile = get_smi_with_pH(smile, pH)
    
    mol = Chem.MolFromSmiles(smile)
    
    if mol is None:
        print(f"エラー: SMILES文字列 '{smile}' から分子を生成できませんでした。")
        return
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    
    if pdb_block is None:
        print("エラー: 分子のPDBブロックを生成できませんでした。")
        return
    
    cmd.read_pdbstr(pdb_block, name)
    cmd.hide(f'({name} and hydro and (elem C extend 1))')
# コマンドをPyMOLに拡張
cmd.extend('smiles', smiles)

from rdkit import Chem
from pymol import cmd
import os

def printsmiles(selection="sele"):
    # 一時SDFファイル名
    tmp = "./tmp_molecule.sdf"
    # PyMOLから選択された分子をSDF形式で保存
    cmd.save(tmp, selection)

    # rdkitで分子をロード、正規化
    suppl = Chem.SDMolSupplier(tmp)
    mol = next(suppl) if suppl and len(suppl) > 0 else None
    if mol is None:
        print("Error: Failed to load the molecule from the SDF file.")
        return
    smiles = Chem.MolToSmiles(mol)

    # SMILES形式に変換し、オブジェクト名を追加して出力する
    object_name = cmd.get_object_list(selection)[0]
    print(f"SMILES:\n{object_name}\t{smiles}")

    # 一時ファイルを削除
    os.remove(tmp)

# コマンドをPyMOLに拡張
cmd.extend("printsmiles", printsmiles)








#     # Open Babelを使用してPDBファイルを読み込み
#     obConversion = openbabel.OBConversion()
#     obConversion.SetInAndOutFormats("pdb", "smi")
#     mol = openbabel.OBMol()
#     obConversion.ReadFile(mol, tmp_pdb)
    
#     # SMILES形式に変換し、オブジェクト名を追加して出力する
#     smiles = obConversion.WriteString(mol).strip()
#     object_name = cmd.get_object_list(selection)[0]
#     smiles = smiles.strip('./tmp_molecule.pdb').replace(".", "\n")
#     print(f"SMILES:\n{object_name}\t{smiles}")
# # コマンドをPyMOLに拡張
# cmd.extend("printsmiles", printsmiles)

def build(selection="sele"):
    # 一時PDBファイル名
    tmp_pdb = f"./{selection}.pdb"
    
    # PyMOLから選択された分子をPDB形式で保存
    cmd.save(tmp_pdb, selection)

    # Open Babelを使用してPDBファイルを読み込み
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, tmp_pdb)
    
    # SMILES形式に変換してオブジェクト名を取得
    smiles = obConversion.WriteString(mol).strip()
    object_name = cmd.get_object_list(selection)[0]
    print(f"SMILES:\n{object_name}\t{smiles}")

    # SMILESから立体構造を生成
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print(f"エラー: SMILES文字列 '{smiles}' から分子を生成できませんでした。")
        return
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    pdb_block = Chem.MolToPDBBlock(mol)
    
    if pdb_block is None:
        print("エラー: 分子のPDBブロックを生成できませんでした。")
        return
    
    # PyMOLに立体構造を読み込み
    cmd.read_pdbstr(pdb_block, object_name)
    cmd.hide(f'({object_name} and hydro and (elem C extend 1))')
    
    print(f"{object_name}の立体構造がPyMOLにロードされました。")

# コマンドをPyMOLに拡張
cmd.extend('build', build)
