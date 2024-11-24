from pymol import cmd
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem

def build(selection="sele"):
    # 一時PDBファイル名
    tmp_pdb = f"/tmp/{selection}.pdb"
    
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
