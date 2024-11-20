from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem

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
    # cmd.hide(f'({name} and hydro and (elem C extend 1))')

# コマンドをPyMOLに拡張
cmd.extend('smiles', smiles)
