from pymol import cmd
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def prep_smi(output_file='ligand.sdf'):
    # 現在のディレクトリ内のすべての.smiファイルを読み込み
    smiles_files = glob.glob('*.smi')
    
    # データフレームを作成
    df = pd.DataFrame(columns=['SMILES'])
    for n in smiles_files:
        smi = pd.read_csv(n, sep='\t', header=None, names=['SMILES'])
        df = pd.concat([df, smi], ignore_index=True)
    
    # SMILESから分子オブジェクトに変換し、三次元化
    df['ROMol'] = df['SMILES'].map(lambda x: Chem.MolFromSmiles(x))
    df['ROMol'] = df['ROMol'].map(lambda x: Chem.AddHs(x))
    df['3d'] = df['ROMol'].map(lambda x: AllChem.EmbedMolecule(x))
    
    # SDFファイルに保存
    writer = Chem.SDWriter(output_file) 
    for mol in df['ROMol']:
        if mol is not None:  # 分子が正常に変換されたか確認 
            writer.write(mol) 
    writer.close()

# PyMOLにコマンドとして登録
cmd.extend("prep_smi", pre_smi)
