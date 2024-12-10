from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd
from openbabel import openbabel
import datetime

# mkdir like unix command
def mkdir(dir:str):
    os.makedirs(dir)
cmd.extend('mkdir', mkdir)

# rm command
def rm(file:str):
    os.remove(file)
cmd.extend('rm', rm)


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

# calling molecule from PubChem database (with pH)
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

# Calling smiles from PubChem database
def pubchem2smiles(compound_name: str):
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
cmd.extend('pubchem2smiles', pubchem2smiles)

# print smiles
def printsmiles(selection="sele"):
    # 一時PDBファイル名
    tmp_pdb = "/tmp/tmp_molecule.pdb"
    # PyMOLから選択された分子をPDB形式で保存
    cmd.save(tmp_pdb, selection)

    # Open Babelを使用してPDBファイルを読み込み
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, tmp_pdb)
    
    # SMILES形式に変換し、オブジェクト名を追加して出力する
    smiles = obConversion.WriteString(mol).strip()
    object_name = cmd.get_object_list(selection)[0]
    smiles = smiles.strip('/tmp/tmp_molecule.pdb').replace(".", "\n")
    print(f"SMILES:\n{object_name}\t{smiles}")
# コマンドをPyMOLに拡張
cmd.extend("printsmiles", printsmiles)

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


# Saving .smi file of selected objects
def savesmiles(output_file=None, selection="sele"):
    # デフォルトのファイル名を現在の日付と時間に基づいて生成
    if output_file is None:
        output_file = datetime.datetime.now().strftime("%Y%m%d%H%M") + ".smi"
    else:
        output_file = output_file + ".smi"
    
    # 一時PDBファイル名
    tmp_pdb = "/tmp/tmp_molecule.pdb"
    # PyMOLから選択された分子をPDB形式で保存
    cmd.save(tmp_pdb, selection)

    # Open Babelを使用してPDBファイルを読み込み
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, tmp_pdb)
    
    # オブジェクト名を取得
    object_names = cmd.get_object_list(selection)
    
    # 各オブジェクトに対してSMILESを生成
    smiles_lines = []
    for obj in object_names:
        mol.Clear()
        cmd.save(tmp_pdb, obj)
        obConversion.ReadFile(mol, tmp_pdb)
        smiles = obConversion.WriteString(mol).strip()
        # 不要なファイルパスを取り除き、オブジェクト名を追加
        smiles = smiles.split()[0]
        smiles_lines.append(f"{smiles}\t{obj}")
    
    # SMILESをファイルに書き込み
    with open(output_file, 'w') as file:
        file.write('\n'.join(smiles_lines) + '\n')
    
    print(f"SMILESが{output_file}に保存されました。")

# コマンドをPyMOLに拡張
cmd.extend("savesmiles", savesmiles)



# Loading the molecules from .smi file and saving them as a single SDF file
def loadsmiles(input_file):
    # 出力する一時SDFファイル名
    tmp_sdf = f"/tmp/{input_file}.sdf"
    
    # Open Babelを使用してSMILESファイルを読み込み、3次元化してSDFファイルに変換
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "sdf")
    
    # OBMolオブジェクトの初期化
    mol = openbabel.OBMol()
    
    # SDFコンテンツを保存するリスト
    sdf_content = []

    # SMILESファイルを1行ずつ読み込み
    with open(input_file, 'r') as file:
        smiles_lines = file.readlines()
    
    for line in smiles_lines:
        line = line.strip()
        if not line:
            continue
        
        smiles, name = line.split()
        
        # SMILES文字列を分子オブジェクトに変換
        mol.Clear()
        obConversion.ReadString(mol, smiles)
        
        if mol.Empty():
            print(f"エラー: SMILES文字列 '{smiles}' から分子を生成できませんでした。")
            continue
        
        # 3次元座標を生成
        builder = openbabel.OBBuilder()
        builder.Build(mol)
        
        # 分子をSDF形式の文字列に変換して追加
        sdf_string = obConversion.WriteString(mol)
        sdf_content.append(sdf_string)
    
    # SDFコンテンツをファイルに書き込み
    with open(tmp_sdf, 'w') as sdf_file:
        sdf_file.writelines(sdf_content)
    
    # PyMOLにSDFファイルをロード
    cmd.load(tmp_sdf)
    print(f"{input_file}の全ての分子が3次元化されてPyMOLにロードされました。")

# コマンドをPyMOLに拡張
cmd.extend("loadsmiles", loadsmiles)

# concat and prep *.smi
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
cmd.extend("prep_smi", prep_smi)

