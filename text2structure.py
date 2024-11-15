from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp
from pymol import cmd
from openbabel import openbabel
import datetime


# https://iwatobipen.wordpress.com/2024/05/14/add-hydrogen-with-user-defined-ph-from-python-openbabel-cheminformatics/
obc = openbabel.OBConversion()
obc.SetInAndOutFormats('smi', 'smi')

def get_smi_with_pH(smi, pH:float=7.4):
    obmol = openbabel.OBMol()
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


def psmi(selection="sele"):
    # 一時PDBファイル名
    tmp_pdb = "/tmp/tmp_molecule.pdb"
    # PyMOLから選択された分子をPDB形式で保存
    cmd.save(tmp_pdb, selection)

    # Open Babelを使用してPDBファイルを読み込み
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, tmp_pdb)
    
    # SMILES形式に変換し、出力する
    smiles = obConversion.WriteString(mol).strip()
    smiles = smiles.strip('/tmp/tmp_molecule.pdb').replace(".", "\n")
    print(f"SMILES:\n{smiles}")
# コマンドをPyMOLに拡張
cmd.extend("psmi", psmi)



from pymol import cmd
from openbabel import openbabel

def psmi(selection="sele"):
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
cmd.extend("psmi", psmi)



def ss(output_file=None, selection="sele"):
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
        # ファイルパスを取り除き、SMILESを前に追加し、オブジェクト名を後に追加
        smiles = smiles.replace('/tmp/tmp_molecule.pdb', '')
        smiles_lines.append(f"{smiles}\t{obj}")
    
    # SMILESをファイルに書き込み
    with open(output_file, 'w') as file:
        file.write('\n'.join(smiles_lines) + '\n')
    
    print(f"SMILESが{output_file}に保存されました。")

    # 一時ファイルを削除
    os.remove(tmp_pdb)

# コマンドをPyMOLに拡張
cmd.extend("ss", ss)


def loadsmi(input_file):
    # 出力する一時SDFファイル名
    output_file = "/tmp/tmp_molecule.sdf"
    
    # Open Babelを使用してSMILESファイルを読み込み、3次元化してSDFファイルに変換
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "sdf")
    
    # SMILESファイルを読み込み
    mol = openbabel.OBMol()
    not_empty = obConversion.ReadFile(mol, input_file)
    
    if not not_empty:
        print(f"エラー: {input_file}から分子を読み込めませんでした。")
        return

    # 3次元座標を生成
    builder = openbabel.OBBuilder()
    builder.Build(mol)
    
    # エネルギー最小化を実行
    forcefield = openbabel.OBForceField.FindForceField("mmff94")
    forcefield.Setup(mol)
    forcefield.ConjugateGradients(500)  # 最大500ステップの最小化を実行
    forcefield.GetCoordinates(mol)
    
    # 一時SDFファイルに書き込み
    obConversion.WriteFile(mol, output_file)

    # 追加の分子がある場合に処理するループ
    while obConversion.Read(mol):
        builder.Build(mol)
        forcefield.Setup(mol)
        forcefield.ConjugateGradients(500)  # 最大500ステップの最小化を実行
        forcefield.GetCoordinates(mol)
        with open(output_file, 'a') as sdf_file:
            sdf_string = obConversion.WriteString(mol)
            sdf_file.write(sdf_string)
    
    # PyMOLにSDFファイルをロード
    cmd.load(output_file)
    print(f"{input_file}の全ての分子が最適化されてPyMOLにロードされました。")

# コマンドをPyMOLに拡張
cmd.extend("loadsmi", loadsmi)
