#%%
from pymol import cmd
import subprocess

def pdbfixer(obj_name=None, ph=7):
    if obj_name is None:
        # デフォルトでアクティブなオブジェクトを指定
        obj_name = cmd.get_names('objects', enabled_only=1)[0]

    new_obj = 'fixed_' + obj_name
    cmd.save(f'{new_obj}.pdb', obj_name)
    # cmd.create('solvent_' + obj_name, obj_name + ' and resn HOH')
    
    # pdbfixerコマンドを実行
    command = 'pdbfixer ' + new_obj + ' --out ' + new_obj + ' --ph=' + str(ph)
    subprocess.run(command.split())
    
    # 修正後の構造を読み込む
    cmd.load(f'{new_obj}.pdb')
    # cmd.remove(f'(hydro) and {new_obj} and resn arg')
    # cmd.h_add(f'{new_obj} and resn arg')

# PyMOLのコマンドラインで実行可能にする
cmd.extend('pdbfixer', pdbfixer)


