#%%
from pymol import cmd
import subprocess

def pdbfixer(obj_name=None):
    if obj_name is None:
        # デフォルトでアクティブなオブジェクトを指定
        obj_name = cmd.get_names('objects', enabled_only=1)[0]

    new_obj = 'fixed_' + obj_name + '.pdb'
    cmd.save(new_obj, obj_name)
    cmd.create('solvent_' + obj_name, obj_name + ' and resn HOH')
    
    # pdbfixerコマンドを実行
    command = 'pdbfixer ' + new_obj + ' --out ' + new_obj
    subprocess.run(command.split())
    
    # 修正後の構造を読み込む
    cmd.load(new_obj)

# PyMOLのコマンドラインで実行可能にする
cmd.extend('pdbfixer', pdbfixer)


