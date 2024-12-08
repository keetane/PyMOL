import os
import subprocess
import pymol.cmd as cmd

def smina(ligand='enabled and organic', receptor='enabled and polymer.protein', grid=None, padding=4.0, n_pose=3, seed=0, comment=None):
    # temporary file of receptor
    rec='/tmp/rec.pdb'
    # save pdb file of protein from enabled object
    cmd.save(rec, receptor)

    # temporary file of ligand
    lig='/tmp/lig.pdb'
    # save pdb file of protein from enabled object
    cmd.save(lig, ligand)
    if grid==None:
        grid=ligand
    else:
        grid = grid + ' and polymer'    

    # defining grid box
    cmd.select("grid", grid)
    [minX, minY, minZ], [maxX, maxY, maxZ] = cmd.get_extent("grid")
    minX -= float(padding)
    minY -= float(padding)
    minZ -= float(padding)
    maxX += float(padding)
    maxY += float(padding)
    maxZ += float(padding)

    size_x = maxX - minX
    size_y = maxY - minY
    size_z = maxZ - minZ

    center_x = (minX + maxX) / 2
    center_y = (minY + maxY) / 2
    center_z = (minZ + maxZ) / 2

    if comment==None:
        comment=''
    
    print("\nExecuting Smina")
    docking = [
        'smina',
        '-r', rec,
        '-l', lig,
        '-o', f'./docked_{comment}.sdf',
        '--center_x', str(center_x),
        '--center_y', str(center_y),
        '--center_z', str(center_z),
        '--size_x', str(size_x),
        '--size_y', str(size_y),
        '--size_z', str(size_z),
        '--seed', str(seed),
        '--num_modes', str(n_pose)
    ]
    
    process = subprocess.Popen(docking, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    print(stdout.decode())
    print(stderr.decode())
    
    cmd.load(f'./docked_{comment}.sdf')

cmd.extend("smina", smina)
