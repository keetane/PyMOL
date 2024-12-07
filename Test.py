from pymol import cmd
import os

def grid(selection="sele", padding=0.0, lw=1.5, r=0.5, g=0.5, b=0.8):
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)
    print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX - minX, maxY - minY, maxZ - minZ))
    center_of_mass = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    print("Center of mass (%.2f, %.2f, %.2f)" % tuple(center_of_mass))

    minX = minX - float(padding)
    minY = minY - float(padding)
    minZ = minZ - float(padding)
    maxX = maxX + float(padding)
    maxY = maxY + float(padding)
    maxZ = maxZ + float(padding)

    box = [
        LINEWIDTH, float(lw),
        BEGIN, LINES,
        COLOR, float(r), float(g), float(b),
        VERTEX, minX, minY, minZ, VERTEX, maxX, minY, minZ,
        VERTEX, minX, maxY, minZ, VERTEX, maxX, maxY, minZ,
        VERTEX, minX, minY, maxZ, VERTEX, maxX, minY, maxZ,
        VERTEX, minX, maxY, maxZ, VERTEX, maxX, maxY, maxZ,
        VERTEX, minX, minY, minZ, VERTEX, minX, maxY, minZ,
        VERTEX, maxX, minY, minZ, VERTEX, maxX, maxY, minZ,
        VERTEX, minX, minY, maxZ, VERTEX, minX, maxY, maxZ,
        VERTEX, maxX, minY, maxZ, VERTEX, maxX, maxY, maxZ,
        VERTEX, minX, minY, minZ, VERTEX, minX, minY, maxZ,
        VERTEX, maxX, minY, minZ, VERTEX, maxX, minY, maxZ,
        VERTEX, minX, maxY, minZ, VERTEX, minX, maxY, maxZ,
        VERTEX, maxX, maxY, minZ, VERTEX, maxX, maxY, maxZ,
        END
    ]

    boxName = "box_" + str(randint(0, 10000))
    while boxName in cmd.get_names():
        boxName = "box_" + str(randint(0, 10000))

    cmd.load_cgo(box, boxName)
    return boxName

cmd.extend("grid", grid)



def search_space(selection="(all)", padding=4.0):
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)
    print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX - minX, maxY - minY, maxZ - minZ))

    center_of_mass = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    print("Center of mass (%.2f, %.2f, %.2f)" % tuple(center_of_mass))

    minX = minX - float(padding)
    minY = minY - float(padding)
    minZ = minZ - float(padding)
    maxX = maxX + float(padding)
    maxY = maxY + float(padding)
    maxZ = maxZ + float(padding)

    size_x = maxX - minX
    size_y = maxY - minY
    size_z = maxZ - minZ

    center_x = (minX + maxX) / 2
    center_y = (minY + maxY) / 2
    center_z = (minZ + maxZ) / 2
    return size_x, size_y, size_z, center_x, center_y, center_z

    print("\nSmina Options:")
    print("--center_x %.3f," % center_x, "--center_y %.3f," % center_y, "--center_z %.3f," % center_z,
          "--size_x %.3f," % size_x, "--size_y %.3f," % size_y, "--size_z %.3f," % size_z)
cmd.extend("search_space", search_space)


# def smina(selection="sele", padding=4.0):
def smina(ligand='enabled and organic', receptor='enabled and polymer.protein', grid='sele', padding=4.0, n_pose=3, seed=0):
    # temporary file of receptor
    rec='/tmp/rec.pdb'
    # save pdb file of protein from enabled object
    cmd.save(rec, receptor)

    # temporary file of ligand
    lig='/tmp/lig.pdb'
    # save pdb file of protein from enabled object
    cmd.save(lig, ligand)
    # defining grid box
    search_space(grid)

    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(grid)
    center_of_mass = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    minX = minX - float(padding)
    minY = minY - float(padding)
    minZ = minZ - float(padding)
    maxX = maxX + float(padding)
    maxY = maxY + float(padding)
    maxZ = maxZ + float(padding)

    size_x = maxX - minX
    size_y = maxY - minY
    size_z = maxZ - minZ

    center_x = (minX + maxX) / 2
    center_y = (minY + maxY) / 2
    center_z = (minZ + maxZ) / 2


    print("\nExcuting Smina")
    # print("--center_x %.3f," % center_x, "--center_y %.3f," % center_y, "--center_z %.3f," % center_z,
    #       "--size_x %.3f," % size_x, "--size_y %.3f," % size_y, "--size_z %.3f," % size_z)
    docking = f'smina -r {rec} -l {lig} -o /tmp/docked.sdf --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} --seed {seed} --num_modes {n_pose}'
    os.system(docking)
    cmd.load('/tmp/dock.sdf')
cmd.extend("smina", smina)
