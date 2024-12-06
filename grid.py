from pymol import cmd
from pymol.cgo import *
import itertools
from random import randint
# This script was described based on drawgridbox@PyMOL wiki. Users must check original script and license especially for commercial use.
# https://pymolwiki.org/index.php/DrawGridBox

"""
DESCRIPTION
    Given selection, draw a box around it.

USAGE:
    drawbox [selection, [padding, [lw, [r, [g, b]]]]]

PARAMETERS:
    selection,    the selection to enboxen
                    defaults to (all)

    padding,      defaults to 0

    lw,           line width
                    defaults to 2.0

    r,            red color component, valid range is [0.0, 1.0]
                    defaults to 1.0

    g,            green color component, valid range is [0.0, 1.0]
                    defaults to 1.0

    b,            blue color component, valid range is [0.0, 1.0]
                    defaults to 1.0

RETURNS
    string, the name of the CGO box

NOTES
    * This function creates a randomly named CGO box. The user can
    specify the width of the lines, the padding and also the color.
"""

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

    print("\nSmina Options:")
    print("--center_x %.3f," % center_x, "--center_y %.3f," % center_y, "--center_z %.3f," % center_z,
          "--size_x %.3f," % size_x, "--size_y %.3f," % size_y, "--size_z %.3f," % size_z)

cmd.extend("search_space", search_space)
