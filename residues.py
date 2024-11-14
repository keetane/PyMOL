from pymol import cmd
# Visualization
### One-letter Labeling
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

cmd.alias('tag', 'label n. CA and sele, "%s%s" % (one_letter[resn], resi)')
cmd.alias('b', 'label n. CA and sele, b')
cmd.alias('untag', 'hide label, sele')

### hide sticks, show residues
cmd.alias('k', 'hide sticks, polymer.protein; show sticks, enabled and resn lys')
cmd.alias('c', 'hide sticks, polymer.protein; show sticks, enabled and resn cys')
cmd.alias('s', 'hide sticks, polymer.protein; show sticks, enabled and resn ser')
cmd.alias('h', 'hide sticks, polymer.protein; show sticks, enabled and resn his')
cmd.alias('t', 'hide sticks, polymer.protein; show sticks, enabled and resn thr')
cmd.alias('r', 'hide sticks, polymer.protein; show sticks, enabled and resn arg')
cmd.alias('d', 'hide sticks, polymer.protein; show sticks, enabled and resn asp')
cmd.alias('e', 'hide sticks, polymer.protein; show sticks, enabled and resn glu')
cmd.alias('n', 'hide sticks, polymer.protein; show sticks, enabled and resn asn')
cmd.alias('q', 'hide sticks, polymer.protein; show sticks, enabled and resn gln')
cmd.alias('y', 'hide sticks, polymer.protein; show sticks, enabled and resn tyr')
cmd.alias('f', 'hide sticks, polymer.protein; show sticks, enabled and resn phe')
cmd.alias('w', 'hide sticks, polymer.protein; show sticks, enabled and resn trp')
cmd.alias('g', 'hide sticks, polymer.protein; show sticks, enabled and resn gly')
cmd.alias('a', 'hide sticks, polymer.protein; show sticks, enabled and resn ala')
cmd.alias('m', 'hide sticks, polymer.protein; show sticks, enabled and resn met')
cmd.alias('p', 'hide sticks, polymer.protein; show sticks, enabled and resn pro')
cmd.alias('i', 'hide sticks, polymer.protein; show sticks, enabled and resn ile')
cmd.alias('l', 'hide sticks, polymer.protein; show sticks, enabled and resn leu')
cmd.alias('v', 'hide sticks, polymer.protein; show sticks, enabled and resn val')
cmd.alias('me', 'hide sticks, polymer.protein; show sticks, enabled and resn val; show sticks, enabled and resn leu; show sticks, enabled and resn ile; show sticks, enabled and resn met')
cmd.alias('z', 'hide sticks, polymer.protein')
cmd.alias('addh', 'cmd.h_add("sele")')
cmd.alias('remh', 'cmd.remove("sele & hydro & not nbr. (don.|acc.)")')

### water and hydrophobic hydrogens manupilations
cmd.alias('del', 'as nb_spheres, solvent; hide nb_spheres')
cmd.alias('sol', 'as nb_spheres, solvent')
cmd.alias('ph', 'hide (hydro) and enabled and (elem C extend 1)')


### coloring
def white():
    util.cba(144, "enabled", _self=cmd)
cmd.extend('white',white)

def cc():
    util.color_chains("(enabled and name CA)", _self=cmd)
cmd.extend('cc',cc)

def dssp():
    cmd.color('gray', 'enabled and ss l')
    cmd.color('yellow', 'enabled and ss s')
    cmd.color('red', 'enabled and ss h')
cmd.extend('dssp',dssp)

def idr():
    cmd.hide('cartoon', 'v. and b < 50')
cmd.extend('idr',idr)

def plddt():
    cmd.color('0x0053D6', 'enabled and b < 100')
    cmd.color('0x65CBF3', 'enabled and b < 90')
    cmd.color('0xFFDB13', 'enabled and b < 70')
    cmd.color('0xFF7D45', 'enabled and b < 50')
    util.cnc()
    cmd.color('white', 'enabled and sc. and elem C')
cmd.extend('plddt',plddt)

## interaction analysis
### zoom into around ligand
def see(distance=8):
    cmd.zoom('byres organic around %s and enabled' % distance)
cmd.extend('see', see)

### show lines around ligand
def line(distance=4):
    cmd.show('lines', 'byres organic around %s and enabled' % distance)
cmd.extend('line', line)

### select residues with distance around ligand
def grab(sele='organic', distance=4):
    cmd.select('sele', 'byres %s around %s and enabled' % (sele, distance))
    cmd.enable('sele')
cmd.extend('grab', grab)

### show polar contact of ligand to residues in the enabled same object
def hb():
    cmd.select('ligand', 'organic and enabled')
    cmd.dist('contact', '(ligand)', '(byobj (ligand)) and (not(ligand))', quiet=1, mode=2, label=0, reset=1)
    cmd.enable('contact')
cmd.extend('hb',hb)

### show polar contact of ligand to residues between the enabled different object
def dock(ligand='enabled and organic'):
    cmd.delete('polar_conts')
    cmd.dist("polar_conts","%s" % ligand, "(not %s)" % ligand ,quiet=1,mode=2,label=0,reset=1)
    cmd.enable("polar_conts")
cmd.extend('dock', dock)

### show polar contact of sele to residues
def at():
    cmd.dist('contact', '(sele)', '(byobj (sele)) and (not(sele))', quiet=1, mode=2, label=0, reset=1)
    cmd.enable('contact')
    see(distance)
cmd.extend('at',at)

### create object of ligand and residues with distance
def byres(distance=8, ligand='organic and enabled'):
    object_name = 'byres' + str(distance)
    cmd.create('Lig', '%s' % ligand)
    cmd.create(object_name, 'byres %s around %s' % (ligand, distance))
    cmd.disable('!' + object_name)
    cmd.enable('Lig')
    line()
    dock('Lig')
    see()
    grab()
cmd.extend('byres',byres)

# def byres(distance=8):
#     object_name = 'byres' + str(distance)
#     cmd.create('Lig', 'enabled and organic')
#     cmd.create(object_name, 'byres organic around %s and enabled' % distance)
#     cmd.disable('!' + object_name)
#     cmd.enable('Lig')
#     line()
#     dock('Lig')
#     see()
#     grab()
# cmd.extend('byres',byres)

### show surface of residues with distance from ligand
def well(distance=5):
    cmd.show('surface', 'byres organic expand %s and enabled' % distance)
cmd.extend('well',well)

### change surface cavity mode
def cav():
    cmd.set('surface_cavity_mode', 2)
cmd.extend('cav',cav)

def wall():
    cmd.set('surface_cavity_mode', 0)
cmd.extend('wall',wall)

### preset ligand cartoon
def pl():
    preset.ligand_cartoon("enabled", _self=cmd)
    cmd.show('spheres', 'solvent')
cmd.extend('pl',pl)


## sequencing
### print fasta sequence of enabled object
def fasta():
    print(cmd.get_fastastr('enabled'))
cmd.extend('fasta',fasta)

### print fasta sequence of sele object
def seq():
    print(cmd.get_fastastr('sele'))
cmd.extend('seq',seq)

