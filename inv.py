from pymol import cmd

def invert_coordinates():
    model = cmd.get_model('visible')
    for atom in model.atom:
        atom.coord = [-atom.coord[0], -atom.coord[1], -atom.coord[2]]
    cmd.load_model(model, 'inv')
    cmd.hide('lines', 'inv')
    cmd.show('cartoon', 'inv')

cmd.extend('invert_coordinates', invert_coordinates)
