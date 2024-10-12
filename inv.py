from pymol import cmd

def invert_coordinates(selection='visible'):
    model = cmd.get_model(selection)
    for atom in model.atom:
        atom.coord = [-atom.coord[0], -atom.coord[1], -atom.coord[2]]
    cmd.load_model(model, 'inv')
    cmd.hide('lines', 'inv')
    cmd.show('cartoon', 'inv')

cmd.extend('invertcoordinates', invert_coordinates)
