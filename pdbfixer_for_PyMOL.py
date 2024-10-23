from pymol import cmd
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os

def run_pdbfixer_on_enabled_objects_with_solvent(ph=7):
    enabled_objects = cmd.get_names("objects", enabled_only=1)
    
    for obj in enabled_objects:
        pdb_file = f"/tmp/{obj}.pdb"
        fixed_pdb_file = f"/tmp/fixed_{obj}.pdb"
        solvent_obj = f"solvent_{obj}"
        
        # Save the current enabled object to a PDB file
        cmd.save(pdb_file, obj)
        
        # Extract water molecules and create a new object
        cmd.select(solvent_obj, f"{obj} and resn HOH")
        cmd.create(solvent_obj, f"{solvent_obj}")
        
        # Run pdbfixer
        fixer = PDBFixer(filename=pdb_file)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
        
        with open(fixed_pdb_file, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        
        # Check if the fixed PDB file was created successfully
        if os.path.exists(fixed_pdb_file):
            # Load the fixed PDB file back into PyMOL
            cmd.load(fixed_pdb_file, f"fixed_{obj}")
        else:
            print(f"Error: failed to open file {fixed_pdb_file}")

cmd.extend("pdbfixer", run_pdbfixer_on_enabled_objects_with_solvent)

