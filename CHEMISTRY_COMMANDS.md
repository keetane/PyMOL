# chemistry.py command reference

## Loading the commands

Run the script from the PyMOL command line:

```pml
run https://raw.githubusercontent.com/keetane/PyMOL/refs/heads/main/chemistry.py
```

After a successful load, PyMOL prints:

```text
chemistry: Loaded all commands, including reinvent and reinvent_priors.
```

Run the GitHub loading command again after updating the file. This is especially
important after `target`, because that command changes PyMOL's current working
directory to the new target directory.

All messages printed by the script are in English. Arguments after the first
argument use PyMOL's comma-separated command syntax.

## Contents

- [Alias](#alias)
- [Visualization](#visualization)
- [Structure](#structure)
- [Chemistry](#chemistry)

## Alias

### Label aliases

| Alias | Description | Example |
|---|---|---|
| `tag` | Label selected CA atoms with one-letter residue code and number. | `tag` |
| `b` | Label selected CA atoms with their B-factor. | `b` |
| `untag` | Hide labels in `sele`. | `untag` |

### Residue sticks aliases

Each alias first hides protein sticks and then shows the named residue type in
enabled objects.

| Alias | Residue | Example |
|---|---|---|
| `a` | Alanine | `a` |
| `c` | Cysteine | `c` |
| `d` | Aspartate | `d` |
| `e` | Glutamate | `e` |
| `f` | Phenylalanine | `f` |
| `g` | Glycine | `g` |
| `h` | Histidine | `h` |
| `i` | Isoleucine | `i` |
| `k` | Lysine | `k` |
| `l` | Leucine | `l` |
| `m` | Methionine | `m` |
| `n` | Asparagine | `n` |
| `p` | Proline | `p` |
| `q` | Glutamine | `q` |
| `r` | Arginine | `r` |
| `s` | Serine | `s` |
| `t` | Threonine | `t` |
| `v` | Valine | `v` |
| `w` | Tryptophan | `w` |
| `y` | Tyrosine | `y` |

### Other display aliases

| Alias | Description | Example |
|---|---|---|
| `me` | Show Val, Leu, Ile, and Met sticks in enabled objects. | `me` |
| `z` | Hide all protein sticks. | `z` |
| `addh` | Add hydrogens to `sele`. | `addh` |
| `remh` | Remove selected neighboring hydrogens except donor/acceptor neighbors. | `remh` |
| `del` | Render solvent as hidden nonbonded spheres. | `del` |
| `sol` | Show solvent as nonbonded spheres. | `sol` |
| `ph` | Hide hydrogens attached to carbon in enabled objects. | `ph` |

## Visualization

### `white`

Color only enabled protein polymers white by atom type using PyMOL utilities.
Ligands, water, nucleic acids, and other non-protein atoms are not changed.

```pml
white
```

### `cc`

Color enabled protein chains by chain identifier.

```pml
cc
```

### `dssp`

Color enabled secondary structure: loop gray, strand yellow, and helix red.

```pml
dssp
```

### `plddt`

Color enabled AlphaFold structures by pLDDT values stored in B-factors.

```pml
plddt
```

### `see`

Zoom to enabled residues around organic ligands. The default distance is 8 Å.

```pml
see
see 10
```

### `line`

Show lines for enabled residues around organic ligands. The default distance
is 4 Å.

```pml
line
line 5
```

### `grab`

Select enabled residues around a ligand selection into `sele`. The default
selection is `organic` and the default distance is 4 Å.

```pml
grab
grab ligand, 6
```

### `hb`

Use visible, non-water HETATM atoms as the ligand and show polar contacts to
enabled protein receptor atoms and enabled water molecules. Both receptor and
water contacts are stored together in the `contact` measurement object.

For multi-state objects, the measurement object is generated with all input
states. Changing the current PyMOL state therefore changes the displayed
interactions to the matching state.

```pml
hb
set state, 2
```

### `dock`

Show polar contacts between a ligand selection and other atoms as
`polar_conts`. The default ligand selection is `enabled and organic`.

```pml
dock
dock Lig
```

### `at`

Show polar contacts between a selection and enabled protein residues or water.
The default selection is `sele`. Protein and water contacts are stored together
in the `contact` measurement object, with separate measurements for each state.
The command does not zoom or call `see`.

```pml
at
at ligand
at visible and hetatm and not solvent
```

### `byres`

Create `Lig` and a nearby-residue object from a ligand selection. The first
argument is the ligand selection and defaults to `enabled and organic`; the
second argument is the neighborhood distance and defaults to 10 Å. The nearby
object is named `byres_<original_object>`, using the first source object in the
ligand selection.

The source objects are disabled. The new neighborhood object contains complete
residues and complete water molecules within the requested distance. Protein
residues within 4 Å of `Lig` are shown as lines. The command then calls `hb` to
create state-aware receptor and water polar contacts, then runs `sol` to show
water molecules as non-bonded spheres. It does not call `grab`.

```pml
byres
byres sele
byres visible and hetatm and not solvent, 8
```

For example, if `sele` belongs to an object named `docked`, `byres sele`
creates `byres_docked`.

### `well`

Show the surface of enabled residues near organic ligands. The default
distance is 5 Å.

```pml
well
well 7
```

### `cav` and `wall`

Enable or disable PyMOL surface cavity mode.

```pml
cav
wall
```

### `pl`

Apply PyMOL's ligand-cartoon preset to enabled objects and show solvent as
spheres.

```pml
pl
```

## Structure

### `target`

Create and enter a target directory, search reviewed human UniProt entries,
fetch all referenced PDB structures, select the largest non-solvent HETATM
molecule in each structure, align the associated protein chains, prepare the
binding sites, and save a session, structure list, and log.

PDB objects use only their uppercase PDB IDs. Temporary `align_*` and `site_*`
selections are deleted; remaining selections are put in the `select` group and
disabled. After the PSE is saved, fetched `.cif` and `.cif.gz` files in the new
target directory are deleted.

```pml
target JAK1
target PFKFB3
```

### `af2`

Search AlphaFold DB for a human protein and load its predicted structure. The
default object name follows `AF_<target>_xx`.

```pml
af2 TP53
af2 P04637, p53_model
```

### `fasta`

Print the FASTA sequence of enabled objects and copy the printed FASTA text to
the clipboard.

```pml
fasta
```

### `seq`

Print the FASTA sequence of `sele` and copy the printed FASTA text to the
clipboard.

```pml
select sele, chain A
seq
```

### `pubchem`

Search PubChem by name, copy the returned SMILES to the clipboard, and load
the PubChem 3D SDF conformer.

```pml
pubchem aspirin
pubchem ruxolitinib
```

### `smiles`

Create a 3D RDKit conformer from a SMILES string and load it as a sticks
object. An optional second argument specifies the object name.

```pml
smiles "CCO"
smiles "CC(=O)Oc1ccccc1C(=O)O", aspirin
```

### `peptide`

Build a peptide from one-letter residue codes, generate a 3D conformer, and
load it as sticks. The syntax is `peptide sequence, mode, object_name`.

Modes are `0` linear, `1` head-to-tail, `2` simple disulfide, and `3` cysteine
cyclization. Uppercase residues are L forms and lowercase residues are D forms.
The special code `g` is Aib; `B/b` and `J/j` are paired bridge markers.

```pml
peptide ACDE
peptide ACDE, 1, cyclic_ACDE
peptide AcDe, 0, mixed_stereo
```

## Chemistry

### `printsmiles`

Convert HETATM records in `sele`, or another selection, to SMILES, print the
SMILES, and copy it to the clipboard.

```pml
printsmiles
printsmiles ligand
```

### `loadsdf`

Load a multi-record SDF as a multi-state PyMOL object while retaining every SD
tag for its corresponding state. Property names are not fixed, so commercial
library fields such as `MW`, `HBD`, `HBA`, `cLogP`, vendor identifiers, SMILES,
and custom fields are retained alongside docking scores.

The optional second argument sets the object name. If omitted, a legal PyMOL
name is derived from the filename. Existing objects are not overwritten.

```pml
loadsdf /path/to/vendor_library.sdf, library
loadsdf /path/to/6ETJ_BWW_redock.sdf, redock
```

### `sdfinfo`

Print retained SD properties for the currently displayed state. The first
argument must resolve to exactly one object. The optional second argument is an
exact, case-insensitive property name; omit it to print all fields.

```pml
sdfinfo library
sdfinfo library, MW
sdfinfo library, Vendor_ID
```

### `printscore`

Print numeric docking-score properties for the currently displayed PyMOL
state. By default, the command searches enabled molecular objects and detects
property names containing `score`, `affinity`, `energy`, `confidence`, or
`binding`, as well as RMSD lower- and upper-bound properties. Rank, pose, and
timing properties are excluded from automatic detection.

This supports state properties from smina, GNINA, and other SDF outputs without
requiring a fixed score-property name. Use the optional second argument to
request an exact property name. Property-name matching is case-insensitive.

On older PyMOL versions without the object-property API, `printscore` first
tries PyMOL's exported SDF and then looks for `<object_name>.sdf` in the current
directory, the directory containing `chemistry.py`, its immediate
subdirectories, and the user's Downloads directory. The optional third argument
can specify the original SDF explicitly. Use `auto` as the second argument when
automatic score-name detection is wanted with an explicit source file.

```pml
printscore
printscore docked
printscore docked, CNNscore
printscore enabled, iScreeM_Score_affinity
printscore 20261101_2220_6etj, auto, /Users/keetane/Downloads/20261101_2220_6etj.sdf
```

Typical output:

```text
printscore: object=docked, state=1, iScreeM_Score_affinity=-11.6
```

### `sdffilter`

Create a new multi-state object containing only records whose retained SD
properties match a safe filter expression. The original object is unchanged,
and metadata is copied to the corresponding states in the new object.

Supported comparisons are `==`, `!=`, `<`, `<=`, `>`, and `>=`. Combine them
with lowercase `and`, `or`, `not`, and parentheses. Simple property names can
be used directly. Use `prop("property name")` for names containing spaces or
other characters. Arbitrary Python calls and expressions are rejected.

```pml
sdffilter library, MW < 500 and HBD <= 5, leadlike
sdffilter library, cLogP <= 4 and HBA < 10, permeable
sdffilter library, 'prop("Vendor ID") == "CMPD-001"', chosen
```

If the third argument is omitted, the output is named
`filtered_<original_object>`.

### `reinvent_priors`

Locate the REINVENT4 runtime and list discovered Mol2Mol `.prior` files. Set
`REINVENT_PYTHON`, `REINVENT_PRIORS_DIR`, or `REINVENT4_HOME` before running
PyMOL when automatic detection is insufficient.

```pml
reinvent_priors
```

### `reinvent`

Run asynchronous REINVENT4 Mol2Mol sampling for exactly one selected ligand
residue. Syntax: `reinvent selection, prior, num_smiles, temperature`. `prior`
may be `all`, one prior name, or names joined with `+`. Device preference is
CUDA, then MPS, then CPU. Results and logs are written below `reinvent/`.

```pml
reinvent sele
reinvent ligand, all, 10000, 1.0
reinvent ligand, mol2mol_medium, 2000, 0.8
```

### `mkdir`

Create a directory, including missing parent directories.

```pml
mkdir results
mkdir results/models
```

### `rm`

Delete one file or files matching a non-recursive glob. An empty argument is
rejected. Directories are not deleted.

```pml
rm obsolete.pdb
rm *.cif
```
