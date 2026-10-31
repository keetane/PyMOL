# chemistry.py command reference

## Loading the commands

Run the script from the PyMOL command line:

```pml
run /Users/keetane/Documents/modeling/PyMOL/chemistry.py
```

After a successful load, PyMOL prints:

```text
chemistry: Loaded all commands, including reinvent and reinvent_priors.
```

Run the absolute path again after updating the file. This is especially
important after `target`, because that command changes PyMOL's current working
directory to the new target directory.

All messages printed by the script are in English. Arguments after the first
argument use PyMOL's comma-separated command syntax.

## Structure and chemistry commands

### `af2`

Search AlphaFold DB for a human protein and load its predicted structure. The
default object name follows `AF_<target>_xx`.

```pml
af2 TP53
af2 P04637, p53_model
```

### `printsmiles`

Convert HETATM records in `sele`, or another selection, to SMILES, print the
SMILES, and copy it to the clipboard.

```pml
printsmiles
printsmiles ligand
```

### `smiles`

Create a 3D RDKit conformer from a SMILES string and load it as a sticks
object. An optional second argument specifies the object name.

```pml
smiles "CCO"
smiles "CC(=O)Oc1ccccc1C(=O)O", aspirin
```

### `pubchem`

Search PubChem by name, copy the returned SMILES to the clipboard, and load
the PubChem 3D SDF conformer.

```pml
pubchem aspirin
pubchem ruxolitinib
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

## Visualization commands

### `white`

Color enabled objects white by atom type using PyMOL utilities.

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

Create `ligand` from enabled organic atoms and show polar contacts within the
same object as `contact`.

```pml
hb
```

### `dock`

Show polar contacts between a ligand selection and other atoms as
`polar_conts`. The default ligand selection is `enabled and organic`.

```pml
dock
dock Lig
```

### `at`

Show polar contacts for `sele` and zoom to the ligand site. The default zoom
distance is 8 Å.

```pml
at
at 10
```

### `byres`

Create `Lig` and a disabled nearby-residue object, then prepare lines, polar
contacts, zoom, and `sele`. The default distance is 8 Å.

```pml
byres
byres 6, organic and enabled
```

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

## Aliases

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
