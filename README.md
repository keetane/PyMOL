![PyMOL icon](https://cdn.worldvectorlogo.com/logos/pymol-1.svg)
# PyMOL Plugin
This plugin is based on PyMOL, which has been developed as OSS, with respect and thanks to the PyMOL wiki and the voluntary contributors who have made it freely available.  
The basic development environment is Mac M1 or later, but these Plusins have confirmed that it also works on windows.
## PyMOL install
For easy install of PyMOL, I recommend to use `conda`.  
If you already installed Anaconda or miniconda, just install with command in terminal as follow,
```
conda create -n your-env -y -c conda-forge python=3.9 pymol-open-source rdkit
```
The latest version of PyMOL can be downloaded from Schrödinger https://pymol.org/. This is a paid version and I have not tried it. I have not checked if the following scripts are valid.

## PyMOL instructions
For basic information, I will recommend the link below,
1. PyMOL book
    https://yoshitakamo.github.io/pymol-book/  
    The owere of this tutorial, Dr. Yoshitaka Moriwaki, is one of great contributor of PyMOL. His support helped me enormously with my start of PyMOL. I appreciate his support and contribution to OSS community.
2. PyMOL wiki
    https://pymolwiki.org/index.php/Main_Page
3. Structural analysis at homt (in Japanese)  
    https://ouchidekaiseki.com/index.php

## PyMOL Plugin
This Plugins were described by Python, used pymol and rdkit.
### residues.py
![residues.py](./img/residues.gif)
Visualisation is a very important element in OSS-SBDD, and while PyMOL has a rich implementation of commands for visualisation, this plugin aims to make it easier. Let's get started.

When developing kinase inhibitors, it is necessary to find gatekeepers located deep in the ATP pocket, e.g. threonine, methionine and lysine.  
Or, it will be necessary to identify from the structural information whether the lysine is in the right position on the surface of your target for which you want to explore a degrader.  
After fetching the PDF file `1ATP`, type `see` to focus the surroundings of ligand within 8A distance.
![see](./img/see.gif)

One-letter command for visualization of respective amino acid residues helps you to find the all residues you want to know.  
Type One-letter residue and press Enter, then the sticks of residues you want to know will be showed for enabled objects.
![OneLetter](./img/OneLetter.gif)


Annoying water molecules? Type `del`. You would like to show them again? Type `sol`. Don't you need the hydrophobic hydrogens? Type `ph`. You can focus more into proteins or nucleic acid oligomers.
![sol](./img/sol.gif)


Lots of PDB files includes the structural information of not only monomer, but also multiple complexes. Type `cc` to identify the respective chain information by coloring.  
I like sipmple color cordinates. `white` command will make the color of cartoon and carbon chains white.  
If you would like to change the cartoon color by secondary structural information, type `dssp` to change the color of helixes into red, sheets into yellow and loops into gray, respectively.  
AF2 predicted structure? Type `plddt` to turned cartoon's color by plddt scores.  
![coloring](./img/coloring.gif)

Some PDB files includes the ligand infomation, which bounds to the receptor molecule to interact and get some functions.  
This grateful and useful information will help you to design or modify the molecules. Let's find the intaraction information with this plugin.
`byres` command will make objects of ligand and residues within 4A distance, then navigate you there. Type `tag` to show the residue's name selected. Polar contact object is also created by `byres` command.
Not organid ligand? Well, you may see the peptide or nucleic acid ligand. Typing `at` after selecting the ligand to help you find the polar contacts you would like to see.  
Would you like to see the surface around the ligand? Type `well` to show the surface of around the ligand within 5A distance.



### text2structure.py

![text2structure](./img/text2structure.gif)



## PyMOL Plugins from other sources
My favarite Plugins were described below. Great thanks to all of contributors for PyMOL

### how to install Plugins
1. install as addin from PyMOL  


2. load Plugins from source

3. load Plugins from local file



### Color By Mutations
https://pymolwiki.org/index.php/Color_By_Mutations  
Please refer original script and document. For easy to use, `color_by_mutation.py` was placed on this repository.  
For loading this plugin, run command as follows in PyMOL command line,
```
run https://raw.githubusercontent.com/keetane/PyMOL/refs/heads/main/color_by_mutation.py
```
The following commands are then executed to find the differences between wt and mt residues.
```
color_by_mutation wt, mt
```
### centroid
https://pymolwiki.org/index.php/Centroid  
For loading this plugin, run command as follows in PyMOL command line,
```
run https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/centroid.py
```
then, run command as follow to return the value or the geometric center of your selection or object.
```
centroid [selection of object]
```

### findseq


### PymolFold
https://github.com/JinyuanSun/PymolFold  
Load extension into PyMOL. in the PyMOL command prompt.
```
run https://raw.githubusercontent.com/JinyuanSun/PymolFold/main/pf_plugin.py
```
The following commands are then executed to predict the folded structure you have brought.
```
esmfold [your sequence]
```

### DrawGridBox
https://pymolwiki.org/index.php/DrawGridBox  
Before docking with vina, gnina or other docking software, you may set the gridbox to search the poses molecules docked.  
`DrawGridBox` draw grid boxes around a selection with padding.
Load extension as follows,
```
run https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/drawgridbox.py
```
The following commands will show the grid box, it is very useful to place searching area.

### gridbox
https://github.com/Pymol-Scripts/Pymol-script-repo/blob/master/gridbox.py  
gridbox extentions is also useful to show the gridbox.  
Load extension as follows,
```
run https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/gridbox.py
```
The following commands will show the grid box.
```
gridbox [centroid], size_x, size_y, size_z
```