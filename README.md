# InterPhon

A Python Package for Ab initio Interface Phonon Calculations within a 3D Electronic Structure Framework

## Installation

Please clone this repository and install using:

```
$ git clone https://github.com/inwonyeu/interphon.git
$ cd interphon/
$ python setup.py install
```

## Basic usage in conjunction with VASP

InterPhon supports a range of options to manage phonon computations and plotting styles.  
In order to see all of the available options and their default values:

```
$ interphon --help
```



1. ### **Pre-process**
In the InterPhon pre-process, a file of supercell (**SUPERCELL**) and files of displaced supercells (**POSCAR-0***) are generated:

```
$ interphon -enlarge "2 2 1" -pbc "1 1 0"
```

-> (2×2×1) supercell and displaced supercells  
-> Periodic boundary conditions (1 or True) along a1, a2 lattice directions, while open (0 or False) along a3 direction



2. ### **Post-process**
After the DFT force calculations for the displaced supercells (**POSCAR-0***) are finished in each **FORCE-0*** folder *(folder names are arbitrary—only the order of the folder numbers is important)*, the evaluation of interfacial phonons can be executed by the following ways.

- ***Density of states (DOS):***
```
$ interphon FORCE-0*/vasprun.xml -kdos KPOINTS_dos
```

- ***Thermal properties:***
```
$ interphon FORCE-0*/vasprun.xml -kdos KPOINTS_dos -thermal
```

- ***Band:***
```
$ interphon FORCE-0*/vasprun.xml -kband KPOINTS_band
```

- ***Phonon mode:***
```
$ interphon FORCE-0*/vasprun.xml -kband KPOINTS_band -mode
```

where ‘KPOINTS_dos’ and ‘KPOINTS_band’ are files for the mesh sampling of k-points supported in VASP format. (<https://www.vasp.at/wiki/index.php/KPOINTS>)

***KPOINTS_dos (file name is arbitrary):***
```
kpoint
0
MP  # Monkhorst-Pack grids, use the first character ‘G’ for Gamma-centered grids.
9 9 1
0.0 0.0 0.0
```

***KPOINTS_band (file name is arbitrary):***
```
kpoint
41
L
0.00 0.00 0.00  # G
0.00 0.50 0.00  # M

0.00 0.50 0.00  # M
0.333333 0.333333 0.00  # K

0.333333 0.333333 0.00  # K
0.00 0.00 0.00  # G
```

## Learn more
Detailed manual will be available soon.
