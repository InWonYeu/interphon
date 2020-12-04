# InterPhon

A Python Package for Ab initio Interface Phonon Calculations within a 3D Electronic Structure Framework

</br>
## Installation

Please clone this repository and install using:

```
$ git clone https://github.com/inwonyeu/interphon.git
$ cd interphon/
$ python setup.py install
```

</br>
## Basic usage in conjunction with VASP

InterPhon supports a range of options to manage phonon computations and plotting styles.  
In order to see all of the available options and their default values:

```
$ interphon --help
```

### 1. Pre-process
In the InterPhon pre-process, a file of supercell (**SUPERCELL**) and files of displaced supercells (**POSCAR-0***) are generated:

```
$ interphon -enlarge "2 2 1" -pbc "1 1 0"
```

-> (2×2×1) supercell and displaced supercells  
-> Periodic boundary conditions (1 or True) along a1, a2 lattice directions, while open (0 or False) along a3 direction

### 2. Post-process
After the DFT force calculations for the displaced supercells (**POSCAR-0***) are finished in each **FORCE-0*** folder, the evaluation of interfacial phonons can be executed by the following ways:

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

</br>
## Important files

### 1. DFT input file
Within InterPhon, the interfacial region is supposed to be defined through the statement of constraints on atom movements (selective dynamics).
Phonon evaluation proceeds only in the selected atoms.  
See below example of Cu(111) surface where the top three layers are selected as the surface region.

**POSCAR (VASP format):**
```
Unknown
1.00000000000000
   2.5712952614000000    0.0000000000000000    0.0000000000000000
   1.2856476307000000    2.2268070170000001    0.0000000000000000
   0.0000000000000000    0.0000000000000000   27.7901687622000004
Cu
7
Selective dynamics
Cartesian
   0.0000000000000000    0.0000000000000000    6.2983610849999998   F   F   F
   2.5712951080000002    1.4845379230000000    8.3978147799999991   F   F   F
   1.2856475540000001    0.7422689609999999   10.5098654109999998   F   F   F
   0.0000000000000000    0.0000000000000000   12.6153545979999997   F   F   F
   2.5712951080000002    1.4845379230000000   14.7267644020000006   T   T   T
   1.2856475540000001    0.7422689609999999   16.8249826169999999   T   T   T
   0.0000000000000000    0.0000000000000000   18.9121107990000006   T   T   T
```

### 2. K-points file
The above **KPOINTS_dos** and **KPOINTS_band** files, which are supported in VASP format (<https://www.vasp.at/wiki/index.php/KPOINTS>), are used for the mesh sampling of k-points.

**KPOINTS_dos (file name is arbitrary):**
```
kpoint
0
MP  # Monkhorst-Pack grids, use the first character ‘G’ for Gamma-centered grids.
9 9 1
0.0 0.0 0.0
```

**KPOINTS_band (file name is arbitrary):**
```
kpoint
41
L  # Line path
0.00 0.00 0.00  # G
0.00 0.50 0.00  # M

0.00 0.50 0.00  # M
0.333333 0.333333 0.00  # K

0.333333 0.333333 0.00  # K
0.00 0.00 0.00  # G
```

## Learn more
Detailed manual will be available soon.
