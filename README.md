# InterPhon

A Python Library to Calculate 2D Interface Phonon within 3D Electronic Structure Framework

## Installation

Please clone this repository and install using:

```
$ git clone https://github.com/inwonyeu/interphon.git
$ cd interphon/
$ python setup.py install
```

## Basic usage in conjunction with VASP

InterPhon supports a lot of options to manage phonon computations and plotting styles. 
In order to see all of the available options and their default values:

```
$ interphon --help
```

- ### Pre-process
By the InterPhon pre-process, a file of supercell (SUPERCELL) accommodating several unit cells and files of displaced supercell (POSCAR-0*) are generated:

```
$ interphon -enlarge "2 2 1" -pbc "1 1 0"
```

-> (2×2×1) supercell and displaced supercells 

-> Periodic boundary conditions (1 or True) along the a1, a2 lattice directions, while open (0 or False) along the a3 direction

- ### Post-process
After DFT force calculations for the displaced supercells (POSCAR-0*) are finished  
in each folder of FORCE-0* (folder names are arbitrary—only the order of the folder numbers is important),  
the evaluation of interfacial phonon can be executed by the following post-process  
where the phonon properties are printed out in the forms of data files and graphics:

#### Density of states (DOS)
```
$ interphon -fc "FORCE-0*/vasprun.xml" -kdos KPOINTS_dos
```
