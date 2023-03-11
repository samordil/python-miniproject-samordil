# PythonMiniProject
### Mini-project for the Python for Biologist class at Pwani University


This repository consists of a command-line script written in Python 3.8.10 programming language.

The script extracts information and executes basic analysis on a PDB file that refers to 3D structures of protein. The script also integrates bash tab completion for directories and file names completion when providing a path to a given PDB file.

#### The script's user interface
The script provides the user with a command-line interface containing five options to choose from:
1. Opening and loading of a PDB file to be analyzed.
2. Information about the loaded file such as the name of the loaded PDB file, the title of the protein, number of chains.
3. Histogram of the amino acids from the loaded PDB file.
4. Secondary structure of the protein from the loaded PDB file.
5. Option for quitting the running script.


#### Repository Organization 
The repository consists of a package directory called pdbpackage which contains two files: pdbmodule.py and pdb_analysis.py. All the functions are placed in the module pdbmodule.py and accessed by importing them to pdb_analysis.py which acts as the main program.

#### Repository Organization
pdbpackage (directory) containing:
- pdbmodule.py
- A module containing all the functions relevant to the extraction of information from a PDB file. 
- pdb_analysis.py - The main program that imports pdbmodule.py as a module and utilizes the functions.

#### How to run the script from the terminal inside the package directory
	python3 pdb_analysis.py

#### Script written by: 
Odoyo Samuel Oduor


## Credits
This Project was designed by [Gustavo Salazar]()
