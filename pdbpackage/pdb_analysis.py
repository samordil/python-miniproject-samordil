#!/usr/bin/env python3

"""Analysis script
Usage: python3 pdb_analysis.py
A scripts that does basic PDB file analysis
"""
# Import the required funcitons from pdbmodule
import pdbmodule


# Global variables defination
opt1 = ['1', '2', '3', '4', '5']
opt2 = ['O', 'I', 'H', 'S', 'Q']
current_pdb_filename = 'None'
file_flag = 0
in_memory_file = []

no_file ='{:^80}'.format(
    'NO FILE!: First load a PDB file using option [1 or O]')

# Display user interface
while True:
    # Takes the name of the loaded file in memory 
    # when that operation is successfull
    pdbmodule.displayUI(current_pdb_filename)
    try:
        usr_opt = input(': ') 
    except:
        # exit cleanly on ctrl+C or ctrl+D
        print()
        if pdbmodule.closeProgram(): 
            # Exit the program
            break
        else:
            # Go back to the main menu
            print()
            continue
        
    if usr_opt.upper() in opt1 + opt2:
        if usr_opt == opt1[0] or usr_opt.upper() == opt2[0]:
            try:
                path_to_file = input('Enter a Valid PATH for a PDB File: ')
            except (KeyboardInterrupt, EOFError):
                # exit cleanly on ctrl+C or ctrl+D
                print()
                if pdbmodule.closeProgram(): 
                    # Exit the program
                    break
                else:
                    # Go back to the main menu
                    print()
                    continue
             
            # Accept the user file name
            # Return the name of the successfully loaded pdb file,
            in_memory_file, current_pdb_filename, msg, file_flag = (
                pdbmodule.openPdbFile(in_memory_file, path_to_file, file_flag,
                                      current_pdb_filename))
            # print the message return from the function
            print(msg)
        elif usr_opt == opt1[1] or usr_opt.upper() == opt2[1]:
            if current_pdb_filename == 'None':
                print(no_file)
            else:
                # pass in the name of the loaded file and its content
                # to display summary information
                pdbmodule.pdbFileInfo(current_pdb_filename, in_memory_file)
                
        elif usr_opt == opt1[2] or usr_opt.upper() == opt2[2]:
            if current_pdb_filename == 'None':
                print(no_file)
            else:
                # pass in the name of the loaded file and its content
                # Show histogram of amino acid
                pdbmodule.aaHistogram(in_memory_file)
                
        elif usr_opt == opt1[3] or usr_opt.upper() == opt2[3]:
            if current_pdb_filename == 'None':
                print(no_file)
            else:
                pdbmodule.pdbSecondaryStr(in_memory_file)
                
        elif usr_opt == opt1[4] or usr_opt.upper() == opt2[4]:
             if pdbmodule.closeProgram(): break
    else:
        # centre aligns the error message when printed
        print( '{:^80}'.format(
                  'ATTENTION: Kindly Choose an Option from the Menu.'))