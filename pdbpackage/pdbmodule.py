#!/usr/bin/env python3

"""Generic PDB file operations.

A module that contains functions to extract information  
and execute basic anaylsis on a PDB file.
"""
# Import the neccessary modules
import os
import sys
import textwrap
import readline
import glob

def pathCompleter(text, state):
    """ 
    This is a tab completer for *nix systems paths.
    """
    # replace ~ with the user's home dir. 
    if '~' in text:
        text = os.path.expanduser('~')

    # autocomplete directories with having a trailing slash
    if os.path.isdir(text):
        text += '/'
               
    return (glob.glob(text+'*')+[None])[state]


readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(pathCompleter)


def displayUI(current_pdb_file='None'):
    """
    Shows a display of the main window which gives 
    several options. The user can choose an option
    by typing either a letter or a number related 
    to the option they want.
    """
    ast = '*'
    horiz_line = ast * 80
    newline = '{}{:>79}'.format(ast, ast)
    print()
    print(horiz_line)
    print('{} PDB FILE ANALYZER {:>60}'.format(ast, ast))
    print(horiz_line)
    print('{:} Select an option from below: {:>49}'.format(ast, ast))
    print(newline)
    print('{} {:>23} {:>25} {:>28}'.format(ast, '1) Open a PDB File', '(O)', ast))
    print('{} {:>19} {:>29} {:>28}'.format(ast, '2) Information', '(I)', ast))
    print('{} {:>36} {:>12} {:>28}'.format(ast, '3) Show histogram of amino acid', '(H)', ast))
    print('{} {:>35} {:>13} {:>28}'.format(ast, '4) Display secondary structure', '(S)', ast))
    print('{} {:>12} {:>36} {:>28}'.format(ast, '5) Exit', '(Q)', ast))
    print(newline)
    print('{} {:>76} {}'.format(ast,  'Current PDB: ' + current_pdb_file, ast))
    print(horiz_line)



def openPdbFile(mem, file_path='None', flag=0, current_loaded_file='None'):
    """
    Accepts PDB file to be analyzed and loads it into the
    memory.
    Returns two [pdb_file_name, msg, flag]object which are both
    strings.
    """
    if current_loaded_file == 'None': # Check to see if any file is loaded
        pdb_file_name = None # change the 'None' string to None
    else:
        pdb_file_name = current_loaded_file # Record the current file in mem
        
    loaded_file_status = flag
    msg = ''
    # Reload previous file if still in memory
    in_memory_file = mem
    
    try:
        #Incase of '~' expand to user's home dir 
        expanded_path = os.path.expanduser(file_path)
        base_filename = os.path.basename(expanded_path)
        
        # Raise an exception if file not a pdb file
        if not pdbFileVerifier(expanded_path):
            raise FileNotPdbError
        
        with open(expanded_path) as fh_handle:
            # If the no pdf file loaded and status is false, assign a file name
            if not pdb_file_name and not loaded_file_status: 
                # Get only the name of the file without its path
                pdb_file_name = base_filename
                # Load file to memory and store in a list
                in_memory_file = fh_handle.readlines()
                loaded_file_status = 1     # update file status to true
                msg = 'The File {} has been successfully loaded'.format(pdb_file_name) 
            else:
                # if the pdb file holder not empty ask whether to overide
                while True:
                    try:
                        respond = input('Replace the currently loaded file? [Y/n]: ')
                    except:
                        print()
                        break
                    # check for a valid response of Y/N
                    if respond.upper() in ['Y', 'N']:
                        if respond.upper() == 'Y':
                            # Overide the current file name in the pdf file name holder
                            pdb_file_name = base_filename
                            # Load file to memory and store in a list
                            in_memory_file = fh_handle.readlines()
                            msg = 'The File {} has been successfully loaded'.format(pdb_file_name) 
                            break   # exit the while the while loop
                        else:
                                # When response is No, continue with the current loaded 
                            break       # exit the while loop
    except:
        # centre align the error message
        msg = '{} {:^80}'.format('\n', 
                    'ATTENTION: No such file or Not a Valid PDB File')
        if loaded_file_status: pass # in memory file? yes. continue using it
        else:
            pdb_file_name = 'None'
            
    return(in_memory_file, pdb_file_name, msg, loaded_file_status)


 
def pdbFileInfo(filename, filecontent):
    """
    Displays a summary of a PDB file which base_filename, title,
    amino acid sequence and chain).
    Takes the name of a loaded file and its content to 
    display summary information
    """ 
    file_title =  extractPdbTitle(filecontent)  
    print('PDB File: {}'.format(filename))
    print('Title: {}'.format(file_title))
    # Returns a list of chain names in a pdb file
    total_chains, chain_names_str = extractPdbChain(filecontent)
    print(total_chains)
    
    # Iterate through the chain names
    for chain_ltr in chain_names_str: 
        # returns a dictionary
        amino_len_dict = chainSeqLen(filecontent, chain_ltr) 
        
        # Get the number of helix for every chain
        helix_count = helixCount(filecontent, chain_ltr)
        
        # Get the number of sheet for every chain
        sheet_count = sheetCount(filecontent, chain_ltr)
        
        # Get amino acid sequence in one letter sequence
        amino_seq = aminoSeqExtractor(filecontent, chain_ltr)
        
        # Displaying the summary Interface for a pdb file.
        print(' - Chain {}'.format(chain_ltr))
        print('    Number of amino acids:{:>6}'.format(
                            amino_len_dict[chain_ltr]))
        print('    Number of helix:{:>12}'.format(helix_count))
        print('    Number of sheet:{:>12}'.format(sheet_count))
        aa_seq = textwrap.fill(amino_seq, width=50)
        print('    Sequence:  {}'.format(
                    textwrap.indent(text=aa_seq,
                    # Ignore the first line while identing 
                    predicate=lambda line_ident: not aa_seq.splitlines()[0] in line_ident,
                    # prefix tab=8 and 7 spaces 
                    prefix='\t       ')))
    print()
    
    
    
def extractPdbTitle(content):
    """
    Scans through the TITLE section of a PDB file to 
    extract protein title.
    """
    file_title = ''
      
    for line in content:
        # stop iteration once no line starts with 'COMPD'
        if line.startswith('COMPND'): break
        if line.startswith('TITLE'):
            # remove the word 'TITLE'
           title1 = line.split(maxsplit=1)[1]
           
           # Checking from title 2 which starts with a digit
           if title1.split(maxsplit=1)[0].isdigit():
               file_title += title1.split(maxsplit=1)[1].strip()
           else:
               file_title += title1.strip()
    
    return(textwrap.fill(file_title, width=72))



def extractPdbChain(content):
    """
    Looks for the available chain names in an open pdb file
    and returns the chain names in a string.
    """
    chain_names = []
    total_chains = ''
    
    for line in content:
        # stop iteration once no line starts with 'SOURCE'
        if line.startswith('SOURCE'): break
        if line.startswith('COMPND'):
            # remove the word 'COMPND'
            c = line.split(maxsplit=1)[1].split()
            # Select the lines in which second word is 'CHAIN' 
            if c[1] == 'CHAIN:':
                # Extract the chain(s) name
                chain_names.extend(c[2:])
    
    # Sort the chain letters
    chain_names.sort()  
    # Convert the list to a string          
    chain_names_str = ''.join(chain_names)
    # Replace any ; and , with nothing
    chain_names_str = chain_names_str.replace(',', '')
    chain_names_str = chain_names_str.replace(';', '')
            
    if len(chain_names_str) == 1:
        total_chains = 'CHAIN: {}'.format(chain_names_str[0])
        
    elif len(chain_names_str) == 2:
        total_chains = 'CHAINS: {} and {}'.format(
                chain_names_str[0], chain_names_str[1])
        
    elif len(chain_names_str) > 2:
        subchain = ','.join(chain_names_str[:-1])
        total_chains = 'CHAINS: {} and {}'.format(
                subchain, chain_names_str[-1])
        
    else:
        total_chains = 'CHAIN: {}'.format('NO CHAIN FOUND!')
        
    return(total_chains, chain_names_str)


    
def chainSeqLen (file_content, chain_letter):
    """
    Takes the content of a pdb file and iterates throught 
    the SEQRES section to get the total amino length (integer)
    of each protein chain in the file.
    Returns a dictionary contains chain name as key and amino
    length as the value: eg {"A" -> 167, 'B' -> 10} 
    """
    # Holds the length of every chain in the file 
    aa_length = {}
    
    for line in file_content:
        # stop iteration once no line starts with 'SEQRES'
        if line.startswith('HET'): break
        if line.startswith('SEQRES'):
            # get the third item of the line -> chain letter
            if chain_letter == line.split()[2]:
                # Check the for existence of a chain_letter
                if aa_length.get(chain_letter):
                    # Ignore the chain_letter if it exists
                    continue
                else:
                    # get the forth item from the string and
                    # add to the aa_length dictionary
                    aa_length[chain_letter] = line.split()[3]
    return(aa_length)



def helixCount(filecontent, chain_letter):
    """
    Counts the number of helix for every chain in a pdb file
    and returns the number (int)
    """
    helix_count = 0
    for line in filecontent:
        # Stop iteration once no line starts with 'SQRES'
        if line.startswith('LINK'): break
        if line.startswith('HELIX'):
            # Check the chain letter and record the count
            if chain_letter == line.split()[4]:
                helix_count += 1
    return(helix_count)



def sheetCount(filecontent, chain_letter):
    """
    Counts the number of sheet for every chain in a pdb file
    and returns the number (int)
    """
    sheet_count = 0
    for line in filecontent:
        # Stop iteration once no line starts with 'SQRES'
        if line.startswith('LINK'): break
        if line.startswith('SHEET'):
            # Check the chain letter and record the count
            if chain_letter == line.split()[5]:
                sheet_count += 1
    return(sheet_count)
    


def aminoSeqExtractor(filecontent, chain_ltr):
    """
    Extracts amino acid sequence from the SEQRE line and convert
    the three letter amino sequence to one amino acid letter
    then return the string.
    """
    code_map = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    amino_seq_lst = []
    amino_seq_str = ''
    
    for line in filecontent:
        # stop iteration once no line starts with 'SEQRES'
        if line.startswith('HET'): break
        if line.startswith('SEQRES'):
            # Remove newlines characters
            line = line.strip()
            # check the chain and extract the amino acid sequence
            if chain_ltr == line.split()[2]:
                if len(line.split()[4]) == 3:
                    amino_seq_lst.extend(line.split()[4:])
                else:
                    amino_seq_lst.append('NO SEQUENCE')

    # Convert the three letter amino names to one letter 
    for aa_3_letter in amino_seq_lst:
        if len(aa_3_letter) == 3:
            amino_seq_str += code_map[aa_3_letter]
        else:
            amino_seq_str = 'NO SEQUENCE'        

    return(amino_seq_str)
         
         
         
def aaHistogram(filecontent):
    """
    Displays the frequency of occurence of a given amino
    acid in a sequence.
    """
    available_opts = ['an', 'dn', 'aa', 'da']
    aa_dict = aminoAcidDict(filecontent)
    ast  = '*'
    
    while True:
        print('{}'.format('Choose an option to order by:'))
        print('   {}'.format('number of amino acids - ascending  (an)'))
        print('   {}'.format('number of amino acids - descending (dn)'))
        print('   {}'.format('alphabetically - ascending         (aa)'))
        print('   {}'.format('alphabetically - descending        (da)'))
        
        try:
            user_opt = input('order by: ')    
        except (KeyboardInterrupt, EOFError):
            print() # print a new line
            print()
            break # go back to the main menu on ctrl+C or ctrl+D
        
        print()
        
        if user_opt.lower() not in available_opts:
            # Give another chance for the user to input correct option
            continue
        else:
            if user_opt.lower() == available_opts[0]:
                # Arranged according to number of amino acid in
                # ascending order
                for key, value in sorted(
                                        aa_dict.items(), 
                                        key=lambda dict_item: dict_item[1]):
                    print('{} ({:>3}) : {}'.format(key, value, ast*value))
                print()

                break
            elif user_opt.lower() == available_opts[1]:
                # Arranged according to number of amino acid in
                # descending order
                for key, value in sorted(
                                        aa_dict.items(), 
                                        key=lambda dict_item: dict_item[1],
                                        reverse=True):
                    print('{} ({:>3}) : {}'.format(key, value, ast*value))
                print()
                break
            elif user_opt.lower() == available_opts[2]:
                # Arranged alphabetically in ascending order
                for key, value in sorted(aa_dict.items()):
                    print('{} ({:>3}) : {}'.format(key, value, ast*value))
                print()
                break
            elif user_opt.lower() == available_opts[3]:
                # Arranged alphabetically in descending order
                for key, value in sorted(aa_dict.items(), reverse=True):
                    print('{} ({:>3}) : {}'.format(key, value, ast*value))
                print()
                break
            
            
    
def aminoAcidDict (filecontent):
    """
    Takes a pdb file and extracts the amino acid sequence.
    Returns a count of each amino acid as a dictionary.
    """
    amino_seq_lst = []
    amino_seq_dct = {}
    
    for line in filecontent:
        # stop iteration once no line starts with 'SEQRES'
        if line.startswith('HET'): break
        if line.startswith('SEQRES'):
            # Remove newlines characters
            line = line.strip()
            # check the chain and extract the amino acid sequence
            if len(line.split()[4]) == 3:
                amino_seq_lst.extend(line.split()[4:])
                
    for aa in amino_seq_lst:
        if amino_seq_dct.get(aa):
            amino_seq_dct[aa] += 1
        else:
            amino_seq_dct[aa] = 1 
            
    return(amino_seq_dct)



def pdbSecondaryStr(filecontent):
    """
    Display a representation of the secondary structure
    for each chain of a loaded PDB file.
    """
    # Get the last item of HEADER line which is the PDB id.
    pdb_id = filecontent[0].split()[-1]
    print('Secondary Structure of the PDB id {}:'.format(pdb_id))
    
     # Get a list of all chain names in the pdb file
    _, chain_names_str = extractPdbChain(filecontent)    
    
    # Iterate through the chain names in the pdb file
    for chain_ltr in chain_names_str:
        # Get the amino acid sequence for the current chain name
        aa_seq = aminoSeqExtractor(filecontent, chain_ltr)
        
        print('Chain {}:'.format(chain_ltr))
        if aa_seq.startswith('NO'):
            print('({})'.format(0))
            print('{}'.format(aa_seq))
            print('({})'.format(0))
            print()
        else:
            # get sysmbol representing helix structure for the current 
            # chain name (letter).
            helix_symbol, chain_id = helixSymbolInserter(
                                        aa_seq, chain_ltr, 
                                        filecontent)
            
            helix_symbol = textwrap.wrap(helix_symbol, width=80)
            
            chain_id = textwrap.wrap(chain_id, width=80)
            
            print('({})'.format(1))
            # Get the amino acid in a list of 80 chars lines 
            aa_lines = textwrap.wrap(text=aa_seq, width=80)
            
            # Loop over the amino sequence of length 80
            for i, line in enumerate(aa_lines):
                print('{}'.format(line))
                print('{}'.format(helix_symbol[i]))
                print('{}'.format(chain_id[i].replace('*', ' ')))
                print()
                
            print('({})'.format(len(aa_seq)))
            print()

    
def helixSymbolInserter (aa_seq, chain_letter, filecontent):
    """
    Identifies all the positions of helix in a sequnce
    and marks those position with a '/' symbol
    """
    # produce a string of length same as amino acid seq
    # containing '-' then converts to list  
    aa_seq_lst = list('-' * len(aa_seq))
    aa_seq_id_lst = list('*' * len(aa_seq))
    
    for line in filecontent:
        # Stop iteration once no line starts with 'SEQRES'
        if line.startswith('LINK') or line.startswith('ATOM'): break
        
        if line.startswith('HELIX'):
            # Check the chain letter and record the count
            if chain_letter == line.split()[4]:
                # Beginning of helix
                start_pos = int(line.split()[5])
                # End of helix 
                end_pos = int(line.split()[8])
                # chain id
                chain_id = line.split()[1]
                                                
                # Get the index of every character in the list
                n = 1
                for i, _ in enumerate(aa_seq_lst, start=n):
                    if i >= (start_pos) and i <= (end_pos):
                        aa_seq_lst[i-n] = '/'
                        if i == int(start_pos): 
                            if len(chain_id) == n:
                                # Replace two elements from the list
                                aa_seq_id_lst[i-n] = chain_id 
                            elif len(chain_id) == 2:
                                aa_seq_id_lst[i-n:i+n] = chain_id
                            elif len(chain_id) == 3:
                                aa_seq_id_lst[i-n:i+2] = chain_id    
                    
    # passing the list for insertion of sheet symbols                
    aa_seq_lst, aa_seq_id_lst = sheetSymbolInserter(
                            aa_seq_lst, chain_letter, 
                            filecontent, aa_seq_id_lst)
    # Chain identifier
    aa_seq_str = ''.join(aa_seq_lst)
    
    aa_seq_id_str = ''.join(aa_seq_id_lst)
    
    return(aa_seq_str, aa_seq_id_str)



def sheetSymbolInserter (aa_seq_lst, chain_letter, filecontent, aa_seq_id_lst):
    """
    Identifies all the positions of sheet in a sequnce
    and marks those position with a '|' symbol
    Depends on a list processed by helixSymbolInserter
    """
    for line in filecontent:
        # Stop iteration once no line starts with 'SEQRES'
        if line.startswith('LINK') or line.startswith('ATOM'): break
        
        if line.startswith('SHEET'):
            # Check the chain letter
            if chain_letter == line.split()[5]:
                # Beginning of sheet
                start_pos = int(line.split()[6])
                # End of sheet 
                end_pos = int(line.split()[9])
                
                # Extract chain identifier
                part_a = line.split()[1]
                part_b = line.split()[2]
                chain_id = part_a + part_b
                                                
                # Get the index of every character in the list
                n = 1
                for i, _ in enumerate(aa_seq_lst, start=n):
                    if i >= (start_pos) and i <= (end_pos):
                        aa_seq_lst[i-n] = '|'
                        if i == int(start_pos):    
                            if len(chain_id) == n:
                                # Replace elements from the list
                                aa_seq_id_lst[i-n] = chain_id 
                            elif len(chain_id) == 2:
                                aa_seq_id_lst[i-n:i+n] = chain_id
                            elif len(chain_id) == 3:
                                aa_seq_id_lst[i-n:i+2] = chain_id
                            elif len(chain_id) == 4:
                                aa_seq_id_lst[i-n:i+3] = chain_id
                            elif len(chain_id) == 5:
                                aa_seq_id_lst[i-n:i+4] = chain_id
                            elif len(chain_id) == 6:
                                aa_seq_id_lst[i-n:i+5] = chain_id

    return(aa_seq_lst, aa_seq_id_lst)
    

def pdbFileVerifier(path_to_file):
    """
    Takes a file and confirms whether it is a pdb file
    or not.
    If it is a pdb file, returns true otherwise returns
    false
    """
    with open(path_to_file) as fh:
        header_line = fh.readline() # Get the first line
        for last_line in fh: pass # get the last line and save in variable end
        
        if header_line.startswith('HEADER') and last_line.startswith('END'):
            return True
        else:
            return False
        
        
        
def closeProgram():
    """
    Asks the user to confirm whether they want to terminate the
    running program or go back to the main menu.
    """
    quiz = 'Do you want to exit (E) or do you want go back to the menu (M):'

    while True:
        try:
            answer = input(quiz)
        except:
            print()
            sys.exit() # halt the program
        if answer.upper() == 'E':
            return(True) 
        elif answer.upper() == 'M':
            return(False) 

# A class to deal with a file that is not pdb        
class FileNotPdbError(Exception): pass