## Script by Christopher Walker and Ryan McDowell
## Digital Biology - CMSC 27610
import math, pdb_parser

class Residue(object):
    # Residue object can be either an n-terminus or a dehydron
    # Both instances have attributes of length and constituent atoms
    def __init__(self, atom_list, length):
        self.atoms = atom_list
        self.length = length

def get_residues(residue_path):
    # Input: String
    # Output: List of lists
    # Function will return a list of residues, represented as a list of residue objects
    chain_num = "0"
    begin_chain_num = 1
    atom_list = []
    residue_list = []
    parser = pdb_parser.PdbParser(residue_path)
    for line in parser.get_lines():
        if line[:4] == "ATOM":
            line_list = line.split()
            if chain_num != line_list[5] and chain_num != "0":
                residue_list.append(Residue(atom_list, int(line_list[1]) - begin_chain_num))
                atom_list = []
                begin_chain_num = int(line_list[1])
            atom_list.append([line_list[6], line_list[7], line_list[8]])
            chain_num = line_list[5]
    return residue_list

def get_center(atom_list, length):
    center = []
    for i in xrange(len(atom_list[0])): # All atoms in atom list must have same number of dimensions
        center.append(sum(map(lambda x: float(x[i]), atom_list))/length)
    return center

def get_distance(n_term, dehydron):
    # Input: (List, List)
    # Output: Float
    distance = 0
    nterm_center = get_center(n_term.atoms, n_term.length)
    dehydron_center = get_center(dehydron.atoms, dehydron.length)
    # n-termini and dehydrons must have same number of dimensions
    for i in xrange(len(nterm_center)):
        distance = (dehydron_center[i] - nterm_center[i])**2
    return math.sqrt(distance)

def get_dehydron_and_nterm_distance(dehydron_file, nterm_file):
    # Input: (String, String)
    # Output: Dictionary of dictionaries
    # Function will compute the distances between a pdb file full of dehydrons
    # and a pdb file full of ntermini. We'll get the distance between each dehydron
    # and n-termini returned in a dictionary of the form
    #  nterm_<n-terminus_num>: {
    #                     n_term_size: size_of_nterminus,
    #                     dehydron_distances: [dehydron_distance_1, dehydron_distance_2, ...]
    #                  }
    # Where <n-terminus_num> is the number the n-terminus occured in the list

    nterm_list = get_residues(nterm_file)
    dehydron_list = get_residues(dehydron_file)

    nterm_info = {}

    for i in xrange(len(nterm_list)):
        nterm_info['nterm_' + str(i)] = {
            'size': nterm_list[i].length,
            'dehydron_distances': []
        }
        for dehydron in dehydron_list:
            nterm_info['nterm_' + str(i)]['dehydron_distances'].append(get_distance(nterm_list[i], dehydron))
        nterm_info['nterm_' + str(i)]['dehydron_distances'].sort()
    return nterm_info
