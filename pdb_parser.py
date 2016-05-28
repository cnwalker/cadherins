import os, csv

class PdbParser(object):
    def __init__(self, filepath):
        self.file = open(filepath, 'r').read()
        self.file_name = filepath.split('/')[-1]

    def get_lines(self):
        return self.file.split('\n')

    def count_atoms(self):
        acc = 0
        for line in self.get_lines():
            if line[:4] == 'ATOM':
                acc += 1
        print(self.file_name + ': ' + str(acc))

    def count_amino(self, acid):
        acid_nums = []
        for line in self.get_lines():
            if line[:4] == 'ATOM':
                line = line.split()
                if line[3] == acid and line[5] not in acid_nums:
                    acid_nums.append(line[5])
        print (self.file_name + ' ' + acid + ': ' + str(len(acid_nums)))

    def get_sidechains(self, amino_acids, start_atom, end_atom):
        acc_atoms = False
        sidechain = []
        chain_list = []
        for line in self.get_lines():
            if 'ATOM' == line[:4] and amino_acid in line:
                if start_atom in line:
                    acc_atoms = True

                if acc_atoms:
                    # Get the coordinates
                    print line.split()[6:9]
                    line_list = list(map(lambda x: float(x), line.split()[6:9]))
                    # Append a tuple of the coordinate
                    sidechain.append((line_list[0], line_list[1], line_list[2]))

                if end_atom in line:
                    acc_atoms = False
                    chain_list.append(sidechain)
                    sidechain = []
        return chain_list

    def get_discrete_sidechains(self, amino_acid, atom_list):
        sidechain_len = len(atom_list)
        cur_atom_index = 0
        sidechain = []
        chain_list = []


        for line in self.get_lines():
            if 'ATOM' == line[:4] and amino_acid in line:
                if atom_list[cur_atom_index] in line:
                    if (len(line.split()) <= 10):
                        line_list = list(map(lambda x: float(x), line.split()[5:8]))
                    else:
                        try:
                            line_list = list(map(lambda x: float(x), line.split()[6:9]))
                        except:
                            return chain_list
                    sidechain.append((line_list[0], line_list[1], line_list[2]))
                    if cur_atom_index < sidechain_len - 1:
                        cur_atom_index += 1
                    else:
                        chain_list.append(sidechain)
                        if len(chain_list) > 150:
                            return chain_list
                        sidechain = []
                        cur_atom_index = 0
        return chain_list

    # Retrieves mainchain hydrogen bonds from i + lower_i to i + upper_i
    def get_mainchain_hbonds(self):
        if 'HELIX' in self.file:
            helix_dict = {}
            start_dex = 0
            all_lines = self.get_lines()
            for line in all_lines:
                if 'HELIX' == line[:5]:
                    curline = line.split()
                    helix_dict[curline[5]] = curline[8]
                    start_dex += 1
                elif 'ATOM' == line[:4]:
                    break
                else:
                    start_dex += 1

            # initialize a dictionary to store the frequency of i_3, i_4, and i_5 hydrogen bonds
            hbond_num_dict = {
                'i_plus_3': 0,
                'i_plus_4': 0,
                'i_plus_5': 0
            }

            hbond_candidates = []
            collect_nitrogens = False
            for line in all_lines[start_dex:]:
                if 'ATOM' == line[:4]:
                    curline = line.split()
                    if helix_dict.get(curline[5]):
                        # Not interested in i > 5 h bonds
                        if len(hbond_candidates) > 2:
                            hbond_candidates = []
                        if curline[-1] == 'N':
                            # push nitrogens on the stack
                            hbond_candidates.push(list(map(lambda x: float(x), curline[6:9])))





# Normal functions not nested in an object are here
def compute_from_list(func, in_list):
    return map(lambda x: func(x[0], x[1], x[2], x[3]), in_list)

def make_readers(pdb_dir='/Users/Christopher/Desktop/digital_bio/psets/pdb_files'):
    reader_list = []
    for file in os.listdir(pdb_dir):
        if file.endswith('.pdb'):
            reader_list.append(PdbParser(pdb_dir + '/' + file))
    return reader_list

def run_function_and_output_results(func, amino_acid, atom_list, target_file):
    all_sidechains = reduce(lambda x,y: x + y, map(lambda x: x.get_discrete_sidechains(amino_acid, atom_list),
                                                                                              make_readers()))
    results = map(lambda x: func(x[0], x[1], x[2], x[3]), all_sidechains)
    result_file = open(target_file, 'wb')
    result_writer = csv.writer(result_file)
    for result in results:
        result_writer.writerow([result])
    result_file.close()
