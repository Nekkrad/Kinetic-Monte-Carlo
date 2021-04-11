import pandas as pd 
import numpy as np
#lengths = []
#num_mol = []
#mol_type = []
#mol_mass = []
#ave_mass = []
#ave_length = []
chains_data = []

def get_data(molecule, lengths, num_mol, mol_type, mol_mass):
    ave_mass = 0
    ave_length = 0
        lengths.append(molecule.chain_length)
        num_mol.append(molecule.nmol)
        mol_type.append(molecule.molecule_type_string)
        mol_mass.append(molecule.mass)
        ave_mass = np.average(mol_mass, weights=num_mol)
        ave_length = np.average(lengths, weights=num_mol)
    return lengths, num_mol, mol_type, mol_mass, ave_mass, ave_length

def get_chain(molecule, chains, chains_data):
    if molecule.chain_length <= chains:
        chains_data[molecule.chain_length - 2] = molecule.nmol
    return chains_data
 
def get_chain_dict(chains_data, molecule, t, nstep, GroupDataFrame):

    chains_data.insert(0, molecule.nmol) 
    chains_data.insert(0, t) 
    chains_data.insert(0, nstep) 
    G_series = pd.Series(chains_data, index=GroupDataFrame.columns)
    GroupDataFrame = GroupDataFrame.append(G_series, ignore_index=True)
    return GroupDataFrame

def get_matrix(molecule,row,column, monomer1, monomer2, polymermatrix):
    m1num = 0
    m2num = 0
    if molecule.chain_length <= (row + column):
        m1num = molecule.molecule_type.count(monomer1)
        m2num = molecule.molecule_type.count(monomer2)
        if m1num <= row and m2num <= column:
            polymermatrix[m1num][m2num] +=molecule.nmol
    return polymermatrix

def save_matrix(polymermatrix, name):
    matrix = pd.DataFrame(polymermatrix) 
    matrix.to_csv(name + '.csv')
    return

 #def get_type_and_length(molecule):
 #       tmp_chainlength = finallength
 #       tmp_chaintype = finalmoltype
 #       tmp_chainnumber = finalmol
 #       tmp_chainlength, tmp_chainnumber, tmp_chaintype = zip(*sorted(zip(tmp_chainlength, tmp_chainnumber, tmp_chaintype)))
 #       tmp_chainlength, tmp_chainnumber, tmp_chaintype = (list(t) for t in zip(*sorted(zip(tmp_chainlength, tmp_chainnumber, tmp_chaintype))))
 #       for i in range(2, 41):
 #           try:
 #               index1 = tmp_chainlength.index(i)
 #               index2 = len(tmp_chainlength) - tmp_chainlength[::-1].index(i)
 #               tmp2_chainnumber = tmp_chainnumber[index1:index2+1]
 #               tmp2_chaintype = tmp_chaintype[index1:index2+1]
 #               tmp2_chainnumber, tmp2_chaintype = zip(*sorted(zip(tmp2_chainnumber, tmp2_chaintype)))
 #               tmp2_chainnumber, tmp2_chaintype = (list(t) for t in zip(*sorted(zip(tmp2_chainnumber, tmp2_chaintype))))
#              tmp2_chainnumber = tmp2_chainnumber[::-1]
#                tmp2_chaintype = tmp2_chaintype[::-1]
#                tmp_dict = {"Type": tmp2_chaintype, "Number": tmp2_chainnumber}
#                df = pd.DataFrame(tmp_dict)
 #               dumpval = int(nstep/dump)
 #               df.to_csv(f"chain{i}_step{dumpval}" + ".csv")
 #           except:
  #              pass
        
        
        

