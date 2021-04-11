import pandas as pd 
import numpy as np

lengths = []
num_mol = []
mol_type = []
mol_mass = []
ave_mass = []
ave_length = []
chains_data = []

def get_data(molecule, lengths, num_mol, mol_type, mol_mass, ave_mass, ave_length):

    
    if molecule.is_chain and molecule.nmol != 0:
        lengths.append(Current_species[key].chain_length)
        num_mol.append(Current_species[key].nmol)
        mol_type.append(Current_species[key].molecule_type_string)
        mol_mass.append(Current_species[key].mass)
        ave_mass = np.average(mol_mass, weights=num_mol)
        ave_length = np.average(lengths, weights=num_mol)
    return lengths, num_mol, mol_type, mol_mass, ave_mass, ave_length


def get_chain(molecule,chains):
    if molecule.is_chain and molecule.nmol != 0:
        if molecule.chain_length <= chains:
            chains_data[molecule.chain_length - 2] = molecule.nmol
    return chains_data
 
def get_chain_dict(chains_data, molecule, t, nstep):

    chains_data.insert(0, molecule.nmol) 
    chains_data.insert(0, t) 
    chains_data.insert(0, nstep) 
    G_series = pd.Series(chains_data, index=GroupDataFrame.columns)
    GroupDataFrame = GroupDataFrame.append(G_series, ignore_index=True)
    return GroupDataFrame

def get_matrix(molecule,length,row,column):
    if Current_species[key].chain_length <= length:
        m1num = Current_species[key].molecule_type.count(GA)
        m2num = Current_species[key].molecule_type.count(LA)
        if m1num <= row and m2num <= column:
            polymermatrix[m1num][m2num] += Current_species[key].nmol
    return polymermatrix

def save_matrix(polymermatrix, name):
    matrix = pd.DataFrame(polymermatrix) 
    matrix.to_csv(name + '.csv')
    return

 def get_type_and_length(molecule)       
        tmp_chainlength = finallength
        tmp_chaintype = finalmoltype
        tmp_chainnumber = finalmol
        tmp_chainlength, tmp_chainnumber, tmp_chaintype = zip(*sorted(zip(tmp_chainlength, tmp_chainnumber, tmp_chaintype)))
        tmp_chainlength, tmp_chainnumber, tmp_chaintype = (list(t) for t in zip(*sorted(zip(tmp_chainlength, tmp_chainnumber, tmp_chaintype))))
        for i in range(2, 41):
            try:
                index1 = tmp_chainlength.index(i)
                index2 = len(tmp_chainlength) - tmp_chainlength[::-1].index(i)
                tmp2_chainnumber = tmp_chainnumber[index1:index2+1]
                tmp2_chaintype = tmp_chaintype[index1:index2+1]
                tmp2_chainnumber, tmp2_chaintype = zip(*sorted(zip(tmp2_chainnumber, tmp2_chaintype)))
                tmp2_chainnumber, tmp2_chaintype = (list(t) for t in zip(*sorted(zip(tmp2_chainnumber, tmp2_chaintype))))
                tmp2_chainnumber = tmp2_chainnumber[::-1]
                tmp2_chaintype = tmp2_chaintype[::-1]
                tmp_dict = {"Type": tmp2_chaintype, "Number": tmp2_chainnumber}
                df = pd.DataFrame(tmp_dict)
                dumpval = int(nstep/dump)
                df.to_csv(f"chain{i}_step{dumpval}" + ".csv")
            except:
                pass
        
        
        

