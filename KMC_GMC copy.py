# This is a Kinetic Monte Carlo program to simulate step-growth polymerization
#              of esters based on Gillespie Algorithm
# JOURNAL OF COMPUTATIONAL PHYSICS 22, 403-434(1976) (Gillespie algorithm)


# ==============================================================================
#                Importing the modules needed for the program
# =============================================================================

import numpy as np
import pandas as pd
import Molecules
import Reactor
import time
from utils import *
reactor = Reactor.Reactor()
# ==============================================================================
#             Setting up System variable such as Volume etc..
# ==============================================================================

fr = open("KMC.log", "a")  # LOG OUTPUT

fr.truncate(0)

N_A = 6*10**(23)  # Avogadro's number

seed = 123456  # random seed

iseed = 0  # Put 1 if you want the same seed (useful to test the code)

t_final = 20  # Final time

Vol = 10**(-17)  # System Volume in L

totstep = 1000000  # Final numbers of reactions

dump = 5000  # Dumps data into final dictionary every dump step.

ngroups = 300  # Longest polymer chain we're interested into

matrixrows = 21  # Number of Monomer  to be included in the polymermatrix

matrixcolumns = 21  # Number of Monomer  to be included in the polymermatrix

# range of matrix values i want to study
frm = 2 
to = 41

conc = {
        "GA": 0.5,
        "LA": 0.5,
        "MA": 0,
        "SA": 0,
        "PA": 0,
        }

GA_rate_constant = {
    "GA": 2*0.5/(Vol*N_A),
    "LA": 100/(Vol*N_A),
    "MA": 1/(Vol*N_A),
    "SA": 1/(Vol*N_A),
    "PA": 1/(Vol*N_A)
    }

LA_rate_constant = {
    "GA": 100/(Vol*N_A),
    "LA": 2*0.5/(Vol*N_A),
    "MA": 1/(Vol*N_A),
    "SA": 1/(Vol*N_A),
    "PA": 1/(Vol*N_A)
    }

MA_rate_constant = {
    "GA": 1/(Vol*N_A),
    "LA": 1/(Vol*N_A),
    "MA": 2*1/(Vol*N_A),
    "SA": 1/(Vol*N_A),
    "PA": 1/(Vol*N_A)
    }

SA_rate_constant = {
    "GA": 1/(Vol*N_A),
    "LA": 1/(Vol*N_A),
    "MA": 1/(Vol*N_A),
    "SA": 2*1/(Vol*N_A),
    "PA": 1/(Vol*N_A)
    }

PA_rate_constant = {
    "GA": 1/(Vol*N_A),
    "LA": 1/(Vol*N_A),
    "MA": 1/(Vol*N_A),
    "SA": 1/(Vol*N_A),
    "PA": 2*1/(Vol*N_A)
    }


# ==============================================================================
# Variables and constants which do not need to be changed
# ==============================================================================
start = time.time() # start time to calculate process time

t = 0  # Starting time

dt = 0  # step time

nstep = 0  # Starting value of reactions step

GA_conversion = 0

LA_conversion = 0

MA_conversion = 0

SA_conversion = 0

PA_conversion = 0
#========================================================================================
#                         ARRAY, DICTIONARY, LIST
#========================================================================================

reaction_i = []  # list containing all the possible couple of reaction

Current_species = {}  # Dictionary containing all the species during simulation

prob_i = np.zeros(25)  # List containing all the reaction probability

mw_list = []
mn_list = []
Total_chain_group_data = {}  # Dict containing values for only ngroups chains at the end of the sim

GroupData = {}  # Dict containing values for only ngroups chains at ith step

Complete_dataset = {}  # Dict containing all the data

Timestep = {}  # Dict containing ith step data

# List which contains information for polymers of up to ngroups length
Grouplength = list(range(2, ngroups+1))

Grouplength.insert(0, "GA")

Grouplength.insert(0, "Time")

Grouplength.insert(0, "nstep")

GroupDataFrame = pd.DataFrame(columns=Grouplength)
# ==============================================================================
#          Setting up the initial number of molecules and creating instances
# ==============================================================================
GA_in_mol = int(conc["GA"]*Vol*N_A)
LA_in_mol = int(conc["LA"]*Vol*N_A)
MA_in_mol = int(conc["MA"]*Vol*N_A)
SA_in_mol = int(conc["SA"]*Vol*N_A)
PA_in_mol = int(conc["PA"]*Vol*N_A)

GA = Molecules.GA(GA_in_mol, GA_rate_constant)
LA = Molecules.LA(LA_in_mol, LA_rate_constant)
MA = Molecules.MA(MA_in_mol, MA_rate_constant)
SA = Molecules.SA(SA_in_mol, SA_rate_constant)
PA = Molecules.PA(PA_in_mol, PA_rate_constant)

# ==============================================================================
# Setting up the species in the species dictionary and list the reaction
# ==============================================================================
Current_species[GA] = GA
Current_species[LA] = LA
Current_species[MA] = MA
Current_species[SA] = SA
Current_species[PA] = PA


reaction_i.append((GA, GA))
reaction_i.append((GA, LA))
reaction_i.append((GA, MA))
reaction_i.append((GA, SA))
reaction_i.append((GA, PA))
reaction_i.append((LA, GA))
reaction_i.append((LA, LA))
reaction_i.append((LA, MA))
reaction_i.append((LA, SA))
reaction_i.append((LA, PA))
reaction_i.append((MA, GA))
reaction_i.append((MA, LA))
reaction_i.append((MA, MA))
reaction_i.append((MA, SA))
reaction_i.append((MA, PA))
reaction_i.append((SA, GA))
reaction_i.append((SA, LA))
reaction_i.append((SA, MA))
reaction_i.append((SA, SA))
reaction_i.append((SA, PA))
reaction_i.append((PA, GA))
reaction_i.append((PA, LA))
reaction_i.append((PA, MA))
reaction_i.append((PA, SA))
reaction_i.append((PA, PA))


# ==============================================================================
# STARTING WITH THE MAIN LOOP AND THE SIMULATIONS
# ==============================================================================
if iseed == 1:
    np.random.seed(seed)


while nstep < totstep:
    
    reactant_array = []  # list of all the species in the simulations
    
    lengths = []  # list of polymers chain lenght
    
    num_mol = []  # number of molecules per polymeric chain
    
    mol_type = []  # list of strings with polymer structure
    
    mol_mass = []  # list of polymers mass
    
    chains_data = [0] * (ngroups-1)  # list of 0 as big as the longest polymer you want to study (from 2)
    
    polymermatrix = np.zeros(shape=(matrixrows, matrixcolumns)) # Matrix containing polymers
    
    for key in Current_species.keys():
        if Current_species[key].nmol > 0:
            reactant_array.append(key)

    # Select the number of chain with X (GA,LA,...etc) as last element
    GA_first_chain_nmol = reactor.reactant_nmol_first(reactant_array, GA)
    LA_first_chain_nmol = reactor.reactant_nmol_first(reactant_array, LA)
    SA_first_chain_nmol = reactor.reactant_nmol_first(reactant_array, SA)
    MA_first_chain_nmol = reactor.reactant_nmol_first(reactant_array, MA)
    PA_first_chain_nmol = reactor.reactant_nmol_first(reactant_array, PA)

    # Select the number of chain with X (GA,LA,...etc) as first element
    GA_second_chain_nmol = reactor.reactant_nmol_second(reactant_array, GA)
    LA_second_chain_nmol = reactor.reactant_nmol_second(reactant_array, LA)
    SA_second_chain_nmol = reactor.reactant_nmol_second(reactant_array, SA)
    MA_second_chain_nmol = reactor.reactant_nmol_second(reactant_array, MA)
    PA_second_chain_nmol = reactor.reactant_nmol_second(reactant_array, PA)

    # CALCULATION OF REACTION PROBABILITY

    # GA + GA -----> GA-GA
    prob_i[0] = GA.rate_constant["GA"] * ((GA.nmol) * (GA.nmol - 1) + (GA_first_chain_nmol + GA_second_chain_nmol))

    # GA + LA -----> GA-LA
    prob_i[1] = GA.rate_constant["LA"]*(GA.nmol + GA_first_chain_nmol)*(LA.nmol + LA_second_chain_nmol)

    # GA + MA -----> GA-MA
    prob_i[2] = GA.rate_constant["MA"]*(GA.nmol + GA_first_chain_nmol)*(MA.nmol + MA_second_chain_nmol)

    # GA + SA -----> GA-SA
    prob_i[3] = GA.rate_constant["SA"]*(GA.nmol + GA_first_chain_nmol)*(SA.nmol + SA_second_chain_nmol)

    # GA + PA -----> GA-PA
    prob_i[4] = GA.rate_constant["PA"]*(GA.nmol + GA_first_chain_nmol)*(PA.nmol + PA_second_chain_nmol)

    # LA + GA -----> LA-GA
    prob_i[5] = LA.rate_constant["GA"]*(LA.nmol + LA_first_chain_nmol)*(GA.nmol + GA_second_chain_nmol)

    # LA + LA -----> LA-LA
    prob_i[6] = LA.rate_constant["LA"]*((LA.nmol)*(LA.nmol - 1) + (LA_first_chain_nmol + LA_second_chain_nmol))

    # LA + MA -----> LA-MA
    prob_i[7] = LA.rate_constant["MA"]*(LA.nmol + LA_first_chain_nmol)*(MA.nmol + MA_second_chain_nmol)

    # LA + SA -----> LA-SA
    prob_i[8] = LA.rate_constant["SA"]*(LA.nmol + LA_first_chain_nmol)*(SA.nmol + SA_second_chain_nmol)

    # LA + PA -----> LA-PA
    prob_i[9] = LA.rate_constant["PA"]*(LA.nmol + LA_first_chain_nmol)*(PA.nmol + PA_second_chain_nmol)

    # MA + GA -----> MA-GA
    prob_i[10] = MA.rate_constant["GA"]*(MA.nmol + MA_first_chain_nmol)*(MA.nmol + MA_second_chain_nmol)

    # MA + LA -----> MA-LA
    prob_i[11] = MA.rate_constant["LA"]*(MA.nmol + MA_first_chain_nmol)*(GA.nmol + GA_second_chain_nmol)

    # MA + MA -----> MA-MA
    prob_i[12] = MA.rate_constant["MA"]*((MA.nmol)*(MA.nmol - 1) + (MA_first_chain_nmol + MA_second_chain_nmol))

    # MA + SA -----> MA-SA
    prob_i[13] = MA.rate_constant["SA"]*(MA.nmol + MA_first_chain_nmol)*(SA.nmol + SA_second_chain_nmol)

    # MA + PA -----> MA-PA
    prob_i[14] = MA.rate_constant["PA"]*(MA.nmol + MA_first_chain_nmol)*(SA.nmol + SA_second_chain_nmol)

    # SA + GA -----> SA-GA
    prob_i[15] = SA.rate_constant["GA"]*(SA.nmol + SA_first_chain_nmol)*(GA.nmol + GA_second_chain_nmol)

    # SA + LA -----> SA-LA
    prob_i[16] = SA.rate_constant["LA"]*(SA.nmol + SA_first_chain_nmol)*(LA.nmol + LA_second_chain_nmol)

    # SA + MA -----> SA-MA
    prob_i[17] = SA.rate_constant["MA"]*(SA.nmol + SA_first_chain_nmol)*(MA.nmol + MA_second_chain_nmol)

    # SA + SA -----> SA-SA
    prob_i[18] = SA.rate_constant["SA"]*((SA.nmol)*(SA.nmol - 1) + (SA_first_chain_nmol + SA_second_chain_nmol))

    # SA + PA -----> SA-PA
    prob_i[19] = MA.rate_constant["PA"]*(SA.nmol + SA_first_chain_nmol)*(PA.nmol + PA_second_chain_nmol)

    # PA + GA -----> PA-GA
    prob_i[20] = PA.rate_constant["GA"]*(PA.nmol + PA_first_chain_nmol)*(GA.nmol + GA_second_chain_nmol)

    # PA + LA -----> PA-LA
    prob_i[21] = PA.rate_constant["LA"]*(PA.nmol + PA_first_chain_nmol)*(LA.nmol + LA_second_chain_nmol)

    # PA + MA -----> PA-MA
    prob_i[22] = PA.rate_constant["MA"]*(PA.nmol + PA_first_chain_nmol)*(MA.nmol + MA_second_chain_nmol)

    # PA + SA -----> PA-SA
    prob_i[23] = PA.rate_constant["SA"]*(PA.nmol + PA_first_chain_nmol)*(SA.nmol + SA_second_chain_nmol)

    # PA + PA -----> PA-PA
    prob_i[24] = PA.rate_constant["PA"]*((PA.nmol)*(PA.nmol - 1) + (PA_first_chain_nmol + PA_second_chain_nmol))

    totalprob = np.sum(prob_i)

    # selects time of the next reaction, index of the reaction that occurs
    dt, next_r = reactor.choose_t_r(totalprob, prob_i)
    
    # Chooses what molecule type is reactiong (GA or LA or SA ..)
    molecule_searched_one, molecule_searched_two = reaction_i[next_r]
    
    # Chooses which molecules react
    reactant_one, reactant_two = reactor.choose_which_molecule(reactant_array, molecule_searched_one, molecule_searched_two)

    # reaction occurs and updates the values
    product = reactor.react(reactant_one, reactant_two)
    
    # updating dictionary of the species
    if product is not None:
        if product not in Current_species:
            Current_species[product] = product
        else:
            Current_species[product].nmol += 1

    # updating iterators
    t += dt
    nstep = nstep + 1

    
    # Updating Final Data Dictionary
    if (nstep % dump) == 0 or nstep == 1:
        total_mass = 0  # to calculate Mw and Mn
        total_mass_squared = 0  # to calculate Mw and Mn
        total_num = 0 # to calculate Mw and Mn
        for key in Current_species:
            if Current_species[key].is_chain and Current_species[key].nmol != 0: 
                # get mass and number of mol
                total_mass += Current_species[key].mass*Current_species[key].nmol
                total_mass_squared += Current_species[key].mass**2 * Current_species[key].nmol
                total_num += Current_species[key].nmol
                lengths1, num_mol1, mol_type1, mol_mass1, ave_mass1, ave_length1 = get_data(Current_species[key], lengths, num_mol, mol_type, mol_mass)
                chains_data1 = get_chain(Current_species[key], ngroups, chains_data)
                polymer_matrix = get_matrix(Current_species[key],matrixrows,matrixcolumns, GA, LA, polymermatrix)
                get_type_and_length(lengths1, num_mol1, mol_type1, frm, to, nstep, dump)
        save_matrix(polymer_matrix, name = str(nstep/dump))
        chain_df = get_chain_dict(chains_data1, GA, t, nstep, GroupDataFrame)
        Mn,Mw = get_Mn_Mw(total_mass,total_num,total_mass_squared)
        mw_list.append(Mw)
        mn_list.append(Mn)
    # Calculating Conversion
        try:
            GA_conversion = GA.nmol/GA_in_mol
            LA_conversion = LA.nmol/LA_in_mol
            MA_conversion = MA.nmol/MA_in_mol
            SA_conversion = SA.nmol/SA_in_mol
            PA_conversion = PA.nmol/PA_in_mol
        except ZeroDivisionError:
            pass

        Timestep = {
           "Nstep": nstep,
           "Time": t,
           # "GA": GA.nmol,
           # "LA": LA.nmol,
           # "MA": MA.nmol,
           # "SA": SA.nmol,
           # "PA": PA.nmol,
           "Chain length": lengths1,
           "Chain number": num_mol1,
           "Chain Type": len(mol_type1),
           "Chain Value": mol_type1,
           "Chain Mass": mol_mass1,
           "Average Mass": ave_mass1,
           "Average Length": ave_length1,
           "GA conversion": GA_conversion,
           # "LA conversion": LA_conversion,
           # "MA conversion": MA_conversion,
           # "SA conversion": SA_conversion,
           # "PA conversion": PA_conversion,
            }
        Complete_dataset[f"Step {nstep}"] = Timestep
        inter_dataframe = pd.DataFrame.from_dict(Timestep)
        inter_dataframe.to_csv(f"Data_{int(nstep/dump)}" + ".csv")


# Creating Pandas Data frame and csv
Complete_datasetframe = pd.DataFrame.from_dict(Complete_dataset, orient="index")
Complete_datasetframe.to_csv("Data.csv")
chain_df.to_csv("Polymersgroup.csv", sep=' ', index=False)
Mw_Mn_dataframe = pd.DataFrame(list(zip(mn_list,mw_list)), columns=["Mn","Mw"])
Mw_Mn_dataframe.to_csv("Mw_Mn.csv",sep = " ")
end = time.time()
delta = (end - start)/60
fr.write("The program run for :")
fr.write(str(delta))
fr.close()
