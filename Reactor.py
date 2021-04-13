import numpy as np
import Molecules


class Reactor:
    """Contains the reactions type that can occur, polymer-polymer, mon-mon
    polymer-monomer and the useful methods"""

    def __init__(self):
        pass

    def react(self, molecule_a, molecule_b):

        if molecule_a.is_chain and molecule_b.is_chain:
            if molecule_a is molecule_b:
                if molecule_a.nmol >= 2:
                    return self.reaction_poly_poly(molecule_a, molecule_b)
                else:
                    return None
            else:
                return self.reaction_poly_poly(molecule_a, molecule_b)

        if molecule_a.is_chain and not molecule_b.is_chain:
            if molecule_a.nmol >= 1 and molecule_b.nmol >= 1:
                return self.reaction_poly_mono(molecule_a, molecule_b)
            else:
                return None

        if not molecule_a.is_chain and molecule_b.is_chain:
            if molecule_a.nmol >= 1 and molecule_b.nmol >= 1:

                return self.reaction_mono_poly(molecule_a, molecule_b)
            else:
                return None

        if not molecule_a.is_chain and not molecule_b.is_chain:
            if molecule_a is molecule_b:
                if molecule_a.nmol >= 2:
                    return self.reaction_mono_mono(molecule_a, molecule_b)

                else:
                    return None
            else:
                return self.reaction_mono_mono(molecule_a, molecule_b)

    def reaction_poly_poly(self, molecule_a, molecule_b):
        """It updates the number of molecules after the reaction and creates
        a new chain which length is the sum of the reactants chain length and
        it updates the molecule_type"""
        molecule_a.nmol -= 1
        molecule_b.nmol -= 1
        total_chain_length = molecule_a.chain_length + molecule_b.chain_length
        total_mass = molecule_a.mass + molecule_b.mass - 18
        new_type = molecule_b.get_new_type(molecule_a.molecule_type, molecule_b.molecule_type)
        new_molecule_type_hash = molecule_a.molecule_type_string + molecule_b.molecule_type_string
        new_molecule_type_string = molecule_a.molecule_type_string + "-" + molecule_b.molecule_type_string
        return Molecules.Polymers(total_chain_length, new_type, total_mass,  new_molecule_type_hash, new_molecule_type_string, molecule_a.get_first_element(), molecule_b.get_last_element())

    def reaction_poly_mono(self, molecule_a, molecule_b):
        """It updates the number of molecules after the reaction and update
        the polymer chain length and molecule_type"""
        molecule_a.nmol -= 1
        molecule_b.nmol -= 1
        total_chain_length = molecule_a.chain_length + 1
        total_mass = molecule_a.mass + molecule_b.mass - 18
        new_type = molecule_b.get_new_type(molecule_a.molecule_type, [molecule_b])
        new_molecule_type_hash = molecule_a.molecule_type_string + molecule_b.molecule_type_string
        new_molecule_type_string = molecule_a.molecule_type_string + "-" + molecule_b.molecule_type_string
        return Molecules.Polymers(total_chain_length, new_type, total_mass, new_molecule_type_hash, new_molecule_type_string, molecule_a.get_first_element(), molecule_b)

    def reaction_mono_poly(self, molecule_a, molecule_b):
        """It updates the number of molecules after the reaction and update
        the polymer chain length and molecule_type"""
        molecule_a.nmol -= 1
        molecule_b.nmol -= 1
        total_chain_length = molecule_b.chain_length + 1
        total_mass = molecule_a.mass + molecule_b.mass - 18
        new_type = molecule_b.get_new_type([molecule_a], molecule_b.molecule_type)
        new_molecule_type_hash = molecule_a.molecule_type_string + molecule_b.molecule_type_string
        new_molecule_type_string = molecule_a.molecule_type_string + "-" + molecule_b.molecule_type_string
        return Molecules.Polymers(total_chain_length, new_type, total_mass, new_molecule_type_hash, new_molecule_type_string, molecule_a, molecule_b.get_last_element())

    def reaction_mono_mono(self, molecule_a, molecule_b):
        """It updates the number of molecules after the reaction and creates
        a new chain which length is the sum of the reactants chain length and
        it updates the molecule_type"""
        molecule_a.nmol -= 1
        molecule_b.nmol -= 1
        total_chain_length = 2
        total_mass = molecule_a.mass + molecule_b.mass - 18
        new_type = [molecule_a, molecule_b]
        new_molecule_type_hash = molecule_a.molecule_type_string + molecule_b.molecule_type_string
        new_molecule_type_string = molecule_a.molecule_type_string + "-" + molecule_b.molecule_type_string
        return Molecules.Polymers(total_chain_length, new_type, total_mass, new_molecule_type_hash, new_molecule_type_string, molecule_a, molecule_b)

    def choose_t_r(self, totalprob, prob_i):
        """kMC step selecting the time step when the next reaction occurs
        and which reaction will occur"""
        r = np.random.random(2)
        r1 = r[0]
        r2 = r[1]
        dt = (1/totalprob)*np.log(1/r1)
        i = 0
        N = r2*totalprob - prob_i[i]
        while N > 0:
            i = i + 1
            N = N - prob_i[i]

        next_r = i

        return dt, next_r

    def filter_first_reactant(self, reactant_array, molecule_searched):
        """It gets a dictionary with allt the monomers + chain in the system
        and return an array of all the possible chain and monomer involved in
        a reaction e.g. if a GA reaction occurs it creates an array with all
        the chain having all GA monomer and all the chains with GA
        as LAST element"""
        final_array_first = []
        for i in range(len(reactant_array)):
            if reactant_array[i].is_chain and reactant_array[i].nmol >= 1:
                if reactant_array[i].get_last_element() is molecule_searched:
                    final_array_first.append(reactant_array[i])
            elif reactant_array[i] is molecule_searched and reactant_array[i].nmol >= 1:
                final_array_first.append(reactant_array[i])
        return final_array_first

    def filter_second_reactant(self, reactant_array, molecule_searched):
        """It gets a dictionary with allt the monomers + chain in the system
        and return an array of all the possible chain and monomer involved in a reaction
        e.g. if a GA reaction occurs it creates an array with all the chain having all GA monomer
        and all the chains with GA as FIRST element. The function also removes the chain or monomer
        which was selected by the filter_first_reactant module"""
        final_array_second = []
        for i in range(len(reactant_array)):
            if reactant_array[i].is_chain and reactant_array[i].nmol >= 1:
                if reactant_array[i].get_first_element() is molecule_searched:
                    final_array_second.append(reactant_array[i])
            elif reactant_array[i] is molecule_searched and reactant_array[i].nmol >= 1:
                final_array_second.append(reactant_array[i])
        return final_array_second

    def choose_which_molecule(self, reactant_array, molecule_searched_one, molecule_searched_two):
        """selects the two molecule reacting based on the list of all possible molecules"""
        sumnmol1 = 0
        probability1 = []
        sumnmol2 = 0
        probability2 = []
        final_array_first = self.filter_first_reactant(reactant_array, molecule_searched_one)
        for molecule in final_array_first:
            sumnmol1 += molecule.nmol
        for molecule in final_array_first:
            probability1.append(molecule.nmol/(sumnmol1*molecule.chain_length))
        norm1 = sum(probability1)
        probability1[:] = [x / norm1 for x in probability1]
        first_molecule_selected = np.random.choice(final_array_first, p=probability1)
        first_molecule_selected.nmol -= 1
        final_array_second = self.filter_second_reactant(reactant_array, molecule_searched_two)
        for molecule in final_array_second:
            sumnmol2 += molecule.nmol
        for molecule in final_array_second:
            probability2.append(molecule.nmol/(sumnmol2*molecule.chain_length))

        norm2 = sum(probability2)
        probability2[:] = [x / norm2 for x in probability2]
        second_molecule_selected = np.random.choice(final_array_second, p=probability2)
        first_molecule_selected.nmol += 1
        return first_molecule_selected, second_molecule_selected

    def reactant_nmol_first(self, reactant_array, molecule_searched):
        totnmol = 0

        for i in range(len(reactant_array)):
            if reactant_array[i].is_chain:
                if reactant_array[i].nmol >= 1:
                    if reactant_array[i].get_last_element() is molecule_searched:
                        totnmol += reactant_array[i].nmol

        return totnmol

    def reactant_nmol_second(self, reactant_array, molecule_searched):
        totnmol = 0
        for i in range(len(reactant_array)):
            if reactant_array[i].is_chain:
                if reactant_array[i].nmol >= 1:
                    if reactant_array[i].get_first_element() is molecule_searched:
                        totnmol += reactant_array[i].nmol
        return totnmol
