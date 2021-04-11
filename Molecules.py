class Molecules:
    """Contains the information regarding all the monomers and polymers in the polymerization process"""

    def __init__(self):
        pass

    def get_new_type(self, first, second):
        """Return the new type of a chain"""
        tmp_type = []
        for element in first:
            tmp_type.append(element)
        for element in second:
            tmp_type.append(element)
        return tmp_type


class Monomers(Molecules):
    """Regarding all that concerns monomers in this code"""

    def __init__(self, nmol, molecule_type, mass):
        self.nmol = nmol  # Number of molecule
        self.molecule_type = molecule_type  # Monomer type --> GA, PA etc
        self.is_chain = False  # True if is a polymer
        self.mass = mass
        self.chain_length = 1

    def __hash__(self):
        return hash(self.molecule_type_string)

    def __eq__(self, other):
        return (self.molecule_type_string) == (other.molecule_type_string)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def get_last_element(self):
        """Get last element"""
        return self

    def get_first_element(self):
        """Get first element"""
        return self


class Polymers(Molecules):
    """Contains information regarding polymers created in the process  """

    def __init__(self, chain_length, molecule_type, mass, molecule_type_hash, molecule_type_string, first_element, last_element):
        self.chain_length = chain_length  # Chain length value
        self.molecule_type = molecule_type  # molecule type GA-PA-SA-GA-GA as a list
        self.nmol = 1  # Number of chain is fixed at 1 and it gets removed if a reaction between two chains occurs (one chain "disappear")
        self.is_chain = True  # True if it's a polymer
        self.mass = mass
        self.molecule_type_string = molecule_type_string
        self.molecule_type_hash = molecule_type_hash
        self.first_element = first_element
        self.last_element = last_element

    def __str__(self):
        return self.molecule_type_string

    def __hash__(self):
        return hash((self.molecule_type_string))

    def __eq__(self, other):
        return (self.molecule_type_string) == (other.molecule_type_string)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def get_first_element(self):
        """"Get the first element of a polymer chain"""
        return self.first_element

    def get_last_element(self):
        """Get the last element of a polymer chain"""
        return self.last_element


class GA(Monomers):
    """Glycolic acid monomers (GA), the value of rate_constants are the kinetic constats
    they change bases on what's the other molecule GA reacts with"""

    def __init__(self, nmol, rate_constant):
        self.nmol = nmol
        self.molecule_type = self
        self.molecule_type_string = "GA"
        self.rate_constant = rate_constant
        self.is_chain = False  # True if is a polymer
        self.mass = 76.05
        self.chain_length = 1

    def __str__(self):
        return str(self.molecule_type)


class LA(Monomers):
    """DL-Lactic Acid  monomers (LA) """

    def __init__(self, nmol, rate_constant):
        self.nmol = nmol
        self.molecule_type = self
        self.molecule_type_string = "LA"
        self.rate_constant = rate_constant
        self.is_chain = False  # True if is a polymer
        self.mass = 90.08
        self.chain_length = 1

    def __str__(self):
        return str(self.molecule_type)


class MA(Monomers):
    """ DL-2-hydroxy-4-methylpentanoic acid monomers (MA)"""

    def __init__(self, nmol, rate_constant):
        self.nmol = nmol
        self.molecule_type = self
        self.molecule_type_string = "MA"
        self.rate_constant = rate_constant
        self.is_chain = False  # True if is a polymer
        self.mass = 132.16
        self.chain_length = 1

    def __str__(self):
        return str(self.molecule_type)


class SA(Monomers):
    """DL-2-hydroxy-4-(methylsulfanyl)butanoic acid monomers (SA)"""

    def __init__(self, nmol, rate_constant):
        self.nmol = nmol
        self.molecule_type = self
        self.molecule_type_string = "SA"
        self.rate_constant = rate_constant
        self.is_chain = False  # True if is a polymer
        self.mass = 150.2
        self.chain_length = 1

    def __str__(self):
        return str(self.molecule_type)


class PA(Monomers):
    """ DL-3-phenyllactic acid (PA)"""

    def __init__(self, nmol, rate_constant):
        self.nmol = nmol
        self.molecule_type = self
        self.molecule_type_string = "PA"
        self.rate_constant = rate_constant
        self.is_chain = False  # True if is a polymer
        self.mass = 166.17
        self.chain_length = 1

    def __str__(self):
        return str(self.molecule_type)
