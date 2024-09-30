import numpy as np
import chemparse as chp

from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString
# Importing internal packages
from utilities import print_model
import cellml

import chebi_fetch as chf

import sys






####################################################################################
####################################################################################
####################################################################################
###############----------------- Main -------------------###########################

def cellml_compound_element_sorter ( components ):
    
    """
    CellML reader now gives a dictionary containing variable names as keys and variable ids as values. Variable id is like 'id_chebicode', e.g., id_18654. so we nned to extract the chebicode
    and give to chemparse to parse it into its constituents and return it as a dictionary.
    Then we put these chemicals into compounds dictionary where keys are chemical formula and values are dictionaries containing its constituents and their corresponding quantity
    Elements dictionary which contains indices of elements is also constructed and returned
    In the latest update (April 16, 2024), variable annotation can include the composition of the compound which is given by its name to CellML: 
    'va_mekppmapk-mek.mapk': This annotation containes two parts separated by an underline(_). 2nd part contains the compound name and its composition. The minus sign(-) separates component name and composition.
    So the first part (mekppmapk) is the component name. The 2nd part (mek.mapk) shows that it is composed of 'mek' and 'mapk'.
    """


    compound_to_composition = {}                # This dictionary will store the compounds and their chemical composition

    reaction_indices = {}                       # This dictionary keeps an index for each reaction that is used in constructing stoichiometric matrix ==> {'r1': 0, 'r2': 1, ...} "It can be unordered numbers for reactions or any alphabets ==> {'r10': 0, ...}"

    compound_indices = {}                       # This dictionary keeps an index for each compound to be used as a reference to construct stoichiometric matrix ==> {'C6H12O6': 0, 'CO2': 1, ...} so we know that row 0 in matrix is for C6H12O6

    element_indices ={}                         # This dictionary contains elements as keys and assigns a number for each in order to construct the elemental matrix as values ==> {'C': 0, 'H': 1, 'O': 2}

    compound_index = 0                          # This variable keeps the ordered number for assigning any compound a number

    element_index = 0                           # This index keeps track of the order each new element will be put in the elemental matrix

    reaction_index = 0                          # This index keeps track of the order for reactions to be used as a reference for the column in stoichiometric matrix

    variables = []                              # This list contains the variables of cellml file [In CellML file, I defined two types of parameters, first type is variable which can change in the model and the second type is coefficient which is constant]

    coefficients = []                           # The difference between variables and coefficients is identified by the first part of their id which is like "va_12345_r1", first part indicates that this parameter is a variable

    symbols_list = []                           # This list contains the list of symbols

    bcvirtual_compound_coefficients = {}        # This is a dictionary to store the compound names for boundary conditions and a value showing if the compound is an external or internal one.

    variables, coefficients, _ , _ , _ , boundary_conditions, _ , _ = variable_sorter( components )

    # Going through all variables: 1-to get their chemical composition (Or any name assigned to it in ID) 2-Assign indices to the compounds to assign a row in stoichiometric matrix by the index or a column in the elemental matrix

    for variable in variables:

        name = variable.name()  # This is the name the user has given to a variable in CellML

        chebi_code =  variable.id().split('_')[1]   # A variables id is like "va_12345" or "va_NO-N.O", so the 1st part after splitting it by '_' will be the chebi code or its name

        # Since the variable sometimes does not have a chebi code, the user should write its name and comosition in the id, so we check to see if the id is all_digits or alphabets

        if all_digits( chebi_code ):
            
            # If the chebi code is all digits, so it will be easy for us to find its chemical name and composition
            compound, composition = chf.chebi_comp_parser( chebi_code )

        else:
            
            # If there is no chebi code for the variable, then we have to split the code more
            # A '-' after the name of the compound shows its composition, and each species is separated by a '.'
            if ( len(chebi_code.split('-')) == 1 ):
                compound = chebi_code
                composition={ chebi_code: 1 }
            
            else:

                compound = chebi_code.split('-')[0]
                formula = chebi_code.split('-')[1].split('.')

                composition={}
                for element in formula:
                    composition[element] = 1

        # Storing compound and its corresponding composition in a dictionary if it hasn't been done previously
        if compound not in compound_to_composition:
                
            compound_to_composition[compound] = composition

        # Assigning an index to the compound if it hasn't been done previously
        if compound not in compound_indices:
                
            compound_indices[compound] = compound_index
            compound_index += 1
        
        # Storing the variable's name in CellML (User's input) in the list of symbols
        if name not in symbols_list:

            symbols_list.append( name )

        # To check the conservation of mass, we need an elemental matrix
        # Elements make the rows, and compounds make columns, so we need indices for elements as weel
        # Going through the elements in the chemical compoistion of the compound and assigning an index to each element
        for element in composition:
                        
            if element not in element_indices:

                element_indices[element] = element_index
                element_index += 1


    # Assigning a number index to reactions
    for coefficient in coefficients:
            
        chebi_code = coefficient.id().split('_')[1]

        reaction_number_parts = coefficient.id().split('_')[2]

        for reaction_number_part in reaction_number_parts.split('-'):

            reaction_number = reaction_number_part.split('.')[0]

            if reaction_number not in reaction_indices:
                reaction_indices[reaction_number] = reaction_index
                reaction_index += 1



    # Now I check the boundary conditions
    # Boundary conditions are considered as reactions with single a participant inside and outside of the membrane
    # Therefore, for each boundary condition for a species, I have to create two virtual species, one inside and one outside the membrane
    # I name them as "species_name" + "_i" for the virtual internal species and "species_name" + "_e" for the virtual external species
    for bc in boundary_conditions:

        # The reaction name should be the same as id
        reaction_number = bc.id()

        # if bc.id() == 'bc_p16.CycDCdk4_o.17':
        #     print( bc.id() )

        # I assign an index to the reaction if it is not assigned previously
        if reaction_number not in reaction_indices:
            reaction_indices[reaction_number] = reaction_index
            reaction_index += 1

        # The boundary condition id is "bc_12345.1" where "bc" means it is a boundary condition variable, "12345" is the chebi code or it can be compound's name, and "1" is the enumeration number for multiple boundary conditons for a species
        chebi_codes =  bc.id().split('_')[1].split('.')

        for chebi_code in chebi_codes:

            # If the species has a chebi code, so we need to find its name and check if it exists in the compounds dictionary
            if all_digits( chebi_code ):

                compound, composition = chf.chebi_comp_parser( chebi_code )

                # Here I check to see if it already exists in the compounds dictionary, if it exists, then it's OK, if not, so there is something wrong in the user's input
                try:
                    compound_to_composition[compound]

                except KeyError:
                    print(f"\nThe species '{compound}' with the chebi code '{chebi_code}' does not match any available species.")
                    sys.exit("\nModify the equations and run the simulations again to see the figures\n")    

            else:
                
                try:
                    compound = chebi_code
                    composition = compound_to_composition[compound]

                except KeyError:
                    print(f"The species '{chebi_code}' does not match any available species.")
                    sys.exit("\nModify the equations and run the simulations again to see the figures")
            
            # Now, we will create the virtual external species for the compound
            compound_bc = compound + '_e'

            # Storing compound and its corresponding composition in a dictionary
            if compound_bc not in compound_to_composition:
                    
                compound_to_composition[compound_bc] = composition

            # Assigning an index to the compound
            if compound_bc not in compound_indices:
                    
                compound_indices[compound_bc] = compound_index
                compound_index += 1

            if compound_bc not in symbols_list:

                symbols_list.append( compound_bc )

            # Since we can have more than one boundary condition for a species, we need to assign the virtual species to the correct boundary condtion reaction
            encoded_coefficient_name = compound_bc + '-' + reaction_number

            bcvirtual_compound_coefficients[ encoded_coefficient_name ] = -1

            # We also need to assign the existing variable as the internal compound for the reaction which represents the flow of the compound through the boundary 
            compound_bc = compound

            # Storing compound and its corresponding composition in a dictionary
            if compound_bc not in compound_to_composition:
                    
                compound_to_composition[compound_bc] = composition

            # Assigning an index to the compound
            if compound_bc not in compound_indices:
                    
                compound_indices[compound_bc] = compound_index
                compound_index += 1

            # Since we can have more than one boundary condition for a species, we need to assign the virtual species to the correct boundary condtion reaction
            encoded_coefficient_name = compound_bc + '-' + reaction_number

            bcvirtual_compound_coefficients[ encoded_coefficient_name ] = +1
        

    return element_indices, compound_indices, reaction_indices, symbols_list, compound_to_composition, bcvirtual_compound_coefficients




# ################################################################
# -------------------- Auxiliary functions -----------------------


# This function determines if the characters used in a variable are all digits or alphabets
def all_digits( word ):

    return all( char.isdigit() for char in word )


# This function sorts the variables of the CellMl file into the corresponding lists
def variable_sorter( components, all_vars = 'Y' ):

    """
    This function takes all the components of the CellML file and sorts all the variables of the CellML file into the corresponding lists
    If "all_vars" is selected, all variables will be sorted, if not, only the main variables which are compound concentration variables and stoichiometric coefficients will be sorted and returned
    """

    variables = []

    coefficients = []

    enzymes = []

    rate_constants = []

    rates = []

    boundary_conditions = []

    equation_variables = []

    boundary_values = []

    for component in components:

        number_of_variables = component.variableCount()

        if ( all_vars == "Y" or all_vars == "y" ): 
            
            for v in range( 0, number_of_variables ):
                
                if component.variable(v).id():      # There are some variables in CellML like time (t) that I don't want to import to python,so I left it without id and here I can import those variables with id that I need here

                    id = component.variable(v).id()

                    identifier = id.split('_')[0]

                    if identifier == 'va': variables.append( component.variable(v) )    # Since we have two different types of parameters in CellML, I put them in different lists
                    elif identifier == 'co': coefficients.append( component.variable(v) )
                    elif identifier == 'en': enzymes.append( component.variable(v) )
                    elif identifier == 'rc': rate_constants.append( component.variable(v) )
                    elif identifier == 'ra': rates.append( component.variable(v) )
                    elif identifier == 'bc': boundary_conditions.append( component.variable(v) )
                    elif identifier == 'ev': equation_variables.append( component.variable(v) )
                    elif identifier == 'bv': boundary_values.append( component.variable(v) )

        else:
                    
            for v in range( 0, number_of_variables ):
                
                if component.variable(v).id():      # There are some variables in CellML like time (t) that I don't want to import to python,so I left it without id and here I can import those variables with id that I need here

                    id = component.variable(v).id()

                    identifier = id.split('_')[0]

                    if identifier == 'va': variables.append( component.variable(v) ) 
                    elif identifier == 'co': coefficients.append( component.variable(v) )
                    

    if ( all_vars == "Y" or all_vars == "y" ): 

        return variables, coefficients, enzymes, rates, rate_constants, boundary_conditions, equation_variables, boundary_values
    
    else:

        return variables, coefficients
    


def variable_name_mapper( components ):

    '''
    This function receives the components of a CellML file and maps the chebi codes to the variable name,
    and returns a dictionary mapping chebi code to the CellML variable name
    chebi_toCellML = { 'C2H6O12': x_glu, ... }
    '''

    chebi_to_CellML = {}

    chebi_initialvalues = {}

    variables, _ = variable_sorter( components , 'n')

    for variable in variables:      # After putting the parameters in their corresponding list, I get their ChEBI code and build the matrices
            
        name = variable.name()

        initial_value = variable.initialValue()

        chebi_code =  variable.id().split('_')[1]

        if all_digits( chebi_code ):

            compound, _ = chf.chebi_comp_parser( chebi_code )

        else:

            compound = chebi_code.split('-')[0]

        chebi_to_CellML[compound] = name

    return chebi_to_CellML


def initial_value_finder( components, general_equations ):

    chebi_initialvalues = {}

    variables, _ = variable_sorter( components , 'n')

    for variable in variables:      # After putting the parameters in their corresponding list, I get their ChEBI code and build the matrices

        v_id = variable.id()

        chebi_code =  v_id.split('_')[1]

        if all_digits( chebi_code ):

            compound, _ = chf.chebi_comp_parser( chebi_code )

        else:

            compound = chebi_code.split('-')[0]

        if  ( not variable.initialValue().strip() ) and ( len(v_id.split('_')) < 4 ):

            continue

        if not variable.initialValue().strip():

            name = variable.name()

            rhs = general_equations[name]

            eq_all_symbols = list( rhs.free_symbols )

            for to_find_symbol in eq_all_symbols:

                for to_find_variable in variables:

                    if str( to_find_symbol ) == to_find_variable.name():

                        rhs = rhs.subs( to_find_symbol, to_find_variable.initialValue() )

                        break

            initial_value = str( rhs ).rstrip('0')

        else:

            initial_value = variable.initialValue()
            name = variable.name()


        chebi_initialvalues[compound] = initial_value

    return chebi_initialvalues