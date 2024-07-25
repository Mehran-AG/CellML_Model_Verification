
# Importing external or built-in packages
import numpy as np
import sys


# Importing internal packages
import compound_element_sorter as ces
import chebi_fetch as chf




def stoichiometric_matrix_builder( reaction_indices, compound_indices, coefficients, bc_coefficients ):

    """
    This function receives a list containing a numpy array of compounds and their corresponding stoichiometric coefficients
    And returns a dictionary containing compounds and their indices showing which row they are assigned in the stoichiometric matrix
    And also returns the stoichiometic matrix (numpy array)
    """

    columns = len( reaction_indices )                                               # Number of columns will be equal to number of reactions in the Stoichiometric matrix

    rows = len( compound_indices )                                                  # Number of rows will be equal to number of compounds in the Stoichiometric matrix

    stoichiometric_matrix = np.zeros(( rows, columns), dtype = int)                 # A numpy matrix is made to construct the stoichiometric matrix


    # to construct the stoichiometric matrix, we need to go through all coeffcients for all compounds
    for coefficient in coefficients:

        chebi_code = coefficient.id().split('_')[1]                                 # At first we get the compound chebi code which this coeffcient belongs to

        reaction_parts = coefficient.id().split('_')[2].split('-')                            # Then we get the reaction which this coeffcient is applied for this compound
        
        already_added = []

        for reaction_part in reaction_parts:

            reaction_number = reaction_part.split('.')[0]

            if reaction_number in already_added:

                continue
            
            already_added.append( reaction_number )

            # We construct the stoichiometric matrix based on the chemical formulas of the compounds not the variable name that the user has given to it
            if ces.all_digits( chebi_code ):

                formula = chf.chebi_formula( chebi_code )                               # If the chebi code part is all digits, it is considered as chebi code, otherwise it is the name of the compound

            else:

                formula = chebi_code

            # Now we will try getting the row index for this compound. We have stored indices for the compounds as we were reading all of them from the CellML file
            try:

                row = compound_indices[formula]

            except KeyError as KE:

                print("TThere is not a variable for compound {v} in CellMl while there is a stoichiometric coefficient defined for it.".format( v = KE ))
                sys.exit("Exiting due to an error\nModify CellML file and check to see if you have this compound or errors in a ChEBI code in reaction {r}".format( r = reaction_number))


            column = reaction_indices[reaction_number]                                  # Now we try to get the column for this compound which is the reaction in which this compound participates with this coefficient

            stoichiometric_matrix[ row, column ] = coefficient.initialValue()           # Then we store the value of the coefficient in the stoichiometric matrix



    # Constructing the rows and columns for the boundary conditions in the stoichiometric matrix
    for bc_coefficient in bc_coefficients.keys():
        
        compound = bc_coefficient.split('-')[0]                                      # bc_coefficient is like 'NO_e-bc_12345.1' and we split it in two by '-' and the one we get is virtual compound name: such as 'NO_e', 'No_i'

        reaction_number = bc_coefficient.split('-')[1]                               # Second part after being split is considered as the reaction's specific name

        # Now we will get the row in the stoichiometric matrix for this virtual compound
        try:

            row = compound_indices[compound]

        except KeyError as KE:

            print("There is not a variable for compound {v} in CellMl while there is a boundary condition defined for it.".format( v = KE ))
            sys.exit("Exiting due to an error\nModify CellML file and check to see if you have this compound or errors in a ChEBI code in reaction {r}".format( r = reaction_number))

        # here we get the column number for this virtual reaction
        column = reaction_indices[reaction_number]
        
        stoichiometric_matrix[ row, column ] = bc_coefficients[bc_coefficient]
        

    return stoichiometric_matrix