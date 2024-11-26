import numpy as np
import sympy as sp


def stoichiometric_matrix_modifier( modified_stoichiometric_matrix, modified_elemental_matrix, modified_compound_indices, modified_element_indices, reaction_indices, species_index, reaction_index, compound_stoichio_coefficient, modified_rate_matrix ):

    for key, value in modified_element_indices.items():
                
        if value == species_index:
                    
            species = key
            break

    if species == "Pi":

        if "ADP" not in modified_compound_indices:

            largest_index = max( modified_compound_indices.values() )

            modified_compound_indices["ADP"] = largest_index + 1

            modified_stoichiometric_matrix = np.vstack ( [ modified_stoichiometric_matrix, np.zeros( modified_stoichiometric_matrix.shape[1] ) ] )

            zero_column = np.zeros((modified_elemental_matrix.shape[0], 1))

            modified_elemental_matrix = np.hstack([modified_elemental_matrix, zero_column])

            new_symbol = sp.symbols('ADP')
            modified_rate_matrix = modified_rate_matrix.row_insert(modified_rate_matrix.rows, sp.Matrix([new_symbol]))

        if "ADP" not in modified_element_indices:

            largest_index = max( modified_element_indices.values() )

            modified_element_indices["ADP"] = largest_index + 1

            modified_elemental_matrix = np.vstack ( [ modified_elemental_matrix, np.zeros( modified_elemental_matrix.shape[1] ) ] )

            row = modified_element_indices["ADP"]

            column = modified_compound_indices["ADP"]

            modified_elemental_matrix[row][column] = 1

        if "ATP" not in modified_compound_indices:

            largest_index = max( modified_compound_indices.values() )

            modified_compound_indices["ATP"] = largest_index + 1

            modified_stoichiometric_matrix = np.vstack ( [ modified_stoichiometric_matrix, np.zeros( modified_stoichiometric_matrix.shape[1] ) ] )

            zero_column = np.zeros((modified_elemental_matrix.shape[0], 1))

            modified_elemental_matrix = np.hstack([modified_elemental_matrix, zero_column])

            new_symbol = sp.symbols('ATP')
            modified_rate_matrix = modified_rate_matrix.row_insert(modified_rate_matrix.rows, sp.Matrix([new_symbol]))

            if "ADP" in modified_compound_indices:

                column = modified_compound_indices["ATP"]

                row = modified_element_indices["ADP"]

                modified_elemental_matrix[row][column] = 1

                row = modified_element_indices["Pi"]

                modified_elemental_matrix[row][column] = 1


        '''
            If the stoichiometric coefficient for a compound like Pi is positive, it means that the compound is a product.
            To balance the mass on both sides of the reaction, we need to add "ADP" to the side that there is "Pi"
            Then we need to add "ATP" to the other side of the reaction
        '''

        if compound_stoichio_coefficient > 0:

            column = reaction_index

            row = modified_compound_indices["ADP"]

            modified_stoichiometric_matrix[row][column] = 1

            row = modified_compound_indices["ATP"]

            modified_stoichiometric_matrix[row][column] = -1

        elif compound_stoichio_coefficient < 0:

            column = reaction_index

            row = modified_compound_indices["ADP"]

            modified_stoichiometric_matrix[row][column] = -1

            row = modified_compound_indices["ATP"]

            modified_stoichiometric_matrix[row][column] = 1


    return modified_stoichiometric_matrix, modified_elemental_matrix, modified_compound_indices, modified_element_indices, modified_rate_matrix