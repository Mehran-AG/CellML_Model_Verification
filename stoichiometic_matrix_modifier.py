import numpy as np


def stoichiometric_matrix_modifier( modified_stoichiometric_matrix, modified_elemental_matrix, modified_compound_indices, modified_element_indices, reaction_indices, species_index, reaction_index, compound_stoichio_coefficient ):

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

        if "ATP" not in modified_compound_indices:

            largest_index = max( modified_compound_indices.values() )

            modified_compound_indices["ATP"] = largest_index + 1

            modified_stoichiometric_matrix = np.vstack ( [ modified_stoichiometric_matrix, np.zeros( modified_stoichiometric_matrix.shape[1] ) ] )

            zero_column = np.zeros((modified_elemental_matrix.shape[0], 1))

            modified_elemental_matrix = np.hstack([modified_elemental_matrix, zero_column])

        if "ADP" not in modified_element_indices:

            largest_index = max( modified_element_indices.values() )

            modified_element_indices["ADP"] = largest_index + 1

            modified_elemental_matrix = np.vstack ( [ modified_elemental_matrix, np.zeros( modified_elemental_matrix.shape[1] ) ] )

        




    return