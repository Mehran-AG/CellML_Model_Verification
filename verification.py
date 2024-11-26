
# importing external or built-in packages
import numpy as np
import sympy as sp
import sys
from colorama import Fore, Back, Style, init
from pprint import pprint

def verification( stoichiometric_array, elemental_array, element_indices, compound_indices, reaction_indices, rate_array = 0, modified = False ):

    """
    This function receives stoichiometric matrix (numpy array), elemental matrix (numpy array), and symbolic parameters as a rates array (sympy array)
    and does some operations to check the model's validity
    """

    init( autoreset=True )

    #np.set_printoptions(threshold=np.inf)
    np.set_printoptions(threshold=np.inf)

    elemental_array = elemental_array.astype(int)
    stoichiometric_array = stoichiometric_array.astype(int)

    elemental_matrix = sp.Matrix(elemental_array)
    stoichiometric_matrix = sp.Matrix(stoichiometric_array)

    if rate_array != 0:

        conservation_matrix = elemental_matrix * stoichiometric_matrix 

        conservation_array = np.array( conservation_matrix )

        if np.all( ( conservation_array == 0 ) ) == True:

            if modified:

                print( Style.BRIGHT + Fore.YELLOW + "\nThere are Phosphates in the reactions which are not linked to any ATP. So, we have added ATP and ADP to the reactions to satisfy mass balance for the reactions\n", )

                input("\nPress Enter to continue...")

                print( Fore.RED + "\nElemental matrix (modified) is:\n", elemental_array )
                print( Fore.YELLOW + "\nStoichiometric matrix (modified) is:\n", stoichiometric_array )

            else:

                print( Fore.RED + "\nElemental matrix is:\n", elemental_array )
                print( Fore.YELLOW + "\nStoichiometric matrix is:\n", stoichiometric_array )

            stoichiometric_matrix_transposed =stoichiometric_matrix.transpose()
            nullspace_transposed = stoichiometric_matrix_transposed.nullspace()
            l = len(nullspace_transposed)
        
            if l == 0:
                print(Fore.CYAN + "There is no Left Null Space for this matrix")
            elif l == 1:
                nullspace = np.transpose(np.array(nullspace_transposed[0]))

                conservation_equations_array = nullspace * rate_array

                if modified:

                    print( Fore.CYAN + "\nThe Left Null Sapce (for modified equations) is:\n", nullspace )
                    
                    print( Style.BRIGHT + Fore.GREEN + "\nCONGRATULATIONS!!! >>>> Your MODIFIED model passed the mass conservation verification test <<<<")
                    print('\nConservation equations (for the modified reactions) are:\n', conservation_equations_array[0], ' = 0\n', conservation_equations_array[1], ' = 0\n' )

                else:
                    
                    print( Fore.CYAN + "\nThe Left Null Sapce is:\n", nullspace )
                    
                    print( Style.BRIGHT + Fore.GREEN + "\nCONGRATULATIONS!!! >>>> Your model passed the mass conservation verification test <<<<")
                    print('\nConservation equations are:\n', conservation_equations_array[0], ' = 0\n', conservation_equations_array[1], ' = 0\n' )
                    
                return { "Pass": "Yes", "items": [conservation_equations_array] }
            else:
                n_t = np.array(nullspace_transposed[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate( ( n_t, nullspace_transposed[counter] ), axis=1 )
                    counter+=1
                nullspace = np.transpose(n_t)

                conservation_equations_array = nullspace * rate_array

                if modified:

                    print( Fore.CYAN + "\nThe Left Null Space (for modified reactions) is:\n", nullspace)
                    
                    print( Style.BRIGHT + Fore.GREEN + "\nCONGRATULATIONS!!! >>>> Your MODIFIED model passed the mass conservation verification test <<<<")
                    print( Fore.MAGENTA + "\nRate of conservation equations (for modified reactions) are:")

                else:

                    print( Fore.CYAN + "\nThe Left Null Space is:\n", nullspace)
                    
                    print( Style.BRIGHT + Fore.GREEN + "\nCONGRATULATIONS!!! >>>> Your model passed the mass conservation verification test <<<<")
                    print( Fore.MAGENTA + "\nRate of conservation equations are:")
                
                count = 0
                while count < l:
                    print( conservation_equations_array[count], ' = 0')
                    count += 1
                print("\n\n")
                return { "Pass": "Yes", "items": [conservation_equations_array] }

        else:
            
            
            length = conservation_array.shape[1]
            for i in  range(0,length):
                if np.all( conservation_array[:, i] == 0 ):
                    pass
                else:

                    for key, value in reaction_indices.items():

                        if value == i:
                            reaction = key
                            reaction_index = value
                            break

                    column = conservation_array[:, i]

                    non_zero_indices = np.nonzero(column)[0]

                    non_zero_values = {}

                    for non_zero_index in non_zero_indices:

                        non_zero_values[non_zero_index] = column[non_zero_index]

                    for species_index in non_zero_indices:

                        for key, value in element_indices.items():

                            if value == species_index:

                                if key == "Pi" or key == "ATP" or key == "ADP":
                                
                                    #print( f"It is ralted to a Phosphate in reaction {reaction}")

                                    compound_stoichio_coefficient = non_zero_values[species_index]

                                    return { "Pass": False, "items": [species_index, reaction_index, compound_stoichio_coefficient] } 

                                else:

                                    print( Style.BRIGHT + Fore.RED + "\nConservation of Mass is violated" )

                                    print( Style.BRIGHT + Fore.CYAN + "\nSpecies", end='' )
                                    print( Style.NORMAL + Fore.YELLOW + f" {key} ", end='')
                                    print( Style.BRIGHT + Fore.CYAN + f"is not conserved in reaction {reaction}" )
                                    break

                    

                    coefficients = stoichiometric_matrix.T [i,:]

                    reactants = []
                    products =[]

                    for j, coeff in enumerate( coefficients ):

                        if coeff != 0:

                            for key, value in compound_indices.items():

                                if value == j:
                                    compound = key
                                    break

                            if coeff < 0:
                                reactants.append( f"{-coeff} * {compound}" )
                            elif coeff > 0:
                                products.append( f"{coeff} * {compound}" )

                    reactants_str = "  +  ".join(reactants)
                    products_str = "  +  ".join(products)

                    reaction_str = f"{reactants_str} --> {products_str}"

                    print( Style.DIM + f"\nReaction {reaction}:    ", end='' )

                    print( Style.BRIGHT + Fore.MAGENTA + "{rr}".format( rr = reaction_str ) )
                    
                    break
            
            sys.exit( Fore.GREEN + "\nModify the equations and run the simulations again to see the figures\n" )

    else:

        conservation_matrix = elemental_matrix * stoichiometric_matrix 

        conservation_array = np.array( conservation_matrix )

        if np.all( ( conservation_array == 0 ) ) == True:

            print( Fore.RED + "\nElemental matrix is:\n", elemental_array )
            print( Fore.YELLOW + "\nStoichiometric matrix is:\n", stoichiometric_array )
        
            stoichiometric_matrix_transposed=stoichiometric_matrix.transpose()
            nullspace_transposed= stoichiometric_matrix_transposed.nullspace()
            l = len(nullspace_transposed)

            if l == 0:
                print('There is no Left Null Space for this matrix')
            elif l == 1:
                nullspace = np.transpose(np.array(nullspace_transposed[0]))
                print('\nThe Left Null Sapce is:\n', nullspace)
                conservation_equations_array = nullspace * rate_array

                return { "Pass": "Yes", "items": [conservation_equations_array] } 
            
            else:
                n_t = np.array(nullspace_transposed[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate((n_t, nullspace_transposed[counter]), axis=1)
                    counter+=1
                nullspace = np.transpose(n_t)
                print('\nThe Left Null Space is:\n', nullspace)
                conservation_equations_array = nullspace * rate_array

                return { "Pass": "Yes", "items": [conservation_equations_array] } 

        else:
            
            length = conservation_array.shape[1]
            for i in  range(0,length):
                if np.all(conservation_array[:, i] == 0):
                    pass
                
                else:

                    for key, value in reaction_indices.items():

                        if value == i:
                            reaction = key
                            reaction_index = value
                            break

                    column = conservation_array[:, i]

                    non_zero_indices = np.nonzero(column)[0]

                    non_zero_values = {}

                    for non_zero_index in non_zero_indices:

                        non_zero_values[non_zero_index] = column[non_zero_index]

                    for species_index in non_zero_indices:

                        for key, value in element_indices.items():

                            if value == species_index:

                                if key == "Pi" or key == "ATP" or key == "ADP":
                                
                                    print( "It is related to a Phosphate")

                                    compound_stoichio_coefficient = non_zero_values[species_index]

                                    return { "Pass": "No", "items": [species_index, reaction_index, compound_stoichio_coefficient] } 

                    print('Conservation of Mass is violated\n')
                    print ('Reaction number %i is not correctly defined' %i)
                    break

            sys.exit("\nModify the equations and run the simulations again to see the figures")