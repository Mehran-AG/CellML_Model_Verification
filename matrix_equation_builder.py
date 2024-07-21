
import numpy as np
from sympy import symbols
import sympy as sp
from colorama import Fore, Back, Style, init

import compound_element_sorter as ces
import chebi_fetch as chf

def matrix_equation_builder ( stoichiometric_matrix, rows, columns, reaction_rate_equations_dict, general_equations, components, imported_bc_equations, printing = 'off' ):

    """
    This function creates the concentration rate equations from the stoichiometric matrix that is constructed from the variables in CellML file
    Inputs are: the stoichiometric matrix, rows and columns which are two dictionaries mapping row and column names to their ordering numbers, reaction rate equations, and components of the CellML file as a list
    And the function returns the concentration rate equations
    """

    _ , _ , _ , reaction_rates, _ , boundary_conditions, equation_variables = ces.variable_sorter( components )

    concentration_rate_equations = {}                                                           # This dictionary will map the compound name to the corresponding equation for it

    chebi_initial_values = ces.initial_value_finder( components, general_equations )

    chebi_to_CellML = ces.variable_name_mapper( components )

    ev_variables = []

    for equation_variable in equation_variables:

        ev_variables.append( equation_variable.name() )

    # ------------ << Since each row shows the reaction a compound participates, we go through each row and find the reactions the compound participates
    # 'rows' is a dictionary mapping compound names to their row position in the stoichiometric matrix, so we get the row number by checking the 'rows' dictionary >> --------------------
    for compound, row_number in rows.items():

        try:

            to_find = chebi_to_CellML[ compound ]

        except:

            to_find = compound

        if ( to_find not in ev_variables ) and ( compound in chebi_initial_values ):

            temporary_reactions = {}                                                                # To construct the right hand side of the equations, I need to store reaction name with their stoichiometric coefficient which shows the rate of consumption or production of a variable in a reaction
                                                                                                    # I will multiply this reaction rate to its coefficient later, so I need to keep both of them for later use as a mapping
            
            # ------------ << Now we want to get the value of the element in the stoichiometric matrix for this specific compound to see if it participates in which reaction >> -----------------
            for column_number, element in enumerate( stoichiometric_matrix[row_number] ):

                if element != 0:                                                                    # If the element value is not zero, then it participates in this reaction. We have to get the reaction name here.

                    reaction_name = next( ( key for key, value in columns.items() if value == column_number ), None )
                    temporary_reactions[reaction_name] = element

            rhs = 0                                                                                 # Here, I construct an empty right hand side for the equation

            for reaction in temporary_reactions.keys():                                             # I nned to look for the reaction in the reactions list to find its name since I only have ids of the reaction which might be different with its variable name in CellML

                for reaction_rate in reaction_rates:

                    if reaction == reaction_rate.id().split('_')[1]:                                # Here I find the CellML reaction component and retrieve its variable name
                        
                        rate_symbol = symbols(reaction_rate.name())
                        
                        rhs = rhs + temporary_reactions[reaction] * rate_symbol                     # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                        rhs = rhs.subs( rate_symbol, reaction_rate_equations_dict[reaction_rate.name()] )

                # ----------- << We will go through all boundary conditions to see which one belongs to this compound that the concentration rate equation being written
                for bc in boundary_conditions:

                    chebi_code = bc.id().split('_')[1]                                              # Chebi code stored in the boundaru condition's id

                    # Getting the compound name for the boundary condition since the stoichionetric is built upon the compound names
                    if ces.all_digits( chebi_code ):

                        bc_compound, _ = chf.chebi_comp_parser( chebi_code )

                    else:

                        bc_compound = chebi_code.split('-')[0]

                    # Checking to see if it matches with the compound that its rate is being written
                    # If it matches with the compound, then it will be added to the equation
                    if bc_compound == compound:

                        bc_name = bc.name()

                        rate_symbol = symbols( bc.name() )

                        rhs = rhs + rate_symbol                                                     # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                        if not bc.initialValue():

                            bc_value = None

                            try:

                                imported_bc_equations

                                for eq in imported_bc_equations:

                                    if str(eq.lhs) == bc_name:

                                        bc_value = eq.rhs
                                        break

                            except:

                                print("There is no initial value for boundary condition {bc_condition}, and no imported equation is found for it.".format( bc_condition = bc_name ))
                                exit()

                        else:

                            bc_value = bc.initialValue()                            # The value of the boundary condition is stored in a variable

                        rhs = rhs.subs( rate_symbol, bc_value )

            # Since there are rows for boundary conditions in the stoichionetric matrix, the concentration rate equation for these species will be zero, so we try not to write these equations
            if rhs != 0:

                concentration_rate_equations[compound] = rhs

    

    if printing == 'on' or printing =='On' or printing == 'ON':
        
        printer( concentration_rate_equations, 'Concentration rate equations generated by Stoichiometric Matrix:' )

    return concentration_rate_equations
        



def printer( equations, description ):

    # Initialize colorama
    init(autoreset=True)

    print(Fore.GREEN + "\n {d} \n                \u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193".format(d=description))

    #print('\n', description, '\n                \u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193')
    for compound in equations.keys():

        lhs = 'd[' + compound + ']/dt'

        rhs = equations[compound]

        print( Style.BRIGHT + Fore.RED + "{c}".format( c = lhs ), end='')
        print( Style.BRIGHT + " = ", end='' )
        print( Style.BRIGHT + Fore.BLUE + "{rh}".format( rh = rhs ) )

    print("**********************************************************************")