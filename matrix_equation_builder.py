
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

    _ , _ , _ , reaction_rates, _ , boundary_conditions, equation_variables, _ = ces.variable_sorter( components )

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

            already_added_BC = []

            for reaction in temporary_reactions.keys():                                             # I nned to look for the reaction in the reactions list to find its name since I only have ids of the reaction which might be different with its variable name in CellML

                reaction_and_parts = reaction.split('-')

                previous_reactions_to_check = []

                for reaction_and_part in reaction_and_parts:

                    should_continue = False

                    if len( reaction_and_part.split('.') ) > 1:

                        reaction_to_check = reaction_and_part.split('.')[0]

                        for previous_reaction_to_check in previous_reactions_to_check:

                            if previous_reaction_to_check == reaction_to_check:

                                should_continue = True

                                break

                        if should_continue:
                            continue

                        for reaction_rate in reaction_rates:

                            if reaction_to_check == reaction_rate.id().split('_')[1]:                                # Here I find the CellML reaction component and retrieve its variable name
                                
                                rate_symbol = symbols(reaction_rate.name())
                                
                                rhs = rhs + temporary_reactions[reaction] * rate_symbol                     # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                                previous_reactions_to_check.append( reaction_to_check )

                                try:

                                    rhs = rhs.subs( rate_symbol, reaction_rate_equations_dict[reaction_rate.name()] )

                                except:

                                    rhs = rhs.subs( rate_symbol, 0 )

                    else:
                        
                        reaction_to_check = reaction_and_part

                        for reaction_rate in reaction_rates:

                            if reaction_to_check == reaction_rate.id().split('_')[1]:                                # Here I find the CellML reaction component and retrieve its variable name
                                
                                rate_symbol = symbols(reaction_rate.name())
                                
                                rhs = rhs + temporary_reactions[reaction] * rate_symbol                     # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                                try:

                                    rhs = rhs.subs( rate_symbol, reaction_rate_equations_dict[reaction_rate.name()] )

                                except:

                                    rhs = rhs.subs( rate_symbol, 0 )

                # ----------- << We will go through all boundary conditions to see which one belongs to this compound that the concentration rate equation being written
                for bc in boundary_conditions:

                    chebi_codes = bc.id().split('_')[1].split('.')                                              # Chebi code stored in the boundary condition's id

                    for chebi_code in chebi_codes:

                        bc_id = bc.id()

                        fake_bc_id = None
                        
                        if len( bc_id.split('_') ) > 2:

                            fake_bc_id = chebi_code + bc.id().split('_')[2]

                        # Getting the compound name for the boundary condition since the stoichionetric is built upon the compound names
                        if ces.all_digits( chebi_code ):

                            bc_compound, _ = chf.chebi_comp_parser( chebi_code )

                        else:

                            bc_compound = chebi_code

                        # Checking to see if it matches with the compound that its rate is being written
                        # If it matches with the compound, then it will be added to the equation
                        if bc_compound == compound and fake_bc_id not in already_added_BC:

                            already_added_BC.append( fake_bc_id )

                            bc_name = bc.name()

                            rate_symbol = symbols( bc.name() )

                            bc_in_out_sign = bc_id.split('_')[2].split('.')[0]

                            if bc_in_out_sign == 'i':

                                rhs = rhs + ( 1 * rate_symbol )     # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                            elif bc_in_out_sign == 'o':

                                rhs = rhs - ( 1 * rate_symbol )

                            else:

                                print("There is no identifier of whether the flow is in or out for boundary condition {bc_condition}".format( bc_condition = bc_name ) )
                                exit()
                            
                            #rhs += rate_symbol                           

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

    temporary_bc_verifier = []

    for i, boundary_condition in enumerate( boundary_conditions ):

        bc_id = boundary_condition.id()

        bc_chebis = bc_id.split('_')[1].split('.')

        for bc_chebi in bc_chebis:

            if ces.all_digits( bc_chebi ):                          # if the chebi code aprt of the id is all digits, so it is considered as the chebi code, otherwise it shold be the compound name and composition given by the user

                compound, _ = chf.chebi_comp_parser( bc_chebi )

            else:

                compound = bc_chebi

            if compound not in concentration_rate_equations.keys() or compound in temporary_bc_verifier:

                temporary_bc_verifier.append( compound )

                boundary_condition_name = boundary_condition.name()

                boundary_condition_symbol = symbols( boundary_condition_name )

                #rhs = boundary_condition_symbol

                bc_in_out_sign = bc_id.split('_')[2].split('.')[0]

                if bc_in_out_sign == 'i':

                    rhs = +1 * boundary_condition_symbol

                elif bc_in_out_sign == 'o':

                    rhs = -1 * boundary_condition_symbol

                else:

                    print("There is no identifier of whether the flow is in or out for boundary condition {bc_condition}".format( bc_condition = bc_name ) )
                    exit()
                

                bc_value = None

                if not boundary_condition.initialValue():

                    try:

                        imported_bc_equations

                        for eq in imported_bc_equations:

                            if str(eq.lhs) == boundary_condition_name:

                                bc_value = eq.rhs
                                break

                    except:

                        print("There is no initial value for boundary condition {bc_condition}, and no imported equation is found for it.".format( bc_condition = bc_name ))
                        exit()

                else:

                    bc_value = abs( boundary_condition.initialValue() )                            # The value of the boundary condition is stored in a variable

                rhs = rhs.subs( boundary_condition_symbol, bc_value )

                for j, another_boundary_condition in enumerate( boundary_conditions ):

                    another_bc_id = another_boundary_condition.id()

                    another_bc_chebi = another_bc_id.split('_')[1]

                    if ces.all_digits( another_bc_chebi ):                          # if the chebi code aprt of the id is all digits, so it is considered as the chebi code, otherwise it shold be the compound name and composition given by the user

                        another_compound, _ = chf.chebi_comp_parser( another_bc_chebi )

                    else:

                        another_compound = another_bc_chebi.split('-')[0]

                    if j != i and another_compound == compound:

                        another_boundary_condition_name = another_boundary_condition.name()

                        another_boundary_condition_symbol = symbols( another_boundary_condition_name )

                        rhs += another_boundary_condition_symbol

                        bc_value = None

                        if not another_boundary_condition.initialValue():

                            try:

                                imported_bc_equations

                                for eq in imported_bc_equations:

                                    if str(eq.lhs) == another_boundary_condition_name:

                                        bc_value = eq.rhs
                                        break

                            except:

                                print("There is no initial value for boundary condition {bc_condition}, and no imported equation is found for it.".format( bc_condition = bc_name ))
                                exit()

                        else:

                            bc_value = another_boundary_condition.initialValue()                            # The value of the boundary condition is stored in a variable

                        rhs = rhs.subs( another_boundary_condition_symbol, bc_value )




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