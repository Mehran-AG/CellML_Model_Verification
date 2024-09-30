# importing external or built-in packages
from sympy import symbols, lambdify, Symbol, Eq, Basic, Float
import numpy as np
import scipy.integrate
import sys
from colorama import Fore, Back, Style, init

# importing internal packages
from compound_element_sorter import variable_name_mapper, initial_value_finder, variable_sorter, all_digits
from matrix_equation_builder import printer
import chebi_fetch as chf

def external_species_concentrations( components, imported_bc_equations, solutions, t_f, delta_t, sympy_to_CellML ):

    '''
    This function gets the conservation equations, which is the change of the concentration of species in the simulated time, and checks if the conservation equations are violated
    '''

    # ------------ << Calling sorter function to sort the variables into their corresponding lists >> --------------
    variables, coefficients, enzymes, reaction_rates, reaction_rate_constants, boundary_conditions, equation_variables, boundary_values = variable_sorter( components )

    external_species_eqs = {}

    # -------------------- << Constructing the virtual compounds and reaction for the boundary conditions >> ----------------------
    for bc in boundary_conditions:

        bc_name = bc.name()                                     # For each boundary condition, I will assign a name to distinguish from each other

        bc_id = bc.id()

        if bc_id.split('_')[2].split('.')[0] == 'i':

            bc_sign = -1

        elif bc_id.split('_')[2].split('.')[0] == 'o':

            bc_sign = +1

        else:

            print( "There is no identifier for boundary {b}!".format( b = bc_name ) )
            exit( "Please add an identifier to the boundary condition {b} so that the direction of flow can be identified".format( b = bc_name ) )

        bc_chebi = bc.id().split('_')[1].split('.')[0]

        # --------- { Now I want to get the name of the compound, it can be a } --------
        if all_digits( bc_chebi ):                          # if the chebi code aprt of the id is all digits, so it is considered as the chebi code, otherwise it shold be the compound name and composition given by the user

            compound, _ = chf.chebi_comp_parser( bc_chebi )

        else:

            compound = bc_chebi.split('-')[0]

        bc_value = None

        if not bc.initialValue():

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


        #bc_display_name = 'v_' + compound                       # The name of the boundary condition to be displayed is composed of 'V' showing flow and the compound's name

        if compound in external_species_eqs.keys():

            rhs = external_species_eqs[ compound ]

            if bc_value is not None:

                try:

                    rhs_new = float( bc_value )

                except:

                    rhs_new = bc_value

            else:

                print("There is no value or equation defined for the boundary condition named {bc_condition}".format( bc_condition = bc_name ) )
                exit("Modify the equations and define a value or an equation for {bc_condition}".format( bc_condition = bc_name ) )

            rhs += bc_sign * rhs_new

            if rhs == 0 or rhs is None:
                print( "Right hand side of {equ} is zero".format( equ = compound ) )

            external_species_eqs[ compound ] = rhs

        else:

            if bc_value is not None:

                try:

                    rhs = bc_sign * float( bc_value )

                except:

                    rhs = bc_sign * bc_value

            else:

                print("There is no value or equation defined for the boundary condition named {bc_condition}".format( bc_condition = bc_name ) )
                exit("Modify the equations and define a value or an equation for {bc_condition}".format( bc_condition = bc_name ) )

            if rhs == 0 or rhs is None:
                print( "Right hand side of {equ} is zero".format( equ = compound ) )

            external_species_eqs[ compound ] = rhs

    chebi_to_CellML = variable_name_mapper( components )

    all_symbols = set().union( *[eq.free_symbols for eq in external_species_eqs.values() if hasattr(eq, 'free_symbols')] )

    compounds_of_equations = [ symbols(chebi_to_CellML[key]) for key in external_species_eqs.keys() ]

    not_state_variables = list( all_symbols - set( compounds_of_equations ) )

    for variable_to_find in not_state_variables:
                
        for key, equation_to_check in external_species_eqs.items():

            if hasattr(equation_to_check, 'free_symbols'):

                if variable_to_find in equation_to_check.free_symbols:

                    for variable_with_value in reaction_rate_constants:

                        if variable_with_value.name() == str( variable_to_find ):

                            initial_value = variable_with_value.initialValue()

                            external_species_eqs[key] = equation_to_check.subs( variable_to_find, initial_value )

                    for variable_with_value in boundary_values:

                        if variable_with_value.name() == str( variable_to_find ):

                            initial_value = variable_with_value.initialValue()

                            external_species_eqs[key] = equation_to_check.subs( variable_to_find, initial_value )

    ext_spcs_chng = {}

    for species, equation in external_species_eqs.items():

        param_values = {}

        if hasattr(equation, 'free_symbols'):

            eq_symbols = equation.free_symbols

            if len( eq_symbols ) != 0:

                for symb in eq_symbols:

                    for key, value in sympy_to_CellML.items():

                        if value == str(symb):

                            digits = ''.join([char for char in str(key) if char.isdigit()])

                            param_values[symb] = solutions[int(digits)].tolist()

                            break

                # Extract the list of values for each symbol
                param_values_lists = list(param_values.values())
                symbolic_vars = list(param_values.keys())

                chng = 0

                for values in zip(*param_values_lists):
                    
                    subs_dict = dict(zip(symbolic_vars, values))
                    evaluated_equation = equation.subs(subs_dict)
                    result = evaluated_equation.evalf()  # Evaluate the equation to get a numerical result
                    delta_chng = delta_t * result

                    chng += delta_chng

                ext_spcs_chng[species] = chng

                print(chng)

            else:

                chng = t_f * float( equation )

                ext_spcs_chng[species] = chng

                print(chng)

        else:

            chng = t_f * float( equation )

            ext_spcs_chng[species] = chng

            print(chng)

    return ext_spcs_chng