# importing external or built-in packages
from sympy import symbols, lambdify, Symbol, Eq
import numpy as np
import scipy.integrate
import sys
from colorama import Fore, Back, Style, init

# importing internal packages
from compound_element_sorter import variable_name_mapper, initial_value_finder, variable_sorter, all_digits
from matrix_equation_builder import printer
import chebi_fetch as chf

def conservation_eqs_verification( conservations_array, solution, sympy_to_CellML, ext_spcs_chng):

    for conservation_eq in conservations_array:

        symbols_list = conservation_eq.free_symbols

        for symb in symbols_list:

            if len( str(symb).split('_') ) == 1:

                spcs_chng = None

                for sp_vrbl, variable in sympy_to_CellML.items():

                    if str(symb) == variable:

                        digits = ''.join([char for char in str(sp_vrbl) if char.isdigit()])

                        spcs_values = solution[int(digits)].tolist()

                        spcs_chng = spcs_values[-1] - spcs_values[0]

                        break

                if spcs_chng != None:

                    conservation_eq = conservation_eq.subs( symb, spcs_chng )

            elif len( str(symb).split('_') ) > 1:

                if str(symb).split('_')[1] == 'e':

                    variable = str(symb).split('_')[0]

                    chng = ext_spcs_chng[variable]

                    conservation_eq = conservation_eq.subs( symb, chng )

                else:

                    spcs_chng = None

                    for sp_vrbl, variable in sympy_to_CellML.items():

                        if str(symb) == variable:

                            digits = ''.join([char for char in str(sp_vrbl) if char.isdigit()])

                            spcs_values = solution[int(digits)].tolist()

                            spcs_chng = spcs_values[-1] - spcs_values[0]

                            break

                    if spcs_chng != None:

                        conservation_eq = conservation_eq.subs( symb, spcs_chng )

        
        evaluated_eqn = conservation_eq.evalf()

        print(evaluated_eqn)
        

    return