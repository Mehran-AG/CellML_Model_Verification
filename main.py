
"""
*******************************************************************************************************************************
*                                                                                                                             *
*  This module reads CellML files containing reaction stoichiometries and equations                                           *
*                                                                                                                             *
*  Then using the modules provided constructs stoichiometric and elemental matrices and verifies the model                    *
*                                                                                                                             *
*  Then the corresponding Concentration Rate Equations are generated using the Stoichiometric matrix and solved               *
*                                                                                                                             *
*  Full details of the code can be found at the end of this file                                                              *
*                                                                                                                             *
*******************************************************************************************************************************
"""

# Importing external or built-in packages
import numpy as np
import sympy as sp
import os

# Importing internal packages
import excel_read as xls
import compounds_extractor as comex
import stochiometric_matrix_builer as smb
import elemental_matrix_builder as emb
import rate_matrix_builder as rmb
import verification as vf
import CellML_reader as cmlr
import compound_element_sorter as ces
import equation_builder as eb
import matrix_equation_builder as meb
import sympy_ode_solver as sos
import equations_from_text as eft
import external_concentrations as ec
import conservation_check as cc
import reversibility_check as rs
import kinetic_thermo_conversion_matrix as ktcm
import kinetic_thermo_consts_convertor as ktcc


command = 'cls' if os.name == 'nt' else 'clear'
os.system(command)

#file_path = './docs/aguda_tang_1999.cellml'
#file_path = './docs/aguda_b_1999.cellml'
#file_path = './docs/huang_ferrell_1996.cellml'
file_path = './docs/modified_huang_ferrell_1996.cellml'
#file_path ='./docs/NitrosylBromide.cellml'
cellml_file_dir = file_path
cellml_file = file_path
cellml_strict_mode = False



bc_file_path = None
#bc_file_path = './docs/boundary_conditions.txt'
#bc_file_path = './docs/boundary_conditions_AT.txt'
#bc_file_path = './docs/boundary_conditions_NB.txt'

#eq_file_path = './docs/equations.txt'
eq_file_path = None

read_eqs_from_file = False
solve_equations = False
delta_t = 0.1
t_f = 500


components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

variables, coefficients, enzymes, reaction_rates, rate_constants, boundary_conditions, equation_variables, boundary_values = ces.variable_sorter( components )

try:
    
    bc_file_path

    imported_bc_equations = eft.read_equations_from_file( bc_file_path )

except:

    imported_bc_equations = None


if solve_equations == True:

    if read_eqs_from_file:

        reaction_rate_equations_dict, bc_equations_dict, general_equations = eb.equation_reader( components, eq_file_path, imported_bc_equations, 'on' ) #print

    else:

        reaction_rate_equations_dict, bc_equations_dict, general_equations = eb.equation_builder( components, imported_bc_equations, 'on' ) #print


element_indices, compound_indices, reaction_indices, symbols_list, compound_to_composition, bcvirtual_compound_coefficients = ces.cellml_compound_element_sorter ( components )

element_matrix = emb.elemental_matrix_builder( compound_indices, element_indices, compound_to_composition )

stoichiometric_matrix = smb.stoichiometric_matrix_builder( reaction_indices, compound_indices, coefficients, bcvirtual_compound_coefficients )

if solve_equations == True:

    concentration_rate_equations = meb.matrix_equation_builder ( stoichiometric_matrix, compound_indices, reaction_indices, reaction_rate_equations_dict, general_equations, components, imported_bc_equations, 'on' ) #print

rate_matrix = rmb.rate_matrix_builder ( symbols_list )

# Calling the function
equations_array = vf.verification( stoichiometric_matrix, element_matrix, element_indices, compound_indices, reaction_indices, rate_matrix )

if solve_equations == True:

    solution, time, x, sympy_to_CellML = sos.sympy_ode_solver( components, concentration_rate_equations, general_equations, t_f, delta_t, 'on' )

    variables_to_plot = []

    sos.plotter(  solution, time, variables_to_plot, x, sympy_to_CellML, show_legends = 'on' ) # Show Legends


rs.reversibility_check( components )

conversion_matrix, kinetic_constants = ktcm.kietic_thermo_convertor( components, reaction_indices, compound_indices, coefficients )

ktcc.kienetic_thermo_consts_convertor( conversion_matrix, kinetic_constants )










'''
*****************************************************************************************************************************************************************************

The code is structured as below:

The main file calls the functions one by one
The functions are written in separate py files

1- The code starts by calling "CellML_reader" function
   CellML reader reads the CellML file and returns a list containing the components of CellML as classes
   The list is named 'components'

2- "variable_sorter" gets the list of components and sorts the variables in components into corresponding lists of 'variables', 'coefficients', 'reaction_rates', 'rate_constants', 'boundary_conditions
   'variables' are the concentrations of the compounds
   'coefficients' are the stoichiometric coefficients of the compounds in reactions
   'reaction_rates' are the variables for the reaction rates
   'rate_constants' are constants for forward and revese reaction rates
   'boundary_conditions' are the flow rates for the compounds flowing in and out of a compartment

3- "equation_builder" gets the list of components and generates the equations of reactions such as reaction rates and virtual reactions for the boundary conditions
    This function returns dictionaries 'reaction_rate_equations_dict' and 'bc_equations_dict'
    These dictionaries map the compound to their corresponding equations. Reaction rate equations will be used in the concentration rate equations generated from the stoichiometric matrix

4- "cellml_compound_element_sorter" gets the components and sorts the compounds and elements. By sorting I mean that the compounds and elements are given an index number to be identified by that number in the code
    This number will be used to generate the Elemental and Stoichiometric matrices. In stoichiometric matrix, each row is representative of a compoun and each column is representative of a reaction
    In Elemental matrix, each row is representative of an element and each column is representative of a compound
    The function returns 'element_indices', 'compound_indices', 'reaction_indices', 'symbols_list', 'compound_to_composition', 'bcvirtual_compound_coefficients'
    'element_indices' maps each element to a specific number which shows its row in the Elemental matrix
    'compound_indices' maps each compound to a specific number which shows its row in stoichiometric and its column in Elemental matrix
    'reaction_indices' maps each reaction to a specific number which shows its column in stoichiometric matrix
    'symbols_list' is a list containing variables of concentration rates for compounds in the same order as the compounds have been enumerated
    'compound_to_composition' maps the compound to its chemical composition. Chemical composition is a dictionary containing the details of the elements
    'bcvirtual_compound_coefficients' maps the virtual external and internal boundary condition compound to its stoichiometric coefficient in the virtual reaction

5- "elemental_matrix_builder" gets the 'element_indices', 'compound_indices', and 'compound_to_composition' dictionaries and generates the Elemental matrix

6- "stoichiometric_matrix_builder" gets 'reaction_indices', 'compound_indices', 'coefficients', 'bcvirtual_compound_coefficients' and generates the stoichiometric matrix

7- "matrix_equation_builder" gets 'stoichiometric_matrix', 'compound_indices', 'reaction_indices', 'reaction_rate_equations_dict', 'components' and generates concentration rate equations from the Stoichiometric Matrix

8- "rate_matrix_builder" gets the symbols list and creates the rate matrix containing Sympy symbols of the variables in symbols_list which can be multiplied to the left null space to get the conservation equations

9- "verification" gets 'stoichiometric_matrix', 'element_matrix', 'element_indices', 'compound_indices', 'reaction_indices', 'rate_matrix' and checks the verification of the model
    The function also prints the species violating the mass balance and the reaction in which this violation occurs

10- "sympy_ode_solver" gets 'components', 'concentration_rate_equations', final_time, time_step and returns the 'solution', 'time', 'x', 'sympy_to_CellML'
    'solution' is the class of sympy containing the solutions for the variables
    'x' is the list of the dependent variables
    'time' is a list containing all time steps taken by the solver. It is to plot the figure
    'sympy_to_CellML' is a dictionary mapping sympy variables which are 'x's to their coresponding CellML variables. This will be used to create the legend for the figure

11- "plotter" gets 'solution', 'time', 'variables_to_plot', 'x', 'sympy_to_CellML' ) and plots the figure
    'variables_to_plot' is a list containing the variables that we want to draw the curves for them. If it is left empty, the plotter will draw all curves

*****************************************************************************************************************************************************************************
'''