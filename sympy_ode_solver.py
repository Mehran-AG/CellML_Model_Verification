
# importing external or built-in packages
from sympy import symbols, lambdify, Symbol
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import sys
from colorama import Fore, Back, Style, init
import signal
import time

# importing internal packages
from compound_element_sorter import variable_name_mapper, initial_value_finder, variable_sorter
from matrix_equation_builder import printer


def sympy_ode_solver( components, concentration_rate_equations, general_equations, t_max = 10, delta_t = 0.001, printing = 'off' ):

    """
    This function solves a set of Ordinary Differential Equations (ODEs) by Sympy internal solver which is Runge Kutta 4th order.
    It gets a dictionary mapping the 'compound', represented by its chemical formula, to the equation created by Sympy variables using the stoichiometric matrix
    'components' is the list containing the imported CellML components
    """


    all_symbols = set().union( *[eq.free_symbols for eq in concentration_rate_equations.values()] )     # I need to know all variables that I have in these equations. This is why I should have replaced all variables that have values and are not meant ot be solved for

    number_of_variables = len( all_symbols )                                                              # I need to count the number of variables that I have so I can create the list of sympy variables as 'x' shown below

    number_of_equations = len( concentration_rate_equations )

    chebi_to_CellML = variable_name_mapper( components )

    # ------------ << Calling sorter function to sort the variables into their corresponding lists >> --------------
    _ , _ , _ , _ , reaction_rate_constants, boundary_conditions, _ , boundary_values = variable_sorter( components )

    # Initialize colorama
    init( autoreset=True )

    if number_of_variables != number_of_equations:

        if number_of_equations > number_of_variables:

            print( Style.BRIGHT + Fore.RED + "The number of equations, {e}, does not match the number of variables, {v}.".format( e = number_of_equations, v = number_of_variables ) )
            print( Fore.BLUE + "\nYour variables are:\n{var}".format( var = all_symbols ) )
            sys.exit( Fore.YELLOW + "\nExiting due to an error" + Style.BRIGHT + Fore.GREEN + "\nModify CellML file and check to see if you have defined the equations correctly\n\n")

        elif number_of_equations < number_of_variables:

            compounds_of_equations = [symbols(chebi_to_CellML[key]) for key in concentration_rate_equations.keys()]

            not_state_variables = list( all_symbols - set( compounds_of_equations ) )

            for variable_to_find in not_state_variables:
                
                for key, equation_to_check in concentration_rate_equations.items():

                    if variable_to_find in equation_to_check.free_symbols:

                        for variable_with_value in reaction_rate_constants:

                            if variable_with_value.name() == str( variable_to_find ):

                                initial_value = variable_with_value.initialValue()

                                concentration_rate_equations[key] = equation_to_check.subs( variable_to_find, initial_value )

                        for variable_with_value in boundary_values:

                            if variable_with_value.name() == str( variable_to_find ):

                                initial_value = variable_with_value.initialValue()

                                concentration_rate_equations[key] = equation_to_check.subs( variable_to_find, initial_value )

            # After replacing the variables which are not state variables and were unknown at first, We need to find the variables again
            all_symbols = set().union( *[eq.free_symbols for eq in concentration_rate_equations.values()] )

            number_of_variables = len( all_symbols )

            if number_of_variables < number_of_equations:

                print( Fore.RED + "The number of equations, {e}, does not match the number of variables, {v}.".format( e = number_of_equations, v = number_of_variables ) )
                print( Fore.BLUE + "\nYour variables are:\n{var}".format( var = all_symbols ) )
                
                print( Style.BRIGHT + Fore.MAGENTA + "\n" + "*"*50 )
                print( Style.BRIGHT + Fore.MAGENTA + "*"*50 )
                print( Back.WHITE + Fore.BLACK + "\nIMPORTANT: User input required." )

                condition = True

                while condition:
                    response = input( Style.BRIGHT + Fore.CYAN + "\nDo you want to continue running the code? (yes/no): ").strip().lower()

                    if response == 'no':
            
                        sys.exit( Fore.YELLOW + "\nExiting due to an error" + Style.BRIGHT + Fore.GREEN + "\nModify CellML file and check to see if you have defined the equations correctly\n\n")

                    elif response == 'yes':
                        
                        print("\nContinuing the code...")
                        break

                    else:

                        print( Style.BRIGHT + Fore.RED + "Invalid response. Please enter 'yes' or 'no'.")

                    print("*"*20 + "\n")

                time.sleep(2)
                print("\nCode continues running...")

            if printing == 'on' or printing =='On' or printing == 'ON':
        
                printer( concentration_rate_equations, 'Concentration rate equations generated by Stoichiometric Matrix after substitution of constants:' )


    sympy_to_CellML = {}                                                                                # I need to map sympy variables to CellML variables written by user. This is because I have ChEBI code and variable name for the variables that I have
                                                                                                        # For construction of the stoichiometric matrix, I used ChEBI codes and here I need variable names

    sympy_equations = {}                                                                                # This dictionary will store the sympy equations for concentration rate equations

    initial_values = []                                                                                 # Initial values should be stored in the same order as variables in 'x'

    xdot = []                                                                                           # Sympy needs the equations in the set as a list containing all the equations

    x = symbols( 'x:' + str(number_of_variables) )                                                           # Creation of the number of Sympy variables that I need

    chebi_initial_values = initial_value_finder( components, general_equations )

    # I need to replace the CellML variables I have used to write the equations with the Sympy variables which are like 'x0', 'x1', 'x3', ...

    for i, variable in enumerate ( all_symbols ):

        # here I will map the sympy variable to the CellML variable and store them in a dictionary for future reference
        sympy_to_CellML[x[i]] = str(variable)

        # Now, I need to replace the CellML variables with the corresponding Sympy variable in all equations of concentration rates
        # So, I need a 'for loop' to check all the equations
        for key in concentration_rate_equations.keys():

            concentration_rate_equations[key] = concentration_rate_equations[key].subs( variable, x[i] )

    for key in concentration_rate_equations.keys():

        # key in the concentration_rate_equations dictionary is the compound name as a string in chemical formula format
        try:

            to_find = chebi_to_CellML[key]  # I get the CellML variable of the chemical formulation. then I need to find its correspondent Sympy variable. Hence, I need a 'for' loop as below

        except:

            to_find = key

        for x_find, value in sympy_to_CellML.items():
            
            # to_find is a CellML variable, x_find is sympy variable, value is again CellMl variable mapped from x_find. So, I'll compare value and to_find
            #If they are the same, I'll map the sympy variable to its corresponding equation in a new dictionary which belongs to it

            if value == to_find:

                sympy_equations[x_find] = concentration_rate_equations[key]
                break

    # Now, I have the equations all in sympy format. A Sympy type variable mapped to its corresponding equation in Sympy formatted variables
    # Now, I have to find the initial values for these Sympy variables to solve the set of ODEs
    for i in range( len(x) ):

        variable = sympy_to_CellML[symbols('x'+str(i))]

        for key, value in chebi_to_CellML.items():

            if variable == value:
                
                try:

                    value_to_replace = chebi_initial_values[key]
                    initial_values.append( float(value_to_replace) )
                    break

                except ValueError as KE:

                    print(Style.BRIGHT + Fore.RED + "Initial Value is not defined for \"{v}\" and the ODEs cannot be solved without having initial values for all variables".format( v = key ))
                    sys.exit(Style.BRIGHT + Fore.MAGENTA +"Exiting due to an error\nModify CellML file and check to see if you have defined the inital value for this variable")
        
        variable = symbols('x'+str(i))

        # xdot is the list of equations that Sympy needs to solve
        for key in sympy_equations.keys():

            if variable == key:

                xdot.append(sympy_equations[key])
                break


    #################################################################
    #################################################################

    # Now, lets solve the obtained equations

    t = symbols( 't' )
    f = lambdify( ( t, x ), xdot )

    n = int(t_max/delta_t)+1

    t_eval = np.linspace( 0, t_max, n )

    # Timeout handler
    def timeout_handler():
        raise TimeoutError("Execution time exceeded")
    
    # Set the timeout duration (in seconds)
    timeout_duration = 15

    # Register the timeout handler
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout_duration)

    try:

        print(Style.BRIGHT + Fore.YELLOW + "\nSolving ODEs ...")

        solution = scipy.integrate.solve_ivp( f, (0, t_max), initial_values, t_eval = t_eval )

        # Cancel the alarm if the execution completes within the timeout duration
        signal.alarm(0)

        print(Style.BRIGHT + Fore.GREEN + "\n**** ODEs Solved ****\n\n")
    
        y = solution.y

        return y, t_eval, x, sympy_to_CellML
    
    except TimeoutError:

        print(Style.BRIGHT + Fore.RED + "\nTimeout, {te} seconds, exceeded and solution didn't converge!".format( te = timeout_duration ))
        exit(Style.BRIGHT + Fore.MAGENTA + "\nModify boundary conditions and run the code again!\n")

    except Exception:

        print(Style.BRIGHT + Fore.RED + "\nTimeout, {te} seconds, exceeded and solution didn't converge!".format( te = timeout_duration ))
        exit(Style.BRIGHT + Fore.MAGENTA + "\nModify boundary conditions and run the code again!\n")




# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# This function plots the solution

def plotter( solutions, time, variables_to_plot, x, sympy_to_CellML, show_legends = 'on' ):

    # Initialize colorama
    init( autoreset=True )

    if not variables_to_plot:

        legends = []

        # This 'for' loop creates the legend based on the CellMl variable names for the independent variables in ODEs
        for variable in x:

            legends.append( sympy_to_CellML[variable] )

        try:


            plt.plot( time, solutions.T )
            plt.title( 'Chemical Kinetics' )
            if show_legends == 'on' or show_legends == 'On' or show_legends == 'ON':
                plt.legend( legends, shadow = True )
            plt.xlabel('time')
            plt.ylabel('concentration')
            plt.grid(True)

            plt.show()

        except:

            print( Style.BRIGHT + Fore.RED + "\nSolution has failed and cannot be plotted")
            exit("\n")

    else:

        sympy_variables_to_plot = []

        legends = []

        for variable in variables_to_plot:

            flag = False

            for to_find in sympy_to_CellML.keys():

                if sympy_to_CellML[to_find] == variable:

                    sympy_variables_to_plot.append( to_find )
                    legends.append( variable )
                    digits = ''.join([char for char in str(to_find) if char.isdigit()])
                    flag = True
                    break
            
            if flag == True:

                try:

                    plt.plot(time, solutions[int(digits)])

                except:

                    print( Style.BRIGHT + Fore.RED + "\nSolution has failed and cannot be plotted." )
                    exit("\n")

            else:

                print( "**************************************************************\n**************************************************************\nThe variable {v} can't be found in the variables, make sure you have entered it correctly".format( v=variable ) )


        plt.title( 'Chemical Kinetics' )
        if show_legends == 'on' or show_legends == 'On' or show_legends == 'ON':
            plt.legend( legends, shadow = True )
        plt.xlabel('time')
        plt.ylabel('concentration')
        plt.grid(True)

        plt.show()