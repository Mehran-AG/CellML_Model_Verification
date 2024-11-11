import numpy as np

from colorama import Fore, Back, Style, init



def kienetic_thermo_consts_convertor( conversion_matrix, kinetic_constants, printing = 'off' ):


    log_kinetic_constants = np.log( kinetic_constants )

    conversion_matrix_psuedo_inv = np.linalg.pinv(conversion_matrix)

    log_thrmo_cnsts = conversion_matrix_psuedo_inv @ log_kinetic_constants

    thrmo_cnsts = np.exp( log_thrmo_cnsts )

    # Initialize colorama
    init(autoreset=True)

    print( Style.BRIGHT + Fore.BLUE + "\nThermodynamic rate constants are successfully obtained" )

    print( Style.BRIGHT + Fore.GREEN + "\n   ******** The model is Thermodynamically CONSISTENT ********\n" )

    if printing == 'on' or printing =='On' or printing == 'ON':

        print( Style.BRIGHT + Fore.YELLOW + "\nThermodynamic constants are successfully obtained as below:\n{arr}".format( arr = thrmo_cnsts ) )

    return thrmo_cnsts