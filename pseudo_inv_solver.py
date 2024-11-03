import numpy as np

k_plus_1 = 0.5
k_plus_2 = 0.1
k_minus_1 = 0.1
k_minus_2 = 0.5
constant = 1

dpndncy_array = np.array([1, 1, -1, -1, -1])

#dpndncy_array = dpndncy_array.reshape( 1, -1 )

kinetic_rt_cnst_array = np.array([k_plus_1, k_plus_2, k_minus_1, k_minus_2, constant])

log_kinetic_rt_cnst_array = np.log(kinetic_rt_cnst_array)

#log_kinetic_rt_cnst_array = log_kinetic_rt_cnst_array.reshape( -1, 1 )

tst_rslt = dpndncy_array @ log_kinetic_rt_cnst_array

print( tst_rslt )

N = np.array([[1, 0, 1, 0, 1, 0], [0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1], [0, 1, 0, 1, 1, 0], [0, 0, 1, -1, 0, 0]])

N_psuedo_inv = np.linalg.pinv(N)

log_thrmo_cnsts = N_psuedo_inv @ log_kinetic_rt_cnst_array

thrmo_cnsts = np.exp( log_thrmo_cnsts )

print( thrmo_cnsts )