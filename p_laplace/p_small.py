from p_laplace import p_laplace
import numpy as np
import matplotlib.pyplot as plt

epsilons = [10.0**i for i in range(-10, -22, -2)]
energies = []
maximums = []

# Converge with regularization
for eps in epsilons:
    result = p_laplace(1.1, eps)
    u = result[0]
    energies.append(result[1:])
    # Maxmimal nodal value (correct maximum only for P1 function)
    maximums.append(u.vector().max())
energies = np.array(energies)

# Report
print(' epsilon  |  energy_reg  |  energy | u max \n',
    np.concatenate( ( np.array([epsilons,]).T,
                      energies,
                      np.array([maximums,]).T,
                    ), axis=1))

# Plot energies
plt.figure()
plt.plot(epsilons, energies[:,0], 'o-', label='energy regularized')
plt.plot(epsilons, energies[:,1], 'o-', label='energy')
plt.xscale('log')
plt.legend(loc='upper left')
plt.show()
