import numpy as np
import matplotlib.pyplot as plt
L=6
T=0.5
deltaElist = []
energylist = []
sim_energy = []

# for a ferromagnetic system the coupling constant is positive:
J=2.0

#repeat simulation 100 times
for a in range(100):
    # create a random spin configuration
    # the possible values of the spins are +1 or -1
    spins = 2*np.random.randint(2, size=(L,L))-1

    # calculate the energy of the system
    energy = 0
    for i in range(L):
        for j in range(L):
            S = spins[i,j]
            # this sums the neighbors with periodic boundary conditions by using the modulo operator % to wrap around the edges
            nb = spins[(i+1)%L, j] + spins[i,(j+1)%L] + spins[(i-1)%L, j] + spins[i,(j-1)%L]
            #instead of summing neighbor times spin ij, we sum all neighbors and multiply by the spin ij
            energy += -nb*J*S

    # Divide by 2 to correct for double-counting the pairwise interactions
    energy /= 2

    print("initial E ",energy)
    for i in range(400):
        #gonna flip a random spin
        i = np.random.randint(L)
        j = np.random.randint(L)
        #copy the spins
        new_spins = np.copy(spins)
        #flip the spin
        new_spins[i,j] *= -1

        #calculate the energy of the new configuration
        new_energy = 0
        for i in range(L):
            for j in range(L):
                S = new_spins[i,j]
                nb = new_spins[(i+1)%L, j] + new_spins[i,(j+1)%L] + new_spins[(i-1)%L, j] + new_spins[i,(j-1)%L]
                new_energy += -nb*J*S

        # Divide by 2 to correct for double-counting the pairwise interactions
        new_energy /= 2
        #calculate the energy difference
        deltaE = new_energy - energy
        deltaElist.append(deltaE)
        print("deltaE ",deltaE)


        #calculate the acceptance probability
        # kb=1 for simplicity
        p = np.exp(-deltaE/T)
        r = np.random.rand()
        #accept or reject the move
        if deltaE < 0:
            spins = new_spins
            energy = new_energy
            print("accepted")
        elif  r < p:
            spins = new_spins
            energy = new_energy
            print("accepted")
        else:
            print("rejected")
        energylist.append(energy)
    sim_energy.append(energy)

print(sim_energy)
plt.scatter(range(100),sim_energy)
plt.show()


print("final E ",energy)
print(spins)


