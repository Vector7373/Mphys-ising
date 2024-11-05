import numpy as np
import matplotlib.pyplot as plt
L=16
T=0.1
iterations = 6000
avgpool = 20
0
kb= 1#.380649e-23
templist = []
deltaElist = []
energylist = []
Magnetisationlist = []
Magnetisationavglist = []
probacceptlist = []
randomlist = []
magerrorlist = []

for t in range(78):
    templist.append(T)
    Magnetisationavglist = []
    for a in range(avgpool):
        # for a ferromagnetic system the coupling constant is positive:
        J=1#.60218e-19

        # create a random spin configuration
        # the possible values of the spins are +1 or -1
        spins = 2*np.random.randint(2, size=(L,L))-1
        #print(spins)

        for i in range(iterations):
                # calculate the energy of the system
            energy = 0
            for i in range(L):
                for j in range(L):
                    S = spins[i,j]
                    # this sums the neighbors with periodic boundary conditions by using the modulo operator % to wrap around the edges
                    neighbours = spins[(i+1)%L, j] + spins[i,(j+1)%L] + spins[(i-1)%L, j] + spins[i,(j-1)%L]
                    #instead of summing neighbor times spin ij, we sum all neighbors and multiply by the spin ij
                    energy += -neighbours*J*S

            # Divide by 2 to correct for double-counting
            energy /= 2
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
                    neighbours = new_spins[(i+1)%L, j] + new_spins[i,(j+1)%L] + new_spins[(i-1)%L, j] + new_spins[i,(j-1)%L]
                    new_energy += -neighbours*J*S

            # Divide by 2 to correct for double-counting
            new_energy /= 2
            #calculate the energy difference
            deltaE = new_energy - energy
            deltaElist.append(deltaE)
            #print("deltaE ",deltaE)

            #accept or reject the move
            if deltaE < 0:
                probacceptlist.append(1)
                spins = new_spins
                energy = new_energy
                #print("accepted <0")
            else:
                #calculate the acceptance probability
                p = np.exp(-deltaE/kb*T)
                probacceptlist.append(p)
                r=0.8
                #r = np.random.rand()
                randomlist.append(r)
                #print("p ",p)
                #print("r ",r)
                if r < p:
                    spins = new_spins
                    energy = new_energy
                    #print("accepted by probability")
            energylist.append(energy)
        magnetisation = np.abs(np.sum(spins)) * (1 / (L * L))
        Magnetisationavglist.append(magnetisation)
        print("magnetisation ",magnetisation," T ",T)
    Magnetisationlist.append(np.mean(Magnetisationavglist))
    magerrorlist.append(np.std(Magnetisationavglist))
    T+=0.5

print(spins)

plt.errorbar(templist,Magnetisationlist,yerr=magerrorlist,fmt='o')
plt.xlabel('Temperature')
plt.ylabel('Magnetisation')
plt.title('Magnetisation vs Temperature')
plt.figtext(0.15, 0.85, f'J = {J}', fontsize=12)
plt.figtext(0.15, 0.80, f'L = {L}', fontsize=12)
plt.figtext(0.15, 0.75, f'iterations = {iterations}', fontsize=12)
plt.figtext(0.15, 0.70, f'avgpool = {avgpool}', fontsize=12)
plt.show()
