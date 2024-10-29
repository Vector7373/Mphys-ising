import numpy as np
import matplotlib.pyplot as plt
L=6
T=1
kb= 1.380649e-23
templist = []
deltaElist = []
energylist = []
Magnetisationlist = []
probacceptlist = []
randomlist = []

# for a ferromagnetic system the coupling constant is positive:
J=1

# create a random spin configuration
# the possible values of the spins are +1 or -1
spins = 2*np.random.randint(2, size=(L,L))-1
#print(spins)

for i in range(10):
        # calculate the energy of the system
    energy = 0
    for i in range(L):
        for j in range(L):
            S = spins[i,j]
            # this sums the neighbors with periodic boundary conditions by using the modulo operator % to wrap around the edges
            nb = spins[(i+1)%L, j] + spins[i,(j+1)%L] + spins[(i-1)%L, j] + spins[i,(j-1)%L]
            #instead of summing neighbor times spin ij, we sum all neighbors and multiply by the spin ij
            energy += -nb*J*S

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
            nb = new_spins[(i+1)%L, j] + new_spins[i,(j+1)%L] + new_spins[(i-1)%L, j] + new_spins[i,(j-1)%L]
            new_energy += -nb*J*S

    # Divide by 2 to correct for double-counting
    new_energy /= 2
    #calculate the energy difference
    deltaE = new_energy - energy
    deltaElist.append(deltaE)
    print("deltaE ",deltaE)

    #accept or reject the move
    if deltaE < 0:
        probacceptlist.append(1)
        spins = new_spins
        energy = new_energy
        print("accepted <0")
    else:
        #calculate the acceptance probability
        p = np.exp(-deltaE/kb*T)
        probacceptlist.append(p)
        r=0.7
        #r = np.random.rand()
        randomlist.append(r)
        print("p ",p)
        print("r ",r)
        if r < p:
            spins = new_spins
            energy = new_energy
            print("accepted by probability")
        if r > p:
            print("rejected by probability")
    energylist.append(energy)
magnetisation = np.abs(np.sum(spins)) * (1 / (L * L))
print("mag ",magnetisation)
print("final E ",energy)
print(spins)

fig, axs = plt.subplots(2, 1, figsize=(10, 8))

axs[0].plot(energylist)
axs[0].set_xlabel('Iteration')
axs[0].set_ylabel('Energy')
axs[0].set_title('Energy vs. Iteration')

axs[1].plot(deltaElist)
axs[1].set_xlabel('Iteration')
axs[1].set_ylabel('Delta E')
axs[1].set_title('Delta E vs. Iteration')

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(randomlist, label='Random List')
ax.plot(probacceptlist, label='Acceptance Probability')
ax.set_xlabel('Iteration')
ax.set_ylabel('Value')
ax.set_title('Random List and Acceptance Probability vs. Iteration')
ax.legend()

#print net magnetisation density
magnetisation = np.abs(np.sum(spins)) * (1 / (L * L))
print("final M ",magnetisation)

plt.tight_layout()
plt.show()
