import numpy as np
import matplotlib.pyplot as plt
L=4
T=0.1
kb= 1#.380649e-23
templist = []
deltaElist = []
energylist = []
Magnetisationlist = []
probacceptlist = []
randomlist = []
np.random.seed(0)
# for a ferromagnetic system the coupling constant is positive:
J=1#.60218e-19

# create a random spin configuration
# the possible values of the spins are +1 or -1
spins = 2*np.random.randint(2, size=(L,L))-1
#print(spins)

for a in range(1):
        # calculate the energy of the system
    energy = 0
    for i in range(L):
        for j in range(L):
            S = spins[i,j]
            # this sums the neighbors with periodic boundary conditions
            if i+1 == L:
                neighbour1 = spins[0, j]
            else:
                neighbour1 = spins[(i+1), j]
            if j+1 == L:
                neighbour2 = spins[i,0]
            else:
                neighbour2 = spins[i,(j+1)]
            if i-1 == -1:
                neighbour3 = spins[L-1, j]
            else:
                neighbour3 = spins[(i-1), j]
            if j-1 == -1:
                neighbour4 = spins[i,L-1]
            else:
                neighbour4 = spins[i,(j-1)]
            neighbours = neighbour1 + neighbour2 + neighbour3 + neighbour4
            print("neighbours bf +=",neighbours*J*S)
            #instead of summing neighbor times spin ij, we sum all neighbors and multiply by the spin ij
            energy += -neighbours*J*S
            print("E ",energy)

    print(spins, "before flip")
    
    # Divide by 2 to correct for double-counting
    energy /= 2
    print("initial E ",energy)

    #gonna flip a random spin
    x = np.random.randint(L)
    y = np.random.randint(L)
    #copy the spins
    new_spins = np.copy(spins)
    #flip the spin
    new_spins[x,y] *= -1
    print(x,y)
    print(new_spins, "after flip")
    #calculate the energy of the new configuration
    new_energy = 0
    i=0
    j=0
    for i in range(L):
        for j in range(L):
            S = new_spins[i,j]
            if i+1 == L:
                neighbour1 = new_spins[0, j]
            else:
                neighbour1 = new_spins[(i+1), j]
            if j+1 == L:
                neighbour2 = new_spins[i,0]
            else:
                neighbour2 = new_spins[i,(j+1)]
            if i-1 == -1:
                neighbour3 = new_spins[L-1, j]
            else:
                neighbour3 = new_spins[(i-1), j]
            if j-1 == -1:
                neighbour4 = new_spins[i,L-1]
            else:
                neighbour4 = new_spins[i,(j-1)]
            neighbours = neighbour1 + neighbour2 + neighbour3 + neighbour4
            print("neighbours af +=",neighbours*J*S)
            new_energy += -neighbours*J*S
            print("new E ",new_energy)

    # Divide by 2 to correct for double-counting
    new_energy /= 2
    print("new E ",new_energy)

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
        p = np.exp(-deltaE/(kb*T))
        probacceptlist.append(p)
        r=0.85
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

# fig, axs = plt.subplots(3, 1, figsize=(10, 8))

# axs[0].plot(energylist)
# axs[0].set_xlabel('Iteration')
# axs[0].set_ylabel('Energy')
# axs[0].set_title('Energy vs. Iteration')

# axs[1].plot(deltaElist)
# axs[1].set_xlabel('Iteration')
# axs[1].set_ylabel('Delta E')
# axs[1].set_title('Delta E vs. Iteration')

# axs[2].plot(randomlist, label='Random List')
# axs[2].plot(probacceptlist, label='Acceptance Probability')
# axs[2].set_xlabel('Iteration')
# axs[2].set_ylabel('Value')
# axs[2].set_title('Random List and Acceptance Probability vs. Iteration')
# axs[2].legend()

# #print net magnetisation density
# magnetisation = np.abs(np.sum(spins)) * (1 / (L * L))
# print("final M ",magnetisation)

# plt.tight_layout()
# plt.show()
