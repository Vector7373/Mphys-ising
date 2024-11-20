import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# If running this program be prepared for it to save animation gifs to the save location
# path = Path(...) at the top of the file

L=10
T0=0.5
kb= 1#.380649e-23
Iterations = 2000
templist = []
deltaElist = []
energylist = []
Magnetisationlist = []
Magnetisationavglist = []
Magerrors = []
probacceptlist = []
randomlist = []
avgpool=10#10
r=0.1
#np.random.seed(0)
# for a ferromagnetic system the coupling constant is positive:
J=1#.60218e-19

#make a directory to save the gifs
custom_dir_name = f'1D_Crit_temp_AVG_Gif_set_L={L}_iters={Iterations}_random avg'
path = Path(r'') / custom_dir_name
path.mkdir(parents=True, exist_ok=True)

def initialise(L):
    """
    create a random spin configuration
    the possible values of the spins are +1 or -1
    returns an L length array of spins
    """
    return 2*np.random.randint(2, size=(1,L))-1

def MCSIM(L, T, J, spins, Iterations, r, g):
    """
    Monte Carlo simulation of the Ising model
    L: size of the 1D lattice
    T: temperature
    J: coupling constant    
    spins: initial spin config
    iterations: number of Monte Carlo steps
    This function finds E, then flips a spin, then finds E again, then finds the difference between them.
    if DeltaE < 0, accept the move
    if DeltaE > 0, accept the move with probability p = exp(-DeltaE/T)
    if the move is accepted then set the old spin configuration to the flipped spin configuration.
    then move onto the next iteration and repeat.
    This function returns the magnetisation of the final spin configuration for each temperature
    returns the magnetisation of the final spin configuration
    """
  
    for a in range(Iterations):
        energy = 0
        for j in range(L):
            S = spins[0,j]
            # this sums the neighbors with periodic boundary conditions
            if j+1 == L:
                neighbour2 = spins[0,0]
            else:
                neighbour2 = spins[0,j+1]
            if j-1 == -1:
                neighbour4 = spins[0,L-1]
            else:
                neighbour4 = spins[0,j-1]
            neighbours = neighbour2 + neighbour4
            #instead of summing neighbor times spin j, we sum all neighbors and multiply by the spin j
            energy += -neighbours*J*S

        # Divide by 2 to correct for double-counting (is this still applicable? check this)
        energy /= 2
        #gonna flip a random spin
        x = np.random.randint(L)
        #copy the spins
        new_spins = np.copy(spins)
        #flip the spin
        new_spins[0,x] *= -1
        #calculate the energy of the new configuration
        new_energy = 0
        j=0
        for j in range(L):
            S = new_spins[0,j]
            # this sums the neighbors with periodic boundary conditions
            if j+1 == L:
                neighbour2 = new_spins[0,0]
            else:
                neighbour2 = new_spins[0,j+1]
            if j-1 == -1:
                neighbour4 = new_spins[0,L-1]
            else:
                neighbour4 = new_spins[0,j-1]
            neighbours = neighbour2 + neighbour4
            #instead of summing neighbor times spin j, we sum all neighbors and multiply by the spin j
            new_energy += -neighbours*J*S

        # Divide by 2 to correct for double-counting
        new_energy /= 2

        #calculate the energy difference
        deltaE = new_energy - energy
        #print("deltaE", deltaE)

        #accept or reject the move
        if deltaE < 0:
            spins = new_spins
            energy = new_energy
            #print("accepted <0")
        else:
            #calculate the acceptance probability
            p = np.exp(-deltaE/(T))
            #print("p ",p)
            #print("r ",r)
            if r < p:
                spins = new_spins
                energy = new_energy
                #print("accepted by probability")
    # Magnetisation Calculation
    magnetisation = np.abs(np.sum(spins)) * (1 / (L))
    #print(spins)
    return magnetisation

# This runs the entire simulation for each val of r, increments of 0.1, and saves the plot
for c in range(9):
    templist=[]
    T=T0
    for t in range(25):#25
        templist.append(T)
        Magnetisationavglist = []
        for g in range(avgpool):
            spins = initialise(L)
            Magnetisationavglist.append(MCSIM(L, T, J, spins, Iterations, r, g))
            #print("Magnetisationavglist ",Magnetisationavglist, g)
        Magnetisationlist.append(np.mean(Magnetisationavglist))
        print("Magnetisationlist ",Magnetisationlist)
        Magerrors.append(np.std(Magnetisationavglist))
        T+=0.5
    # plotting and saving for each r
    fig=plt.figure()
    plt.plot(templist,Magnetisationlist)
    plt.errorbar(templist, Magnetisationlist, yerr=Magerrors)
    plt.xlabel('Temperature')
    plt.ylabel('Magnetisation')
    output_plot_path = path / f'Magnetisation_vs_Temperature_r={r}.png'
    plt.title(f'Magnetisation vs Temperature, r={r}')
    plt.savefig(output_plot_path)
    plt.close(fig)
    Magnetisationlist=[]
    Magerrors=[]
    r+=0.1
