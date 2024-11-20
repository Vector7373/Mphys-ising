import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pathlib import Path
from matplotlib.animation import PillowWriter

# If running this program be prepared for it to save animation gifs to the save location
# path = Path(...) at the top of the file

L=10
T=1
kb= 1#.380649e-23
Iterations = 10000
templist = []
deltaElist = []
energylist = []
Magnetisationlist = []
Magnetisationavglist = []
Magerrors = []
probacceptlist = []
randomlist = []
avgpool=10
r=0.25
#np.random.seed(0)
# for a ferromagnetic system the coupling constant is positive:
J=1#.60218e-19

#make a directory to save the gifs
custom_dir_name = f'1DSingle_TESTMK1_Gif_set_L={L}_iters={Iterations}_r={r}_random avg'
path = Path(r'') / custom_dir_name
path.mkdir(parents=True, exist_ok=True)

def initialise(L):
    """
    create a random spin configuration
    the possible values of the spins are +1 or -1
    returns an L array of spins
    """
    return 2*np.random.randint(2, size=(1,L))-1

#print(initialise(L))

def MCSIM(L, T, J, spins, Iterations, r):
    """
    Monte Carlo simulation of the Ising model
    L: size of the lattice
    T: temperature
    J: coupling constant    
    spins: initial spin config
    iterations: number of Monte Carlo steps
    This function finds E, then flips a spin, then finds E again, then finds the difference between them.
    if DeltaE < 0, accept the move
    if DeltaE > 0, accept the move with probability p = exp(-DeltaE/T)
    if the move is accepted then set the old spin configuration to the flipped spin configuration.
    then move onto the next iteration and repeat.
    This function also animates the spins at every 200th iteration.
    This function returns the magnetisation of the final spin configuration for each temperature
    returns the magnetisation of the final spin configuration
    """

    fig = plt.figure()
    images = []   
    for a in range(Iterations):
        print("old spins", spins)
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
        print("energy", energy)
        #gonna flip a random spin
        x = np.random.randint(L)
        print("X=",x)
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
        print("new spins", new_spins)
        print("new energy", energy)

        #calculate the energy difference
        deltaE = new_energy - energy
        print("deltaE ",deltaE)

        #accept or reject the move
        if deltaE < 0:
            spins = new_spins
            energy = new_energy
            print("accepted <0")
        else:
            #calculate the acceptance probability
            p = np.exp(-deltaE/(T))
            #r = np.random.rand()
            print("p ",p)
            #print("r ",r)
            if r < p:
                spins = new_spins
                energy = new_energy
                print("accepted by probability")
    # ANIMATION -----------
    # only plot if its the end of the avgpool
        if a%10 == 0:
            im = plt.imshow( spins , cmap = 'viridis', animated=True)
            images.append([im])
    animate = ani.ArtistAnimation(fig, images, interval=25, blit=True)
    output_path = path / f'1D_ising_model_L{L}_T{T:.1f}.gif'
    animate.save(output_path, writer=PillowWriter())
    print("animation saved")
    plt.close(fig)
    # ANIMATION -----------
    # Magnetisation Calculation
    magnetisation = np.abs(np.sum(spins)) * (1 / (L))
    print("magnetisation ",magnetisation," T ",T)
    return magnetisation




spins = initialise(L)
Magnetisationlist.append(MCSIM(L, T, J, spins, Iterations, r))
print("Magnetisationlist ",Magnetisationlist)
