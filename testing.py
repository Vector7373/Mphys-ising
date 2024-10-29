import numpy as np

def apply_periodic_boundary_conditions(spins, i, j):
    L = spins.shape[0]
    nb = spins[(i+1)%L, j] + spins[i, (j+1)%L] + spins[(i-1)%L, j] + spins[i, (j-1)%L]
    return nb

def main():
    # Example spin configuration (3x3 grid)
    spins = np.array([[1, -1, 1],
                      [-10, 5, -3],
                      [1, -1, 1]])
    
    i, j = 1, 1  # Position to check neighbors
    nb = apply_periodic_boundary_conditions(spins, i, j)

    #hamiltonian
    J = 1.0
    S = spins[i, j]
    
    
    print(f"Spin configuration:\n{spins}")
    print(f"Sum of neighbors for spin at ({i}, {j}): {nb}")

if __name__ == "__main__":
    main()