# Code to calculate the stiffness matrix and displacements in a single 2D truss
import numpy as np
# Data
E = 200e6 # Young's modulus
A = 0.07071 # Section area
L = 1.41421 # Length
θ = 45*np.pi/180 # Slope
# Stiffness matrix for a 1D truss element
K1D = E*A/L*np.array([ [1, -1],
[-1, 1] ])
# Transformation matrix
T = np.array([
    [np.cos(θ), np.sin(θ), 0, 0],
    [0, 0, np.cos(θ), np.sin(θ) ]
])
# Stiffness matrix for a 2D truss element
K = np.transpose(T)@K1D@T
print("K = ", K)
