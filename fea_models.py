import numpy as np
class truss_2d:
    def __init__(self, E, A, nodes_list, connectivity):
        self.E = E # Young's modulus
        self.A = A # Section area
        self.nodes_list = nodes_list
        self.connectivity = connectivity
        self.dof = 2*len(nodes_list)      # Degrees of freedom
    def stiffness_matrix(self):
        self.Kg = np.zeros((self.dof, self.dof))  # global stiffness matrix
        for truss in self.connectivity:
            first_node = truss[0]
            last_node = truss[1]
            vector = self.nodes_list[last_node] - self.nodes_list[first_node]
            L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
            theta = np.arctan2(vector[1],vector[0])
            
            # Stiffness matrix for a 1D truss element
            K1D = (self.E*self.A/L)*np.array([ [1, -1],
            [-1, 1] ])
            
            # Transformation matrix
            T = np.array([
                [np.cos(theta), np.sin(theta), 0, 0],
                [0, 0, np.cos(theta), np.sin(theta) ]
            ])
            # Stiffness matrix for a 2D truss element
            matrix_index = [2*first_node, 2*first_node+1, 2*last_node, 2*last_node+1]
            K2D = np.transpose(T)@K1D@T
            for i in range(4):
                for j in range(4):
                    self.Kg[matrix_index[i], matrix_index[j]] += K2D[i,j]
        return self.Kg
        
    def free_displacements(self, forces, displacements):
        reduced_index = []  # Selecting the rows and columns with known boundary conditions
        Kg = self.stiffness_matrix()
        for i in range(len(displacements)):
            if displacements[i] != 0:
               reduced_index += [i]
        reduced_Kg = np.zeros((len(reduced_index),len(reduced_index)))
        reduced_forces = np.zeros(len(reduced_index))
        for i in range(len(reduced_index)):
            reduced_forces[i] +=  forces[reduced_index[i]]
            for j in range(len(reduced_index)):
                reduced_Kg[reduced_index[i], reduced_index[j]] += Kg[i,j]
        free_displacements = np.linalg.solve(reduced_Kg, reduced_forces)   # Calculating free displacements associated with the boundary conditions
        for i in range(len(reduced_index)):
            displacements[reduced_index[i]] = free_displacements[i]
        support_forces = Kg @ displacements  # Calculating the reaction forces considering the displacements in each node
        self.reactions = np.array([round(x,10) for x in support_forces])
        displacements = np.array(displacements)
        return [displacements]
