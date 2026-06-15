import numpy as np
class Truss2d:
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
        self.displacements = np.array(displacements)
        
        return self.displacements

    def stress_strain(self):
        forces_list = []
        strain = []
        for truss in self.connectivity:
            i, j = truss[0], truss[1]
            u_element = self.displacements[[2*i, 2*i+1, 2*j, 2*j+1]]
            vector = self.nodes_list[j] - self.nodes_list[i]
            L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
            theta = np.arctan2(vector[1],vector[0])
            
            # Transformation matrix
            T = np.array([-np.cos(theta), -np.sin(theta), np.cos(theta), np.sin(theta)])
            
            # delta_L is the truss deformation value (Li - L0)
            delta_L = T @ u_element
            f_normal = (self.E * self.A / L) * delta_L
            forces_list.append(round(f_normal, 3))
            strain.append(delta_L/L)
            
        self.strain = strain
        self.normal = np.array(forces_list)
        self.stress = self.E*self.strain
        
        return self.stress

    def interpolate(self, bar_num, x):
        
        truss = self.connectivity[bar_num]
        i, j = truss[0], truss[1]
        u_element = self.displacements[[2*i, 2*i+1, 2*j, 2*j+1]]
        vector = self.nodes_list[j] - self.nodes_list[i]
        L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
        theta = np.arctan2(vector[1],vector[0])
        for i in range(len())
        # shape functions
        N_1 = (L-x)/L
        N_2 = x/L

        # interpolated displacement vector
        u = [0,0]
        u[0] = np.round(N_1 * u_element[0] + N_2 * u_element[2], 6)
        u[1] = np.round(N_1 * u_element[1] + N_2 * u_element[3], 6) 
        
        return np.array(u)

class Beam2d:
    def __init__(self, E, A, I, nodes_list, connectivity):
        self.E = E # Young's modulus
        self.I = I # Moment of Inercia
        self.A = A
        self.nodes_list = nodes_list
        self.connectivity = connectivity
        self.dof = 3*len(nodes_list)      # Degrees of freedom (u,v,θ)
        self.axial = Truss2d(E, A, nodes_list, connectivity)
        
    def stiffness_matrix(self):
        self.Kg = np.zeros((self.dof, self.dof))  # global stiffness matrix
        for truss in self.connectivity:
            i, j = truss[0], truss[1]
            matrix_index = [3*i+1, 3*i+2, 3*j+1, 3*j+2]
            axial_index = [3*i, 3*i +1, 3*j, 3*j+1]
            vector = self.nodes_list[j] - self.nodes_list[i]
            L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
            theta = np.arctan2(vector[1],vector[0])
            
            # Stiffness matrix for a beam element
            Kb = (self.E * self.I / L**3)*np.array([ 
                [12, 6*L, -12, 6*L],
                [6*L, 4*L**2, -6*L, 2*L**2],
                [-12, -6*L, 12, -6*L],
                [6*L, 2*L**2, -6*L, 4*L**2]
            ])
            K1D = (self.E*self.A/L)*np.array([ [1, 0, -1, 0],
            [-1, 0, 1, 0] ])
            # Stiffness matrix for a 2D truss element
            for m in range(4):
                for n in range(4):
                    self.Kg[matrix_index[m], matrix_index[n]] += Kb[m,n]
                    self.Kg[axial_index[m], axial_index[n]] += K1D[m,n]
            

        return self.Kg


class Frame2d:
    def __init__(self, nodes_list, connectivity, E = 1.0, A = 1.0, I = 1.0,):
        self.E = E # Young's modulus
        self.I = I # Moment of Inercia
        self.A = A
        self.nodes_list = nodes_list
        self.connectivity = connectivity
        self.dof = 3*len(nodes_list)      # Degrees of freedom (u,v,θ)

    def calculate_load_type(self, load_type, params, bar):
        truss = self.connectivity[bar]
        i, j = truss[0], truss[1]
        vector = self.nodes_list[j] - self.nodes_list[i]
        L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
        theta = np.arctan2(vector[1],vector[0])
        
        local_load = np.zeros(6)
        
        if load_type == "point":
            P = params["P"]  # load value
            a = params["Pa"] # distance from the last node
            b = L - a        # distance  from the first node

            # reactions vector
            local_load[1] = (P * b**2 * (3*a + b)) / (L**3)
            local_load[2] = (P * a * b**2) / (L**2)
            local_load[4] = (P * a**2 * (3*b + a)) / (L**3)
            local_load[5] = -(P * a**2 * b) / (L**2)

        elif load_type == "moment":
            M = params["M"]
            a = params["Ma"]
            b = L-a
            
            # reactions vector
            local_load[1] += 6*M*a*b/(L**3)
            local_load[2] += M*b*(2*a-b)/(L**2)
            local_load[4] += -6*M*a*b/(L**3)
            local_load[5] += M*a*(2*b-a)/(L**2)
        
        elif load_type == "uniform":
            q = params["q"]

            local_load[1] = q * L/2
            local_load[2] = q * (L**2)/12
            local_load[4] = q * L/2
            local_load[5] = - q * (L**2)/12
        
        elif load_type == "trapezoidal":
            q1 = params["q1"]
            q2 = params["q2"]
            
            #reactions vector
            local_load[1] = 7*q1*L/20 + 3*q2*L/20
            local_load[2] = (q1*L**2)/20 + (q2*L**2)/30
            local_load[4] = (3*q1*L)/20 + (7*q2*L)/20
            local_load[5] = -(q1*L**2)/30 - (q2*L**2)/20
            
        else:
            print('''Please, try again with one of the avaiable load conditions
                    "point" for concentrated load
                    "moment" for concentrated moment
                    "uniform" for uniform distributed load
                    "trapezoidal" for trapezoidal distributed load''')

        # transformation matrix
        c = np.cos(theta)
        s = np.sin(theta)

        T = np.array([
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1]
            ])
        local_load = T.T @ local_load
        
        return local_load
        
    def stiffness_matrix(self, hinges = None):
        self.Kg = np.zeros((self.dof, self.dof)) # global stiffness matrix
        if hinges == None:
            hinges = [0] * len(self.connectivity)
    
        for k, frame in enumerate(self.connectivity):
            i, j = frame[0], frame[1]
            matrix_index = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
            vector = self.nodes_list[j] - self.nodes_list[i]
            L = (vector[0]**2 + vector[1]**2)**(0.5)  #length
            theta = np.arctan2(vector[1],vector[0])
            
            # Stiffness matrix for a beam element
            rig_v = (self.E * self.I / L**3)
            rig_u = (self.E*self.A/L)

            hinge_type = hinges[k]

            if hinge_type == 0:
                frame_stiff = np.array([ 
                    [   rig_u,          0,             0,   -1*rig_u,           0,               0],
                    [       0,   12*rig_v,     6*L*rig_v,          0,   -12*rig_v,       6*L*rig_v],
                    [       0,  6*L*rig_v,  4*L**2*rig_v,          0,   -6*L*rig_v,   2*L**2*rig_v],
                    [-1*rig_u,          0,             0,    1*rig_u,            0,              0],
                    [       0,  -12*rig_v,    -6*L*rig_v,          0,     12*rig_v,     -6*L*rig_v],
                    [       0,  6*L*rig_v,  2*L**2*rig_v,          0,   -6*L*rig_v,   4*L**2*rig_v]
    
                ])

            elif hinge_type == 1: 
                #hinge in the initial node
                frame_stiff = np.array([
                    [   rig_u,          0,   0,   -1*rig_u,           0,               0],
                    [       0,    3*rig_v,   0,          0,    -3*rig_v,       3*L*rig_v],
                    [       0,          0,   0,          0,           0,               0],
                    [-1*rig_u,          0,   0,    1*rig_u,            0,              0],
                    [       0,   -3*rig_v,   0,          0,     3*rig_v,      -3*L*rig_v],
                    [       0,  3*L*rig_v,   0,          0,  -3*L*rig_v,    3*L**2*rig_v]
                ])

            elif hinge_type == 2:
                #hinge in the final node
                frame_stiff = np.array([
                    [   rig_u,          0,             0,   -1*rig_u,           0,   0],
                    [       0,    3*rig_v,     3*L*rig_v,          0,    -3*rig_v,   0],
                    [       0,  3*L*rig_v,  3*L**2*rig_v,          0,  -3*L*rig_v,   0],
                    [-1*rig_u,          0,             0,    1*rig_u,            0,   0],
                    [       0,   -3*rig_v,    -3*L*rig_v,          0,     3*rig_v,   0],
                    [       0,          0,             0,          0,           0,   0]
                ])
                
            # transformation matrix
            c = np.cos(theta)
            s = np.sin(theta)

            T = np.array([
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1]
            ])

            # Stiffness matrix for a 2D frame element
            frame_stiff = T.T @ frame_stiff @ T
            
            for m in range(len(matrix_index)):
                for n in range(len(matrix_index)):
                   self.Kg[matrix_index[m], matrix_index[n]] += frame_stiff[m,n]
        
        return self.Kg
        
    def reactions_solver(self, local_load, displacements, hinges=None): # displacements is a list with ones if there are free displacements in node and 0 if there are not 
        reduced_index = []
        Kg = self.stiffness_matrix(hinges)
        for i in range(len(displacements)):
            if displacements[i] != 0: 
                reduced_index += [i]
            else:
                continue
        reduced_Kg = Kg[np.ix_(reduced_index, reduced_index)]
        reduced_load = local_load[reduced_index]
        
        free_displacements = np.linalg.solve(reduced_Kg, reduced_load)
        for i in range(len(reduced_index)):
            displacements[reduced_index[i]] = free_displacements[i]
        support_forces = Kg @ displacements - local_load
        self.reactions = np.array([round(x,10) for x in support_forces])
        self.displacements = np.array(displacements)
    
        return [self.displacements, self.reactions]
