from fea_models import Frame2d

# Coordinates [X, Y] of each important point
nodes = np.array([
    [0.0, 0.0],  # Node 0 (Support A)
    [5.0, 0.0],  # Node 1 (Support B)
    [0.0, 3.0],  # Node 2 (Point C - Hinge/Intersection)
    [5.0, 3.0],  # Node 3 (Point D)
    [0.0, 6.0],  # Node 4 (Point E)
    [2.5, 6.0],  # Node 5 (Point F - Central hinge)
    [5.0, 6.0],  # Node 6 (Point G)
    [5.0, 4.5],  # Node 7 (Point H - Hinge on the column)
    [2.5, 3.0]   # Node 8 (Moment 40 kNm)
])

# Connectivity (which nodes form each member/bar)
bar = np.array([
    [0, 2],  # Member 0: Column AC (with distributed load q=20)
    [2, 4],  # Member 1: Column CE
    [4, 5],  # Member 2: Beam EF (with distributed load q=30 and moment of 60 at E)
    [5, 6],  # Member 3: Beam FG (with distributed load q=30)
    [6, 7],  # Member 4: Column GH
    [7, 3],  # Member 5: Column HD
    [2, 8],  # Member 6: Hinge-beam C (with moment of 40 at the center and force of 50 at node D)
    [8, 3],  # Member 7: Hinge-D
    [3, 1]   # Member 8
])

frame_structure = Frame2d(nodes, bar)
x = np.zeros(3 * 9)
params1 = {"q": -20, "P": F}
x[[0, 1, 2, 6, 7, 8]] += truss.calculate_load_type('uniform', params1, 0)
params2 = {"q": -30, "P": F}
x[[12, 13, 14, 15, 16, 17]] += frame_structure.calculate_load_type('uniform', params2, 2)
x[[15, 16, 17, 18, 19, 20]] += frame_structure.calculate_load_type('uniform', params2, 3)

display = [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
hinges = [0, 1, 0, 1, 0, 1, 1, 0, 0]  # for each bar where the hinge is positioned
x = x + np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40])

displacements, reactions = frame_structure.reactions_solver(x, display, hinges)

