# asynchronous implicit variational integrators for 2D-elastodynamics
# @Prashanth Setty (LTD, FAU Erlangen)
# last update 24/07/2024

# imports
from fenics import *
import numpy as np
import heapq
import mesh_converter as mc
import differentiator as df
import os
import matplotlib.pyplot as plt

# function to calculate linear strain
def epsilon(u):
    return sym(grad(u))

# function to calculate strain energy density
def psi(u):
	return 0.5*lmbda*(tr(epsilon(u)))**2+mu*inner(epsilon(u),epsilon(u))

def AVI_explicit(mesh):
    # definition of function spaces for velocity, displacement, nodal time
    V = VectorFunctionSpace(mesh, 'CG', 1)
    dg= FunctionSpace(mesh, 'DG', 0)
    position = Function(V)     # position vector
    velocity = Function(V)     # velocity vector
    ta = Function(V)           # nodal time vector
    mass = Function(V)         # mass vector: just another way to enforce mass lumping
    updates = Function(dg)
    t_step = Function(dg)

    # Initialize dictionaries to store displacement vectors
    displacements = {
        't_minus_2': Function(V),
        't_minus_1': Function(V),
        't_current': Function(V)
    }

    # splitting vector v in two components and making dof map
    x_comp, y_comp = V.split()
    dof_coordinate = V.tabulate_dof_coordinates() #each node has 2 dof in 2D,so it contains 2*No. of nodes amount of terms representing co.ord

    # priority que
    dts = []
    # nodal time vector
    ta.vector()[:] = 0.0 #size = dof of domain = 2*no. of nodes
    # elemental time
    tk=[]
    #mass vector
    mass.vector()[:] = 0.0 #size = dof of domain = 2*no. of nodes

    # boolean function for fixed bottom
    def bottom(x, on_boundary):
        return near(x[1], 0.0) and on_boundary

    # initial condition and boundary conditions
    a = 0.005
    disp = Expression(("0", "x[1]* a"), t=0, a=a, degree=2)
    disp2 = Expression(("0", "0"), t=0, degree=2)
    velo = Expression(("0", "0"), t=0, degree=2)
    bc1 = DirichletBC(V, disp2, bottom)
    bcdisp1 = []
    bcdisp1 = bc1.get_boundary_values().keys()

    # Project the initial displacement expression to the first two time steps
    displacements['t_minus_2'].assign(project(disp, V))
    displacements['t_minus_1'].assign(project(disp, V))

    # post-processing .txt file
    filespath = "resultsavi/"
    filename = "2d_AVI_implicit_elastodynamics"
    fl = open(filename + '.txt', 'w')
    if not os.path.exists(filespath):
        os.makedirs(filespath, exist_ok=True)

    all_xnodes=[]; all_ynodes=[]
    # initialization of the problem
    for cell in range(mesh.num_cells()):
        # making a priority que from the time steps
        element = Cell(mesh, cell)
        dt = (f* element.inradius())/c
        tk.append(0.0)
        dts.append((dt, cell))
        
        # mass of an element
        m = rho * element.volume()
        
        x_nodes = x_comp.dofmap().cell_dofs(cell)
        y_nodes = y_comp.dofmap().cell_dofs(cell)

        all_xnodes.append(x_nodes) 
        all_ynodes.append(y_nodes)

        # definition of initial position here
        for node in all_xnodes[cell]:
            position.vector()[node] = dof_coordinate[node, 0] + displacements['t_current'].vector()[node]
            mass.vector()[node] = mass.vector()[node] + m / 3.0
        for node in all_ynodes[cell]:
            position.vector()[node] = dof_coordinate[node, 1] + displacements['t_current'].vector()[node]
            mass.vector()[node] = mass.vector()[node] + m / 3.0
        updates.vector()[cell] = 0.0
        t_step.vector()[cell] = (f* element.inradius())/c
    # initial time
    t=0.0
    # discrete action
    Sd = 0.0

    # making heap from the previous array
    heapq.heapify(dts)
    if not dts:  # Check if the list is empty
        raise ValueError("Priority queue is empty before simulation starts.")

    while(len(dts)!=0):
        # ascending order of time steps by sorting
        dts.sort()
        # pop an element with highest priority from queue
        el = heapq.heappop(dts)
        cell=el[1]

        # updating time
        t = el[0]
        print(f"elemental time = {t}")

        # coordinates  of element under consideration
        dofx_def=[]; dofy_def=[]

        # displacement of element under consideration at t-2 time step
        ux_old_def = []; uy_old_def = []
        
        # displacement of element under consideration at t-1 time step
        ux_cur_def=[]; uy_cur_def=[]

        # mass of nodes of element under consideration
        nodal_mass = []

        # Integration time of nodes of element under consideration
        current_time = []

        # Energy lists
        ke_list = []
        pe_list = []
        te_list = []
        time_list = []

        for node in all_xnodes[cell]:            
            ta.vector()[node]= t
            ux_old_def.append(displacements['t_minus_2'].vector()[node])
            ux_cur_def.append(displacements['t_minus_1'].vector()[node])
            dofx_def.append(dof_coordinate[node, 0])
            nodal_mass.append(mass.vector()[node])
            current_time.append(ta.vector()[node])

        for node in all_ynodes[cell]:           
            ta.vector()[node] = t
            uy_old_def.append(displacements['t_minus_2'].vector()[node])
            uy_cur_def.append(displacements['t_minus_1'].vector()[node])
            dofy_def.append(dof_coordinate[node, 1])
            nodal_mass.append(mass.vector()[node])
            current_time.append(ta.vector()[node])
   
        # nodal displacement of current element based on disp at previous 2 timesteps
        u_old_def = [ux_old_def[0], uy_old_def[0], ux_old_def[1], uy_old_def[1], ux_old_def[2], uy_old_def[2]]
        u_def = [ux_cur_def[0], uy_cur_def[0], ux_cur_def[1], uy_cur_def[1], ux_cur_def[2], uy_cur_def[2]]
        dof_def = [dofx_def[0], dofy_def[0], dofx_def[1], dofy_def[1], dofx_def[2],
                        dofy_def[2]]
        
        print(f"u_old_def = {u_old_def} ")
        print(f"u_def = {u_def} ")
        print(f"dof_def = {dof_def} ")
        #input("Press Enter to continue...")
        
        # calling the differentiator
        u_newest = df.differentiate(u_old_def,u_def,dof_def, lmbda, mu, E, nu, nodal_mass, current_time) # a elemental displacement vector
        print(f"u_newest = {u_newest} ")

        for idx, node in enumerate(all_xnodes[cell]):
            if node not in bcdisp1:  # Check if node is not a boundary node
                displacements['t_current'].vector()[node] = u_newest[2 * idx]  # ux1, ux2, ux3

        for idx, node in enumerate(all_ynodes[cell]):
            if node not in bcdisp1:  # Check if node is not a boundary node
                displacements['t_current'].vector()[node] = u_newest[2 * idx + 1]  # uy1, uy2, uy3

        # Re-apply Dirichlet boundary conditions to ensure boundary nodes are fixed
        bc1.apply(displacements['t_current'].vector())

        velocity.assign(project((displacements['t_current'] - displacements['t_minus_1']) / t, V))
        KE = assemble(0.5 * rho * inner(velocity, velocity) * dx)
        PE = assemble(psi(displacements['t_current']) * dx)
        TE = PE + KE

        # Store energies
        ke_list.append(KE)
        pe_list.append(PE)
        te_list.append(TE)
        time_list.append(t)

        ## postprocessing
        fl.write(str(t) + "\t")
        fl.write(str(PE) + "\t")
        fl.write(str(KE) + "\t")
        fl.write(str(TE) + "\n")

        dt = t + el[0]
        tk[cell]= t
        updates.vector()[cell] = updates.vector()[cell] + 1.0

        print(f"dt = {dt}")
        if (dt<tf):
            dts.append((dt, cell))
            displacements['t_minus_2'].assign(displacements['t_minus_1'])
            displacements['t_minus_1'].assign(displacements['t_current'])

    # calculation of total number of updates
    total_updates=0.0
    for cell in range(mesh.num_cells()):
        total_updates=total_updates+updates.vector()[cell]
    print("total number of cell updates are",total_updates)  

    return ke_list, pe_list, te_list, time_list
    
def main():
    # converting the file to .xml file
    file = "meshes/square_lamina_smaller.msh"
    mesh = Mesh(mc.convertMesh(file, "xml"))
    ke_list, pe_list, te_list, time_list = AVI_explicit(mesh)

    plt.figure()
    plt.plot(time_list, ke_list, label='Kinetic Energy')
    plt.plot(time_list, pe_list, label='Potential Energy')
    plt.plot(time_list, te_list, label='Total Energy')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.title('Energy vs. Time')
    plt.legend()
    plt.show()
    

if __name__ == "__main__":
    # Material properties (All in SI units)
    E = 1.0e6                                               # Young's modulus (N/mm2)
    nu = 0.3                                                # Poission's ratio
    rho = 2.0                                               # density (g/mm3)
    mu = (E/(2.0*(1.0 + nu)))                               # mu-lames constant 
    lmbda = (E*nu/((1.0 + nu)*(1.0-2.0*nu)))                # lambda-lames constant
    f = 0.1                                                # time constant
    #f = 0.5                                                 # time constant
    c = np.sqrt(E/rho)                                      # speed volumetric waves in solids
    t0 = 0.0                                                # start time of the simulation
    tf = 0.004 + 1.0e-8                                     # end time of the simulation
    #tf = 0.00009                                   # end time of the simulation
    main()