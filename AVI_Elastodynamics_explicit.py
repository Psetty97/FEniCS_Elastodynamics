# asynchronous variational integrators for elastodynamics

from fenics import *
import numpy as np
import heapq
import mesh_converter as mc
import differentiator as df
import os

# function to calculate strain
def epsilon(u):
    return sym(grad(u))

# function to calculate total strain energy
def psi(u):
    return 0.5*lmbda*(tr(epsilon(u)))**2+mu*inner(epsilon(u),epsilon(u))

def AVI_explicit(mesh):
    # definition of function spaces for velocity, displacement, nodal time
    V = VectorFunctionSpace(mesh, 'CG', 1)
    position = Function(V)     # position vector
    velocity = Function(V)     # velocity vector
    ta = Function(V)           # nodal time vector
    u_new = Function(V)        # displacement vector
    u_old = Function(V)        # displacement vector
    v = TestFunction(V)        # test function
    u = TrialFunction(V)       # trial function

    # mass lumping based on Gauas quadrature at the vertices
    M_lumped2 = assemble(rho * inner(v, u) * dx(scheme="vertex", metadata={"degree": 1, "representation": "quadrature"}))
    # splitting vector v in two components and making dof map
    x_comp, y_comp = V.split()
    dof_coordinate = V.tabulate_dof_coordinates()

    # priority que
    np.dts = []
    # nodal time vector
    ta.vector()[:] = 0.0
    # elemental time
    np.tk=[]

    # boolean function for fixed bottom
    def bottom(x, on_boundary):
        return near(x[1], 0.0) and on_boundary

        # boolean function for fixed bottom
    def top(x, on_boundary):
        return near(x[1], 1.0) and on_boundary

    # initial condition and boundary conditions
    a = 0.005
    disp = Expression(("0", "x[1]* a"), t=0, a=a, degree=2)
    disp2 = Expression(("0", "0"), t=0, a=a, lmbda=lmbda, mu=mu, rho=rho, degree=2)
    velo =  Expression(("0", "0"), t=0, a=a, lmbda=lmbda, mu=mu, rho=rho, degree=2)
    u_old = project(disp, V)
    velocity = project(velo, V)
    u_new.assign(u_old)
    bc1 = DirichletBC(V, disp2, bottom)
    np.bcdisp1 = []
    bcdisp1 = bc1.get_boundary_values().keys()

    # post-processing with .xdmf file
    xdmf_file = XDMFFile("elastodynamics-AVI.xdmf")
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True
    xdmf_file.parameters["rewrite_function_mesh"] = False

    # post-processing .txt file
    filespath = "/home/fenics/shared/Hamilton_AVI/Elastodynamics/"
    filename = "2d_AVI_elastodynamics"
    fl = open(filename + '.txt', 'w')
    if not os.path.exists(filespath):
        os.makedirs(filespath, exist_ok=True)

    # initialization of the problem
    for cell in range(mesh.num_cells()):
        # making a priority que from the time steps
        element = Cell(mesh, cell)
        dt = (f* element.inradius())/c
        np.tk.append(0.0)
        np.dts.append([dt, cell])

        xnodes = x_comp.dofmap().cell_dofs(cell)
        ynodes = y_comp.dofmap().cell_dofs(cell)

        # definition of initial position here
        for node in xnodes:
            position.vector()[node] = dof_coordinate[node, 0] + u_new.vector()[node]
        for node in ynodes:
            position.vector()[node] = dof_coordinate[node, 1] + u_new.vector()[node]

    # initial time
    t=0.0
    # discrete action
    Sd = 0.0

    ## postprocessing
    xdmf_file.write(u_new, t)
    xdmf_file.write(velocity, t)
    PE = assemble(psi(u_new) * dx)
    KE = assemble(rho * inner(velocity, velocity) * dx)
    TE = PE + KE
    fl.write(str(t) + "\t")
    fl.write(str(PE) + "\t")
    fl.write(str(KE) + "\t")
    fl.write(str(TE) + "\n")

    while(len(np.dts)>0):
        # making heap from the previous array
        heapq.heapify(np.dts)
        # ascending order of time steps by sorting
        np.dts.sort()
        # pop an element with highest priority from queue
        el = heapq.heappop(np.dts)
        cell=el[1]
        element = Cell(mesh, cell)

        # updating time if t<=tf || at the end just interpolate to time tf
        if(t<=tf):
            t = el[0]
        else:
            t=t
        print(t)

        # updating positions
        xnodes = x_comp.dofmap().cell_dofs(cell)
        ynodes = y_comp.dofmap().cell_dofs(cell)

        # displacement of element under consideration
        np.ux_def=[]; np.uy_def=[]
        # coordinates  of element under consideration
        np.dofx_def=[]; np.dofy_def=[]

        for node in xnodes:
            position.vector()[node] = position.vector()[node] + velocity.vector()[node] * (t-ta.vector()[node])
            u_new.vector()[node] = u_new.vector()[node] + velocity.vector()[node] * (t-ta.vector()[node])
            ta.vector()[node]= t
            np.ux_def.append(u_new.vector()[node])
            np.dofx_def.append(dof_coordinate[node, 0])

        for node in ynodes:
            position.vector()[node] = position.vector()[node] + velocity.vector()[node] * (t-ta.vector()[node])
            u_new.vector()[node] = u_new.vector()[node] + velocity.vector()[node] * (t-ta.vector()[node])
            ta.vector()[node] = t
            np.uy_def.append(u_new.vector()[node])
            np.dofy_def.append(dof_coordinate[node, 1])

        if(t<tf):
            # post-processing
            xdmf_file.write(u_new, t)
            xdmf_file.write(velocity, t)
            PE = assemble(psi(u_new) * dx)
            KE = assemble(rho*inner(velocity, velocity)*dx)
            TE = PE + KE
            fl.write(str(t) + "\t")
            fl.write(str(PE) + "\t")
            fl.write(str(KE) + "\t")
            fl.write(str(TE) + "\n")

        if((t>tf) and len(np.dts)==0):
            xdmf_file.write(u_new, t)
            xdmf_file.write(velocity, t)
            PE = assemble(psi(u_new) * dx)
            KE = assemble(rho*inner(velocity, velocity)*dx)
            TE = PE + KE
            fl.write(str(t) + "\t")
            fl.write(str(PE) + "\t")
            fl.write(str(KE) + "\t")
            fl.write(str(TE) + "\n")

        # updating velocities if t <= tf : calculate velocity for the next update
        if(t<tf):
            # nodal displacement of changed element to calculate the increment
            np.u_def = [np.ux_def[0], np.uy_def[0], np.ux_def[1], np.uy_def[1], np.ux_def[2], np.uy_def[2]]
            np.dof_def = [np.dofx_def[0], np.dofy_def[0], np.dofx_def[1], np.dofy_def[1], np.dofx_def[2],
                          np.dofy_def[2]]

            # calling the differentiator
            psi_elem, np.f_int = df.diffrentiate(np.u_def, np.dof_def, lmbda, mu, E, nu)
            PE_elem = psi_elem*(1.0/(4.0 *element.volume()))                      # integral of psi*dx over the element
            f_elem = np.f_int                                                     # actually a nodal force vector

            # empty array nodal forces
            f_nodal = np.zeros(1000000)
            i=0; j=1
            for node in xnodes:
                f_nodal[node]= f_elem[i]*(1.0/(4.0 *element.volume()))
                i=i+2
            for node in ynodes:
                f_nodal[node] = f_elem[j]*(1.0/(4.0 *element.volume()))
                j=j+2

            for node in xnodes:
                if (node in bcdisp1):
                    velocity.vector()[node] = velocity.vector()[node]
                else:
                    velocity.vector()[node] = velocity.vector()[node]\
                                              -(1.0/(M_lumped2.array()[node, node]))*\
                                              (t-np.tk[cell])*f_nodal[node]
            for node in ynodes:
                if (node in bcdisp1):
                    velocity.vector()[node] = velocity.vector()[node]
                else:
                    velocity.vector()[node] = velocity.vector()[node]\
                                              -(1.0/(M_lumped2.array()[node, node]))*\
                                              (t-np.tk[cell])*f_nodal[node]

            dt = t + el[0]
            np.tk[cell]= t
            # appending instead of pushing
            np.dts.append([dt, cell])
    PE = assemble(psi(u_new) * dx)
    print(PE)

def main():
    # converting the file to .xml file
    file = "/home/fenics/shared/Hamilton_AVI/Elastodynamics/meshes/square_lamina_0.1.msh"
    mesh = Mesh(mc.convertMesh(file, "xml"))
    AVI_explicit(mesh)

if __name__ == "__main__":
    # Material properties (All in SI units)
    E = 210.0e6                                             # Young's modulus (N/m2)
    nu = 0.3                                                # Poission's ratio
    rho = 10000.0                                           # density (Kg/m3)
    mu = (E/(2.0*(1.0 + nu)))                               # mu-lames constant
    lmbda = (E*nu/((1.0 + nu)*(1.0-2.0*nu)))                # lambda-lames constant
    f = 0.1                                                 # time constant
    c = np.sqrt((lmbda+2*mu)/rho)                           # speed volumetric waves in solids (m/s)
    t0 = 0.0                                                # start time of the simulation
    tf = 0.1 + 1.0e-8                                       # end time of the simulation
    main()