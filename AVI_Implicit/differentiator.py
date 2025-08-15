# differentiation function to calculate the increment of velocity and elemental potential energy
import numpy as np
import equation_1 as eq1
import equation_2 as eq2
import equation_3 as eq3
import equation_4 as eq4
import equation_5 as eq5
import equation_6 as eq6

def solve_with_regularization(C, D, alpha=1e-10):
    """
    Attempts to solve the linear system C * E = D.
    If solving fails, adds a small value `alpha` to the diagonal elements of C to attempt to stabilize the matrix and try solving again.
    :param C: Coefficient matrix
    :param D: Right-hand side vector
    :param alpha: Regularization parameter
    :return: Solution vector E
    """
    try:
        # First attempt to solve the system normally
        E = np.linalg.solve(C, D)
        return E
    except np.linalg.LinAlgError:
        # If normal solving fails, apply regularization and try again
        regularized_C = C + alpha * np.eye(C.shape[0])
        try:
            E = np.linalg.solve(regularized_C, D)
            return E
        except np.linalg.LinAlgError:
            # If regularization also fails, raise an exception
            raise ValueError("Regularized system is still singular")

def differentiate(u_old_def,u_def,dof_def, lmbda, mu, E, nu, nodal_mass,current_time):

    n_array = np.array([[1.0, dof_def[0], dof_def[1]], [1.0, dof_def[2], dof_def[3]],
                        [1.0, dof_def[4], dof_def[5]]])

    A = 0.5 * np.linalg.det(n_array)

    # derivatives of shape function  x=0,2,4 y=1,3,5 ----> Coefficients of B matrix
    b1 = (dof_def[3] - dof_def[5])*(1.0/(2.0*A)) #(dN1/dx)
    b2 = (dof_def[5] - dof_def[1])*(1.0/(2.0*A)) # (dN2/dx)
    b3 = (dof_def[1] - dof_def[3])*(1.0/(2.0*A)) #(dN3/dx)
    d1 = (dof_def[4] - dof_def[2])*(1.0/(2.0*A)) #(dN1/dy)
    d2 = (dof_def[0] - dof_def[4])*(1.0/(2.0*A)) #(dN2/dy)
    d3 = (dof_def[2] - dof_def[0])*(1.0/(2.0*A)) #(dN3/dy)

    # strains
    ex = (b1 * u_def[0] + b2 * u_def[2] + b3 * u_def[4])
    ey = (d1 * u_def[1] + d2 * u_def[3] + d3 * u_def[5])
    exy = 0.5*(d1 * u_def[0] + d2 * u_def[2] + d3 * u_def[4] + b1 * u_def[1] + b2 * u_def[3] + b3 * u_def[5])
    tr_eps = ex + ey

    # strain tensor
    #eps = np.array([[ex, exy], [exy, ey]])

    # strain energy of an element
    #psi_elem = (0.5 * lmbda * (tr_eps) ** 2 + mu * np.tensordot(eps, eps, axes=2))*np.abs(A)
    m1 = nodal_mass[0]
    h1 = current_time[0]
    m2 = nodal_mass[1] 
    h2 = current_time[1]
    m3 = nodal_mass[2] 
    h3 = current_time[2]
    m4 = nodal_mass[3] 
    h4 = current_time[3]
    m5 = nodal_mass[4] 
    h5 = current_time[4]
    m6 = nodal_mass[5] 
    h6 = current_time[5]
    
    a11, a12, a13, a14, a15, a16, rhs_1 = eq1.equation_1(A, m1, h1, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)
    a21, a22, a23, a24, a25, a26, rhs_2 = eq2.equation_2(A, m2, h2, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)
    a31, a32, a33, a34, a35, a36, rhs_3 = eq3.equation_3(A, m3, h3, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)
    a41, a42, a43, a44, a45, a46, rhs_4 = eq4.equation_4(A, m4, h4, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)
    a51, a52, a53, a54, a55, a56, rhs_5 = eq5.equation_5(A, m5, h5, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)
    a61, a62, a63, a64, a65, a66, rhs_6 = eq6.equation_6(A, m6, h6, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def)

    C = np.array([
    [a11, a12, a13, a14, a15, a16],
    [a21, a22, a23, a24, a25, a26],
    [a31, a32, a33, a34, a35, a36],
    [a41, a42, a43, a44, a45, a46],
    [a51, a52, a53, a54, a55, a56],
    [a61, a62, a63, a64, a65, a66]])
                              
    D = np.array([rhs_1, rhs_2, rhs_3, rhs_4, rhs_5, rhs_6])

    #input("reached solver")
    E = solve_with_regularization(C, D)

    return E