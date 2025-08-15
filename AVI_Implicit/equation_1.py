# Derivative of strain energy density wrt (u1)^x,i
def equation_1(A, m1, h1, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def):

    constant_1_u1_x_ip1 = (
        (m1 / (h1)) +
        (A * (h1) * lmbda / 8) * (b1) * (b1)+
        (A * (h1) * mu / 2) * (b1) * (b1) +
        (A * (h1) * mu / 4) * (d1) * (d1)
        )
    
    constant_1_u2_x_ip1 = (
        (A * (h1) * lmbda/ 4) * (b2) * (b1)+
        (A * (h1) * mu / 2) * (b2) * (b1) +
        (A * (h1) * mu / 4) * (d2) * (d1)
        )
    
    constant_1_u3_x_ip1 = (
        (A * (h1) * lmbda/ 4) * (b3) * (b1)+
        (A * (h1) * mu / 2) * (b3) * (b1) +
        (A * (h1) * mu / 4) * (d3) * (d1)
        )

    constant_1_u1_y_ip1 = (
        (A * (h1) * lmbda/ 4) * (d1) * (b1)+
        (A * (h1) * mu / 4) * (b1) * (d1)
        )

    constant_1_u2_y_ip1 = (
        (A * (h1) * lmbda/ 4) * (d2) * (b1)+
        (A * (h1) * mu / 4) * (b2) * (d1)
        )

    constant_1_u3_y_ip1 = (
        (A * (h1) * lmbda/ 4) * (d3) * (b1)+
        (A * (h1) * mu / 4) * (b3) * (d1)
        )

    rhs_1 = ((2*u_def[0])*((-constant_1_u1_x_ip1) + (3 * m1 / (h1)))) + \
            ((2*u_old_def[0])*(-constant_1_u1_x_ip1)) + \
            ((2*u_def[1])*(-constant_1_u1_y_ip1)) + \
            ((2*u_old_def[1])*(-constant_1_u1_y_ip1)) + \
            ((2*u_def[2])*(-constant_1_u2_x_ip1)) + \
            ((2*u_old_def[2])*(-constant_1_u2_x_ip1)) + \
            ((2*u_def[3])*(-constant_1_u2_y_ip1)) + \
            ((2*u_old_def[3])*(-constant_1_u2_y_ip1)) + \
            ((2*u_def[4])*(-constant_1_u3_x_ip1)) + \
            ((2*u_old_def[4])*(-constant_1_u3_x_ip1)) + \
            ((2*u_def[5])*(-constant_1_u3_y_ip1)) + \
            ((2*u_old_def[5])*(-constant_1_u3_y_ip1))
    
    return constant_1_u1_x_ip1, constant_1_u1_y_ip1, constant_1_u2_x_ip1, constant_1_u2_y_ip1, constant_1_u3_x_ip1, constant_1_u3_y_ip1, rhs_1