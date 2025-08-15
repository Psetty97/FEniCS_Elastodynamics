# Derivative of strain energy density wrt (u1)^y,i

def equation_2(A, m2, h2, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def):

    constant_2_u1_x_ip1 = (
        (A * (h2) * lmbda / 8) * (b1) * (d1) +
        (A * (h2) * mu / 4) * (d1) * (b1)
        )

    constant_2_u2_x_ip1 = (
        (A * (h2) * lmbda/ 4) * (b2) * (d1) +
        (A * (h2) * mu / 4) * (d2) * (b1)
        )

    constant_2_u3_x_ip1 = (
        (A * (h2) * lmbda/ 4) * (b3) * (d1) +
        (A * (h2) * mu / 4) * (d3) * (b1)
        )

    constant_2_u1_y_ip1 = (
        (m2 / (h2)) +
        (A * (h2) * lmbda/ 4) * (d1) * (d1) +
        (A * (h2) * mu / 2) * (b1) * (d1) +
        (A * (h2) * mu / 4) * (b1) * (b1)
        )

    constant_2_u2_y_ip1 = (
        (A * (h2) * lmbda/ 4) * (d2) * (d1) +
        (A * (h2) * mu / 2) * (b2) * (d1) +
        (A * (h2) * mu / 4) * (b2) * (b1)
        )

    constant_2_u3_y_ip1 = (
        (A * (h2) * lmbda/ 4) * (d3) * (d1) +
        (A * (h2) * mu / 2) * (b3) * (d1) +
        (A * (h2) * mu / 4) * (b3) * (b1)
        )

    rhs_2 = ((2*u_def[0])*(-constant_2_u1_x_ip1)) + \
            ((2*u_old_def[0])*(-constant_2_u1_x_ip1)) + \
            ((2*u_def[1])*((-constant_2_u1_y_ip1) + (3 * m2 / (h2)))) + \
            ((2*u_old_def[1])*(-constant_2_u1_y_ip1)) + \
            ((2*u_def[2])*(-constant_2_u2_x_ip1)) + \
            ((2*u_old_def[2])*(-constant_2_u2_x_ip1)) + \
            ((2*u_def[3])*(-constant_2_u2_y_ip1)) + \
            ((2*u_old_def[3])*(-constant_2_u2_y_ip1)) + \
            ((2*u_def[4])*(-constant_2_u3_x_ip1)) + \
            ((2*u_old_def[4])*(-constant_2_u3_x_ip1)) + \
            ((2*u_def[5])*(-constant_2_u3_y_ip1)) + \
            ((2*u_old_def[5])*(-constant_2_u3_y_ip1))
    
    return constant_2_u1_x_ip1, constant_2_u1_y_ip1, constant_2_u2_x_ip1, constant_2_u2_y_ip1, constant_2_u3_x_ip1, constant_2_u3_y_ip1, rhs_2