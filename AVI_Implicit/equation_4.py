# Derivative of strain energy density wrt (u2)^y,i

def equation_4(A, m4, h4, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def):
    
    constant_4_u1_x_ip1 = (
        (A * (h4) * lmbda / 8) * (b1) * (d2) +
        (A * (h4) * mu / 4) * (d1) * (b2)
        )

    constant_4_u2_x_ip1 = (
        (A * (h4) * lmbda/ 4) * (b2) * (d2) +
        (A * (h4) * mu / 4) * (d2) * (b2)
        )

    constant_4_u3_x_ip1 = (
        (A * (h4) * lmbda/ 4) * (b3) * (d2) +
        (A * (h4) * mu / 4) * (d3) * (b2)
        )

    constant_4_u1_y_ip1 = (
        (m4 / (h4)) +
        (A * (h4) * lmbda/ 4) * (d1) * (d2) +
        (A * (h4) * mu / 2) * (b1) * (d2) +
        (A * (h4) * mu / 4) * (b1) * (b2)
        )

    constant_4_u2_y_ip1 = (
        (A * (h4) * lmbda/ 4) * (d2) * (d2) +
        (A * (h4) * mu / 2) * (b2) * (d2) +
        (A * (h4) * mu / 4) * (b2) * (b2)
        )

    constant_4_u3_y_ip1 = (
        (A * (h4) * lmbda/ 4) * (d3) * (d2) +
        (A * (h4) * mu / 2) * (b3) * (d2) +
        (A * (h4) * mu / 4) * (b3) * (b2)
        )


    rhs_4 = ((2*u_def[0])*(-constant_4_u1_x_ip1)) + \
            ((2*u_old_def[0])*(-constant_4_u1_x_ip1)) + \
            ((2*u_def[1])*((-constant_4_u1_y_ip1) + (3 * m4 / (h4)))) + \
            ((2*u_old_def[1])*(-constant_4_u1_y_ip1)) + \
            ((2*u_def[2])*(-constant_4_u2_x_ip1)) + \
            ((2*u_old_def[2])*(-constant_4_u2_x_ip1)) + \
            ((2*u_def[3])*(-constant_4_u2_y_ip1)) + \
            ((2*u_old_def[3])*(-constant_4_u2_y_ip1)) + \
            ((2*u_def[4])*(-constant_4_u3_x_ip1)) + \
            ((2*u_old_def[4])*(-constant_4_u3_x_ip1)) + \
            ((2*u_def[5])*(-constant_4_u3_y_ip1)) + \
            ((2*u_old_def[5])*(-constant_4_u3_y_ip1))
    
    return constant_4_u1_x_ip1, constant_4_u1_y_ip1, constant_4_u2_x_ip1, constant_4_u2_y_ip1, constant_4_u3_x_ip1, constant_4_u3_y_ip1, rhs_4