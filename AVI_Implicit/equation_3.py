# Derivative of strain energy density wrt (u2)^x,i

def equation_3(A, m3, h3, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def):

    constant_3_u1_x_ip1 = (
        (m3 / (h3)) +
        (A * (h3) * lmbda / 8) * (b1) * (b2)+
        (A * (h3) * mu / 2) * (b1) * (b2) +
        (A * (h3) * mu / 4) * (d1) * (d2)
        )
    
    constant_3_u2_x_ip1 = (
        (A * (h3) * lmbda/ 4) * (b2) * (b2)+
        (A * (h3) * mu / 2) * (b2) * (b2) +
        (A * (h3) * mu / 4) * (d2) * (d2)
        )
    
    constant_3_u3_x_ip1 = (
        (A * (h3) * lmbda/ 4) * (b3) * (b2)+
        (A * (h3) * mu / 2) * (b3) * (b2) +
        (A * (h3) * mu / 4) * (d3) * (d2)
        )
    
    constant_3_u1_y_ip1 = (
        (A * (h3) * lmbda/ 4) * (d1) * (b2)+
        (A * (h3) * mu / 4) * (b1) * (d2)
        )

    constant_3_u2_y_ip1 = (
        (A * (h3) * lmbda/ 4) * (d2) * (b2)+
        (A * (h3) * mu / 4) * (b2) * (d2)
        )

    constant_3_u3_y_ip1 = (
        (A * (h3) * lmbda/ 4) * (d3) * (b2)+
        (A * (h3) * mu / 4) * (b3) * (d2)
        )


    rhs_3 = ((2*u_def[0])*((-constant_3_u1_x_ip1) + (3 * m3 / (h3)))) + \
            ((2*u_old_def[0])*(-constant_3_u1_x_ip1)) + \
            ((2*u_def[1])*(-constant_3_u1_y_ip1)) + \
            ((2*u_old_def[1])*(-constant_3_u1_y_ip1)) + \
            ((2*u_def[2])*(-constant_3_u2_x_ip1)) + \
            ((2*u_old_def[2])*(-constant_3_u2_x_ip1)) + \
            ((2*u_def[3])*(-constant_3_u2_y_ip1)) + \
            ((2*u_old_def[3])*(-constant_3_u2_y_ip1)) + \
            ((2*u_def[4])*(-constant_3_u3_x_ip1)) + \
            ((2*u_old_def[4])*(-constant_3_u3_x_ip1)) + \
            ((2*u_def[5])*(-constant_3_u3_y_ip1)) + \
            ((2*u_old_def[5])*(-constant_3_u3_y_ip1))
    
    return constant_3_u1_x_ip1, constant_3_u1_y_ip1, constant_3_u2_x_ip1, constant_3_u2_y_ip1, constant_3_u3_x_ip1, constant_3_u3_y_ip1, rhs_3