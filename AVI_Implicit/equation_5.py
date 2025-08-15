# Derivative of strain energy density wrt (u3)^x,i

def equation_5(A, m5, h5, lmbda, mu, b1, b2, b3, d1, d2, d3, u_old_def, u_def):
    
    constant_5_u1_x_ip1 = (
        (m5 / (h5)) +
        (A * (h5) * lmbda / 8) * (b1) * (b3)+
        (A * (h5) * mu / 2) * (b1) * (b3) +
        (A * (h5) * mu / 4) * (d1) * (d3)
        )
    
    constant_5_u2_x_ip1 = (
        (A * (h5) * lmbda/ 4) * (b2) * (b3)+
        (A * (h5) * mu / 2) * (b2) * (b3) +
        (A * (h5) * mu / 4) * (d2) * (d3)
        )
    
    constant_5_u3_x_ip1 = (
        (A * (h5) * lmbda/ 4) * (b3) * (b3)+
        (A * (h5) * mu / 2) * (b3) * (b3) +
        (A * (h5) * mu / 4) * (d3) * (d3)
        )
    
    constant_5_u1_y_ip1 = (
        (A * (h5) * lmbda/ 4) * (d1) * (b3)+
        (A * (h5) * mu / 4) * (b1) * (d3)
        )

    constant_5_u2_y_ip1 = (
        (A * (h5) * lmbda/ 4) * (d2) * (b3)+
        (A * (h5) * mu / 4) * (b2) * (d3)
        )

    constant_5_u3_y_ip1 = (
        (A * (h5) * lmbda/ 4) * (d3) * (b3)+
        (A * (h5) * mu / 4) * (b3) * (d3)
        )


    rhs_5 = ((2*u_def[0])*((-constant_5_u1_x_ip1) + (3 * m5 / (h5)))) + \
            ((2*u_old_def[0])*(-constant_5_u1_x_ip1)) + \
            ((2*u_def[1])*(-constant_5_u1_y_ip1)) + \
            ((2*u_old_def[1])*(-constant_5_u1_y_ip1)) + \
            ((2*u_def[2])*(-constant_5_u2_x_ip1)) + \
            ((2*u_old_def[2])*(-constant_5_u2_x_ip1)) + \
            ((2*u_def[3])*(-constant_5_u2_y_ip1)) + \
            ((2*u_old_def[3])*(-constant_5_u2_y_ip1)) + \
            ((2*u_def[4])*(-constant_5_u3_x_ip1)) + \
            ((2*u_old_def[4])*(-constant_5_u3_x_ip1)) + \
            ((2*u_def[5])*(-constant_5_u3_y_ip1)) + \
            ((2*u_old_def[5])*(-constant_5_u3_y_ip1))
    
    return constant_5_u1_x_ip1, constant_5_u1_y_ip1, constant_5_u2_x_ip1, constant_5_u2_y_ip1, constant_5_u3_x_ip1, constant_5_u3_y_ip1, rhs_5