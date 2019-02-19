""" @package examples/shocktube_1d/Riemann.py
exact solution to Riemann problem
    
created by Rainer Weinberger, last modified: 19.02.2019
"""

import numpy as np


def f_K(pres, W_K, gamma):
    """
    Flux functions for left or right state of a Riemann problem;
    Auxiliary function for NewtonRaphson.
    
    \param[in] pres float pressure of central state
    \param[in] W_K array, shape (3,), left/right hand side state of Riemann problem
    \param[in] gamma float, adiabatic index of gas
    """
    integy_K = W_K[2] / ( gamma - 1.0 ) / W_K[0]
    a_K = np.sqrt(  gamma * W_K[2] / W_K[0] )
    A = 2.0 / ( gamma + 1.0 ) / W_K[0]
    B = ( gamma - 1.0 ) / ( gamma + 1.0 ) * W_K[2]
    
    if pres > W_K[2] :    ##shock
        return ( pres - W_K[2] ) * np.sqrt( A / ( pres + B ) )
    else:    ## rarefaction
        exponent = ( gamma - 1.0) / 2.0 / gamma
        return ( 2.0 * a_K / ( gamma - 1.0 ) ) * ( ( pres / W_K[2] )**exponent - 1.0 )

def f_K_prime(pres, W_K, gamma_K):
    """
    Derivative of flux functions for left or right state of a Riemann problem;
    Auxiliary function for NewtonRaphson.
    
    \param[in] pres float pressure of central state
    \param[in] W_K array, shape (3,), left/right hand side state of Riemann problem
    \param[in] gamma float, adiabatic index of gas
    """
    integy_K = W_K[2] / ( gamma_K - 1.0 ) / W_K[0]
    a_K = np.sqrt( gamma_K *  W_K[2] / W_K[0] )
    A = 2.0 / ( gamma_K + 1.0 ) / W_K[0]
    B = ( gamma_K - 1.0 ) / ( gamma_K + 1.0 ) * W_K[2]
    
    if pres > W_K[2] :    ##shock
        return np.sqrt( A / ( B + pres ) ) * ( 1.0 - ( pres - W_K[2] ) / ( 2.0 * ( B + pres ) ) )
    else:    ## rarefaction
        exponent = -( gamma_K + 1.0 ) / 2.0 / gamma_K
        return ( pres / W_K[2] )**exponent / W_K[0] / a_K

def NewtonRaphson(pres, W_L, W_R, gamma):
    """
    NewtonRaphson iteration to determine central pressure of a Riemann problem 

    \param[in] W_L array, shape (3,) left hand side density, velocity and pressure
    \param[in] W_R array, shape (3,) right hand side density, velocity and pressure
    \param[in] gamma float, adiabatic index of gas
    
    \return pres pressure of central region
    """
    change = 1.0
    iter = 0
    while np.abs(change) > 1e-6:
        iter += 1
        ## root finding for function
        function = f_K(pres, W_L, gamma) + f_K(pres, W_R, gamma) + W_R[1] - W_L[1]
        function_prime = f_K_prime(pres, W_L, gamma) + f_K_prime(pres, W_R, gamma)
        pres_new = pres - function / function_prime
        change = 2.0 * ( pres_new - pres ) / ( pres_new + pres )
        pres = pres_new
        if(iter > 100):
            print("ERROR: NewtonRaphson: maximum number of iterations reached! Could not converge!")
            break
    return pres
    
def rarefaction_fan_left(W_K, gamma, xx, time):
    gm1 = (gamma - 1.0)
    gp1 = (gamma + 1.0)
    gamma_ratio = gm1 / gp1
    a = np.sqrt( gamma * W_K[2] / W_K[0] )
    rho = W_K[0] * (2.0 / gp1 + gamma_ratio / a * ( W_K[1] - xx / time ) )**( 2.0 / gm1 )## wrong?
    vel = 2.0 / gp1 * ( a + gm1 / 2.0 * W_K[1] + xx / time )
    pres = W_K[2] * ( 2.0 / gp1 + gamma_ratio / a * ( W_K[1] - xx / time ) )**( 2.0 * gamma / gm1 )
    return np.array([rho, vel, pres]).T
        
def rarefaction_fan_right(W_K, gamma, xx, time):
    gm1 = (gamma - 1.0)
    gp1 = (gamma + 1.0)
    gamma_ratio = gm1 / gp1
    a = np.sqrt( gamma * W_K[2] / W_K[0] )
    rho = W_K[0] * (2.0 / gp1 - gamma_ratio / a * ( W_K[1] - xx / time ) )**( 2.0 / gm1 )
    vel = 2.0 / gp1 * ( -a + gm1 / 2.0 * W_K[1] + xx / time )
    pres = W_K[2] * ( 2.0 / gp1 - gamma_ratio / a * ( W_K[1] - xx / time ) )**( 2.0 * gamma / gm1 )
    return np.array([rho, vel, pres]).T
        
def RiemannProblem(xx, x0, W_L, W_R, gamma, time):
    """
    RiemannProblem: Samples the solution of a Riemann problem at a given time at positions xx
    
    \param[in] xx, array, shape (n,), 1d coordinate of grid
    \param[in] x0, float initial coordinate of discontinuety
    \param[in] W_L, array, shape(3,), density, velocity and pressure of left state
    \param[in] W_R, array, shape(3,), density, velocity and pressure of right state
    \param[in] gamma, float, adiabatic index of gas
    \param[in] time, time since initial conditions
    
    \return xx, W, x_shock; xx: as input; W: array, shape(n,3) densiy, velocity and pressure at pos xx
            x_shock: postition of charactaeristic features
    """
    
    ## some auxiliary quantitites
    a_L = np.sqrt( gamma *  W_L[2] / W_L[0] )    ## sound speeds
    a_R = np.sqrt( gamma *  W_R[2] / W_R[0] )
    gm1 = gamma - 1.0    ## gamma +- 1
    gp1 = gamma + 1.0
    gamma_ratio = gm1 / gp1
    PosOfCharacteristics = np.zeros(5) ## 0: left most to 4: right most; 2 is rarefaction; if 0 and 1 (3 and 4) have the same value, there is a shock, otherwise rarefaction wave

    ## calculate central pressure; initial guess: arithmetic mean
    pstar = 0.5 * ( W_L[2] + W_R[2] )
    pstar = NewtonRaphson( pstar, W_L, W_R, gamma)
    
    ## calculate central velocity and sound speed
    ustar = 0.5 * ( W_L[1] + W_R[1] ) + 0.5 * ( f_K(pstar, W_R, gamma) - f_K(pstar, W_L, gamma) )
    astar_L = a_L * ( pstar / W_L[2] )**( (gamma - 1.0) / (2.0 * gamma) )
    astar_R = a_R * ( pstar / W_R[2] )**( (gamma - 1.0) / (2.0 * gamma) )
    
    ## decide on shock or rarefaction on left and right; get velocities of features
    PosOfCharacteristics[2] = ustar 
    if pstar > W_L[2]:    ## shock on left
        PosOfCharacteristics[0] = W_L[1] - a_L * np.sqrt( gp1 / 2.0 / gamma * pstar / W_L[2] + gm1 / 2.0 / gamma )
        PosOfCharacteristics[1] = PosOfCharacteristics[0]
        leftShock = True
    else:    ## rarefaction wave on left
        PosOfCharacteristics[0] = W_L[1] - a_L ## outer
        PosOfCharacteristics[1] = ustar - astar_L ## inner   
    if pstar > W_R[2]:    ## shock on right
        PosOfCharacteristics[4] = W_R[1] + a_R * np.sqrt( gp1 / 2.0 / gamma * pstar / W_R[2] + gm1 / 2.0 / gamma )
        PosOfCharacteristics[3] = PosOfCharacteristics[4]
        rightShock = True
    else:    ## rarefaction wave on right
        PosOfCharacteristics[4] = W_R[1] + a_R ## outer
        PosOfCharacteristics[3] = ustar + astar_R ## inner
    
    ## self-similar, i.e. pos = velocity * time
    PosOfCharacteristics[:] *= time
    
    ## select different regions
    i_L, = np.where(xx-x0 < PosOfCharacteristics[0])
    i_star_L, = np.where( (xx-x0 >= PosOfCharacteristics[1]) & (xx-x0 < PosOfCharacteristics[2]) )
    i_star_R, = np.where( (xx-x0 >= PosOfCharacteristics[2]) & (xx-x0 < PosOfCharacteristics[3]) )
    i_R, = np.where(xx-x0 >= PosOfCharacteristics[4])
    
    ## sample solution
    W = np.zeros( [len(xx), 3], dtype=np.float64 )
    ## left most and right most state: same as initial conditions
    W[i_L,:] = W_L[:]
    W[i_R, :] = W_R[:]

    if PosOfCharacteristics[0] != PosOfCharacteristics[1]: ## rarefaction wave on left
        i_fan_L, = np.where( (xx-x0 >= PosOfCharacteristics[0]) & (xx-x0 < PosOfCharacteristics[1]) )
        W[i_fan_L,:] = rarefaction_fan_left(W_L, gamma, xx[i_fan_L]-x0, time)    
        W[i_star_L, 0] = (pstar / W_L[2])**(1./gamma) * W_L[0]
    else: ## shock on left
        W[i_star_L, 0] = W_L[0] * ( ( pstar / W_L[2] + gamma_ratio ) / ( gamma_ratio * pstar / W_L[2] + 1.0 ) )
    
    W[i_star_L, 1] = ustar
    W[i_star_L, 2] = pstar

    if PosOfCharacteristics[3] != PosOfCharacteristics[4]:  ## rarefaction wave on right
        i_fan_R, = np.where( (xx-x0 >= PosOfCharacteristics[3]) & (xx-x0 < PosOfCharacteristics[4]) )
        W[i_fan_R, :] = rarefaction_fan_right(W_R, gamma, xx[i_fan_R]-x0, time)
        W[i_star_R, 0] = (pstar / W_R[2])**(1./gamma) * W_R[0]
    else: ## shock on right
        W[i_star_R, 0] = W_R[0] * ( ( pstar / W_R[2] + gamma_ratio ) / ( gamma_ratio * pstar / W_R[2] + 1.0 ) )

    W[i_star_R, 1] = ustar
    W[i_star_R, 2] = pstar
    
    return xx, W, PosOfCharacteristics
