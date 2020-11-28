"""
    Python script that solves the Blasius Boundary Value Problem using Euler Method 
    (first order upwind differencing), a third-order nonlinear ODE
        2*f"' + f*f" = 0
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor

import numpy as np 
import math as m 
import matplotlib.pyplot as plt



def blasius_eqn_sovlver(acc):   
    """
        Function that solves the Blasius Boundary Value Problem

        Args:
            acc: float; increment of the self similar varaible
         
        Returns:
            f: ndarry(eta_limit/acc + 1, ); func in Blasius Equation
            f_dash: ndarry(eta_limit/acc + 1, ); first derivative of func in Blasius Equation 
            f_double_dash: ndarry(eta_limit/acc + 1, ); second derivative of func in Blasius 
                                                            Equation 
    """

    # Setting the limit for eta 
    eta_limit = 10

    # Creating required arrays for the 3 First Order ODEs
    f = np.zeros(int(eta_limit/acc) + 1)
    f_dash = np.zeros(int(eta_limit/acc) + 1) 
    f_double_dash= np.zeros(int(eta_limit/acc) + 1) 

    # Setting iteration count 
    count = 1

    # Loop that iterate until convergence of the shooting point method
    while ( abs(f_dash[-1] - 1.0) > 1e-5 ):
    
        # Random guesses required for first 2 iterations of the shooting method 
        if ( count == 1):
            f_double_dash[0] = 0.5
        elif ( count == 2):
            f_double_dash[0] = 0.4
        else:
            # Using Newton's Method for approximating value for shooting to get faster convergence
            f_double_dash[0] = f_double_dash[0] - (f_dash[-1] - 1.0) \
                                * (f_double_dash_for_shooting_1 - f_double_dash_for_shooting_2 ) \
                                / ( f_dash_for_shooting_1 - f_dash_for_shooting_2 ) 
        # Obtaining the right value required for shooting
        if ( count >= 2):
            f_dash_for_shooting_2 = f_dash_for_shooting_1 
            f_double_dash_for_shooting_2 = f_double_dash_for_shooting_1

        f_dash_for_shooting_1 = f_dash[-1] 
        f_double_dash_for_shooting_1 = f_double_dash[0]

        # Loop that marches along increasing the self similar variable, eta 
        for i, ignore_val in enumerate(f[1:]):
            f[i+1] = f[i] + f_dash[i] * acc
            f_dash[i+1] = f_dash[i] + f_double_dash[i] * acc
            f_double_dash[i+1] = f_double_dash[i] - f[i] * f_double_dash[i] * acc / 2

        # Calculating residuals
        print("Residual at ", count, ": ", round(abs(f_dash[-1] - 1.0),6))

        # Incrementing iteration count 
        count += 1

        if ( count >= 10 ):
            print("\n\nITERATION LIMIT EXCEEDED")
            break
    

    # Return statement 
    return f, f_dash, f_double_dash


def post_processing(f, f_dash, f_double_dash):
    """
        Function that is used to obtain the desired results in the form of graph
        or print statements

        Args:
            f: ndarry(eta_limit/acc + 1, ); func in Blasius Equation
            f_dash: ndarry(eta_limit/acc + 1, ); first derivative of func in Blasius Equation 
            f_double_dash: ndarry(eta_limit/acc + 1, ); second derivative of func in Blasius 
                                                            Equation 
            
        Returns:
    """
    # Calculating the wall shear stress

    # Plotting the graphs
    eta_vec = np.linspace(0, 10, num=101)
    plt.plot(eta_vec, f, eta_vec, f_dash, eta_vec, f_double_dash)
    plt.show()
     


# Function call
f, f_dash, f_double_dash = blasius_eqn_sovlver(0.1)

# Sending the array for post processing 
post_processing(f, f_dash, f_double_dash)
	
