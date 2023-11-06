import matplotlib as plt
import numpy as np
import sys
import time


# Bisection method for functions f: R -> R
def bisectionUniVariable(f, a, b, N, precision):
    f_a = f(a)
    f_b = f(b)
    prod = f_a*f_b
    i = 0
    if prod > precision:
        print(f"Bad input, f({a})f({b}) > 0.")
        
    elif np.absolute(np.linalg.norm(prod)) < precision:
        if np.absolute(np.linalg.norm(f_a)) < precision:
            print(f"Input was a zero, {a}")
        else:
            print(f"Input was a zero, {b}")
    
    else:
        for i in range(N):
            c = (a+b)/2
            f_c = f(c)

            if np.absolute(a-b) < precision:
                print(f"Bisection method converges to {c} after {i} iterations")
                return c
            if np.absolute(f_c) < precision:
                print(f"Found a zero for the given function at {c} after {i} iterations")
                return c

            elif f_a*f_c <= precision:
                b = c
                f_b = f_c
            else:
                a = c
                f_a = f_c

        print(f"No root was found after {i} iterations")


# Bisection method for functions f: R^n -> R
def bisectionMultiVariable(f, a, b, N, precision):
    f_a = f(*a)
    f_b = f(*b)
    prod = f_a*f_b
    i = 0
    if prod > precision:
        print(f"Bad input, f({a})f({b}) > 0.")
        
    elif np.absolute(np.linalg.norm(prod)) < precision:
        if np.absolute(np.linalg.norm(f_a)) < precision:
            print(f"Input was a zero, {a}")
        else:
            print(f"Input was a zero, {b}")
    
    else:
        for i in range(N):
            c = (a+b)/2
            f_c = f(*c)

            if np.absolute(np.linalg.norm(a-b)) < precision:
                print(f"Bisection method converges to {c} after {i} iterations")
                return c
            if np.absolute(np.linalg.norm(f_c)) < precision:
                print(f"Found a zero for the given function at {c} after {i} iterations")
                return c

            elif f_a*f_c <= precision:
                b = c
                f_b = f_c
            else:
                a = c
                f_a = f_c

        print(f"No root was found after {i} iterations")

#Newtons method for functions f: R -> R
def newtonsUniVariable(f, f_prime, x_0, N, precision):
    if np.absolute(np.linalg.norm(f(x_0))) < precision:
        print(f"Given input was a zero {f(x_0)}")
        return x_0
    

    new_x = x_0
    i = 0
    for i in range(N):
        if np.absolute(np.linalg.norm(f(new_x))) < precision:
            print(f"Found a zero for the given function at {new_x} after {i} iterations.")
            return x_0
        if f_prime(x_0) == 0:
            print(f"Error: the given function's derivative's value at {new_x} is 0, no root will be found.")
            break 
        else:
            new_x = new_x - (f(new_x)/f_prime(new_x)) #i love my girlfriend 
    
    print(f"No root was found after {i} iterations")

#Newtons method for functions f: R^n -> R
def newtonsMultiVariable(f, f_prime, x_0, N, precision):
    if np.absolute(np.linalg.norm(f(*x_0))) < precision:
        print(f"Given input was a zero {f(*x_0)}")
        return x_0
    

    new_x = x_0
    i = 0
    for i in range(N):
        if np.absolute(np.linalg.norm(f(*new_x))) < precision:
            print(f"Found a zero for the given function at {new_x} after {i} iterations.")
            return new_x
        if f_prime(*new_x) == 0:
            print(f"Error: the given function's derivative's value at {new_x} is 0, no root will be found.")
            break 
        if np.absolute(f(*new_x)) > precision**(-1):
            print(f"Guesses grew too large, no root found.")
            break
        else:
            new_x = new_x - (f(*new_x)/f_prime(*new_x)) #i love my girlfriend 
    
    print(f"No root was found after {i} iterations")

#Fixed point method for functions f: R -> R
def fixedPointUniVariable(f, x_0, N, precision):
    g = lambda x: f(x) + x
    new_x = x_0
    i = 0
    for i in range(N):
        evaluate = g(new_x)
        if np.absolute(np.linalg.norm(evaluate - new_x)) < precision:
            print(f"Found a zero for the given function at {new_x} after {i} iterations.")
            return x_0
        if np.absolute(np.linalg.norm(evaluate)) > precision**(-1):
            print(f"Guesses grew too large, no root found.")
            break
        else:
            new_x = evaluate #i love my girlfriend
        
    print(f"No root was found after {i} iterations")

#Fixed point method for functions f: R^n -> R
def fixedPointMultiVariable(f, x_0, N, precision):
    g = lambda x: f(*x) + x
    new_x = x_0
    i = 0
    for i in range(N):
        evaluate = g(new_x)
        if np.absolute(np.linalg.norm(evaluate - new_x)) < precision:
            print(f"Found a zero for the given function at {new_x} after {i} iterations.")
            return x_0
        if np.absolute(np.linalg.norm(evaluate)) > precision**(-1):
            print(f"Guesses grew too large, no root found.")
            break
        else:
            new_x = evaluate #i love my girlfriend
        
    print(f"No root was found after {i} iterations")


if __name__ == '__main__':
    plt.rcParams.update({'font.size': 16})

    
    #Set some parameters:
    precision = 10**(-6)
    maxIterations = 10**6
    x_0 = 1             #sets an initial guess


    #functions to test
    funcList = []       # holds the functions
    dfuncList = []      # holds the derivatives, if it exists

    #when creating functions, please declare the function, the derivative/gradient, name, if it is multi, and create a tuple from the information

    #a contraction mapping to show power of fixed point
    contraction = lambda x : 1/2 * (np.cos(x)**2)
    contraction_prime = lambda x: -1*np.sin(x)
    contraction_name = "1/2 * cos^2(x)"
    contraction_is_multi = False
    contraction_tuple = (contraction,contraction_prime,contraction_name, x_0, contraction_is_multi)


    #a simple quadratic function
    quad = lambda x,y: x**2 - 16 
    quad_prime = lambda x,y: 2*x
    quad_name = "x^2 - 16"
    quad_is_multi = False
    quad_tuple = (quad, quad_prime, quad_name, x_0 , quad_is_multi)


    #a simple quadratic function
    quad_multi = lambda x,y: x**2 - 16 + y
    quad_multi_prime = lambda x,y: [2*x, 1]
    quad_multi_name = "x^2 - 16 + y"
    quad_multi_is_multi = True
    quad_tuple = (quad_multi, quad_multi_prime, quad_multi_name, np.array([x_0, x_0]), quad_multi_is_multi)

    #add all functions to the list
    funcList.append(contraction_tuple)
    funcList.append(quad_tuple)
    

    for f, f_prime, f_name, initGuess, isMulti in funcList:
        print(f"\n\n\nRunning function: {f_name}\n------------------------------------------------")

        if isMulti:
             # run bisection method
            startTimeNewtons = time.time()
            bisectionMultiVariable(f, initGuess*(-1), initGuess*5, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Bisection's took {endTimeNewtons - startTimeNewtons}s \n")

            # run newton's method
            startTimeNewtons = time.time()
            newtonsMultiVariable(f,f_prime, initGuess, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Newton's took {endTimeNewtons - startTimeNewtons}s \n")

            # run fixed point method
            startTimeNewtons = time.time()
            fixedPointMultiVariable(f, initGuess, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Fixed point's took {endTimeNewtons - startTimeNewtons}s \n")
        else:
            # run bisection method
            startTimeNewtons = time.time()
            bisectionUniVariable(f, initGuess*(-1), initGuess*4, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Bisection's took {endTimeNewtons - startTimeNewtons}s \n")

            # run newton's method
            startTimeNewtons = time.time()
            newtonsUniVariable(f,f_prime, initGuess, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Newton's took {endTimeNewtons - startTimeNewtons}s \n")

            # run fixed point method
            startTimeNewtons = time.time()
            fixedPointUniVariable(f, initGuess, maxIterations, precision)
            endTimeNewtons = time.time()
            print(f"Fixed point's took {endTimeNewtons - startTimeNewtons}s \n")