# Computation Methods in Nonlinear Optimization
# PS2
# David Harar


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#### DEFINE GRADIENT HESSIAN FOR FOR THE Rosenbrock function

# define the function and gradients:
def fxy(x,y): return 100 * np.square(y-np.square(x)) + np.square(1-x)
# gradients
def dx(x,y): return -400 * x * (y-np.square(x)) - 2 *(1 - x)
def dy(x,y): return 200 * (y-np.square(x))
# hessian
def dxx(x,y): return 1200 * np.square(x) - 400 * y + 2
def dxy(x,y): return -400 * x
def dyy(x,y): return 200
def dyx(x,y): return -400 * x

# 1. Investigate the performance of the steepest descent algorithm (gradient method)
# for the function f(x,y) = 100 * (y-x^2)^2 + (1-x)^2
# (a) Write a code for the steepest decent algorithm

def grad(a):
    """
    a - a column vector of shape (1,2)
    """
    x = a[0][0]
    y = a[1][0]
    return np.array([[dx(x,y)],[dy(x,y)]])
def Hess(a):
    """
    a - a column vector of shape (1,2)
    """
    x = a[0][0]
    y = a[1][0]
    return np.array([[dxx(x,y) , dxy(x,y)],[dyx(x,y), dyy(x,y)]])

# We have that phi(t) is a function of only one variable and that it is convex.
# Therefore we will use the newthon method for a function with one input
# and find the root, t

def phi(t,val): return fxy((val - t*grad(val))[0][0],(val - t*grad(val))[1][0]) # Given values, we need to search for t = argmin
def phi_tag(t, val): return np.sum((-1) * grad((val - t * grad(val))) * grad(val))

def phi_tag2(t, val):
    HESS = Hess(val).T.sum(axis = 0).reshape((2,1))
    phi_tt = grad(val) * grad(val) * HESS
    return -np.sum(phi_tt)

def phi_bisec(val, a = 0, b = 7, eps = 0.00001): # works only localy. drop this method and turn to Newthon.
    fa = phi_tag(a, val)
    fb = phi_tag(b, val)
    d = (a+b)/2
    fd = phi_tag(d, val)

    while abs(b-a) > eps:

        if fa*fd < 0:
            b = d
            d = (a+b)/2
            fd = phi_tag(d,val)

            interval = abs(a-b) # corrent interval

        if fd*fb < 0:
            a = d
            d = (a+b)/2
            fd = phi_tag(d,val)

            interval = abs(a-b) # corrent interval
    return d


def sdescent(tol, maxiter, a0 = np.array([0.8,0.5]).reshape((2,1)), prnt = False):
    old_val = a0
    values_fxy = {}
    values_params = {}
    values_fxy['0'] = fxy(old_val[0][0], old_val[1][0])
    values_params['0'] = a0.reshape((1,2))

    for i in range(0,maxiter):
        # minimize phi to find t

        t = phi_bisec(old_val)

        # update a's
        new_val = old_val - t * grad(old_val)

        # stopping rule
        #if np.linalg.norm(new_val-old_val) <= tol: break
        if fxy(old_val[0][0], old_val[1][0]) <= tol: break # this worked better for me

        # update values
        old_val = new_val

        # store values
        values_fxy["{0}".format(i+1)]=fxy(old_val[0][0], old_val[1][0])
        values_params["{0}".format(i+1)] = old_val.reshape((1,2))

        # print
        if prnt == True:
            if i % 10 == 0: print(i)
    return new_val, min(i, maxiter), values_fxy, values_params

roots, iters, values_fxy, values_params = sdescent(tol=0.000001, maxiter=10000)

# turn values into an array
values_list = list(values_fxy.values())

print("x = ", roots[0][0])
print("y = ", roots[1][0])
print("number of iterations until convergence:", iters)
print("function value:", values_list[-1])

# (b) Plot the function and the points produced by the algorithm.
fig = plt.figure()
ax = fig.gca(projection='3d')

PARAMS = np.array(list(values_params.values()))
PARAMS = PARAMS.reshape((PARAMS.shape[0], PARAMS.shape[2]))
X = PARAMS[:,0]
Y = PARAMS[:,1]
X, Y = np.meshgrid(X, Y)
Z = fxy(x = X, y = Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title("Steepest Descent Algotithm")
plt.savefig("PS2_203670872_1b")

# (c) What can you say about the speed of convergence?
'''
It seems that the algorithm converged in 911 iterations. The steepest descent
algorithm is considered relatively slow and therefore it used seldom in the
literature. The reason for that is that in every iteration we have to search
for the steepest direction, alpha. When the problem has a quadratic form there
is a closed solution for alpha which spare us from numericly search for it.
In any other case, which we see more often, the problem is not quadratic then
in each iteration we have to search numericly for the alpha that minimizes the
next step. Consequentaly, even though the graph in (b) shows that it took about
911 iterations to converge, the real number of iterations was much greater.
'''


# 2. Newthon method
def Newthon_grad(a): return a - (np.linalg.inv(Hess(a))).dot(grad(a))

def newthon(tol, maxiter, a0 = np.array([0.8,0.5]).reshape((2,1))):
    """
    newthon impements the Newthon method for a function of two inputs.
    inputs:
    - a0 - initial values, a column vector includes x and y
    - tol - tolerance, tells newthon to stop if ||new_val - old_val||<tol
    - maxiter - number of iterations
    output:
    - Yn - a column vector of the roots of the function.
    """
    old_val = a0
    values_fxy = {}
    values_params = {}
    values_fxy['0'] = fxy(old_val[0][0], old_val[1][0])
    values_params['0'] = a0.reshape((1,2))

    for i in range(0,maxiter):
        new_val = Newthon_grad(old_val)

        # store values
        values_fxy["{0}".format(i+1)]=fxy(old_val[0][0], old_val[1][0])
        values_params["{0}".format(i+1)] = old_val.reshape((1,2))

        # update values
        old_val = new_val

        # stopping rule
        #if np.linalg.norm(new_val-old_val) <= tol: break
        if fxy(old_val[0][0], old_val[1][0]) <= tol: break

    return new_val, min(i, maxiter), values_fxy, values_params

roots, iters, values_fxy, values_params = newthon(tol=0.000001, maxiter=10000)

print("x = ", roots[0][0])
print("y = ", roots[1][0])
print("number of iterations until convergence:", iters)

# b. Plot the function and the points produced by the algorithm
fig = plt.figure()
ax = fig.gca(projection='3d')

PARAMS = np.array(list(values_params.values()))
PARAMS = PARAMS.reshape((PARAMS.shape[0], PARAMS.shape[2]))
X = PARAMS[:,0]
Y = PARAMS[:,1]
X, Y = np.meshgrid(X, Y)
Z = fxy(x = X, y = Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title("Newthon Algotithm")
plt.savefig("PS2_203670872_2b")

# c.
'''
While the steepest descent method took at least 911 iterations to converge
(possibly a lot more, see 1(b) for more details), it took only 4 iterations
until the Newthon method converged.
It is not suprising that the Newthon method converged faster, it has a quadratic
rate of convergence. The steepest descent has liniear rate of convergence and
when the problem is not of a quadratic form we have to iterate inside of every
step in order to numericly look for the steepest directeion of update.
Moreover, the Newthon has a simpler code and computationaly more efficient -
matrix multiplication could done parallelly, while the search for alpha has to
be done sequentially.
I started both algorithms in the same point, (x,y) = (0.8,0.5)
'''
