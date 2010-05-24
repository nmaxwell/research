
from numpy import *
from math import *
import itertools as itts
import random as rand

from numpy.linalg import svd
from numpy import sum,where

def matrixrank(A,tol=1e-8):
    s = svd(A,compute_uv=0)
    return sum( where( s>tol, 1, 0 ) )




nu = 4
N = 2**nu

m = 3

A = zeros((N,m))
B = zeros((m,N))

for i in range(N):
    for j in range(m):
        A[i][j] = rand.uniform(-1.0,1.0)
        B[j][i] = rand.uniform(-1.0,1.0)

f = zeros((N,1))

for i in range(N):
    f[i] = rand.uniform(-1.0,1.0)

K = dot(A,B)

u = dot(K,f)


print u

print matrixrank(K)












"""




def cart_prod2(X,Y ):
    prod=[]
    for x in X:
        for y in Y:
            prod.append((x,y))
    return prod


def FIO(f, X, P):
    u=[]
    for x in X:
        sum = 0.0
        for p in P:
            sum += K(x,p)*f(p)
        u.append(sum)
    
    return u






import waveprop as wp



if __name__ == "__main__":
    
    M = 64
    
    grid = wp.grid2d( n1=M, n2=M, a1=0.0, a2=0.0, b1=1.0, b2=1.0 )
    
    index = cart_prod2(range(M),range(M))
    
    X = grid.X1
    Y = grid.X2
    X = cart_prod2(X,Y)
    
    P = arange(0.5,-0.5, 1.0/4 )
    P = cart_prod2(P,P)
    
    K_function = lambda x1,x2: lambda p1,p2: exp( (-x1*p1-x2*p2)*2.0*pi )
    
    y1 = 1.0
    y2 = 1.0
    f_function = lambda x1,x2: exp( (x1*y1+x2*y2)*2.0*pi )
    
    K=[]
    for i in range(M):
        L = []
        for j in range(M):
            x1,x2 = grid(i,j)
            L.append( grid.evaluate(K_function(x1,x2)) )
        K.append(L)
    
    f = grid.evaluate(f_function)
    #f = grid.zeros()
    #f[2][2] = 1.0
    
    u = grid.zeros()
    for (i,j) in index:
        u += dot(K[i][j], f)
    
    print f 
    
    
    wp.write_png( f, "out.png", major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ), ordering='rm' )
    
    
"""

"""
    def K(x,y):
        m = dot(x,y)*2.0*pi
        return cos(m) + 1.0j*sin(m)
    
    P = arange(-1.0,1.0, 1.0/4 )
    P = cart_prod2(P,P)
    
    X = arange(0.0,1.0, 1.0/4 )
    X = cart_prod2(X,X)
    
    
    
    f = lambda (x,y): exp(-x*x-y*y)
    
    u = FIO(f, X, P)
    
    for x in u:
        print x
"""











