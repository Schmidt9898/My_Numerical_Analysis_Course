
from enum import Enum
import math
import numpy as np
import matplotlib.pyplot as plt

class BoundaryCondition(Enum):
	Dirichlet = 0
	Neumann = 1


def Errornorm( u_e , u_h ):
	pass

def trapazoid(a,b,h):
	return (a+b) / (2*h) 

def LoadVector(N, f):
	h = 1 / (N-1)
	b = np.zeros(N)

	b[0] = f[0]/2.0*h
	for i in range(1,N-1):
		b[i] = f[i]*h
	n=N-1
	b[n] = f[n]/2.0*h
	return b



def StiffnessMatrix(N, c):
	h = 1 / (N-1)
	print(h)
	A = np.zeros((N,N))

	A[0,0] = trapazoid(c[0],c[1],h)
	A[0,1] = -trapazoid(c[0],c[1],h)

	A[1,0] = -trapazoid(c[0],c[1],h)
	for i in range(1,N-1):
		A[i,i] = trapazoid(c[i-1],c[i],h) + trapazoid(c[i],c[i+1],h)
		A[i+1,i] = -trapazoid(c[i],c[i+1],h)
		A[i,i+1] = -trapazoid(c[i],c[i+1],h)

	n=N-1
	A[n,n] = trapazoid(c[n-1],c[n],h)
	return A

def ApplyBoundary(A,b,b_left,b_right,h):
	

	#applying left side
	val,type=b_left
	if type == BoundaryCondition.Neumann:
		b[0] -= val
	else: #BoundaryCondition.Dirichlet
		b[0] = val
		A[0,0] = 1
		A[0,1] = 0

	#applying right side
	val,type=b_right
	if type == BoundaryCondition.Neumann:
		b[-1] += val
	else: #BoundaryCondition.Dirichlet
		b[-1] = val
		A[-1,-1] = 1
		A[-1,-2] = 0
	
	
	return A,b



U_2_der = lambda x : -6*x # -u'' = 6x 
U = lambda x : - math.pow(x,3) + 3*x # u = -x^3 + 3x hand calculated

input_function = U_2_der 
#input_function = lambda x : 1 



N = 100
X = np.linspace(0,1,N)
h = 1 / (N-1)

solution = [U(x) for x in X]
c = [input_function(x) for x in X]


print(N)
print(c)
print(solution)

bound_a = (0,BoundaryCondition.Dirichlet)
bound_b = (0,BoundaryCondition.Neumann)

A_matrix = StiffnessMatrix(N, c)
b_vector = LoadVector(N, c)

A_matrix, b_vector = ApplyBoundary(A_matrix,b_vector,bound_a,bound_b,h)

print(A_matrix)
print(b_vector)


ans = np.linalg.solve(A_matrix, b_vector)
print(np.allclose(np.dot(A_matrix, ans), b_vector))

plt.figure("ans")
plt.plot(X,ans)
plt.plot(X,solution)

plt.show()

#a_val,type = bound_a
#if type == "dir":
#	h*c[-1]/2 + bound_a(0) 