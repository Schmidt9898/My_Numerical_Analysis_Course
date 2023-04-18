
from enum import Enum
import math
import numpy as np
import matplotlib.pyplot as plt

class BoundaryCondition(Enum):
	Dirichlet = 0
	Neumann = 1


def Errornorm( u_e , u_h ):
	u_h = np.asarray(u_h)
	# ERRORNORM This function computes L2 error
	# Computes an estimated distance between the
	# exact solution and the numerical solution .
	# The integral is approximated with three - point
	# Gauss - Legendre quadrature rule
	# u_e : Exact solution ( function )
	# u_h : Numerical solution ( vector of nodal values )
	n = len ( u_h ) # number of subintervals
	h = 1.0 / ( n - 1) # length of subintervals
	# rule for reference interval [0 , 1]
	t = np.asarray([ (1 - math.sqrt (3/5)) / 2 , 1 / 2 , (1 + math.sqrt (3/5)) / 2 ])
	w = np.asarray([ 5/18 , 8/18 , 5/18 ])
	e2 = 0
	for i  in range(0,n-1):
		# evaluate at quadrature points
		points = (i) * h + h * t
		v_e = [u_e(x) for x in points]
		v_h = (1 - t ) * u_h[i] + t * u_h[i+1]
		# square difference and sum with weights
		e2 += sum(np.power(v_e - v_h,2) * w)
	return  math.sqrt ( h * abs( e2 ) )

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
	#print(h)
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


N = 5
X = np.linspace(0,1,N)
h = 1 / (N-1)

# -(cu')' = f

print("Exercise 1")

U_2_der = lambda x : 6*x # -u'' = 6x 
U = lambda x : - math.pow(x,3) + 3*x # u = -x^3 + 3x

f_function = U_2_der 
c_function = lambda x : 1 


solution = [U(x) for x in X]
f = [f_function(x) for x in X]
c = [c_function(x) for x in X]


print("N:",N)
print("c:",c)
print("f:",f)
print("Exact solution:", np.asarray(solution))

bound_a = (0,BoundaryCondition.Dirichlet)
bound_b = (0,BoundaryCondition.Neumann)

A_matrix = StiffnessMatrix(N, c)
b_vector = LoadVector(N, f)

A_matrix, b_vector = ApplyBoundary(A_matrix,b_vector,bound_a,bound_b,h)

print("A_matrix:",A_matrix)
print("b_vector:",b_vector)


ans = np.linalg.solve(A_matrix, b_vector)
#print(np.allclose(np.dot(A_matrix, ans), b_vector))


print("u_h:",ans)
print("Errornorm:",Errornorm(U,ans))



plt.figure("Exercise 1")
plt.plot(X,ans)
plt.plot(X,solution)

print("-----------------------------")

Ns = [4,8,16,32,64]

# building table:
last_diff = None
r = None
for n in Ns:
	n+=1 # because
	XX = np.linspace(0,1,n)
	h = 1 / (n-1)

	solution = [U(x) for x in X]
	f = [f_function(x) for x in XX]
	c = [c_function(x) for x in XX]

	A_matrix = StiffnessMatrix(n, c)
	b_vector = LoadVector(n, f)
	A_matrix, b_vector = ApplyBoundary(A_matrix,b_vector,bound_a,bound_b,h)
	ans = np.linalg.solve(A_matrix, b_vector)
	error = Errornorm(U,ans)
	if last_diff is None:
		last_diff = error
	else:
		r = math.log2(error/last_diff)*-1 ## hehe
		last_diff = error
	print("n-1:",n-1,"Errornorm:",error,"Rate:",r)
print("-----------------------------")


# Exercise 2
print("Exercise 2")


h = 1 / (N-1)
f_function = lambda x : 4*x
c_function = lambda x : 1+x 

U = lambda x: 2*x - math.pow(x,2) + math.log(1+x) - math.log(2)

solution = [U(x) for x in X]
f = [f_function(x) for x in X]
c = [c_function(x) for x in X]




print("N:",N)
print("c:",c)
print("f:",f)
print("Exact solution:", np.asarray(solution))

bound_a = (3,BoundaryCondition.Neumann)
bound_b = (1,BoundaryCondition.Dirichlet)

A_matrix = StiffnessMatrix(N, c)
b_vector = LoadVector(N, f)

A_matrix, b_vector = ApplyBoundary(A_matrix,b_vector,bound_a,bound_b,h)

print("A_matrix:",A_matrix)
print("b_vector:",b_vector)


ans = np.linalg.solve(A_matrix, b_vector)
#print(np.allclose(np.dot(A_matrix, ans), b_vector))

print("u_h:",ans)
print("Errornorm:",Errornorm(U,ans))

plt.figure("Exercise 2")
plt.plot(X,ans)
plt.plot(X,solution)

plt.show()