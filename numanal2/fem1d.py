print("FEM 1D")

from enum import Enum
import math
import numpy as np
import matplotlib.pyplot as plt


#Simple itaration from last semester assigment updated
def simple_iteration(g,v,err=4,MAX_ITER=1000):
	if err > 1:
		err = pow(10,-err)
	
	for i in range(1,MAX_ITER):
		new_v = np.asarray( g(v) )

		if np.linalg.norm(v - new_v) < err:
			#print("iteration =",i, "x =",x)
			return new_v,i
			break
		v = new_v
	print(v )
	raise(Exception("Out of iteration!"))




class BoundaryCondition(Enum):
	Dirichlet = 0
	Neumann = 1

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



def FEM_HeatEq_1D(Theta,M,N,c_f,f_f,ux0,a_,b_,t0,tm,bound_a,bound_b):
	n = N+1
	m = M+1
	h = (b_ - a_) / n
	k = (tm - t0) / m
	X = np.linspace(a_,b_,n,True)
	T = np.linspace(t0,tm,m,True)



	print("M: {}, N: {}, Theta: {}".format(M,N,Theta))
	print("X",X)
	print("k",k)
	#T0 = np.asarray([ux0(x) for x in X]) 
	#print("T0",T0)


	Mass = MassMatrix(n,h)
	#print(Mass)
	y = ApplyIC(ux0,X)


	U = [y]
	t=t0
	f = [f_f(x,t) for x in X]
	c = [c_f(x,t) for x in X]
	A = StiffnessMatrix(n, c)
	b = LoadVector(n, f)
	first = True
	for t in T:
		if first:
			first = False
			continue
		f = [f_f(x,t) for x in X]
		c = [c_f(x,t) for x in X]
		Atn = StiffnessMatrix(n, c)
		btn = LoadVector(n, f)

		y = np.transpose(np.asarray(y))

		#A_hat = Mass + Theta * k * A
		#M_term = Mass - (1 - Theta) * k * A
		#b_hat = np.matmul(M_term , y + k * b)
		A_hat = Mass + (1 - Theta)*Atn
		b_hat = np.matmul((Mass - (Theta)*k*A),y) + k*((Theta)*b+(1 - Theta)*btn)


		A_hat, b_hat = ApplyBoundary(A_hat,b_hat,bound_a,bound_b,h)



		y = np.linalg.solve(A_hat, b_hat)
		#print(y)
		
		U.append(y)
		#print(y)

		A=Atn
		b=btn

	U = np.asarray(U)
	print(U)
	print(U.shape)
	return U
	pass	

def MassMatrix(N,h):
	n = N + 1
	M = np.zeros((N,N)) 
	for i in range(1,N-1): # 
		M[i,i] = 1
	M[0,0] = 0.5
	M[N-1,N-1] = 0.5
	return M * h


def ApplyIC(g,X):
	y = [ g(x) for x in X ] 
	#U = [] # n x m

	return y


		#Theta,M,N
cases =[
		(1,10,10),
		(1,10,200),
		(1,200,10),
		(0.5,10,10),
		(0.5,10,200),
		(0.5,200,10),
		(0,10,10),
		(0,10,200),
		(0,200,10),]


fig = plt.figure("Exercise 1")

idx = 0
for case in cases:
	idx+=1
	

	Theta,M,N = case#


	a,b = 0,1
	t0,tm = 0,1


	bound_a = (0,BoundaryCondition.Dirichlet) #TODO
	bound_b = (0,BoundaryCondition.Dirichlet) #TODO

	ux0 = lambda x : math.sin(math.pi*x)

	c = lambda x,t : 1 # 0 function # TODO
	f = lambda x,t : 0 #3*x*math.exp(-3*t) # 0 function # TODO





	U = FEM_HeatEq_1D(Theta,M,N,c,f,ux0,a,b,t0,tm,bound_a,bound_b)

	n = N+1
	m = M+1

	U_exact = np.zeros((m,n))
	ue = lambda x,t : math.exp(-math.pi*math.pi*t) * math.sin(math.pi*x)

	ti = t0
	for t in range(m):
		xi=a
		for x in range(n):
			U_exact[t,x] = ue(xi,ti)
			xi+= 1/N
		ti+= 1/M
	print(U_exact)

	T = np.linspace(t0,tm,m,True)
	X = np.linspace(a,b,n,True)



	X, Y = np.meshgrid(X,T)

	# Create a figure and axis
	ax = fig.add_subplot(330+idx, projection='3d')
	ax.set_title("Theta_{}_M_{}_N_{}".format(Theta,M,N))

	# Plot the mesh using plot_trisurf
	#axs[int(idx/3),idx%3]
	ax.plot_surface(X,Y,U,label = "estimate")

	#plt.legend()
	
fig = plt.figure("Exact solution.")
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y,U_exact,label = "Exact")


plt.show()

quit()



