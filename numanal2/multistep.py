import numpy as np
import matplotlib.pyplot as plt
import math
import os

#Simple itaration from last assigment
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


def Euler(functions,v0,T,N):
	h = T/N
	v = np.asarray(v0)
	F = [v]
	for i in range(N):
		f_vec = np.asarray([f(v) for f in functions]) 
		v_new = v + h * f_vec
		v = v_new
		F.append(v)
	return F

#backward Euler
def reluE(functions,v0,T,N):
	h = T/N
	v = np.asarray(v0)
	F = [v]
	for i in range(N):
		#f_vec = np.asarray([f(v) for f in functions]) 
		#v_new = v + h * f_vec
		g = lambda vt : v + h * np.asarray([f(vt) for f in functions])
		#g = lambda x : ylast + h*yt(x)
		#do it, nested iteration
		v_new,_ = simple_iteration(g,v)
		#break
		v = v_new
		F.append(v)
	return F

def Cranky_Nicolson(functions,v0,T,N):
	h = T/N
	v = np.asarray(v0)
	F = [v]
	for i in range(N):
		hf = 0.5 * h * np.asarray([f(v) for f in functions]) 
		g = lambda vt : v + hf + 0.5 * h * np.asarray([f(vt) for f in functions])
		#do it, nested iteration
		v_new,_ = simple_iteration(g,v)
		v = v_new
		F.append(v)
	return F

# Adams-Bashford method
def Adam_Smash(functions,v0,T,N):
	h = T/N
	v = np.asarray(v0)
	hf = 0.5 * h * np.asarray([f(v) for f in functions]) 
	g = lambda vt : v + hf + 0.5 * h * np.asarray([f(vt) for f in functions])
	v1,_ = simple_iteration(g,v)
	F = [v,v1]
	for i in range(1,N):
		f0 = np.asarray([f(v0) for f in functions])
		f1 = np.asarray([f(v1) for f in functions])
		v2 = v1 + h * ( (3/2) * f1 - (1/2) * f0 )
		v0 = v1
		v1 = v2
		F.append(v2)
	return F

def find_peeks(array,N,window = 50):
	array = list(array)
	pos_ = []
	for n in range(N):
		max_value = max(array)
		max_index = array.index(max_value)
		for i in range(max(max_index-window,0),min(max_index+window,len(array))):
			array[i] = 0
		pos_.append(max_index)

	return pos_
	pass



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


print("multistep.py")

a,b,c,d = (1,1,1,1)
x0 = 5.0
y0 = 1.0
v = [x0,y0]
T = 50
N = 2**10
f = [
    lambda v : a * v[0] - b * v[0] * v[1] ,
    lambda v : c * v[0] * v[1] - d * v[1] 
]

exact_sol = None

if os.path.exists("./Lotka-Volterra.npy"):
	print("Solution found. Loading..")
	exact_sol = np.load("./Lotka-Volterra.npy")
else:
	print("Give it time, it is calculating...")
	exact_sol = np.asarray(Adam_Smash(f,v,T,10**7))
	np.save("./Lotka-Volterra.npy",exact_sol)
	print("./Lotka-Volterra.npy hs been saved forfaster rerun :)")
plt_x_exact = np.linspace(0,50,exact_sol.shape[0])

print(find_peeks(exact_sol[:,0],6))
#quit()


#solution_x = lambda x : np.interp(x, plt_x_exact, exact_sol[:,0])
#solution_y = lambda x : np.interp(x, plt_x_exact, exact_sol[:,1])

val = Euler(f,v,T,N)
#val = reluE(f,v,T,N)
#val = Cranky_Nicolson(f,v,T,N)
#val = Adam_Smash(f,v,T,N)
val = np.asarray(val)
plt_x = np.linspace(0,50,val.shape[0])

#print([solution_x(x) for x in plt_x])
#print([solution_y(x) for x in plt_x])



#print("X Error:",Errornorm( solution_x , val[:,0]))
#print("Y Error:",Errornorm( solution_y , val[:,1]))
#print(exact_sol.shape)

plt.figure("Lotka-Volterra")
plt.plot(plt_x_exact,exact_sol[:,0])
plt.plot(plt_x_exact,exact_sol[:,1])

plt.plot(plt_x,val[:,0])
plt.plot(plt_x,val[:,1])

#val2 = np.asarray(Euler(f,v,T,N))


plt.show()

#print(val)