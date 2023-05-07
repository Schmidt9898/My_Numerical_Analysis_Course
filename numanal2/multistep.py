import numpy as np
import matplotlib.pyplot as plt
import math
import os

# I tried to write all function so that they can work with any number of dimension

#Simple itaration from last semester assigment
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
def Adam_Bash(functions,v0,T,N):
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

def find_peeks(array):
	array = list(array)
	pos_ = []
	#for n in range(N):
	#	max_value = max(array)
	#	max_index = array.index(max_value)
	#	for i in range(max(max_index-window,0),min(max_index+window,len(array))):
	#		array[i] = 0
	#	pos_.append(max_index)

	last_val = array[0] 
	is_down = False
	for i,val in enumerate(array):
		if is_down:
			if last_val < val:
				is_down = False
			
		else:
			if last_val > val:
				is_down = True 
				pos_.append(i - 1)
		last_val = val

	return pos_
	pass


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
	exact_sol = np.asarray(Adam_Bash(f,v,T,10**7))
	np.save("./Lotka-Volterra.npy",exact_sol)
	print("./Lotka-Volterra.npy has been saved for faster rerun :)")
plt_x_exact = np.linspace(0,T,exact_sol.shape[0])


peeks = find_peeks(exact_sol[:,0])
print(peeks)
last_peek_pos_t = plt_x_exact[peeks[-1]]
last_peek_exact_value = np.interp(last_peek_pos_t, plt_x_exact, exact_sol[:,0]) #we use the X function
print("last peek pos t:",last_peek_pos_t,"value",last_peek_exact_value)


#Solver = Euler
#Solver = reluE
#Solver = Cranky_Nicolson
#Solver = Adam_Bash

Solvers = [	Euler,
			reluE, 			 # backward euler and Crank Nicolson can be removed 
			Cranky_Nicolson, # It was easy to write so i keept them here
			Adam_Bash]

print("a,b,c,d",a,b,c,d)
print("x0,y0",x0,y0)
print("T,N",T,N)


for Solver in Solvers:
	print("Solver is:",Solver.__name__)

	val = Solver(f,v,T,N)
	val = np.asarray(val)
	plt_x = np.linspace(0,T,val.shape[0])

	estimate_peeks_x = find_peeks(val[:,0])
	estimate_peeks_y = find_peeks(val[:,1])


	estimate_vals_at_peeks_x = [val[p,0] for p in estimate_peeks_x]
	estimate_vals_at_peeks_y = [val[p,1] for p in estimate_peeks_y]

	print(Solver.__name__,"Estimate peeks:",estimate_peeks_x )
	print(Solver.__name__,"Estimate values:",estimate_vals_at_peeks_x )
	print(Solver.__name__,"Estimate peeks:",estimate_peeks_y )
	print(Solver.__name__,"Estimate values:",estimate_vals_at_peeks_y )


	print("Max population of X:",max(estimate_vals_at_peeks_x))
	print("Min population of X:",min(estimate_vals_at_peeks_x))
	print("Max population of Y:",max(estimate_vals_at_peeks_y))
	print("Min population of Y:",min(estimate_vals_at_peeks_y))


	plt.figure("Lotka-Volterra with "+Solver.__name__)
	plt.plot(plt_x_exact,exact_sol[:,0],label = "exact X")
	plt.plot(plt_x_exact,exact_sol[:,1],label = "exact Y")

	plt.plot(plt_x,val[:,0],label = "Estimated X")
	plt.plot(plt_x,val[:,1],label = "Estimated Y")
	plt.legend()

	plt.show(block = False)

	print("Calculating error rate 2^10 -> 2^18")

	error_rate = []

	for n in range(10,18):
		val = Solver(f,v,T,2**n)
		val = np.asarray(val)
		plt_x = np.linspace(0,T,val.shape[0])
		calculated_value = np.interp(last_peek_pos_t, plt_x, val[:,0])
		error = ( calculated_value - last_peek_exact_value ) / last_peek_exact_value 
		error_rate.append(error)

	plt.figure("Error rate of "+Solver.__name__)
	plt.plot(range(10,18),error_rate)

	print(Solver.__name__,"Done.")
	print("-"*100)



plt.show()

#print(val)