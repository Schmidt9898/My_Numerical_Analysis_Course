from functools import partial
import time
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def Chebyshev_space(a,b,N,endpoint=True): #in this we get the n-2 degree chebyshev and add the endpoints
	if endpoint:
		N=N-2
	multipliyer = b # in this script i will not consider diferent than [-a,a] form
	x = []
	if N>0:
		for j in reversed(range(1,N+1)):
			x.append(multipliyer * math.cos(( (2*j-1)*math.pi ) / ( 2*N )))
	if endpoint:
		x = [a] + x + [b]
	return x

class Polinome:
	def __init__(self,X,Y):
		self.X = X
		self.Y = Y
	def __call__(self, x):
		raise Exception("Base class does not have implementation for y(x)")


class Linear(Polinome):
	def __init__(self,X,Y):
		self.X = X
		self.Y = Y
	def __call__(self, x):
		xi_1 = -1
		xi = -1
		Yi = 0
		Yi_1 = 0
		for i in range(len(self.X)):
			if self.X[i]>=x:
				xi=self.X[i]
				xi_1=self.X[i-1]
				Yi = self.Y[i]
				Yi_1 = self.Y[i-1]
				break
		sl = (xi-x) / (xi-xi_1) * Yi_1 + (x-xi_1) / (xi-xi_1) * Yi  
		return sl


class Hermite_cubic_spline(Polinome):
	def __init__(self,hi,X,Y,Z):
		self.X = X
		self.Y = Y
		self.Z = Z


		self.H0 = lambda x,xi,xi_1 : (pow((x-xi),2)/pow(hi,2))*(1+(2/hi)*(x-xi_1))
		self.H1 = lambda x,xi,xi_1 : (pow((x-xi_1),2)/pow(hi,2))*(1-(2/hi)*(x-xi))
		self.K0 = lambda x,xi,xi_1 : (pow((x-xi),2)/pow(hi,2))*(x-xi_1)
		self.K1 = lambda x,xi,xi_1 : (pow((x-xi_1),2)/pow(hi,2))*(x-xi)
	def __call__(self, x):
		xi_1 = -1
		xi = -1
		Yi = 0
		Yi_1 = 0
		idx = -1
		for i in range(len(self.X)):
			if self.X[i]>=x:
				idx = i
				xi=self.X[i]
				xi_1=self.X[i-1]
				Yi = self.Y[i]
				Yi_1 = self.Y[i-1]
				Zi = self.Z[i]
				Zi_1 = self.Z[i-1]
				break
		val =  self.H0(x,xi,xi_1)*Yi_1 + self.H1(x,xi,xi_1)*Yi + self.K0(x,xi,xi_1)*Zi_1 + self.K1(x,xi,xi_1)*Zi
		return val

#Lord forgive me for what i am about to code
class Lagrange(Polinome):
	def __init__(self,X,Y):
		self.X = X
		self.Y = Y

	def Lk(self,k,x):
		value = 1
		for i in range(len(self.X)):
			if i != k:
				value *= ((x-self.X[i]) / (self.X[k]-self.X[i]))  
		return value
	def __call__(self, x):
		value = 0
		for i in range(len(self.X)):
			value += self.Lk(i,x)*self.Y[i]
		return value


def do_exercise1():
	print("solving a.")
	a,b=(-5,5) 					#interval
	f=lambda x : 1/(1+pow(x,2)) # exercise function
	N = 4 						# Nth digree polinome


								#n+1 point
	X = list(np.linspace(-5, 5, N+1, endpoint=True)) #samples
	Y =[f(x) for x in X]
	print("Samples:",Y)


	plt.figure("Exercise 1 / a") # raise figure for exercise 1
	plt.subplot(121)
	plt.title("Lagrange interpolation")

	N=1025
	plot_x = list(np.linspace(-5, 5, N, endpoint=True))
	plot_y =[f(x) for x in plot_x]
	plt.plot(plot_x,plot_y,label='f(x)')
	plt.ylabel("y")
	plt.xlabel("x")

	lagrange = Lagrange(X,Y)

	plot_y =[lagrange(x) for x in plot_x]
	plt.plot(plot_x,plot_y,label='Lagrange')

	plt.legend()
	plt.show(block = False)


	errors = []
	samples = []
	#2,4,6,....,24
	for n in range(2,26,2):
		X = list(np.linspace(-5, 5, n+1, endpoint=True))
		Y =[f(x) for x in X]
		pn = Lagrange(X,Y)

		test_x = list(np.linspace(-5, 5, 100, endpoint=True))
		error = 0
		for x in test_x:
			error = max(error, abs( f(x) - pn(x) ))
		errors.append(error)
		samples.append(n)

	print("errors of Lagrange with equidistant points: \n",errors)

	plt.subplot(122)
	plt.title("MAX error")
	plt.plot(samples,errors,label="inf error")
	plt.ylabel("error")
	plt.xlabel("degree")
	plt.show(block = False)

##################################################################################################################
##################################################################################################################
##################################################################################################################
	print("solving b.")

									#n+1 point
	X = list(Chebyshev_space(-5, 5, 10+1,False)) #samples
	Y =[f(x) for x in X]
	print("Samples:",Y)


	plt.figure("Exercise 1 / b") # raise figure for exercise 1
	plt.subplot(121)
	plt.title("Lagrange interpolation")

	N=1025
	plot_x = list(np.linspace(-5, 5, N, endpoint=True))
	plot_y =[f(x) for x in plot_x]
	plt.plot(plot_x,plot_y,label='f(x)')
	plt.ylabel("y")
	plt.xlabel("x")


	lagrange = Lagrange(X,Y)

	plot_y =[lagrange(x) for x in plot_x]
	plt.plot(plot_x,plot_y,label='Lagrange')

	plt.legend()
	plt.show(block = False)


	errors = []
	samples = []
	#2,4,6,....,24
	for n in range(2,26,2):
		X = list(Chebyshev_space(-5, 5, n+1,endpoint=False)) #samples
		#X = list(np.linspace(-5, 5, n+1, endpoint=True))
		Y =[f(x) for x in X]
		pn = Lagrange(X,Y)

		test_x = list(np.linspace(-5, 5, 100, endpoint=True))
		error = 0
		for x in test_x:
			error = max(error, abs( f(x) - pn(x) ))
		errors.append(error)
		samples.append(n)

	print("errors of Lagrange with equidistant points: \n",errors)

	plt.subplot(122)
	plt.title("MAX error")
	plt.plot(samples,errors,label="inf error")
	plt.ylabel("error")
	plt.xlabel("degree")
	plt.show(block = False)

	print("We can observe that the for linear mapping the Runge's phenomenon appears, but for the nonlinear mapping it does not, so we must choose our interpolation point cleverly.")

def plot(ax,f,a=0,b=1):
	x_data = list(np.linspace(a, b, 1000, endpoint=True))
	y_data =[f(x) for x in x_data]
	#ax.clear()
	ax.plot(x_data,y_data)
	ax.set_xlim(0,1)
	ax.set_ylim(-5, 5)


def do_exercise2():

	#plt.figure("Animation exercise 2")

	f_t = lambda x,t : math.sin(5*math.pi*x)*math.cos(10*math.pi*t) + 2*math.sin(7*math.pi*x)*math.cos(14*math.pi*t)
	fder_t = lambda x,t : math.pi*( 5*math.cos(10*math.pi*t)*math.cos(5*math.pi*x) + 14*math.cos(14*math.pi*t)*math.cos(7*math.pi*x) )

	fig, ax = plt.subplots()

	for i in range(6):
		for t in np.linspace(0,1,1000,endpoint=True):
			ax.clear()
			X = np.linspace(0,1,51,endpoint=True)
			Y =[f_t(x,t) for x in X]
			Z =[fder_t(x,t) for x in X]
			
			f = Linear(X,Y)
			plot(ax,f)

			#f = Hermite_cubic_spline(1/50,X,Y,Z)
			#plot(ax,f)

			f = lambda x : f_t(x,t)
			plot(ax,f)

			ax.set_title(f"t = {t}")
			# Note that using time.sleep does *not* work here!
			plt.pause(0.05)

	plt.show()

	#ani = animation.FuncAnimation(
	#	fig,animate,init_func=init_anim,fargs=[f_t],frames = np.linspace(0, 1, 100,endpoint=True),blit = True)

	#plot_y =[Linear(x,X,Y)  for x in plot_x]
	#plt.plot(plot_x,plot_y,label='Linear')




print("assigment3")	

do_exercise1()
do_exercise2()

plt.show()

















