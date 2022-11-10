import numpy as np
import math
import matplotlib.pyplot as plt


def Linear(x,X,Y):
	xi_1 = -1
	xi = -1
	Yi = 0
	Yi_1 = 0
	for i in range(len(X)):
		if X[i]>=x:
			xi=X[i]
			xi_1=X[i-1]
			Yi = Y[i]
			Yi_1 = Y[i-1]
			break
	sl = (xi-x) / (xi-xi_1) * Yi_1 + (x-xi_1) / (xi-xi_1) * Yi  
	return sl


#Lord forgive me for what i am about to code
def Lk(k,x,X):
	value = 1
	for i in range(len(X)):
		if i != k:
			value *= ((x-X[i]) / (X[k]-X[i]))  
	return value
def Lagrange(x,X,Y):
	value = 0
	for i in range(len(X)):
		value += Lk(i,x,X)*Y[i]
	return value
	#print(L0(0))
	pass






print("assigment3")	


a=-5
b=5

f=lambda x : 1/(1+pow(x,2))


N=5


X = list(np.linspace(-5, 5, N, endpoint=True))
Y =[f(x) for x in X]
#print(Y)


N=2048
plot_x = list(np.linspace(-5, 5, N, endpoint=True))
plot_y =[f(x) for x in plot_x]
plt.plot(plot_x,plot_y,label='F')
#plt.title("Theta methods")
plt.ylabel("y")
plt.xlabel("t")
#plt.legend()

plot_y =[Linear(x,X,Y)  for x in plot_x]
plt.plot(plot_x,plot_y,label='Linear')

plot_y =[Lagrange(x,X,Y)  for x in plot_x]
plt.plot(plot_x,plot_y,label='Lagrange')



def inf_error(f,pn,samples,a,b,N = 2048):
	
	X = list(np.linspace(-5, 5, samples, endpoint=True))
	Y =[f(x) for x in X]


	test_x = list(np.linspace(-5, 5, N, endpoint=True))
	error = 0
	for x in test_x:
		error = max(error, abs( f(x) - pn(x,X,Y) ))
		#print(error)
	return error


for n in range(7):
	samples = pow(2,n)
	print(samples)
	#error = inf_error(f,Lagrange,samples,a,b)
	error = inf_error(f,Linear,samples,a,b)
	print(error)
	

















plt.show()
