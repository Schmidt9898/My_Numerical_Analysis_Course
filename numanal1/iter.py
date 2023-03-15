import numpy as np
import math

print("Iteration...")


print("Exercise 1 - look at row 9 to 59")

# return a tuple with (isgood,root,iteration count)
def Newton(x0,f,fder,precision = 4,max_iter = 1000):
	precision = pow(10,-precision)
	x = x0
	try:
		for i in range(max_iter):
			x = x - f(x)/fder(x)
			#print(x)
			if(abs(f(x))<precision):
				return (True,x,i)
		return (x,max_iter)
	except OverflowError:
		print("OverflowError accured.")
		return (False,-1,i)
		

def Secant(x0,x1,f,precision = 4,max_iter = 1000):
	precision = pow(10,-precision)
	x = x1
	xmin1 = x0
	try:
		for i in range(max_iter):
			x_back = x
			x = x - f(x)*((x-xmin1)/(f(x)-f(xmin1)))
			xmin1 = x_back
			if(abs(f(x))<precision):
				return (True,x,i)
	except Exception as e:
		print(e)
		return (False,0,max_iter)	
	print("Max iteration count is reached")
	return (False,x,max_iter)	



# This function can raise exception if parameters are wrong 
def Bisec(a,b,f,precision = 4,max_iter = 1000):
	precision = pow(10,-precision)
	if b-a < 0:
		#raise Exception(" error: a > b ",a,b)
		print("Error: a > b ",a,b)
		return (False,a,0)

	try:

		for i in range(max_iter):
			x = a + (b-a)/2

			if(abs(f(x))<precision):
				return (True,x,i)

			if f(a) * f(x) < 0:
				#left part go
				b = x
				pass
			elif f(x)*f(b) < 0:
				#right part go
				a = x
				pass
			else:
				#raise Exception("Mean value theorem is not met")
				print("Mean value theorem is not met")
				return (False,x,i)
	except Exception:
		return (False,0,i)

	return (False,x,max_iter)	


print("Exercise 2:")

f = lambda x : pow(x,3)+2*pow(x,2)+10*x-20
fder = lambda x : 3*pow(x,2)+4*x+10

isgood ,x , iteration = Newton(1,f,fder,4)
print("	Solution Newton: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Newton")
isgood ,x , iteration = Secant(1,2,f,4)
print("	Solution Secant: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Secant")
isgood ,x , iteration = Bisec(1,2,f,4)
print("	Solution  Bisec: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Bisec")


f = lambda x : math.tanh(x)
fder = lambda x : 1/( math.pow(math.cosh(x),2)) 


print("Exercise 3:")
isgood ,x , iteration = Newton(-5,f,fder,4)
print("	Solution Newton: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Newton")
isgood ,x , iteration = Secant(-5,-4,f,4)
print("	Solution Secant: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Secant")
isgood ,x , iteration = Bisec(5,10,f,4)
print("	Solution  Bisec: x={} f(x)={} iteration count={}".format(x,f(x),iteration)) if isgood else print("	Cannot be solved with Bisec")


print("Exercise 4:")
#calculation of L is in the pdf 

g = lambda x :  20/(pow(x,2)+2*x+10)

L = 80/180

x0 = 1
x1 = 1.5384615384615385

eps = pow(10,-4)

k = (np.log(abs(x0-x1))-np.log(eps*(1-L)))/(np.log(1/L)) + 1

print("K:",k)

#lets do a ceil on this k
k = math.ceil(k)


x = x1
x_last = x1
max_iter = 100
iteration_took = 0
max_iter_is_reached = True
for i in range(1,max_iter):
	x = g(x)

	if abs(x - x_last) < eps:
		print("iteration =",i, "x =",x)
		max_iter_is_reached = False
		iteration_took = i
		break
	x_last = x
	print(x)

if max_iter_is_reached:
	print("Max itaration is reached.")
else:
	print("k =",k,">= iteration =",iteration_took)


print("End")










