
import numpy as np
import time

# quadrature rules 
#		M,p,quadrature weights
rule0 = (2,2,[1/2,1/2])
rule1 = (3,4,[1/6,2/3,1/6])
rule2 = (4,4,[1/8,3/8,3/8,1/8])
rule3 = (5,6,[7/90,16/45,4/30,16/45,7/90])


def equidistant(a,b,f,N=1000):
	print("Equidistant:")
	print(f"	a{a}, b{b}, N{N}")

	X = list(np.linspace(a, b, N, endpoint=True))
	Y = [f(x) for x in X]

	w = float(b-a) / N

	sum_=0
	for i,y in enumerate(Y):
		sum_+=y*w
	#print(sum_)

	
	return sum_


def CalculateQ(a,b,f,rule):
	Q = 0
	M,p,W = rule

	global called_count
	called_count+=M

	w = abs(b-a) # so this needed 
	x = list(np.linspace(a, b, M, endpoint=True))
	for i in range(M):
		Q += W[i] * f(x[i]) * w
	return Q

called_count = 0 # number of point operation
def adaptive(a,b,f,rule,e,Q = None):

	M,p,W = rule
	if Q is None: #so we do not recalculate
		Q = CalculateQ(a,b,f,rule)
	Q1 = CalculateQ(a,(a+b)/2,f,rule)
	Q2 = CalculateQ((a+b)/2,b,f,rule)


	error = abs(Q1 + Q2 - Q) / (pow(2,p)-1) 

	if error < e:
		#accept
		return Q1+Q2
	else:
		return adaptive(a,(a+b)/2,f,rule,e/2,Q1) + adaptive((a+b)/2,b,f,rule,e/2,Q2)


f = lambda x : 1/(0.01+pow(x-0.3,2)) + 1/(0.04 + pow(x-0.9,2)) - 6
#f = lambda x : 1

print("Helo")


a,b = 0,1

e = 0.0001

true_value = 29.858325395498671

start = time.time()
val = equidistant(a,b,f,10000)
duration = time.time() - start
print("equid",val - true_value,"duration",duration )

start = time.time()
val = adaptive(a,b,f,rule1,e)
duration = time.time() - start
print("adapt",val - true_value,"duration",duration )

print(called_count)







