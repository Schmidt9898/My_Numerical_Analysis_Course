
import math
import numpy as np
import time

# quadrature rules 
#		M,p,quadrature weights
rule0 = (2,2,[1/2,1/2])
rule1 = (3,4,[1/6,2/3,1/6])
rule2 = (4,4,[1/8,3/8,3/8,1/8])
rule3 = (5,6,[7/90,16/45,4/30,16/45,7/90])

ruleGaussian = (2,2,[-1/math.sqrt(3),1/math.sqrt(3)])


# N will be multiplied with the M number of the rule
def equidistant(a,b,f,rule,N=1000):
	print("Equidistant:")
	print(f"	a{a}, b{b}, Intervals{N} ")

	intervals = list(np.linspace(a, b, N, endpoint=True))

	sum_=0
	for i in range(len(intervals)-1):
		sum_ += CalculateQ(intervals[i],intervals[i+1],f,rule)
	#print(sum_)
	return sum_


operation_count = 0 # number of point operation
interval_count = 0 # number of point operation

def CalculateQ(a,b,f,rule):
	Q = 0
	M,p,W = rule

	global operation_count
	global interval_count
	operation_count+=M
	interval_count+=1

	w = abs(b-a) # so this needed for good scale
	x = list(np.linspace(a, b, M, endpoint=True))
	for i in range(M):
		Q += W[i] * f(x[i]) * w
	return Q


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


def Gaussian_calc(a,b,f):

	global operation_count
	global interval_count
	operation_count+=2
	interval_count+=1

	c = (b+a)/2.0
	length = b-a
	w = 1.0 * abs(b-a)/2 # universal for 2 point

	x_p = [c-length/(2*math.sqrt(3)),c+length/(2*math.sqrt(3))]

	return f(x_p[0])*w  + f(x_p[1])*w 



def Gaussian_inegrate(a,b,f,N = 100):

	print("Gaussian:")
	print(f"	a{a}, b{b}, Intervals{N} ")

	intervals = list(np.linspace(a, b, N, endpoint=True))

	sum_=0
	for i in range(len(intervals)-1):
		sum_ += Gaussian_calc(intervals[i],intervals[i+1],f)
	#print(sum_)
	return sum_





f = lambda x : 1/(0.01+pow(x-0.3,2)) + 1/(0.04 + pow(x-0.9,2)) - 6
#f = lambda x : 1

print("Helo")


a,b = 0,1
rule = rule0
e = 0.0001

true_value = 29.858325395498671




start = time.time()
val = equidistant(a,b,f,rule,1000)
duration = time.time() - start
print(val," equid",val - true_value,"duration",duration,"intervals:",interval_count,"operation:",operation_count)

operation_count = 0 
interval_count = 0 

start = time.time()
val = adaptive(a,b,f,rule,e)
duration = time.time() - start
print(val," adapt",val - true_value,"duration",duration,"intervals:",interval_count,"operation:",operation_count)

operation_count = 0 
interval_count = 0 

start = time.time()
val = Gaussian_inegrate(a,b,f,1000)
duration = time.time() - start
print(val," Gaussian",val - true_value,"duration",duration,"intervals:",interval_count,"operation:",operation_count)


interval_count = 0 
operation_count = 0 







