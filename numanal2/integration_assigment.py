
import math
import numpy as np
import time

# quadrature rules 
#		M,p,quadrature weights
rule0 = (2,2,[1/2,1/2])
rule1 = (3,4,[1/6,2/3,1/6])
rule2 = (4,4,[1/8,3/8,3/8,1/8])
rule3 = (5,6,[7/90,16/45,4/30,16/45,7/90])

def Error_rate(estimate,actual):
	return (estimate - actual) / actual * 100


# N will be multiplied with the M number of the rule
def equidistant(a,b,f,rule,Calculator,N=1000):
	#print("Equidistant:")
	#print(f"	a{a}, b{b}, Intervals{N} ")

	intervals = list(np.linspace(a, b, N, endpoint=True))

	sum_=0
	for i in range(len(intervals)-1):
		sum_ += Calculator(intervals[i],intervals[i+1],f,rule)
	#print(sum_)
	return sum_


operation_count = 0 # number of point operation
interval_count = 0 # number of point operation

def clear_count():
	global operation_count
	global interval_count
	operation_count=0
	interval_count=0

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


def adaptive(a,b,f,rule,e,Calculator,Q = None):

	M,p,W = rule
	if Q is None: #so we do not recalculate
		Q = Calculator(a,b,f,rule)
	Q1 = Calculator(a,(a+b)/2,f,rule)
	Q2 = Calculator((a+b)/2,b,f,rule)


	error = abs(Q1 + Q2 - Q) / (pow(2,p)-1) 

	if error < e:
		#accept
		return Q1+Q2
	else:
		return adaptive(a,(a+b)/2,f,rule,e/2,Calculator,Q1) + adaptive((a+b)/2,b,f,rule,e/2,Calculator,Q2)


def Gaussian_calc(a,b,f,_):

	global operation_count
	global interval_count
	operation_count+=2
	interval_count+=1

	c = (b+a)/2.0
	length = b-a
	w = 1.0 * abs(b-a)/2 # universal for 2 point

	x_p = [c-length/(2*math.sqrt(3)),c+length/(2*math.sqrt(3))]

	return f(x_p[0])*w  + f(x_p[1])*w 




f = lambda x : 1/(0.01+pow(x-0.3,2)) + 1/(0.04 + pow(x-0.9,2)) - 6
#f = lambda x : 1

print("Helo")


a,b = 0,1
#rule = rule0
e = 0.0001

true_value = 29.858325395498671

print("e:",e)
print("[a,b]",a,b)


#############################################################
#Exercise 1-2
#############################################################
#Note to Mitya
# Miért marad el az adaptív ennyire?


print("Exercise 1-2:")
print("Rule 0 :",rule0)

clear_count()
start = time.time()
val = adaptive(a,b,f,rule0,e,CalculateQ)
duration = time.time() - start
print("Adaptive","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

clear_count()
start = time.time()
val = equidistant(a,b,f,rule0,CalculateQ,1000)
duration = time.time() - start
print("Equidistant","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

print("------------------------------------------------------------------")
print("Rule 1 :",rule1)

clear_count()
start = time.time()
val = adaptive(a,b,f,rule1,e,CalculateQ)
duration = time.time() - start
print("Adaptive","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

clear_count()
start = time.time()
val = equidistant(a,b,f,rule1,CalculateQ,50)
duration = time.time() - start
print("Equidistant","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

print("------------------------------------------------------------------")
print("Rule 2 :",rule2)

clear_count()
start = time.time()
val = adaptive(a,b,f,rule2,e,CalculateQ)
duration = time.time() - start
print("Adaptive","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

clear_count()
start = time.time()
val = equidistant(a,b,f,rule2,CalculateQ,70)
duration = time.time() - start
print("Equidistant","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

print("------------------------------------------------------------------")
print("Rule 3 :",rule3)

clear_count()
start = time.time()
val = adaptive(a,b,f,rule3,e,CalculateQ)
duration = time.time() - start
print("Adaptive","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

clear_count()
start = time.time()
val = equidistant(a,b,f,rule3,CalculateQ,20)
duration = time.time() - start
print("Equidistant","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################
#Exercise 3-4
#############################################################
#Note: I tried to find parameters for every method to behave with similar error rate
# therefore all method have similar error rates.

print("------------------------------------------------------------------\n")

print("Exercise 3-4:")



print("Gaussian Rule")

clear_count()
start = time.time()
val = adaptive(a,b,f,rule0,e,Gaussian_calc)
duration = time.time() - start
print("Adaptive","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

clear_count()
start = time.time()
val = equidistant(a,b,f,rule0,Gaussian_calc,100)
duration = time.time() - start
print("Equidistant","	Integrate: {} Error rate: {} Duration: {}s Intervals: {} Operation: {}".format(val,Error_rate(val,true_value),duration,interval_count,operation_count))

print("------------------------------------------------------------------")

quit()









