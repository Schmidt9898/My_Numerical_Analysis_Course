
import numpy as np


# quadrature rules 
#		M,p,quadrature weights
rule0 = (2,2,[1/2,1/2])
rule1 = (3,4,[1/6,2/3,1/6])
rule2 = (4,4,[1/8,3/8,3/8,1/8])
rule3 = (5,6,[7/90,16/45,4/30,16/45,7/90])


def equidistant(a,b,f,N=1000):
	print(f"a{a}, b{b}, N{N}")
	print("equidistant")

	X = list(np.linspace(a, b, N, endpoint=True))
	Y = [f(x) for x in X]

	w = float(b-a) / N

	sum_=0
	for i,y in enumerate(Y):
		sum_+=y*w
	print(sum_)

	
	return sum_

def adaptive(a,b,f,rule):
	return None

f = lambda x : 1/(0.01+pow(x-0.3,2)) + 1/(0.04 + pow(x-0.9,2)) - 6

print("Helo")


a,b = 0,1

equidistant(a,b,f,10000)



true_value = 29.858325395498671






