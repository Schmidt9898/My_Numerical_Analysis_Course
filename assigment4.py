import math
import numpy as np;


def create_matrix(n):
	print(n)
	mat = np.diag([-1,-1,0])
	#np.concatenate((([-1]*(n-1)), [0]), axis=0)
	#np.concatenate(([0],([-1]*(n-1))), axis=0)
#	mat = np.matrix(mat)
	mat = np.roll(np.diag(np.concatenate((([-1]*(n-1)), [0]), axis=0)),1)+ np.diag([2]*n) +np.roll(np.diag(np.concatenate(([0],([-1]*(n-1))), axis=0)),-1)
	#mat = pow(n+1,2)*mat 
	print(mat)
	return mat

def sign(a):
	return 1 if a >=0 else -1

def inflates_matrix(A,n):
	A = np.asmatrix(A)
	size = A.shape[0]
	if size >= n:
		return A
	if size+1 < n:
		A =  np.asmatrix(inflates_matrix(A,n-1))
		size = A.shape[0]
	A = np.concatenate((np.asmatrix([0]*(n-1)),A), axis=0)
	A = np.concatenate(( np.transpose(np.asmatrix(np.identity(n)[0])) , A ), axis=1)
	return A

def QRdecomp(A):
	A = np.asmatrix(A)
	m = A.shape[0] #row number
	n = A.shape[1] #collum number
	#print(A.shape)
	#A_T = np.transpose(A) # for the collum vectors

	H_matrices = []


	for i in range(n): 
		ni = m-i
		
		a = np.transpose( np.array( A[i:m,i] ) )
		
		v = a + sign(a[0,0]) * np.linalg.norm(a) * np.identity(ni)[0]	
		v = np.asmatrix(v)
		vvt = np.transpose(v)*v
		vtv = np.dot(v[0],np.transpose(v))
		
		H = np.identity(ni) - 2 * (vvt/vtv)

		H = inflates_matrix(H,m)

		H_matrices.append(H)
		
		#print(H)
		#print("*****")
	
	#get R
	H_T = np.identity(m)
	for H in reversed(H_matrices):
		H_T = H_T*H
	R = H_T*A

	Q = np.identity(m)
	for H in H_matrices:
		Q = Q*H

	R = np.asarray(R)
	R = R[0:n,0:n]
	Q = Q[0:m,0:n]

	return Q,R


def get_eigen_values(A0):
	n = A0.shape[1]

	A = A0
	for i in range(2):
		Q,R = QRdecomp(A)
		Q_T = np.transpose(Q)
		A = Q_T  * A * Q

	print(A)


#def QRdecomp(A):
#	Q=R=A
#	n = A.shape[1] #collum number
#	print(A.shape)
#	A = np.transpose(A) # for the collum vectors
#	V = []
#	for i in range(n): 
#		v = np.array(A[i])
#		for e in V:
#			alpha = sum(np.dot(np.array(A[i]),np.transpose(e)))
#			#print(alpha)
#			v = v - (e*alpha) 
#		v= v / np.linalg.norm(v)
#		#print(v)
#		V.append(v)
#	print(V)
#	return Q,R
#


#alpha = np.dot([0,1,2],[0,1,1])
#print(alpha)

#quit()


print("assigment4")
n = 10
A = create_matrix(n)

#A = np.asmatrix([[3,0],[0,1]])

#m = inflates_matrix(np.asmatrix([3]),2)
#m = inflates_matrix(m,2)
#m = inflates_matrix(np.asmatrix([3]),4)



res_eigen = get_eigen_values(A)


gen_eig = lambda k,x : 4*pow(n+1,2)*pow(math.sin((k*math.pi)/(2*(n+1))),2)


gen_eig(0,10)

#Q,R = QRdecomp(A)
#print(Q)
#print("----------")
#print(R)

R.shape[0]

for i in range(R.shape[0]):
	for j in range(R.shape[1]):
		if i == j:
			print(R[i,j])



#Q,R = np.linalg.qr(A, mode='reduced')
#print(Q)
#print(R)







