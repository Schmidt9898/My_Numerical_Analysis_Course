import numpy as np


def create_matrix(n):
	print(n)
	mat = np.diag([-1,-1,0])
	np.concatenate((([-1]*(n-1)), [0]), axis=0)
	np.concatenate(([0],([-1]*(n-1))), axis=0)
#	mat = np.matrix(mat)
	mat = np.roll(np.diag(np.concatenate((([-1]*(n-1)), [0]), axis=0)),1)+ np.diag([2]*n) +np.roll(np.diag(np.concatenate(([0],([-1]*(n-1))), axis=0)),-1)
	#mat = pow(n+1,2)*mat 
	print(mat)
	return mat




def QRdecomp(A):
	Q=R=A
	n = A.shape[1] #collum number
	print(A.shape)
	A = np.transpose(A) # for the collum vectors
	V = []
	for i in range(n): 
		v = np.array(A[i])
		for e in V:
			alpha = sum(np.dot(np.array(A[i]),np.transpose(e)))
			#print(alpha)
			v = v - (e*alpha) 
		v= v / np.linalg.norm(v)
		#print(v)
		V.append(v)
	print(V)
		



	return Q,R



#alpha = np.dot([0,1,2],[0,1,1])
#print(alpha)

#quit()


print("assigment4")
A = create_matrix(3)

A = np.asmatrix([[1,1,1],[0,2,0],[0,0,3]])

Q,R = QRdecomp(A)
#Q,R = np.linalg.qr(A, mode='reduced')
#print(Q)
#print(R)







