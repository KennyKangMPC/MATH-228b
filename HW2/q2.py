import numpy as np

# number of elements
num_elem = 2;

# coordinates of the nodes
p = np.array([[0.0, 0.5], [0.5, 1.0]])

for j in range(num_elem):
		for i in range(4):
				a = p[j, 0]
				b = p[j, 1]
				A = np.array([[a**3, a**2, a, 1], [b**3, b**2, b, 1], \
						[3*(a**2), 2*a, 1, 0], [3*(b**2), 2*b, 1, 0]])
		
				B = np.empty([4, 1])
				B[i] = 1
				x = np.linalg.solve(A, B)
				print("For element {0}, H{1}:".format(j+1, i+1))
				print(x)

				# verify correct calculation - zeroth derivatives
				left = x[0,0]*(a**3) + x[1,0]*(a**2) + x[2,0]*a + x[3,0]
				right = x[0,0]*(b**3) + x[1,0]*(b**2) + x[2,0]*b + x[3,0]
				
				# verify correct calculation - first derivatives
				left1 = 3.0 * x[0,0]*(a**2) + 2.0 * x[1,0]*a + x[2,0]
				right1 = 3.0 * x[0,0]*(b**2) + 2.0 * x[1,0]*b + x[2,0]

				print("{0} {1} {2} {3}".format(left, right, left1, right1))
