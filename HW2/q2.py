import numpy as np

# coordinates for the elements - different for each element
a = 0.0
b = 0.5

# different for each H
B = np.array([[1], [0], [0], [0]])

num_elem = 2;
p = np.array([[0.0, 0.5], [0.5, 1.0]])

for j in range(num_elem):
		for i in range(4):
				a = p[j, 0]
				b = p[j, 1]
				A = np.array([[a**3, a**2, a, 1], [b**3, b**2, b, 1], [3*(a**2), 2*a, 1, 0], [3*(b**2), 2*b, 1, 0]])
		
				B = np.empty([4, 1])
				B[i] = 1
				x = np.linalg.solve(A, B)
				print("For element {0}, H{1}:".format(j+1, i+1))
				print(x)
