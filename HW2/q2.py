import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

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
				#print("For element {0}, H{1}:".format(j+1, i+1))
				#print(x)

				# verify correct calculation - zeroth derivatives
				left = x[0,0]*(a**3) + x[1,0]*(a**2) + x[2,0]*a + x[3,0]
				right = x[0,0]*(b**3) + x[1,0]*(b**2) + x[2,0]*b + x[3,0]
				
				# verify correct calculation - first derivatives
				left1 = 3.0 * x[0,0]*(a**2) + 2.0 * x[1,0]*a + x[2,0]
				right1 = 3.0 * x[0,0]*(b**2) + 2.0 * x[1,0]*b + x[2,0]

				#print("{0} {1} {2} {3}".format(left, right, left1, right1))
				
# plot the shape functions
xx = np.linspace(0.0, 0.5, 100)
xxx = np.linspace(0.5, 1.0, 100)

H1_1 = 16*(xx**3) - 12*(xx**2) + 1
H2_1 = -16*(xx**3) + 12*(xx**2)
H3_1 = 4*(xx**3) - 4*(xx**2) + xx
H4_1 = 4*(xx**3) - 2*(xx**2)

H1_2 = 16*(xxx**3) - 36*(xxx**2) + 24*xxx - 4
H2_2 = -16*(xxx**3) + 36*(xxx**2) - 24*xxx + 5
H3_2 = 4*(xxx**3) - 10*(xxx**2) + 8*xxx - 2
H4_2 = 4*(xxx**3) - 8*(xxx**2) + 5*xxx - 1

#plt.plot(xx, H1_1, linestyle='--', label='H1, e=1')
#plt.plot(xx, H2_1, linestyle='--', label='H2, e=1')
#plt.plot(xx, H3_1, linestyle='--',label='H3, e=1')
#plt.plot(xx, H4_1, linestyle='--',label='H4, e=1')
#plt.plot(xxx, H1_2, label='H1, e=2')
#plt.plot(xxx, H2_2, label='H2, e=2')
#plt.plot(xxx, H3_2, label='H3, e=2')
#plt.plot(xxx, H4_2, label='H4, e=2')
#plt.legend()

#plt.savefig('q2_Hermite_functions.png')

# compute elemntal stiffness matrix K
x = sp.Symbol('x')

H1_1 = 16*(x**3) - 12*(x**2) + 1
H2_1 = -16*(x**3) + 12*(x**2)
H3_1 = 4*(x**3) - 4*(x**2) + x
H4_1 = 4*(x**3) - 2*(x**2)

H1_2 = 16*(x**3) - 36*(x**2) + 24*x - 4
H2_2 = -16*(x**3) + 36*(x**2) - 24*x + 5
H3_2 = 4*(x**3) - 10*(x**2) + 8*x - 2
H4_2 = 4*(x**3) - 8*(x**2) + 5*x - 1

K = np.empty([4, 4])
for e in range(num_elem):
		for i in range(4):
				for j in range(4):
						integrand = sp.diff(eval("H{0}_{1}".format(i+1, e+1)), x, 2) * sp.diff(eval("H{0}_{1}".format(j+1, e+1)), x, 2)
						K[i, j] = sp.integrate(integrand, (x, p[e, 0], p[e, 1]))

