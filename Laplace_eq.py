# Numerical Laplace Equation Solution 
import numpy as np
import matplotlib.pyplot as plt


maxIter = 500

# Set Dimension and delta
lenX = lenY = 11 #we set it rectangular
delta = 1

# Boundary condition
Phi_top = 0 #y=1
Phi_bottom = 1#y=0
Phi_left = 0 #x=1
Phi_right = 0 #x=0

# Initial guess of interior grid
Phi_guess = 0.5

#for plot
colorinterpolation = 50
colourMap = plt.cm.jet #you can try: colourMap = plt.cm.coolwarm

# set grid
X, Y = np.meshgrid(np.arange(0, lenX), np.arange(0, lenY))

# Set array size and set the interior value with Tguess
Phi = np.empty((lenX, lenY))
Phi.fill(Phi_guess)

# Set Boundary condition
Phi[(lenY-1):, :] = Phi_top
Phi[:1, :] = Phi_bottom
Phi[:, (lenX-1):] = Phi_right
Phi[:, :1] = Phi_left


for iteration in range(0, maxIter):
    for i in range(1, lenX-1, delta):
        for j in range(1, lenY-1, delta):
            Phi[i, j] = 0.25 * (Phi[i+1][j] + Phi[i-1][j] + Phi[i][j+1] + Phi[i][j-1])

print("Iteration finished")


plt.title("Voltage contour ")
plt.contourf(X/10, Y/10, Phi, colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()

print("")