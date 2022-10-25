import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa



A = pa.mat() #Create pa.mat object (just as arma::mat in C++)
A.load("./res/out/double.bin") #Load the content of the matrix you saved into your Python program.
A = np.array(A) # Convert to numpy array



steps = int(len(A))
particles = int(len(A[0])/7)
values = 7

M = np.empty((steps, particles, values))

for i in range(7):
	M[:,:,i] = A[:,i::7]
	
fig = plt.figure()
plt.plot(M[:,1,1],M[:,1,2], '-', color='#8a1629', linewidth=2.0, alpha=0.8)
plt.axis("equal")
plt.title("Fig. X")
fig.suptitle("test title", fontsize=20)
plt.xlabel(r"x [$\mu$m]", fontsize=8)
plt.ylabel(r"y [$\mu$m]", fontsize=8)

plt.grid()
plt.show()
