import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa



A = pa.mat() #Create pa.mat object (just as arma::mat in C++)
A.load("./res/out/double_interactions.bin") #Load the content of the matrix you saved into your Python program.
A = np.array(A) # Convert to numpy array



steps = int(len(A))
particles = int(len(A[0])/7)
values = 7

M = np.empty((steps, particles, values))

for i in range(7):
	M[:,:,i] = A[:,i::7]
	
fig = plt.figure()
plt.plot(M[:,0,1],M[:,0,2], '-', color='#8a1629', linewidth=2.0, alpha=0.8, label = "Particle. 1")
plt.plot(M[:,1,1],M[:,1,2], '-', color='#8a8229', linewidth=2.0, alpha=0.8, label = "Particle. 2")
plt.axis("equal")
plt.title("Fig. 3")

plt.legend()
plt.xlabel(r"x [$\mu$m]", fontsize=12)
plt.ylabel(r"y [$\mu$m]", fontsize=12)
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/double_interactions.pdf")
plt.show()