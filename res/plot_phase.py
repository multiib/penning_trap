import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa



A1 = pa.mat() #Create pa.mat object (just as arma::mat in C++)
A1.load("./res/out/double.bin") #Load the content of the matrix you saved into your Python program.
A1 = np.array(A1) # Convert to numpy array

A2 = pa.mat() #Create pa.mat object (just as arma::mat in C++)
A2.load("./res/out/double.bin") #Load the content of the matrix you saved into your Python program.
A2 = np.array(A2) # Convert to numpy array



steps = int(len(A1))
particles = int(len(A1[0])/7)
values = 7

M1 = np.empty((steps, particles, values))
M2 = np.empty((steps, particles, values))

for i in range(7):
	M1[:,:,i] = A1[:,i::7]
	M2[:,:,i] = A2[:,i::7]
fig = plt.figure()
plt.plot(M1[:,0,1],M1[:,0,4], '-', color='#8a1629', linewidth=2.0, alpha=0.8)

plt.axis("equal")


plt.xlabel(r"$x$ [$\mu$m]", fontsize=12)
plt.ylabel(r"$v_x$ [$\mu$m]", fontsize=12)
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/phase_x.pdf")

fig = plt.figure()
plt.plot(M2[:,0,1],M2[:,0,4], '-', color='#8a1629', linewidth=2.0, alpha=0.8)

plt.axis("equal")


plt.xlabel(r"$x$ [$\mu$m]", fontsize=12)
plt.ylabel(r"$v_x$ [$\mu$m]", fontsize=12)
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/phase_z.pdf")




fig = plt.figure()
plt.plot(M1[:,0,3],M1[:,0,6], '-', color='#8a1629', linewidth=2.0, alpha=0.8)

plt.axis("equal")

plt.xlabel(r"$z$ [$\mu$m]", fontsize=12)
plt.ylabel(r"$v_z$ [$\mu$m]", fontsize=12)
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/phase_x_i.pdf")

fig = plt.figure()
plt.plot(M2[:,0,3],M2[:,0,6], '-', color='#8a1629', linewidth=2.0, alpha=0.8)

plt.axis("equal")


plt.xlabel(r"$z$ [$\mu$m]", fontsize=12)
plt.ylabel(r"$v_z$ [$\mu$m]", fontsize=12)
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/phase_z_i.pdf")