import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarma as pa

"""Load data"""
A = pa.mat() #Create pa.mat object (just as arma::mat in C++)
A.load("./res/out/double.bin") #Load the content of the matrix you saved into your Python program.
A = np.array(A) # Convert to numpy array

B = pa.mat() #Create pa.mat object (just as arma::mat in C++)
B.load("./res/out/double_interactions.bin") #Load the content of the matrix you saved into your Python program.
B = np.array(B) # Convert to numpy array


"""Reshape data"""
steps = int(len(A))
particles = int(len(A[0])/7)
values = 7

M = np.empty((steps, particles, values))
I = np.empty((steps, particles, values)) # with interactions
for i in range(7):
	M[:,:,i] = A[:,i::7]
	I[:,:,i] = B[:,i::7]

"""Plots"""
# Data is structured like this ** [steps, particle, axis] **

# Axis
t, x, y, z, vx, vy, vz = 0, 1, 2, 3, 4, 5, 6


# Plot t, z    | particle 1     | interactions off
fig = plt.figure(figsize=(6,4.5))
plt.plot(M[:,0,t],M[:,0,z], "-", color="#8a1629", linewidth=2.0, alpha=0.8)
plt.axvline((200*np.pi)/69, -25, 25, ls="--", color="black")

plt.xlabel(r"t [$\mu s$]", fontsize=16)
plt.ylabel(r"z [$\mu m$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.savefig("./res/figs/plots/tx-1-F.pdf")



# Plot x, y    | particle 1 & 2 | interactions off
fig = plt.figure(figsize=(6,4.5))
plt.plot(M[:,0,x],M[:,0,y], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(M[:,1,x],M[:,1,y], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"x [$\mu m$]", fontsize=16)
plt.ylabel(r"y [$\mu m$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/xy-12-F.pdf")

# Plot x, y    | particle 1 & 2 | interactions on
fig = plt.figure(figsize=(6,4.5))
plt.plot(I[:,0,x],I[:,0,y], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(I[:,1,x],I[:,1,y], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"x [$\mu m$]", fontsize=16)
plt.ylabel(r"y [$\mu m$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/xy-12-T.pdf")

# Plot x, vx   | particle 1 & 2 | interactions off
fig = plt.figure(figsize=(6,4.5))
plt.plot(M[:,0,x],M[:,0,vx], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(M[:,1,x],M[:,1,vx], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"x [$\mu m$]", fontsize=16)
plt.ylabel(r"vx [$\frac{\mu m}{\mu s}$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/xvx-12-F.pdf")


# Plot x, vx   | particle 1 & 2 | interactions on
fig = plt.figure(figsize=(6,4.5))
plt.plot(I[:,0,x],I[:,0,vx], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(I[:,1,x],I[:,1,vx], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"x [$\mu m$]", fontsize=16)
plt.ylabel(r"vx [$\frac{\mu m}{\mu s}$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/xvx-12-T.pdf")

# Plot z, vz   | particle 1 & 2 | interactions off
fig = plt.figure(figsize=(6,4.5))
plt.plot(M[:,0,z],M[:,0,vz], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(M[:,1,z],M[:,1,vz], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"z [$\mu m$]", fontsize=16)
plt.ylabel(r"vz [$\frac{\mu m}{\mu s}$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/zvz-12-F.pdf")

# Plot z, vz   | particle 1 & 2 | interactions on
fig = plt.figure(figsize=(6,4.5))
plt.plot(I[:,0,z],I[:,0,vz], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
plt.plot(I[:,1,z],I[:,1,vz], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
plt.xlabel(r"z [$\mu m$]", fontsize=16)
plt.ylabel(r"vz [$\frac{\mu m}{\mu s}$]", fontsize=16)
plt.axis("equal")
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/zvz-12-T.pdf")

# Plot x, y, x | particle 1 & 2 | interactions off
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot3D(M[:,0,x],M[:,0,y],M[:,0,z], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
ax.plot3D(M[:,1,x],M[:,1,y],M[:,1,z], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
ax.set_xlabel(r"x [$\mu m$]", fontsize=16)
ax.set_ylabel(r"y [$\mu m$]", fontsize=16)
ax.set_zlabel(r"z [$\mu m$]", fontsize=16)
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/3D-12-F.pdf")

# Plot x, y, x | particle 1 & 2 | interactions on
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot3D(I[:,0,x],I[:,0,y],I[:,0,z], "-", color="#8a1629", linewidth=2.0, alpha=0.8, label = "Particle 1")
ax.plot3D(I[:,1,x],I[:,1,y],I[:,1,z], "-", color="#8b8229", linewidth=2.0, alpha=0.8, label = "Particle 2")
ax.set_xlabel(r"x [$\mu m$]", fontsize=16)
ax.set_ylabel(r"y [$\mu m$]", fontsize=16)
ax.set_zlabel(r"z [$\mu m$]", fontsize=16)
plt.tight_layout()
plt.grid()
plt.legend()
plt.savefig("./res/figs/plots/3D-12-T.pdf")
