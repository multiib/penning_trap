import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

# Main frequency scan
A = pa.mat()
A.load("./res/out/rand.bin")
M = np.array(A)
x = np.linspace(0.2, 2.5, 116)

fig = plt.figure(figsize=(6,4.5))
plt.plot(x,M[0,:], '-', color='#8a1629', linewidth=2.0, alpha=0.8, label = "f = 0.1")
plt.plot(x,M[1,:], '-', color='#8b8229', linewidth=2.0, alpha=0.8, label = "f = 0.4")
plt.plot(x,M[2,:], '-', color='#2c1629', linewidth=2.0, alpha=0.8, label = "f = 0.7")

plt.xlabel(r"$\omega_V$ [MHz]", fontsize=16)
plt.ylabel(r"N", fontsize=16)
plt.legend()
plt.grid()
plt.savefig("./res/figs/plots/freq-F.pdf")


# Zoomed frequancy scan without interactions
A = pa.mat()
A.load("./res/out/rand_Z_F.bin")
M = np.array(A)
x = np.linspace(1.1, 1.7, 121)

fig = plt.figure(figsize=(6,4.5))
plt.plot(x,M[0,:], '-', color='#8a1629', linewidth=2.0, alpha=0.8, label = "f = 0.1")
plt.plot(x,M[1,:], '-', color='#8b8229', linewidth=2.0, alpha=0.8, label = "f = 0.4")
plt.plot(x,M[2,:], '-', color='#2c1629', linewidth=2.0, alpha=0.8, label = "f = 0.7")

plt.xlabel(r"$\omega_V$ [MHz]", fontsize=16)
plt.ylabel(r"N", fontsize=16)
plt.legend()
plt.grid()
plt.savefig("./res/figs/plots/freq-ZOOM-F.pdf")


# Zoomed frequancy scan with interactions
A = pa.mat()
A.load("./res/out/rand_Z_T.bin")
M = np.array(A)
x = np.linspace(1.1, 1.7, 121)

fig = plt.figure(figsize=(6,4.5))
plt.plot(x,M[0,:], '-', color='#8a1629', linewidth=2.0, alpha=0.8, label = "f = 0.1")
plt.plot(x,M[1,:], '-', color='#8b8229', linewidth=2.0, alpha=0.8, label = "f = 0.4")
plt.plot(x,M[2,:], '-', color='#2c1629', linewidth=2.0, alpha=0.8, label = "f = 0.7")

plt.xlabel(r"$\omega_V$ [MHz]", fontsize=16)
plt.ylabel(r"N", fontsize=16)
plt.legend()
plt.grid()
plt.savefig("./res/figs/plots/freq-ZOOM-T.pdf")

