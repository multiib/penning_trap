import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

np.seterr(divide='ignore')



A = pa.mat()
A.load("./res/out/re_AN/4000.bin")
a4 = np.array(A).T

A = pa.mat()
A.load("./res/out/re_AN/8000.bin")
a8 = np.array(A).T

A = pa.mat()
A.load("./res/out/re_AN/16000.bin")
a16 = np.array(A).T

A = pa.mat()
A.load("./res/out/re_AN/32000.bin")
a32 = np.array(A).T

A = pa.mat()
A.load("./res/out/re_FE/4000.bin")
n4 = np.array(A)

A = pa.mat()
A.load("./res/out/re_FE/8000.bin")
n8 = np.array(A)

A = pa.mat()
A.load("./res/out/re_FE/16000.bin")
n16 = np.array(A)

A = pa.mat()
A.load("./res/out/re_FE/32000.bin")
n32 = np.array(A)


def relative_error(a, n):
	return abs(np.linalg.norm(a-n)/np.linalg.norm(a))


fig = plt.figure(figsize=(6,4.5))

# 4000
a = np.column_stack((a4[:,1],a4[:,2],a4[:,3]))
n = np.column_stack((n4[:,1],n4[:,2],n4[:,3]))

delta_max4 = 0
for i in range(len(a)):
	delta = np.linalg.norm(a[i]-n[i])
	if delta > delta_max4:
		delta_max4 = delta

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a4[:,0],(re), '-', color='#8a1629', linewidth=2.5, alpha=0.8, label = "n = 4000")

# 8000
a = np.column_stack((a8[:,1],a8[:,2],a8[:,3]))
n = np.column_stack((n8[:,1],n8[:,2],n8[:,3]))

delta_max8 = 0
for i in range(len(a)):
	delta = np.linalg.norm(a[i]-n[i])
	if delta > delta_max8:
		delta_max8 = delta

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a8[:,0],(re), '-', color='#8b8229', linewidth=2.5, alpha=0.8, label = "n = 8000")

# 16000
a = np.column_stack((a16[:,1],a16[:,2],a16[:,3]))
n = np.column_stack((n16[:,1],n16[:,2],n16[:,3]))

delta_max16 = 0
for i in range(len(a)):
	delta = np.linalg.norm(a[i]-n[i])
	if delta > delta_max16:
		delta_max16 = delta


re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a16[:,0],(re), '-', color='#2c1629', linewidth=2.5, alpha=0.8, label = "n = 16000")

# 32000
a = np.column_stack((a32[:,1],a32[:,2],a32[:,3]))
n = np.column_stack((n32[:,1],n32[:,2],n32[:,3]))

delta_max32 = 0
for i in range(len(a)):
	delta = np.linalg.norm(a[i]-n[i])
	if delta > delta_max32:
		delta_max32 = delta

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a32[:,0],(re), '-', color='#4d8229', linewidth=2.5, alpha=0.8, label = "n = 32000")
plt.legend()
plt.xlabel(r"$t$ [$\mu s$]", fontsize=16)
plt.ylabel(r"$E_r$", fontsize=16)
plt.grid()
plt.savefig("./res/figs/plots/rel_FE.pdf")


delta_max = [delta_max4, delta_max8, delta_max16, delta_max32]
h = [50/4000, 50/8000, 50/16000, 50/32000]
rel_err = 0

for k in range(1, 4):
	rel_err += np.log(delta_max[k]/delta_max[k-1])/(np.log(h[k]/h[k-1]))
print((1/3)*rel_err)