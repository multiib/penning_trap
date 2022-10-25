import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa



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




fig = plt.figure()

# 4000
a = np.column_stack((a4[:,1],a4[:,2],a4[:,3]))
n = np.column_stack((n4[:,1],n4[:,2],n4[:,3]))

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a4[:,0],np.log(re), '-', color='#2A5C5E', linewidth=2.5, alpha=0.8, label = "n = 4000")

# 8000
a = np.column_stack((a8[:,1],a8[:,2],a8[:,3]))
n = np.column_stack((n8[:,1],n8[:,2],n8[:,3]))

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a8[:,0],np.log(re), '-', color='#40455E', linewidth=2.5, alpha=0.8, label = "n = 8000")

# 16000
a = np.column_stack((a16[:,1],a16[:,2],a16[:,3]))
n = np.column_stack((n16[:,1],n16[:,2],n16[:,3]))

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a16[:,0],np.log(re), '-', color='#B36E5D', linewidth=2.5, alpha=0.8, label = "n = 16000")

# 32000
a = np.column_stack((a32[:,1],a32[:,2],a32[:,3]))
n = np.column_stack((n32[:,1],n32[:,2],n32[:,3]))

re = np.zeros(len(a))
for i in range(len(a)):
	re[i] = relative_error(a[i],n[i])

plt.plot(a32[:,0],np.log(re), '-', color='#E0AB6C', linewidth=2.5, alpha=0.8, label = "n = 32000")

plt.legend()
plt.axis("equal")
plt.title("Fig. X")
fig.suptitle("test title", fontsize=20)
plt.xlabel(r"time [$\mu$s]", fontsize=8)
plt.ylabel(r"relative error (logn)", fontsize=8)

plt.grid()
plt.show()

