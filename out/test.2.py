import matplotlib.pyplot as plt
import numpy as np

data0 = [90.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data1 = [60.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data2 = [90.000000, 30.000000, 60.000000, 90.000000, 30.000000]

recoveredData0 = [91.063171, 98.753677, 70.324493, 73.691025, 35.162247]
recoveredData1 = [58.785641, 80.001495, 48.207329, 44.362068, 24.103664]
recoveredData2 = [89.716568, 27.666317, 57.247536, 86.350037, 28.623768]

fig = plt.figure()

plt.scatter(data0, data1, label="Points")
plt.scatter(recoveredData0, recoveredData1, marker="+", label="PCA Points")

plt.title("2 of 3 PCA Components")
plt.xlim([26.562317, 102.191360])
plt.ylim([18.565441, 93.401642])
plt.legend()
plt.tight_layout()

fig.savefig("test.2.png")
