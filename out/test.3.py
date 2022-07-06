import matplotlib.pyplot as plt
import numpy as np

data0 = [90.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data1 = [60.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data2 = [90.000000, 30.000000, 60.000000, 90.000000, 30.000000]

recoveredData0 = [90.000000, 89.999985, 59.999996, 60.000000, 29.999998]
recoveredData1 = [60.000000, 89.999985, 60.000004, 60.000000, 30.000002]
recoveredData2 = [90.000008, 30.000011, 60.000000, 90.000000, 30.000000]

fig = plt.figure()

plt.scatter(data0, data1, label="Points")
plt.scatter(recoveredData0, recoveredData1, marker="+", label="PCA Points")

plt.title("3 of 3 PCA Components")
plt.xlim([26.562317, 102.191360])
plt.ylim([18.565441, 93.401642])
plt.legend()
plt.tight_layout()

fig.savefig("test.3.png")
