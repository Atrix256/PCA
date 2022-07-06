import matplotlib.pyplot as plt
import numpy as np

data0 = [90.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data1 = [60.000000, 90.000000, 60.000000, 60.000000, 30.000000]
data2 = [90.000000, 30.000000, 60.000000, 90.000000, 30.000000]

recoveredData0 = [92.251221, 76.257820, 67.130188, 79.348923, 33.565094]
recoveredData1 = [60.374939, 49.907860, 43.934174, 51.930870, 21.967087]
recoveredData2 = [87.363731, 72.217667, 63.573616, 75.144997, 31.786808]

fig = plt.figure()

plt.scatter(data0, data1, label="Points")
plt.scatter(recoveredData0, recoveredData1, marker="+", label="PCA Points")

plt.title("1 of 3 PCA Components")
plt.xlim([26.562317, 102.191360])
plt.ylim([18.565441, 93.401642])
plt.legend()
plt.tight_layout()

fig.savefig("test.1.png")
