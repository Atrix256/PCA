import matplotlib.pyplot as plt
import numpy as np

# Eigen Values & Eigen (Principle Component) Vectors:
# [992.216829] : [0.769509, 0.638636]
# [76.533171] : [0.638636, -0.769509]

data0 = [90.000000, 90.000000, 60.000000, 30.000000]
data1 = [60.000000, 90.000000, 60.000000, 30.000000]

recoveredData0 = [90.000000, 90.000000, 60.000000, 30.000000]
recoveredData1 = [60.000000, 90.000000, 60.000000, 30.000000]

fig = plt.figure()

plt.scatter(data0, data1, label="Points")
plt.scatter(recoveredData0, recoveredData1, marker="+", label="PCA Points")

plt.title("2 of 2 PCA Components\nRMSE = 0.000000")
plt.xlim([26.623888, 100.898345])
plt.ylim([23.827692, 93.151063])
plt.legend()
plt.tight_layout()

fig.savefig("test.2.png")
