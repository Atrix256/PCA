import matplotlib.pyplot as plt
import numpy as np

# Eigen Values & Eigen (Principle Component) Vectors:
# [992.216829] : [0.769509, 0.638636]

data0 = [90.000000, 90.000000, 60.000000, 30.000000]
data1 = [60.000000, 90.000000, 60.000000, 30.000000]

recoveredData0 = [82.779150, 97.522232, 65.014822, 32.507411]
recoveredData1 = [68.700592, 80.936264, 53.957510, 26.978755]

fig = plt.figure()

plt.scatter(data0, data1, label="Points")
plt.scatter(recoveredData0, recoveredData1, marker="+", label="PCA Points")

plt.title("1 of 2 PCA Components\nRMSE = 9.268919")
plt.xlim([26.623888, 100.898345])
plt.ylim([23.827692, 93.151063])
plt.legend()
plt.tight_layout()

fig.savefig("test.1.png")
