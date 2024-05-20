# This example loads pre-calculated PTO force histories for a reduced size joint
# probability distribution. Then uses the WDRT fatigue function to calculate
# equivalent fatigue loads for a 1 hour timeframe and a 1 year timeframe.

import numpy as np
import WDRT.fatigue as fatigue
import os

# Reduced size joint probability distribution
# Wave energy periods
Te = [[6.0, 10.0, 14.0], [6.0, 10.0, 14.0], [6.0, 10.0, 14.0]]
# Significant wave heights
Hs = [[1.25, 1.25, 1.25], [2.75, 2.75, 2.75], [4.25, 4.25, 4.25]]
# Probability
P = np.multiply([[23.45, 24.78, 1.64], [9.18, 28.21, 4.11],
                 [0.05, 5.00, 2.34]], 0.01)

N1h = 0
N1y = 0
[h, t] = np.shape(P)
for i in range(h):
    for j in range(t):
        # Average N in 1 hour (Tavg = 0.82476*Te)
        N1h = N1h + 1 * 60 * 60 * P[i][j] / (Te[i][j] * .82476)
        # Average N in 1 year
        N1y = N1y + 1 * 365 * 24 * 60 * 60 * P[i][j] / (Te[i][j] * .82476)

# Assume an S-N curve slope of 6 (representative of cast iron)
m = float(6)
Feq_1h = np.zeros((h, t))
Feq_1y = 0
for i in range(h):
    for j in range(t):
        # Read pre-calculated PTO force histories for each sea state
        Fpto = np.loadtxt(os.path.join(r'data','FptoH') +
                          str(int(Hs[i][j])) +
                          'T' + str(int(Te[i][j])) + '.txt')
        # Equivalent fatigue load for a 1 hour timeframe
        Feq_1h[i][j] = fatigue.EqLoad(Fpto, N1h, m)
        Feq_1y = Feq_1y + (Feq_1h[i][j]**m) * N1y * P[i][j]
# Equivalent fatigue load for a 1 year timeframe
Feq_1y = (Feq_1y / N1y)**(1 / m)

print('1 hour equivalent fatigue loads: (in Newtons)')
print(Feq_1h)
print('1 year equivalent fatigue load: (in Newtons)')
print(Feq_1y)
