import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
import pandas as pd

x = np.array([-0.219, -0.073, 0.031, 0.1]).reshape(4,1)
y = np.array([0.017, 0.021, 0.025, 0.03]).reshape(4,1)
lr = linear_model.LinearRegression()
lr.fit(x, y)
print('Coefficients: \n', lr.coef_)
plt.scatter(x, y, color = 'black')
plt.plot(x, lr.predict(x), color = 'blue', linewidth = 3)
plt.ylabel('Nz (10^20)')
plt.xlabel('Gamma_z (10^20 s^-1)')
plt.title('tau_z =' + str(np.round(lr.coef_[0][0], 4)*1000.)+'ms')
plt.show()
