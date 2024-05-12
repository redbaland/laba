from array import *
import matplotlib.pyplot as plt
import numpy as np

H_array = array('i',[278,271,259,241,218,192,159])
Q_array = array('i',[1000,2000,3000,4000,5000,6000,7000])
H_Results = array('f',[0,0,0,0,0,0,0])

SumH = 0
SumQ2 = 0
SumHQ2 = 0
SumQ4 = 0

for i in range(0, 7):
    SumH = SumH + H_array[i]
    SumQ2 = SumQ2 + Q_array[i]*Q_array[i]
    SumHQ2 = SumHQ2 + H_array[i]*Q_array[i]*Q_array[i]
    SumQ4 = SumQ4 + Q_array[i]*Q_array[i]*Q_array[i]*Q_array[i]

a = (SumQ2*SumHQ2-SumH*SumQ4)/(SumQ2*SumQ2-7*SumQ4)
b = (7*SumHQ2 - SumH*SumQ2)/(SumQ2*SumQ2-7*SumQ4)

print(a)
print(b)

x = np.linspace(1000, 7000, 100)
plt.scatter(Q_array, H_array)
plt.plot(x, a-b*x*x)

plt.show()
