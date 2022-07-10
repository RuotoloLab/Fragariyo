"""
How to correlate two dimensional arrays
Author: CRR
02/29/20
"""

import numpy as np
import scipy
from scipy import signal
from matplotlib import pyplot as plt

array_one = [0,1,2,3,4,5,4,3,2,1,0]
array_two = [0,1,2,3,4,5,4,3,2,1,0]
array_three = [1,2,3,4,5,5,4,3,2,1,0]
array_four = [1,2,30,4,5,5,90,3,2,1,100]

# print(np.correlate(array_one,array_two))
# print(np.correlate(array_one,array_three))
# print(np.correlate(array_one,array_four))

print(np.corrcoef(array_one,array_two))
print(np.corrcoef(array_one,array_three))
print(np.corrcoef(array_one,array_four))


plt.plot(array_two, label= "Two")
plt.plot(array_three, label= "Three")
plt.plot(array_four, label= "Four")
plt.legend(loc='best')

plt.show()


# print(scipy.signal.correlate(array_one, array_two))
# print(scipy.signal.correlate(array_one, array_four, mode="full"))