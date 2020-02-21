import integration1d
import interpolation1d
import numpy as np
import matplotlib.pyplot as plt
from loadfunction import coordinates
inter = interpolation1d.interpolation()
resultants, torques = integration1d.resultandtorques()
print(torques)
y = input()
def torque(x):

    y = inter.interpolate(x,0,'chord',torques)
    return y
def resultant(x):

    y = inter.interpolate(x,0,'chord',resultants)
    return y
x1 = np.linspace(0,2.76,400)
torquefinal = []
resultfinal = []
for element in x1:
    torquefinal.append(torque(element))
    resultfinal.append(resultant(element))

plt.plot(x1,torquefinal)
plt.plot(x1,resultfinal)
plt.show()



