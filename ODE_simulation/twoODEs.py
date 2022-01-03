import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#def model(z,t):
#    k = z[0]
#    m = z[1]
#    r = 0.75
#    e = 1
#    k1 = 0.5
#    dk = k1 * (1 - k) * m - r * k
#    dm = e * (1 - m) - r * k
#    return [dk,dm]
def model(z,t):
    k = z[0]
    m = z[1]
    r = 0.75
    e = 1
    k1 = 0.5
    KM1 = 0.2
    KM2 = 0.75
    n = 3
    dk = k1 * (1 - k)**n/(KM1**n + (1 - k)**n) * m - r * k/(KM2 + k)
    dm = e * (1 - m) - r * k/(KM2 + k)
    return [dk,dm]

z0=[0.3, 0.3]
t=np.linspace(0,10)

z = odeint(model,z0,t)

k = z[:,0]
m = z[:,1]

#plot results
fig = plt.figure(figsize=(8,4))
plt.plot(t,k,'r-', label = 'k(t)')
plt.plot(t,m,'b--', label = 'm(t)')
#plt.plot(t,np.zeros(len(x)),'g--')
plt.xlabel('time')
plt.ylabel('concentration')
plt.ylim([0,1])
plt.legend(loc = 'upper right', shadow = True)
plt.savefig("figures/ode_result.jpg")
plt.show()
