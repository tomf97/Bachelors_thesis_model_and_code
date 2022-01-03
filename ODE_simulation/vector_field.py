#vector fields
import numpy as np
import matplotlib.pyplot as plt

k_range = m_range = np.linspace(0,1, 10) # from, to, points,
k,m = np.meshgrid(k_range,m_range)

#r = 0.4
#e = 0.5
#k1 = 0.7

#fk = k1*(1-k)*m-r*k
#fm = e*(1-m)-r*k

# with MM-kinetics:
r = 0.25
e = 0.5
k1 = 0.7
KM1 = 0.25
KM2 = 0.3
n = 3
fk = k1 * (1 - k)**n/(KM1**n + (1 - k)**n) * m - r * k/(KM2 + k)
fm = e * (1 - m) - r * k/(KM2 + k)

#calculate equilibrium
mBar = -k1 - e + k1*e/r + np.sqrt((-k1 - e + k1*e/r)**2 + 4*k1/r * e**2)
mBar = mBar / (2*k1*e/r)
kBar = (e/r) * (1-mBar)

fig = plt.figure(figsize=(8,4))

#another way
plt.streamplot(k,m,fk,fm, density=1.5)
#plt.plot(kBar,mBar, 'ro', label = 'Equilibrium')
#plt.legend(loc = 'upper right', shadow = True)
plt.xlabel("m(t)")
plt.ylabel("k(t)")
plt.xlim([0,1])
plt.ylim([0,1])
plt.savefig("figures/vector_field.jpg")
plt.show()

#basic construction
#middle vectors equally spaced instead of the beginning of the vector: pivot
#Q=plt.quiver(k,m,fk,fm, pivot="middle")
#plt.quiverkey(Q, .85,1.02,1,'=1m/s', labelpos="E")

#plt.title("Vector Field")
#plt.axis("scaled")
#plt.show()