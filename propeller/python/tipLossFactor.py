import numpy as n
import matplotlib.pyplot as p

R = 0.25
Nb=2
r=n.linspace(0,0.25,100)
Vn=30
Vt=5000*n.pi/30*r

F_tip = 2/n.pi * n.arccos(n.exp(-Nb/2*(R-r)/R*n.sqrt(Vn**2+Vt**2)/Vn))

p.plot(r,F_tip)
p.show()