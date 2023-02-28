import numpy as np;
import matplotlib.pyplot as plot;
import csv



AR=10
alpha_stall = 15 *np.pi/180
alpha=np.linspace(0,alpha_stall,100)
alpha_s = np.linspace(alpha_stall,90* np.pi/180) 
CDmax = 1.11 + 0.018*AR
CL_stall = 1.2
CD_stall = 0.019
A1 = CDmax/2
A2 = (CL_stall -CDmax*np.sin(alpha_stall)*np.cos(alpha_stall))*np.sin(alpha_stall)/(np.cos(alpha_stall)**2)

B1=CDmax
B2 =(CD_stall- CDmax*np.sin(alpha_stall)**2)/(np.cos(alpha_stall))

cl=CL_stall*(alpha)/(alpha_stall)

cls = A1*np.sin(2*alpha_s)+A2*(np.cos(alpha_s)**2)/np.sin(alpha_s)
cd =B1*np.sin(alpha)**2 + B2*np.cos(alpha)
cds= B1*np.sin(alpha_s)**2 + B2*np.cos(alpha_s)

alpha_s2 = 180*np.pi/180 - alpha_s
cls2 = -cls
cl3 = -cl
alpha3 = 180*np.pi/180 - alpha
plot.plot(alpha*180/np.pi,cl)
plot.plot(alpha_s*180/np.pi,cls)
plot.plot(alpha_s2*180/np.pi,cls2)
plot.plot(alpha3*180/np.pi,cl3)
plot.show()

plot.plot(alpha*180/np.pi,cd)
plot.plot(alpha_s*180/np.pi,cds)
plot.plot(alpha_s2*180/np.pi,cds)
plot.plot(alpha3*180/np.pi,cd)
plot.show()