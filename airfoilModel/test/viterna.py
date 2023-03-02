import numpy as np;
import matplotlib.pyplot as plot;
import csv



AR=20
alpha_stall = 15 *np.pi/180
alphas=np.linspace(-180,180,300)

alpha_stall_n = -15
alpha_stall_p = 15
CL_stall_n = -0.7
CL_stall_p = 1.2
CDmax = 1.11 + 0.018*AR
CD_stall_p = 0.05
CD_stall_n = 0.05

#compute coeficients A1 A2 B1 B2 positive alpha
alpha_stall = alpha_stall_p * np.pi/180
CL_stall = CL_stall_p
CD_stall = CD_stall_p
A1_p = CDmax/2
A2_p = (CL_stall -CDmax*np.sin(alpha_stall)*np.cos(alpha_stall))*np.sin(alpha_stall)/(np.cos(alpha_stall)**2)
B1_p = CDmax
B2_p = (CD_stall- CDmax*np.sin(alpha_stall)**2)/(np.cos(alpha_stall))

#compute coeficients A1 A2 B1 B2 negative alpha
alpha_stall = -alpha_stall_n* np.pi/180
CL_stall = -CL_stall_n 
CD_stall = CD_stall_p
A1_n = CDmax/2
A2_n = (CL_stall -CDmax*np.sin(alpha_stall)*np.cos(alpha_stall))*np.sin(alpha_stall)/(np.cos(alpha_stall)**2)
B1_n = CDmax
B2_n = (CD_stall- CDmax*np.sin(alpha_stall)**2)/(np.cos(alpha_stall))

#init arrays
cl=np.zeros(len(alphas))
cd=np.zeros(len(alphas))

for i in range(len(alphas)):
    alpha = alphas[i]
    #check alpha between (-180 , 180]
    if(alpha>=alpha_stall_n and alpha<=alpha_stall_p):
        #datos
        da = (alpha_stall_p-alpha_stall_n)
        dcl = CL_stall_p - CL_stall_n
        #linear interpolatio between endpoints (just for testing)
        cl[i] = CL_stall_n + dcl/da * (alpha - alpha_stall_n)
        cd[i] = 0.01 + (CD_stall-0.01)/(alpha_stall_p)**2 * alpha**2
    elif(alpha>alpha_stall_p and alpha <=90):
        #viterna
        alpha = alpha * np.pi/180
        cl[i]=A1_p*np.sin(2*alpha)+A2_p*(np.cos(alpha)**2)/np.sin(alpha)
        cd[i] = B1_p*np.sin(alpha)**2 + B2_p*np.cos(alpha)
    elif(alpha>90 and alpha<=180) or (alpha>=-180 and alpha<-90):
        #flat plate
        alpha = alpha * np.pi/180
        cl[i]=2*np.abs(CL_stall_n)*np.sin(alpha)*np.cos(alpha)
        cd[i]= B1_n*np.sin(alpha)**2 + B2_n*np.cos(alpha)
    else:
        #viterna negative
        alpha = -alpha
        alpha = alpha * np.pi/180
        cl[i]= -(A1_n*np.sin(2*alpha)+A2_n*(np.cos(alpha)**2)/np.sin(alpha))
        cd[i] = B1_n*np.sin(alpha)**2 + B2_n*np.cos(alpha)


plot.plot(alphas,cl*np.cos(alphas* np.pi/180)+cd*np.sin(alphas* np.pi/180))
plot.plot(alphas,cd*np.cos(alphas* np.pi/180)-cl*np.sin(alphas* np.pi/180))
plot.grid()
plot.show()




#cl=(CL_stall+0.2)*(alpha)/(alpha_stall)
#cl[-1]=CL_stall
#ls = A1*np.sin(2*alpha_s)+A2*(np.cos(alpha_s)**2)/np.sin(alpha_s)
#ccd =B1*np.sin(alpha)**2 + B2*np.cos(alpha)
#cds= B1*np.sin(alpha_s)**2 + B2*np.cos(alpha_s)
#
#alpha_s2 = 180*np.pi/180 - alpha_s
#cls2 = -cls
#cl3 = -cl
#alpha3 = 180*np.pi/180 - alpha
#
#alpha_all = np.concatenate((alpha, alpha_s, np.flip(alpha_s2), np.flip(alpha3)))
#cl_all = np.concatenate((cl,cls,np.flip(cls2),np.flip(cl3)))
#cd_all = np.concatenate((cd,cds,np.flip(cds),np.flip(cd)))
#
#cn = cl_all*np.cos(alpha_all)+cd_all*np.sin(alpha_all)