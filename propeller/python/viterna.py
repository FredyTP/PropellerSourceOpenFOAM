import numpy as np;
import matplotlib as mpl;
import matplotlib.cm as cm;
import matplotlib.pyplot as plot;
import csv



AR=10
alphas=np.linspace(-180,180,360)

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
        cl[i]= np.abs(CDmax)*np.sin(alpha)*np.cos(alpha)
        cd[i]= CDmax*np.sin(alpha)**2
    else:
        #viterna negative
        alpha = -alpha
        alpha = alpha * np.pi/180
        cl[i]= -(A1_n*np.sin(2*alpha)+A2_n*(np.cos(alpha)**2)/np.sin(alpha))
        cd[i] = B1_n*np.sin(alpha)**2 + B2_n*np.cos(alpha)





alpha_lineal = alphas[180+alpha_stall_n-1:180+alpha_stall_p+1]
cl_lineal=cl[180+alpha_stall_n-1:180+alpha_stall_p+1]
cd_lineal=cd[180+alpha_stall_n-1:180+alpha_stall_p+1]

alpha_viterna=alphas[180-90:180+alpha_stall_n]
cl_viterna=cl[180-90:180+alpha_stall_n]
cd_viterna=cd[180-90:180+alpha_stall_n]

alpha_viterna_pos=alphas[180+alpha_stall_p:180+90]
cl_viterna_pos=cl[180+alpha_stall_p:180+90]
cd_viterna_pos=cd[180+alpha_stall_p:180+90]

alpha_plate_neg = alphas[0:91]
cl_plate_neg = cl[0:91]
cd_plate_neg = cd[0:91]

alpha_plate_pos = alphas[180+89:360]
cl_plate_pos = cl[180+89:360]
cd_plate_pos = cd[180+89:360]


colormap = mpl.colormaps['inferno']

plot.plot(alpha_lineal,cl_lineal,color=colormap(1/4),linewidth = 2)
plot.plot(alpha_viterna,cl_viterna,color=colormap(2/4),linewidth = 2)
plot.plot(alpha_plate_neg,cl_plate_neg,color=colormap(3/4),linewidth = 2)
plot.plot(alpha_viterna_pos,cl_viterna_pos,color=colormap(2/4),linewidth = 2)
plot.plot(alpha_plate_pos,cl_plate_pos,color=colormap(3/4),linewidth = 2)

plot.legend(["datos","viterna","placa plana"])
plot.title("Extrapolación Viterna + placa plana")
plot.xlabel(r'$\alpha$ [º]')
plot.ylabel(r'$cl$ [-]')
plot.grid()
plot.savefig("extrapolacion_viterna_flatplate.pdf", format="pdf", bbox_inches="tight")
#plot.show()


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