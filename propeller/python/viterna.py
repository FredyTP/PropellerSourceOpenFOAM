import numpy as np;
import matplotlib as mpl;
import matplotlib.cm as cm;
import matplotlib.pyplot as plot;
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm
import csv



factor=3
alphas=np.linspace(-180,180,int(3600/factor))

AR=10
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



lim_inf = 1800+alpha_stall_n*10
lim_sup=1800+alpha_stall_p*10
lim_inf=int(lim_inf/factor)
lim_sup=int(lim_sup/factor)
alpha_lineal = alphas[lim_inf:lim_sup]
cl_lineal=cl[lim_inf:lim_sup]
cd_lineal=cd[lim_inf:lim_sup]

lim_inf=1800-900
lim_sup=1800+alpha_stall_n*10
lim_inf=int(lim_inf/factor)
lim_sup=int(lim_sup/factor)
alpha_viterna=alphas[lim_inf:lim_sup]
cl_viterna=cl[lim_inf:lim_sup]
cd_viterna=cd[lim_inf:lim_sup]

lim_inf=1800+alpha_stall_p*10
lim_sup=1800+900
lim_inf=int(lim_inf/factor)
lim_sup=int(lim_sup/factor)
alpha_viterna_pos=alphas[lim_inf:lim_sup]
cl_viterna_pos=cl[lim_inf:lim_sup]
cd_viterna_pos=cd[lim_inf:lim_sup]

lim_inf=0
lim_sup=900
lim_inf=int(lim_inf/factor)
lim_sup=int(lim_sup/factor)
alpha_plate_neg = alphas[lim_inf:lim_sup]
cl_plate_neg = cl[lim_inf:lim_sup]
cd_plate_neg = cd[lim_inf:lim_sup]

lim_inf=1800+900
lim_sup=3600
lim_inf=int(lim_inf/factor)
lim_sup=int(lim_sup/factor)
alpha_plate_pos = alphas[lim_inf:lim_sup]
cl_plate_pos = cl[lim_inf:lim_sup]
cd_plate_pos = cd[lim_inf:lim_sup]


colormap = mpl.colormaps['inferno']
mpl.rc('axes', labelsize=12)



plot.plot(alpha_lineal,cl_lineal,color=colormap(1/4),linewidth = 3)
plot.plot(alpha_viterna,cl_viterna,color=colormap(2/4),linewidth = 3)
plot.plot(alpha_plate_neg,cl_plate_neg,color=colormap(3/4),linewidth = 3)
plot.plot(alpha_viterna_pos,cl_viterna_pos,color=colormap(2/4),linewidth = 3)
plot.plot(alpha_plate_pos,cl_plate_pos,color=colormap(3/4),linewidth = 3)

plot.legend(["datos","viterna","placa plana"])
plot.title(r"Extrapolación $c_L$: Viterna y placa plana")
plot.xlabel(r'$\alpha$ ['+u"\u00b0" +"]")
plot.ylabel(r'$c_L$ [-]')
plot.xticks(np.linspace(-180,180,9))
plot.grid()
plot.savefig("cl_extrapolacion_viterna_flatplate.pdf", format="pdf", bbox_inches="tight")
plot.cla()

plot.plot(alpha_lineal,cd_lineal,color=colormap(1/4),linewidth = 3)
plot.plot(alpha_viterna,cd_viterna,color=colormap(2/4),linewidth = 3)
plot.plot(alpha_plate_neg,cd_plate_neg,color=colormap(3/4),linewidth = 3)
plot.plot(alpha_viterna_pos,cd_viterna_pos,color=colormap(2/4),linewidth = 3)
plot.plot(alpha_plate_pos,cd_plate_pos,color=colormap(3/4),linewidth = 3)

plot.legend(["datos","viterna","placa plana"])
plot.title(r"Extrapolación $c_D$: Viterna y placa plana")
plot.xlabel(r'$\alpha$ ['+u"\u00b0" +"]")
plot.ylabel(r'$c_D$ [-]')
plot.xticks(np.linspace(-180,180,9))

plot.grid()
plot.savefig("cd_extrapolacion_viterna_flatplate.pdf", format="pdf", bbox_inches="tight")
plot.show()
