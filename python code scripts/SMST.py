# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:06:34 2017

@author: Mai Nguyen Van
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 22:53:36 2017

@author: Mai Nguyen Van
"""

#import numpy as np
#from bokeh.plotting import output_file, show, figure
#from bokeh.models import Title,Label,Legend
#import bokeh.palettes as color
#import xlwings as xw
#from scipy import interpolate
#from math import sqrt,degrees,pi,cos,sin,asin,tan,atan,radians
#from numpy import roots
#from sympy.solvers import solve
#from sympy import Symbol, simplify, expand, factor
'######################## GRAPSHS SETUP ############################'
p2=figure(title="OPTIMUM SOLIDITY",
         x_range=(0, 12),
         y_range=(0, 1),
         plot_width=1000) 
p2.title.align="center"
p2.title.text_font_size = "25px"
p2.xaxis[0].axis_label = 'TSR'
p2.yaxis[0].axis_label = 'SOLIDITY'



'###############################################################'
'######################### INNITALIZING THETA ##########################'

#colors=color.d3['Category20'][(len(soli_range))]
#solidity=[];CPmax=[]
Ttaup=np.linspace(-pi/2, pi/2, 22).tolist()
Ttadw=np.linspace(pi/2, 3*pi/2, 22).tolist()
#Ttadw.reverse();Ttadw.pop(0);Ttadw.pop(-1)

#AR=10
#cla=2*pi/(1+(2*pi/(pi*AR)))
'###############################################################'
def runup(tta,TSR,solidity):
    a1=0.1
    a2=0
    loop=0     
    while abs(a2-a1)>0.0001:
        a1=a2
        w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
        aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))
        
        re=w7uinf*uinf*c/v        
        
        cl=noisuyCL(re,degrees(abs(aoa)))
        if aoa>0: cl=cl[0]
        else: cl=-cl[0]
        cd=noisuyCD(re,degrees(abs(aoa)))
        cd=cd[0]
        
#        cl=Cla*sin(aoa)
#        cd=abs(K*cl)+((cl**2)/(pi*1*AR))
        
#        cl=vite_lift(aoa) 
#        cd=vite_drag(aoa)
        
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        cti=cn*cos(tta)-ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(2*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        a2=gpt(3,-5,4,-d)
        loop+=1
        if loop>100:
            a2=cuctieuSMST(tta,0)
            a1=a2
            w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
            aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))
            
            re=w7uinf*uinf*c/v        
            
            cl=noisuyCL(re,degrees(abs(aoa)))
            if aoa>0: cl=cl[0]
            else: cl=-cl[0]
            cd=noisuyCD(re,degrees(abs(aoa)))
            cd=cd[0]
            cn=cl*cos(aoa)+cd*sin(aoa)
            ct=cl*sin(aoa)-cd*cos(aoa)
            cti=cn*cos(tta)-ct*sin(tta)
            cpi=ct*TSR/w7uinf
            break
        
        
        
    cp=cpi*w7uinf**3
    cf=cti*w7uinf**2  
    
    AOA.append(degrees(aoa))
    CF.append(cf)
    CT.append(ct)
    CN.append(cn)
    Loop.append(loop)
    A.append(a2)
    CP.append(cp)
    D.append(d)
    Re.append(re)
    
    return
'###############################################################'
#def optimum(tta,a1,TSR):
#    w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
#    aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))
#    
##    re=160000
##    cl=noisuyCL(re,abs(degrees(aoa)))
##    if aoa>0:cl=cl[0]
##    else:cl=cl[0]*-1
##    cd=noisuyCD(re,abs(degrees(aoa)))
##    cd=cd[0]
#
##    cl=Cla*sin(aoa)
##    cd=abs(K*cl)+((cl**2)/(pi*1*AR))
#    cl=vite_lift(aoa)   
##    cd=abs(K*cl)+((cl**2)/(pi*1*AR))
#    cd=vite_drag(aoa)
#    cn=cl*cos(aoa)+cd*sin(aoa)
#    ct=cl*sin(aoa)-cd*cos(aoa)
#    
#    cti=cn*cos(tta)+ct*sin(tta)
#    cpi=ct*TSR/w7uinf 
#    dti=cti*w7uinf**2
#    dpi=cpi*w7uinf**3
#    
#    Dti.append(dti)
#    Dpi.append(dpi)
#    
#    
#    return

'######################## SMST ################################'
v=0.00001331
uinf=13.5
c=1.84
beta=0
TSR_range=np.linspace(0, 8, 40).tolist()
soli_range=np.linspace(0.105, 0.8, 1).tolist();soli_range = [ round(k,2) for k in soli_range]
a0=0;
AR=100000
'########################################################'
for soli in range(len(soli_range)):
    solidity=soli_range[soli]
    CPn=[]
    for tsr in range(len(TSR_range)):
        TSR=TSR_range[tsr]
        CN=[];CT=[];Loop=[];A=[];CP=[];CF=[];AOA=[];D=[];Re=[]
        for k in range(len(Ttaup)):
            runup(Ttaup[k],TSR,solidity)
            #Cp+=solidity*W7uinf[k]**3*Cpi[k]/len(Ttaup)
            Cf=solidity*(sum(CF)/len(CF))
            Cp=solidity*(sum(CP)/len(CP))
        CPn.append(round(Cp,5))
#        if CPn[-1]> 0.6:CPn.pop(-1);break
#        if TSR_range[tsr]>2:
#            if  CPn[tsr]>CPn[tsr-1] and CPn[tsr-1]<CPn[tsr-2]:
#                CPn.pop(-1);break
#    p1.square(TSR_range,CPn,legend='σ = '+str(soli_range[soli]), line_color=None,fill_color='DarkMagenta',alpha=3,size=10)
#    p1.line(TSR_range,CPn,legend='σ = '+str(soli_range[soli]), line_color='DarkMagenta',alpha=3,line_width=1)
    point=CPn.index(max(CPn))
    label=Label(x=TSR_range[point],y=CPn[point],
                text_font_size='15px',
                x_offset=-25,
                y_offset=-15,
                render_mode='canvas',
                #text=str((TSR_range[point],round(CP[point],2))))
                text='σ = '+str(soli_range[soli]))
#    p1.add_layout(label)
'############################### OPTIMUM ##############################'
#CPop=[];Soli_op=[];aop=1/3
#for tsr in range(len(TSR_range)):
#        TSR=TSR_range[tsr]
#        Dpi=[];Dti=[]
#        for l in range(len(Ttaup)):
#            optimum(Ttaup[l],aop,TSR)
#            Clt=4*aop*(1-aop)
#            soli_op=Clt/(sum(Dti)/len(Dti))
#            cp_op=soli_op*(sum(Dpi)/len(Dpi))
#        CPop.append(round(cp_op,5))
#        Soli_op.append(round(soli_op,5))
#p1.square(TSR_range,CPop,legend='Optimum Line',line_color=None,fill_color='red',alpha=3)
#p1.line(TSR_range,CPop,legend='Optimum Line',line_color='red',alpha=3,line_width=2)
#Soli_op.pop(0)
#p2.square(TSR_range,Soli_op,line_color=None,fill_color='red',alpha=3)
#p2.line(TSR_range,Soli_op,line_color='red',alpha=3,line_width=2) 
##show(p2)  
#show(p1) 
'############################### plots ##############################'    
Tta=Ttaup+Ttadw

Tta=[degrees(i) for i in Tta]

Adw=A[:]
Adw.reverse()

A=A+Adw

CNdw=CN[:]
CNdw.reverse()
CNdw=[-i for i in CNdw]
CN=CN+CNdw
CN=[-i for i in CN]

CTdw=CT[:]
CTdw.reverse()

CT=CT+CTdw

AOAdw=AOA[:]
AOAdw.reverse()
AOAdw=[-i for i in AOAdw]
AOA=AOA+AOAdw


#output_file("SMST.html")
#p1=figure(title="CN VS θ",
#         x_range=(0, 8),
##         y_range=(0, 0.5),
#         plot_width=1000) 
#p1.title.align="center"
#p1.title.text_font_size = "25px"
#
#p1.xaxis[0].axis_label = 'Azimuthal angle, θ'
#p1.yaxis[0].axis_label = 'Normal force coefficient'
#p1.xaxis[0].axis_label_text_font_size ="15px"
#p1.yaxis[0].axis_label_text_font_size ="15px"
#p1.xaxis[0].major_label_text_font_size ="15px"
#p1.yaxis[0].major_label_text_font_size ="15px"



p1.triangle(TSR_range,CPn,legend='Multiple Streamtube',line_color='DarkMagenta',fill_color=None,alpha=3,size=13)
p1.line(TSR_range,CPn,legend='Multiple Streamtube', line_color='DarkMagenta',alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)

show(p1)
