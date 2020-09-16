# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 21:17:50 2017

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

'###############################################################'

#colors=color.d3['Category20'][(len(soli_range))]
Ttaup=np.linspace(-pi/2, pi/2, 22).tolist()
Ttadw=np.linspace(pi/2, 3*pi/2, 22).tolist()
Ttadw.reverse();Ttadw.pop(0);Ttadw.pop(-1)
dtta=Ttaup[1]-Ttaup[0]





'###############################################################'
def run(tta,a1,TSR):
        w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
        aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))+radians(beta)
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
        
        cp=cpi*w7uinf**3
        cf=cti*w7uinf**2
    
        AOA.append(degrees(aoa))
        Cti.append(cti)
#        Cpi.append(cpi)
        Re.append(re)
        CN.append(cn)
        CT.append(ct)
        W7uinf.append(w7uinf)
        CP.append(cp)
        CF.append(cf)
        return

'############################################################' 
Tta=np.linspace(-pi/2, 3*pi/2, 50).tolist()
v=0.00001331
uinf=13.5
c=1.84
beta=0
TSR_range=np.linspace(0, 8, 40).tolist()
soli_range=np.linspace(0.105, 0.8, 1).tolist();soli_range = [ round(k,2) for k in soli_range]
solidity=soli_range[0]
'############################################################'  
#TSR=5
#
#a1=0.1 
#a2=0
#loop=0
#TSR=TSR_range[j]     
#while abs(a2-a1)>0.0001:
#    a1=a2
#    CF=[];CP=[];AOA=[];Cti=[];W7uinf=[];Re=[];CN=[];CT=[]
#    for i in range(len(Tta)):
#        run(Tta[i],a1,TSR)
#    CPn=solidity*sum(CP)/len(CP)
#    CFn=solidity*sum(CF)/len(CF)
#    a2=gpt(3,-5,4,-CFn)
#    loop+=1
#    if loop>100:
#        break
'############################ ITERATION #########################'
power=[];Loop=[];CFN=[]
for j in range(len(TSR_range)):    
    a1=0.1 
    a2=0
    loop=0
    TSR=TSR_range[j]     
    while abs(a2-a1)>0.0001:
        a1=a2
        CF=[];CP=[];AOA=[];Cti=[];W7uinf=[];Re=[];CN=[];CT=[]
        for i in range(len(Tta)):
            run(Tta[i],a1,TSR)
        CPn=solidity*sum(CP)/len(CP)
        CFn=solidity*sum(CF)/len(CF)
        a2=gpt(3,-5,4,-CFn)
        loop+=1
        if loop>100:
            a_range=np.linspace(0,1,200).tolist()
            Devi=[]
            for a1 in a_range:
                CF=[];CP=[]#;AOA=[];Cti=[];W7uinf=[];Re=[];CN=[];CT=[]
                for i in range(len(Tta)):
                    run(Tta[i],a1,TSR)    
                CFn=solidity*sum(CF)/len(CF)
                d1=3*a1**3-5*a1**2+4*a1
                devi=abs(d1-CFn)
                Devi.append(devi)
            point=Devi.index(min(Devi))
            a1=a_range[point]
            CF=[];CP=[];AOA=[];Cti=[];W7uinf=[];Re=[];CN=[];CT=[]
            for i in range(len(Tta)):
                run(Tta[i],a1,TSR)
            CPn=solidity*sum(CP)/len(CP)
            CFn=solidity*sum(CF)/len(CF)
            break
  
    power.append(CPn)
    CFN.append(CFn)    
    Loop.append(loop)
'############################ input a #########################'
#CF=[];CP=[];AOA=[];Cti=[];W7uinf=[];Re=[]
#a1=1/3
#for i in range(len(Tta)):
#    run(Tta[i],a1,TSR)
#CPn=solidity*sum(CP)/len(CP)
#CFn=solidity*sum(CF)/len(CF)
'############################ graphs parameter #########################'    
#output_file("SST.html")
#p1=figure(title="Power curve",
#          x_range=(0, 8),
##         y_range=(-90, 270),
#         plot_width=1000) 
#p1.title.align="center"
#p1.title.text_font_size = "25px"
#p1.xaxis[0].axis_label = 'TSR'
#p1.yaxis[0].axis_label = 'Power coefficient'
#
#p1.xaxis[0].axis_label_text_font_size ="15px"
#p1.yaxis[0].axis_label_text_font_size ="15px"
#p1.xaxis[0].major_label_text_font_size ="15px"
#p1.yaxis[0].major_label_text_font_size ="15px"



#Tta=[degrees(i) for i in Tta]
#CN=[-i for i in CN]
p1.square(TSR_range,power,legend='Single Streamtube', line_color='green',fill_color=None,alpha=3,size=9)
p1.line(TSR_range,power,legend='Single Streamtube', line_color='green',alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)
show(p1)

#p_trim=power[:19]+power[-4:]
#TSR_trim=TSR_range[:19]+TSR_range[-4:]
#z=np.polyfit(TSR_range,p_trim,20)
#f = np.poly1d(z)
#ptrim=[f(i) for i in TSR_range]
#p1.square(TSR_range,ptrim,legend=' σ = 0.1', line_color='green',fill_color=None,alpha=3,size=10)
#p1.line(TSR_range,ptrim,legend=' σ = 0.1', line_color='green',alpha=3,line_width=1)
#show(p1)