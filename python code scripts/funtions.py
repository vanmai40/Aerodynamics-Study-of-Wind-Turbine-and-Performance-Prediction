# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 16:10:02 2017

@author: Mai Nguyen Van
"""
import numpy as np
from bokeh.plotting import output_file, show, figure
from bokeh.models import Title,Label,Legend
import bokeh.palettes as color
import xlwings as xw
from scipy import interpolate
from math import sqrt,degrees,pi,cos,sin,asin,tan,atan,radians,acos,e
from numpy import roots
from sympy.solvers import solve
from sympy import Symbol, simplify, expand, factor

'#################### VITERNIA METHOD #############################'
AR=8
K=0.05
alfas=radians(10)
cla=1*pi
Cla=cla/(1+(cla/(pi*AR)))

Cla=cla/(1+(cla/(pi*AR)))
CLs=Cla*sin(alfas)
CDmax=1.11+0.018*AR
CDs=0.01


A1=CDmax/2
B1=CDmax
A2=(CLs-CDmax*sin(alfas)*cos(alfas))*sin(alfas)/(cos(alfas)**2)
B2=(CDs-CDmax*(sin(alfas)**2))/cos(alfas)

vite_cl=lambda a:A1*sin(2*a)+A2*cos(a)**2/sin(a)
vite_cd=lambda a:B1*sin(a)**2+B2*cos(a)

def vite_lift(aoa):
    if abs(aoa) <=alfas:
        return Cla*sin(aoa)
    elif abs(aoa)>alfas and abs(aoa)<=pi/2:
        return vite_cl(aoa)
    if abs(aoa)>pi/2:
        if aoa>0:return -vite_lift(pi-aoa)
        elif aoa <0: return -vite_lift(-pi-aoa)
        
def vite_drag(aoa):
    if abs(aoa) <=2*alfas:
        return CDs
    elif abs(aoa)>alfas and abs(aoa)<=pi/2:
        return vite_cd(aoa)
    if abs(aoa)>pi/2:
        if aoa>0:return vite_drag(pi-aoa)
        elif aoa <0: return vite_drag(-pi-aoa)
'######################### DRAG POLARS ######################'
wb1=xw.Book(r'E:\Google Drive\AAA VAWT\Luận văn VAWT\airfoil polars database.xlsx')
sht1=wb1.sheets['0021Vite']
#sht1=wb1.sheets['CLCD0021']


#Re=sht1.range('D3:N3').value
#alfa=sht1.range('C4:C62').value  #CLCD0015
#CL=sht1.range('D4:N62').value
#CD=sht1.range('Q4:AA62').value

#Re=sht1.range('D3:N3').value
#alfa=sht1.range('C4:C62').value  #CLCD0018
#CL=sht1.range('D4:N62').value
#CD=sht1.range('Q4:AA62').value

Re=sht1.range('H3:O3').value
alfa=sht1.range('G5:G125').value  #0021Vite
CL=sht1.range('H5:O125').value
CD=sht1.range('S5:Z125').value

#Re=sht1.range('H3:O3').value
#alfa=sht1.range('G5:G124').value  #0018Vite
#CL=sht1.range('H5:O124').value
#CD=sht1.range('S5:Z124').value
#
#Re=sht1.range('D3:N3').value
#alfa=sht1.range('C4:C62').value  #CLCD0012
#CL=sht1.range('D4:N62').value
#CD=sht1.range('Q4:AA62').value

#Re=sht1.range('D3:M3').value
#alfa=sht1.range('C4:C52').value  #CLCD0021
#CL=sht1.range('D4:M52').value
#CD=sht1.range('P4:Y52').value

noisuyCL=interpolate.interp2d(Re,alfa,CL)
noisuyCD=interpolate.interp2d(Re,alfa,CD)
'###########################################################'  
def gpt(a,b,c,d):
    r=roots([a,b,c,d])
    real_root = r[(abs(r - r.real)).argmin()]
    r=real_root.real
    r=float(r)
#    if d<-2 :
#        r=0.5
    if d<-1000 or d>1000:
        r=0  #cho trường hợp tta = 270 và 90  

    return r
'########################################################'
'########################################################################'  
def newton(a,b,c,d):
    r1=0.1
    r2=0
    i=0
    Fa=lambda r:a*r**3+b*r**2+c*r+d
    fa=lambda r:3*a*r**2+2*b*r+c
    if d<-4:return 0
    elif d<-2: return 0.5
    elif d>0: return 0
    while abs(r2-r1)>0.001:
        r1=r2
        r2=r1-(Fa(r1)/fa(r1))
        #print(r2)
        i+=1
        if i>100:
            break
    return r2
'########################################################################'
def bisec(a,b,c,d):
    i=0
    low=0
    high=1
    r1=0;r2=1
    while abs(r2-r1)>0.0001:
        r1=r2
        r2=(low+high)/2
        value=lambda x:a*x**3+b*x**2+c*x+d
        if value(r2)*value(low)<0:
            high=r2
        elif value(r2)*value(high)<0:
            low=r2
        else: 
            if abs(value(low))<abs(value(high)):
                r2=low;break
            else: r2=high;break
    return r2    
'#######################################################################'
def cuctieuSMST(tta,z):
    a_range=np.linspace(-0.5,1.5,200).tolist()
    Devi=[]
    for a1 in a_range:

        aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))+radians(beta)
        
        if (z*sin(abs(aoa)))==0  : F=1
        else:
            f=(N/2)*(1-z)/(z*sin(abs(aoa)))
            F=2/pi*acos(e**(-f))

        
        w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*F*cos(tta))**2)
        aoa=atan(((1-a1)*F*cos(tta))/(TSR-(1-a1)*sin(tta)))+radians(beta)
        
        re=w7uinf*uinf*c/v
        cl=noisuyCL(re,degrees(abs(aoa)))
        if aoa>0: cl=cl[0]
        else: cl=-cl[0]
        cl=cl/(1+a0/(pi*AR))
        cd=noisuyCD(re,degrees(abs(aoa)))
        cd=cd[0]+cl**2/(pi*AR)
        
        aoa=aoa-radians(beta)
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        cti=cn*cos(tta)-ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(2*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        d1=(3*a1**3-5*a1**2+4*a1)*F
        devi=abs(d-d1)
        Devi.append(devi)
    point=Devi.index(min(Devi))
    return a_range[point]
'#######################################################################'
def cuctieuDMSTup(tta,z):
    a_range=np.linspace(-0.5,1.5,200).tolist()
    Devi=[]
    for a1 in a_range:
        aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))+radians(beta)
        
        if (z*sin(abs(aoa)))==0  : F=1
        else:
            f=(N/2)*(1-z)/(z*sin(abs(aoa)))
            F=2/pi*acos(e**(-f))

        
        w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*F*cos(tta))**2)
        aoa=atan(((1-a1)*F*cos(tta))/(TSR-(1-a1)*sin(tta)))+radians(beta)
        
        re=w7uinf*uinf*c/v
        cl=noisuyCL(re,degrees(abs(aoa)))
        if aoa>0: cl=cl[0]
        else: cl=-cl[0]
        cl=cl/(1+a0/(pi*AR))
        cd=noisuyCD(re,degrees(abs(aoa)))
        cd=cd[0]+cl**2/(pi*AR)
        
        aoa=aoa-radians(beta)
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        cti=cn*cos(tta)-ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(1*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        d1=(3*a1**3-5*a1**2+4*a1)*F
        devi=abs(d-d1)
        Devi.append(devi)
    point=Devi.index(min(Devi))
    return a_range[point]  
'#######################################################################'
def cuctieuDMSTdw(tta,z):
    a_range=np.linspace(-0.5,1.5,200).tolist()
    Devi=[]
    for a1 in a_range:
        TSRdw=TSR/(1-2*a)
        aoa=atan(((1-a1)*cos(tta))/(TSRdw-(1-a1)*sin(tta)))+radians(beta)
        
        if (z*sin(abs(aoa)))==0  : F=1
        else:
            f=(N/2)*(1-z)/(z*sin(abs(aoa)))
            F=2/pi*acos(e**(-f))  
            
        w7uw=sqrt((TSRdw-(1-a1)*sin(tta))**2+((1-a1)*F*cos(tta))**2)
        aoa=atan(((1-a1)*F*cos(tta))/(TSRdw-(1-a1)*sin(tta)))+radians(beta)
        
        uw=uinf*(1-2*a)
        re=w7uw*uw*c/v        
        
        cl=noisuyCL(re,degrees(abs(aoa)))
        if aoa>0: cl=cl[0]
        else: cl=-cl[0]
        cl=cl/(1+a0/(pi*AR))
        cd=noisuyCD(re,degrees(abs(aoa)))
        cd=cd[0]+cl**2/(pi*AR)
        
        aoa=aoa-radians(beta)
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        
        cti=cn*cos(tta)+ct*sin(tta)
        cpi=ct*TSRdw/w7uw
        
        d=(solidity/pi)*(cti*w7uw**2)/(abs(cos(tta)))
        
        d1=(3*a1**3-5*a1**2+4*a1)*F
        devi=abs(d-d1)
        Devi.append(devi)
    point=Devi.index(min(Devi))
    return a_range[point]
'#######################################################################'
def addlist(a,b):
    if len(a)==0 :return b
    c=[];
    for ind in range(len(a)):
        c.append(a[ind]+b[ind])
    return c
'#######################################################################'
'#######################################################################'
def devlist(a,b):
    c=[]
    for ind in range(len(a)):
        c.append(a[ind]/b)
    return c
'#######################################################################'
def atube(a,b):
    c=[];
    for j in range(len(a)):
        ci=0.5-((1-2*a[j])*(1-2*b[j])/2)
        c.append(ci)
    return c
'#######################################################################'
def ghep(a,b):
    b[0]=a[0]
    b[-1]=a[-1]
    b.reverse()
    c=a+b

    return c
        