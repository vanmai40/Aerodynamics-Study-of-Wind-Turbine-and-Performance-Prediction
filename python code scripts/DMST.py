# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:10:04 2017

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
'###############################################################'


solidity=[];CPmax=[]
Ttaup=np.linspace(-pi/2, pi/2, 22).tolist()
Ttadw=np.linspace(pi/2, 3*pi/2, 22).tolist()
Ttadw.reverse();#Ttadw.pop(0);Ttadw.pop(-1)


#AR=10
#cla=2*pi/(1+(2*pi/(pi*AR)))
'###############################################################'
'##########function tính upwind factor ##########################'
def upfactor(tta,TSR,solidity):
    i=0 #tinh số lần lặp
    a1=0.1
    a2=0 # giá trị ban đầu
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
        
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        
        cti=cn*cos(tta)+ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(1*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        a2=gpt(3,-5,4,-d)
        i+=1
        if i>100:
            a2=cuctieuDMSTup(tta,0)
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
        
    dti=cti*w7uinf**2
    dpi=cpi*w7uinf**3
    '# tính force coefficient tại mỗi tta upwind ####'
    
    '# tính torque tại mỗi tta upwind ##  '
    
    AOA.append(aoa)
    CN.append(cn)
    CT.append(ct)
    Aup.append(a2)
    loop.append(i)
    Dti.append(dti)
    Dpi.append(dpi)        
    return 
'############################################'

'##########        DOWN  WIND      ########################'

'##########function tính downwind factor ##########################'

def dwfactor(tta,a,TSR,solidity):
    w7uw=0;aoa=0;cl=0;cd=0;cn=0;ct=0
    cti=0;cpi=0;d=0;i=0
    a1=0
    a2=a
    while abs(a2-a1)>0.0001:
        a1=a2
        TSRdw=TSR/(1-2*a)
        w7uw=sqrt((TSRdw-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
        aoa=atan(((1-a1)*cos(tta))/(TSRdw-(1-a1)*sin(tta)))
        uw=uinf*(1-2*a)
        re=w7uw*uw*c/v        
        
        cl=noisuyCL(re,degrees(abs(aoa)))
        if aoa>0: cl=cl[0]
        else: cl=-cl[0]
        cd=noisuyCD(re,degrees(abs(aoa)))
        cd=cd[0]
        
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        
        cti=cn*cos(tta)+ct*sin(tta)
        cpi=ct*TSRdw/w7uw
        
        d=(solidity/pi)*(cti*w7uw**2)/(abs(cos(tta)))
        a2=gpt(3,-5,4,-d)
        i+=1
        if i>100:
            a2=cuctieuDMSTdw(tta,0)
            a1=a2
            TSRdw=TSR/(1-2*a)
            w7uw=sqrt((TSRdw-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
            aoa=atan(((1-a1)*cos(tta))/(TSRdw-(1-a1)*sin(tta)))
            uw=uinf*(1-2*a)
            re=w7uw*uw*c/v        
            
            cl=noisuyCL(re,degrees(abs(aoa)))
            if aoa>0: cl=cl[0]
            else: cl=-cl[0]
            cd=noisuyCD(re,degrees(abs(aoa)))
            cd=cd[0]
            
            cn=cl*cos(aoa)+cd*sin(aoa)
            ct=cl*sin(aoa)-cd*cos(aoa)
            
            cti=cn*cos(tta)+ct*sin(tta)
            cpi=ct*TSRdw/w7uw
            break
        
    dti=cti*(w7uw*(1-2*a))**2
    dpi=cpi*(w7uw*(1-2*a))**3        
    
    '######  tính force tại mỗi tta dw ####'
    
    '# tính torque tại mỗi tta upwind ##'
    
    '# tích lũy kết quả các giá trị upwind'

    AOAdw.append(aoa)
    CNdw.append(cn)
    CTdw.append(ct)
    Adw.append(a2)
    loopdw.append(i)
    Dtidw.append(dti)
    Dpidw.append(dpi)    
    return 
'############################################'
#def opup(tta,a1,TSR):
#    w7uinf=sqrt((TSR-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
#    aoa=atan(((1-a1)*cos(tta))/(TSR-(1-a1)*sin(tta)))
#    
#    cl=vite_lift(aoa)   
##    cd=abs(K*cl)+((cl**2)/(pi*1*AR))
#    cd=vite_drag(aoa)
#    cn=cl*cos(aoa)+cd*sin(aoa)
#    ct=cl*sin(aoa)-cd*cos(aoa)
#    
#    cti=cn*cos(tta)+ct*sin(tta)
#    cpi=ct*TSR/w7uinf 
#    
#    dti=cti*w7uinf**2
#    dpi=cpi*w7uinf**3
#
#    
#    Dtiopup.append(dti)
#    Dpiopup.append(dpi)
#    
#    
#    return
#'#################################################################' 
'############################################'
#def opdw(tta,a,a1,TSR):
#    
#    TSRdw=TSR/(1-2*a)
#    w7uw=sqrt((TSRdw-(1-a1)*sin(tta))**2+((1-a1)*cos(tta))**2)
#    aoa=atan(((1-a1)*cos(tta))/(TSRdw-(1-a1)*sin(tta)))
#
#    cl=vite_lift(aoa)   
##    cd=abs(K*cl)+((cl**2)/(pi*1*AR))
#    cd=vite_drag(aoa)
#    cn=cl*cos(aoa)+cd*sin(aoa)
#    ct=cl*sin(aoa)-cd*cos(aoa)
#    
#    cti=cn*cos(tta)+ct*sin(tta)
#    cpi=ct*TSRdw/w7uw
#    
#    dti=cti*(w7uw*(1-2*a))**2
#    dpi=cpi*(w7uw*(1-2*a))**3 
#    
#    Dtiopdw.append(dti)
#    Dpiopdw.append(dpi)
#    
#    
#    return
'########################## SETTING UP GRAPHS ##########################' 
#output_file("DMST.html")
#p1=figure(title="POWER CURVES",
#         x_range=(0, 8),
#         y_range=(0, 0.6),
#         plot_width=1000) 
#p1.title.align="center"
#p1.title.text_font_size = "25px"
#p1.xaxis[0].axis_label = 'TSR'
#p1.yaxis[0].axis_label = 'CP'

#p2=figure(title="OPTIMUM SOLIDITY",
#         x_range=(0, 12),
#         y_range=(0, 1),
#         plot_width=1000) 
#p2.title.align="center"
#p2.title.text_font_size = "25px"
#p2.xaxis[0].axis_label = 'TSR'
#p2.yaxis[0].axis_label = 'SOLIDITY'

#colors=color.d3['Category20'][(len(soli_range))]

'######################## INPUT ################################'
v=0.00001331
uinf=13.5
c=1.84
beta=0
TSR_range=np.linspace(0, 8, 40).tolist()
soli_range=np.linspace(0.105, 0.8, 1).tolist();soli_range = [ round(k,2) for k in soli_range]
a0=0;
AR=100000
'########################################################'
'###### tính toán upwind ###########'
for soli in range(len(soli_range)):
    solidity=soli_range[soli]
    CPtt=[]
    for tsr in range(len(TSR_range)):
        TSR=TSR_range[tsr]
        AOA=[];CN=[];CT=[];Aup=[];loop=[];Dti=[];Dpi=[]
        for i in range(len(Ttaup)):
            upfactor(Ttaup[i],TSR,solidity)
            CF1=solidity*(sum(Dti)/len(Dti))
            CP1=solidity*(sum(Dpi)/len(Dpi))
        '###### tính toán downwind ###########'
        AOAdw=[];CNdw=[];CTdw=[];Adw=[];loopdw=[];Dtidw=[];Dpidw=[]
        for j in range(len(Ttadw)):
            tta=Ttadw[j]
            a=Aup[j]
            dwfactor(tta,a,TSR,solidity)
            CF2=solidity*(sum(Dtidw)/len(Dtidw))
            CP2=solidity*(sum(Dpidw)/len(Dpidw))
        CFn=(CF1+CF2)/2
        CPn=(CP1+CP2)/2
        CPtt.append(CPn)
#    p1.square(TSR_range,CPtt,legend='σ = '+str(soli_range[soli]), line_color=None,fill_color=colors[soli],alpha=3)
#    p1.line(TSR_range,CPtt,legend='σ = '+str(soli_range[soli]), line_color=colors[soli],alpha=3,line_width=1)
    
    point=CPtt.index(max(CPtt))
    label=Label(x=TSR_range[point],y=CPtt[point],
                text_font_size='15px',
                x_offset=-25,
                y_offset=-15,
                render_mode='canvas',
                #text=str((TSR_range[point],round(CP[point],2))))
                text='σ = '+str(soli_range[soli]))
#    p1.add_layout(label)
'############################### OPTIMUM ##############################'
#CPop=[];Soli_op=[];aup=0.25;adw=0.5
#for tsr in range(len(TSR_range)):
#        TSR=TSR_range[tsr]
#        Dpiopup=[];Dtiopup=[];Dpiopdw=[];Dtiopdw=[]
#        for l in range(len(Ttaup)):
#            opup(Ttaup[l],aup,TSR)
#        for l in range(len(Ttadw)):
#            opdw(Ttadw[l],aup,adw,TSR)
#            Clt=4*aup*(1-aup)+(1-2*aup)**2*4*adw*(1*adw)
#            Dtiop=Dtiopup+Dtiopdw
#            Dpiop=Dpiopup+Dpiopdw
#            soli_op=Clt/(sum(Dtiop)/len(Dtiop))
#            cp_op=soli_op*(sum(Dpiop)/len(Dpiop))
#        CPop.append(round(cp_op,5))
#        Soli_op.append(round(soli_op,5))
#p1.square(TSR_range,CPop,legend='Optimum Line',line_color=None,fill_color='red',alpha=3)
#p1.line(TSR_range,CPop,legend='Optimum Line',line_color='red',alpha=3,line_width=2)
#Soli_op.pop(0)
#p2.square(TSR_range,Soli_op,line_color=None,fill_color='red',alpha=3)
#p2.line(TSR_range,Soli_op,line_color='red',alpha=3,line_width=2) 
#show(p2)  
#show(p1) 
'############################### plots ##############################'    
Ttadw.reverse()
Tta=Ttaup+Ttadw
Tta=[degrees(i) for i in Tta]

Adw.reverse()
A=Aup+Adw


CNdw.reverse()
CNdw.pop(0)
CN=CN+CNdw
CN=[-i for i in CN]


CTdw.reverse()
CTdw[0]=CT[-1]
CT=CT+CTdw


AOAdw.reverse()
AOA=AOA+AOAdw
AOA=[degrees(i) for i in AOA]

#output_file("DMST.html")
#p1=figure(title="Power curve",
#         x_range=(0, 9),
#         y_range=(0, 0.5),
#         plot_width=1000) 
#p1.title.align="center"
#p1.title.text_font_size = "25px"
#
#p1.xaxis[0].axis_label = 'TSR '
#p1.yaxis[0].axis_label = 'Power coefficient'
#p1.xaxis[0].axis_label_text_font_size ="15px"
#p1.yaxis[0].axis_label_text_font_size ="15px"
#p1.xaxis[0].major_label_text_font_size ="15px"
#p1.yaxis[0].major_label_text_font_size ="15px"




p1.asterisk(TSR_range,CPtt,legend='D-Multiple Streamtube',line_color='orangered',fill_color=None,alpha=3,size=13)
p1.line(TSR_range,CPtt,legend='D-Multiple Streamtube', line_color='orangered',alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)

show(p1)
