# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:09:14 2017

@author: Mai Nguyen Van
"""

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

#TSR_range.pop(0)

#colors=color.d3['Category20'][(len(soli_range))]
#solidity=[];CPmax=[]
Ttaup=np.linspace(-pi/2, pi/2, 22).tolist()
Ttadw=np.linspace(pi/2, 3*pi/2, 22).tolist()
#Ttadw.reverse();Ttadw.pop(0);Ttadw.pop(-1)


'###############################################################'
def runup(tta,TSR,solidity,z):
    a1=0.1
    a2=0
    loop=0     
    while abs(a2-a1)>0.0001:
        a1=a2
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
        
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        cti=cn*cos(tta)-ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(2*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        a2=gpt(3*F,-5*F,4*F,-d)
        loop+=1
        if loop>100:
            a2=cuctieuSMST(tta,z)
            a1=a2
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
beta=0
uinf=13.5
c=1.84
D=35
AR=13.2
N=2
soli_range=np.linspace(N*c/D, 0.4, 1).tolist();soli_range = [ round(k,2) for k in soli_range]
TSR_range=np.linspace(0, 8, 20).tolist()
z_range=np.linspace(0, 0.9, 10).tolist()
t7c=0.18
a0=1.8*pi*(1+0.8*t7c)
'########################################################'
for soli in range(len(soli_range)):
    solidity=soli_range[soli]
    CPn=[]
    for tsr in range(len(TSR_range)):
        TSR=TSR_range[tsr]
        print('TSR =',TSR)
        CP_tl=[];CN_tl=[]
        for z in z_range:
            CN=[];CT=[];Loop=[];A=[];CP=[];CF=[];AOA=[];D=[];Re=[]
            for k in range(len(Ttaup)):
                runup(Ttaup[k],TSR,solidity,z)
                #Cp+=solidity*W7uinf[k]**3*Cpi[k]/len(Ttaup)
                Cf=solidity*(sum(CF)/len(CF))
                Cp=solidity*(sum(CP)/len(CP))
            CP_tl.append(Cp)
            CN_tl=addlist(CN_tl,CN)
        
        Cp=sum(CP_tl)/len(CP_tl)
        CN_tl=devlist(CN_tl,len(z_range))
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
Ttaup=[degrees(i) for i in Ttaup]

#Adw=A[:]
#Adw.reverse()
#
#A=A+Adw

CNdw=CN_tl[:]
CNdw.reverse()
CNdw=[-i for i in CNdw]
CN=CN_tl+CNdw
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
#         x_range=(-90, 270),
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

colors=['blue','green','darkmagenta','darkcyan','olive','darkred','lime','cyan','magenta','black','brown','gold']
mau=colors[2]

p1.triangle(TSR_range,CPn,legend='Multiple Streamtube',line_color=mau,fill_color=None,alpha=3,size=13)
p1.line(TSR_range,CPn,legend='Multiple Streamtube', line_color=mau,alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)
#p3.legend.location='top_right'
show(p1)

