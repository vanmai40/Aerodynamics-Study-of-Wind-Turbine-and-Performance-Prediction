# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 14:32:30 2017

@author: Mai Nguyen Van
"""

'######################### DRAG POLARS ######################'

'###############################################################'
'###############################################################'

#TSR_range.pop(0)


solidity=[];CPmax=[]
Ttaup=np.linspace(-pi/2, pi/2, 30).tolist()

Ttadw=np.linspace(pi/2, 3*pi/2, 30).tolist()
Ttadw.reverse();#Ttadw.pop(0);Ttadw.pop(-1)


'###############################################################'
'##########function tính upwind factor ##########################'
def upfactor(tta,TSR,solidity):
    i=0 #tinh số lần lặp
    a1=0.1
    a2=0 # giá trị ban đầu
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
        
        cti=cn*cos(tta)+ct*sin(tta)
        cpi=ct*TSR/w7uinf
        
        d=(1*solidity/(1*pi))*w7uinf**2*(cti/abs(cos(tta)))
        a2=gpt(3*F,-5*F,4*F,-d)
        i+=1
        if i>100:
            a2=cuctieuDMSTup(tta,z)
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
    CD.append(cd)
    W7uinf.append(w7uinf)        
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
        aoa=atan(((1-a1)*cos(tta))/(TSRdw-(1-a1)*sin(tta)))+radians(beta)
        
        if (z*sin(abs(aoa)))==0: F=1
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
        
        cn=cl*cos(aoa)+cd*sin(aoa)
        ct=cl*sin(aoa)-cd*cos(aoa)
        
        cti=cn*cos(tta)+ct*sin(tta)
        cpi=ct*TSRdw/w7uw
        
        d=(solidity/pi)*(cti*w7uw**2)/(abs(cos(tta)))
        a2=gpt(3*F,-5*F,4*F,-d)
        i+=1
        if i>100:
            a2=cuctieuDMSTdw(tta,z)
            a1=a2
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
    CDdw.append(cd)
    W7uw.append(w7uw)           
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
Tta=np.linspace(-pi/2, 3*pi/2, 50).tolist()

v=0.00001331
beta=-2
uinf=4
c=0.13
D=2
AR=2.5/0.13
N=3
soli_range=np.linspace(N*c/D, 0.4, 1).tolist();soli_range = [ round(k,2) for k in soli_range]
TSR_range=np.linspace(2, 4, 15).tolist()
z_range=np.linspace(0, 0.999, 10).tolist()
t7c=0.21
a0=1.8*pi*(1+0.8*t7c)
'########################################################'
'###### tính toán upwind ###########'
for soli in range(len(soli_range)):
    solidity=soli_range[soli]
    CPtt=[]
    for tsr in range(len(TSR_range)):
        TSR=TSR_range[tsr]
        print('TSR =',TSR)
        CP_tl=[];CT_tl=[];CTdw_tl=[];
        CN_tl=[];CNdw_tl=[]
        
        CD_tl=[];W7uinf_tl=[];Aup_tl=[]
        CDdw_tl=[];W7uw_tl=[];Adw_tl=[]
        
        
        for z in z_range:
            AOA=[];CN=[];CT=[];Aup=[];loop=[];Dti=[];Dpi=[];
            W7uinf=[];CD=[];
            for i in range(len(Ttaup)):
                upfactor(Ttaup[i],TSR,solidity)
            CF1=solidity*(sum(Dti)/len(Dti))
            CP1=solidity*(sum(Dpi)/len(Dpi))
            
            CT_tl=addlist(CT_tl,CT)
            CN_tl=addlist(CN_tl,CN)
            W7uinf_tl=addlist(W7uinf_tl,W7uinf)
            CD_tl=addlist(CD_tl,CD)
            Aup_tl=addlist(Aup_tl,Aup)
            
            
            
            '###### tính toán downwind ###########'
            AOAdw=[];CNdw=[];CTdw=[];Adw=[];loopdw=[];Dtidw=[];Dpidw=[]
            W7uw=[];CDdw=[];
            for j in range(len(Ttadw)):
                tta=Ttadw[j]
                a=Aup[j]
                dwfactor(tta,a,TSR,solidity)
            CF2=solidity*(sum(Dtidw)/len(Dtidw))
            CP2=solidity*(sum(Dpidw)/len(Dpidw))
            
            CFn=(CF1+CF2)/2
            CPn=(CP1+CP2)/2
            
            CP_tl.append(CPn)
            CTdw_tl=addlist(CTdw_tl,CTdw)
            CNdw_tl=addlist(CNdw_tl,CNdw)
            W7uw_tl=addlist(W7uw_tl,W7uw)
            CDdw_tl=addlist(CDdw_tl,CDdw)
            Adw_tl=addlist(Adw_tl,Adw)
            
        CT_tl=devlist(CT_tl,len(z_range))
        CTdw_tl=devlist(CTdw_tl,len(z_range))
        CN_tl=devlist(CN_tl,len(z_range))
        CNdw_tl=devlist(CNdw_tl,len(z_range))
        
        W7uinf_tl=devlist(W7uinf_tl,len(z_range))
        CD_tl=devlist(CD_tl,len(z_range))
        Aup_tl=devlist(Aup_tl,len(z_range))
        
        W7uw_tl=devlist(W7uw_tl,len(z_range))
        CDdw_tl=devlist(CDdw_tl,len(z_range))
        Adw_tl=devlist(Adw_tl,len(z_range)) 
        
        
        
        CPn=sum(CP_tl)/len(CP_tl)
        CPtt.append(CPn)
#    p1.square(TSR_range,CPtt,legend='σ = '+str(soli_range[soli]), line_color=None,fill_color=colors[soli],alpha=3)
#    p1.line(TSR_range,CPtt,legend='σ = '+str(soli_range[soli]), line_color=colors[soli],alpha=3,line_width=1)
    
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
#A1=[];A2=[]
#k=0.5*1.225*(D/2)*H/len(Ttaup)
#for cs in range(len(Ttaup)):
#    a=Aup_tl[cs]
#    a1=Adw_tl[cs]
#    tta=Ttaup[cs]
#    cd=CD_tl[cs]
#    cd1=CDdw_tl[cs]
#    w7uinf=W7uinf[cs]
#    w7uw=W7uw[cs]
#    a1=k*pi*uinf**3*abs(cos(tta))*(1-a)*(1-(1-2*a)**2*(1-2*a1)**2)
#    a2=k*solidity*((w7uinf*uinf)**3*cd+(w7uw*uinf*(1-2*a))**3*cd1)
#    A1.append(a1)
#    A2.append(a2)
#
#KE=sum(A1)
#dp=sum(A2)
#C_KE=KE/(0.5*1.225*uinf**3*D*H)
#C_dp=dp/(0.5*1.225*uinf**3*D*H)

'############################### plots ##############################'    
Ttadw.reverse()
Tta=Ttaup+Ttadw
Tta=[degrees(i) for i in Tta]
Ttaup=[degrees(i) for i in Ttaup]

Atube=atube(Aup,Adw)
#Adw.reverse()
#A=Aup+Adw
A=ghep(Aup,Adw)


CN=ghep(CN_tl,CNdw_tl)
CN=[-i for i in CN]




CT=ghep(CT_tl,CTdw_tl)



AOA=ghep(AOA,AOAdw)
AOA=[degrees(i) for i in AOA]

#output_file("DMST.html")
#p1=figure(title="Power curve",
#         x_range=(-90, 270),
##         y_range=(0, 0.5),
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

colors=['blue','green','darkmagenta','orangered','darkcyan','olive','darkred','lime','cyan','magenta','black','brown','gold']
mau=colors[int(beta)]

#p1.asterisk(Tta,AOA,legend='D-Multiple Streamtube',line_color=mau,fill_color=None,alpha=3,size=13)
#p1.line(Tta,AOA,legend='D-Multiple Streamtube', line_color=mau,alpha=3,line_width=1)
#p1.asterisk(Ttaup,Atube,legend='At each tube',line_color='darkorange',fill_color=None,alpha=3,size=13)
#p1.line(Ttaup,Atube,legend='At each tube', line_color='darkorange',alpha=3,line_width=1)
#p3.asterisk(Ttaup,Atube,legend='D-Multiple Streamtube',line_color='darkorange',fill_color=None,alpha=3,size=13)
#p3.line(Ttaup,Atube,legend='D-Multiple Streamtube', line_color='darkorange',alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)
p1.asterisk(TSR_range,CPtt,legend='β='+str(beta),line_color=mau,fill_color=None,alpha=3,size=13)
p1.line(TSR_range,CPtt,legend='β='+str(beta), line_color=mau,alpha=3,line_width=1)
#p1.legend.location='top_right'
show(p1)
#print('CPtt=',CPtt[0])

