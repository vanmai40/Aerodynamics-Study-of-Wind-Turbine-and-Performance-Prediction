# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 21:29:35 2017

@author: Mai Nguyen Van
"""

output_file("comparison.html")
p1=figure(#x_range=(1.8,4.2),
          #y_range=(-0.02, 0.56),
         plot_width=1000) 
p1.title.align="center"
p1.title.text_font_size = "25px"
p1.xaxis[0].axis_label = 'TSR'
p1.yaxis[0].axis_label = 'Power coefficient'

p1.xaxis[0].axis_label_text_font_size ="15px"
p1.yaxis[0].axis_label_text_font_size ="15px"
p1.xaxis[0].major_label_text_font_size ="15px"
p1.yaxis[0].major_label_text_font_size ="15px"


p3=figure(x_range=(0,8),
          y_range=(0, 1),
         plot_width=1000) 
p3.title.align="center"
p3.title.text_font_size = "25px"
p3.xaxis[0].axis_label = 'TSR'
p3.yaxis[0].axis_label = 'solidity'

p3.xaxis[0].axis_label_text_font_size ="15px"
p3.yaxis[0].axis_label_text_font_size ="15px"
p3.xaxis[0].major_label_text_font_size ="15px"
p3.yaxis[0].major_label_text_font_size ="15px"

wb=xw.Book(r'E:\Google Drive\AAA VAWT\Luận văn VAWT\H type materials\experimental data\experimental data.xlsx')
sht=wb.sheets['12kW']
#sht1=wb.sheets['predict']

#TSR1=sht.range('J8:J39').value
#CP1=sht.range('K8:K39').value

#TSR2=sht.range('I5:I15').value
#CP2=sht.range('J5:J15').value
#
#TSR3=sht.range('T9:T32').value
#CP3=sht.range('U9:U32').value
#
tta=sht.range('F94:F179').value
tta=[i -90 for i in tta]
a=sht.range('H94:H179').value
#a=[-i for i in a]

p1.diamond(tta,a,legend='Experimental data', line_color='blue',fill_color=None,alpha=3,size=15)
#p1.line(TSR_range,power,legend='SST', line_color='blue',alpha=3,line_width=1)
#p1.line(Tta,0, line_color='black',alpha=3,line_width=1)

#p1.legend.location='bottom_right'
show(p1)
#sht1.range('J27').options(transpose=True).value=power
