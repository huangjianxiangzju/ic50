# -*- coding: utf-8 -*-
"""
huangjianxiang. 2017/12/4 20:34
This should be the most powerful and acurate version so far.
No error is not guaranteed. Use this script critically...
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import newton
import os
import math  
from timeit import default_timer as timer
import textwrap

global N; N=10  #the number of data to be fitted
#Change the number if have different values. 
#Most of the time the default 10 is large enough. 
text="The script is designed to fit the IC50 values for MTT assays\
 and work in batch mode. The max concentration data it can handle\
 is NINE! Use one csv file (data are seperated by comma) for one\
 MTT aasay. So you need the same number of csv files with your MTT\
 assays.Feel free to change the code.Enjoy! "
 
# print textwrap.fill(text,width=80) 
# print ""
#Read MTT data from csv files and plot the graph and output the result files
def readfromfile(csvfile,name):
    read= np.genfromtxt(csvfile,delimiter=',')
    a = open(csvfile, "r")
    yer=[]
	#读取行的数目
    lines=len(a.readlines())
	#假如numpy数组的列数为4，实际上重复试验为3，以此类推
    if read.shape[1]==4:
        ydata=(read[1:N,1]+read[1:N,2]+read[1:N,3])/3.0
        for i in np.arange(1,lines):
            yer.append(read[i,1:4].std())
    elif read.shape[1]==5:
        ydata=(read[1:N,1]+read[1:N,2]+read[1:N,3]+read[1:N,4])/4.0
        for i in np.arange(1,lines):
            yer.append(read[i,1:5].std())
    elif read.shape[1]==6:
        ydata=(read[1:N,1]+read[1:N,2]+read[1:N,3]+read[1:N,4]+read[1:N,5])/5.0
        for i in np.arange(1,lines):
            yer.append(read[i,1:6].std())
    else:
        pass
        # print "ERROR! The number of repeated assays should be between 3 to 5"
#error bar 数组
    yerr=np.array(yer)

    xdata=np.log10(read[1:N,0])
    A=ydata.max()
    D=ydata.min()
#Comment:The function sigmoid is defined to fit the data, 
#the origin program use 10 as base, so is here. change 10 
#to e if you want. Both will work fine. If you want to 
#find out the difference yourself. replace the line 
#"y = (A-D) / (1 + 10**(k*(x-x0))) + D" with
#"y = (A-D) / (1 + np.exp(k*(x-x0))) + D"
    def sigmoid(x, x0, k):
        y = (A-D) / (1 + 10**(k*(x-x0))) + D
        return y
    # sigmoid.func_closure
    popt, pcov = curve_fit(sigmoid, xdata, ydata)
    cmin=math.floor(np.log10(read[1:N,0].min()))
    cmax=math.ceil(np.log10(read[1:N,0].max()))
    x = np.linspace(cmin,cmax, 1000)
    y = sigmoid(x, *popt)
#输出拟合文件  
##############################################################
    xout = np.linspace(cmin,cmax, 50)
    yout = sigmoid(xout, *popt)
    for i in range(50):
        outfile=open("fit_data_"+name+".dat","a")
        outfile.write('{0:<18}'.format(str(10**xout[i]))+str(yout[i])+"\n")
    outfile.close()
###############################################################  
    #Optimal values for the parameters so that the sum of the squared
    #residuals of f(xdata, *popt) - ydata is minimized
    plt.plot(xdata, ydata, 'o', color="red", label='data',markersize=3)
    plt.errorbar(xdata, ydata,yerr=yerr,linestyle='none')
    plt.plot(x,y,color="darkviolet",linewidth=2, label='fit' )
########################################################    
# 修改X轴的label,修改这里,替换里面的QN即可             #
    plt.xlabel(r'Conc.($\mu M)$',fontsize=13)          #
########################################################    
    plt.ylabel("Cell viability(%)",fontsize=13)
    plt.ylim(-5, 110)
    ConcRange=['0.0001','0.001', '0.01', '0.1', '1', '10','100','1000','10000']
    slice1=int(cmin+4)
    slice2=int(cmax+5)
    plt.xticks(np.arange(cmin, cmax+1), ConcRange[slice1:slice2])
    plt.yticks([0,10,20,30,40,50,60,70,80,90,100,110])   
    n=100#default value
    for i in np.arange(cmin,cmax,0.0001):
        if abs(sigmoid(i, *popt)-50)<0.01:
            n=i
    string="The IC50 for "+name+" is " + str('{:.2f}'.format(10**(n)))    
	#在图上标上注释
    plt.annotate(string, xy=(n,52),fontsize=8,color='fuchsia')
    #输出结果文件output.txt
    dfile=open("output.txt","a")
    dfile.write(string+"\n")
    dfile.close()
    #plot the blue horizontal line
    x1 = np.linspace(cmin,cmax, 10)
    y1=np.linspace(50,50,10)
    plt.plot(x1,y1,'--',c='b',linewidth=1.5)
    #plot the vertical line
    y2=np.linspace(0, 100, 10)
    x2=np.linspace(n,n,10)
    plt.plot(x2,y2,'--',c='g',linewidth=1.5)
    plt.legend(loc='best')
    plt.title(name)
    plt.savefig(name+".png",dpi=300)
    plt.clf()
    # readfromfile.func_closure
#find all the csv files under the curent directory for plot and fit IC50 values
#remove all "output.txt", "dat" file and "png" file under the current directory
list1=os.listdir(os.getcwd())
for mm in list1:
    if (".dat" in mm) or (".png" in mm) or ("output.txt" in mm) :
        os.remove(mm)
#start the timer
s = timer()
for m in list1:
    if ".csv" in m:
        name=m.split('.',1)[0]
        readfromfile(m,name)
        # print "\t\tProcessing "+ name
#end the timer
e = timer()
# print "\nThe time spent was " +'{:.1f}'.format(e - s)+" seconds!"