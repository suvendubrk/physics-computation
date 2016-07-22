# -*- coding: utf-8 -*--
"""
Created on Tue Apr 19 13:02:32 2016
@authors: Suvendu Barik and Pranjal Agrawal

This program plots the voltage fluctuations around the capacitors of each LCR connected, with
consequent circuit connected to previous capacitor as voltage source

This code is also available in my github website : Do refer it
"""
#----------------------

#I am working in numpy
import numpy as np

#This is the list for X as time, Y as Voltage
X=[]
Y=[]
voltage=10 #This one is the initial constant given (The max voltage in AC)

#---------------------------------------------
#Modified Runga Kutta method for this purpose
def run_kut4(F,x,y,c,h,m=0):
    K0 = h*F(x,y,c,m)
    K1 = h*F(x+h/2,y+K0/2,c,m)
    K2 = h*F(x+h/2,y+K1/2,c,m)
    K3 = h*F(x+h,y+K2,c,m)
    return (K0+(K1+K2)*2+K3)/6

def RK4_diff(F,x,y,xStop,h,c,M="~"):        
    X=[]
    Y=[]
    X.append(x)
    Y.append(y)
    i=0
    m=0
    elemental=True
    if(type(M)==str):
         m=voltage
         elemental=False
    while x<xStop:
        h = min(h,xStop-x)
        y = y + run_kut4(F,x,y,c,h,m)
        x=x+h
        X.append(x)
        Y.append(y)
        i+=1        
        if(elemental):
         m=M[i] 
    return np.array(X),np.array(Y)
#----------------------------------------------
    
#Function for the differential equation we are working with    
def f(x,y,c,m):
    F=np.zeros(2)
    F[0]=y[1]
    F[1]=-(c[0]/c[1])*y[1]-(1/(c[1]*c[2]))*y[0]+m/c[1]
    return F

#This is the plotter function, to get the results here
def plot(xStop,const):
    
    #For the first circuit, connected with inital AC source
    X,Y=RK4_diff(f,0.0,np.array([0.0,0.0]),xStop,0.1,const[0])
    #I import matplotlib here...    
    import matplotlib.pyplot as plt
    #This is the plotter function for the first one
    plt.plot(X,Y[:,0]) 
    plt.xlabel("Time (in unit)")
    plt.ylabel("Voltage (in unit)")
    #This iteration enables to take the voltage fluctuation of the first circuit
    #and for that discritised time, indicate the voltage fluctuation for the
    #other circuit. This goes on for the next circuit. 
    index=[]
    for i in range(1,len(const)):
        M=Y[:,0]/const[i-1][2]
        X,Y=RK4_diff(f,0.0,np.array([0.0,0.0]),xStop,0.1,const[i],M)
        plt.plot(X,Y[:,0])
        index.append(i)
    plt.legend(index)    
        

#The initiation, 
#first parameter is for what time I want to see the fluctuation, 
#and the second is the information of [Resistance, Inductance, Capacitance] of 
#each circuit from the sequence


#Just a sample output
def demo():
   print("The circuit 1. R=2.8, L=1, C=1\nThe circuit --parallel-- 2. R=1, L=1, C=1")
   print("The circuit --parallel-- 3. R=0.01, L=1, C=1\nThe circuit --parallel-- 4. R=9, L=1, C=1")
   print("The circuit --parallel-- 5. R=1, L=1, C=0.7")
   plot(1000,[[2.8,1,1],[1,1,1],[1,1,1],[0.01,1,1],[9,1,1],[1,1,0.7]]) 
   
demo()   