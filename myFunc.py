# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 03:33:31 2015
@author: Suvendu Kumar Barik
Welcome to myFunc : Version 0.05 dated 11/1/2015
This is the one-stop class where all the necessary functions will be gathered

Till now, it has functions for integration (more added), differentiation, interpolations, 
graphing, fits, roots, linear equations, eigenvalues/eigenvectors and partial derivatives

log - 14-04-2016
Methods on solving the linear equations, finding eigenvalues and eigenvectors and
solving partial derivatives is included after extensive bugging.

Net method errors = 1
"""
import matplotlib.pyplot as plt
import numpy as np

#Some fancy code, but never in use.. till next assignment ;)
def welcome(BusyText=""):
    print("\nHi, I am Virtual Suvendu.",BusyText)

def endTag(Text=""):
    print(Text,"\nEnd of Program")
#Yo, fancy things done
#-------------------------------------Serious Stuff------------------------------------
 
#The averaging integration function
def integrate_average(f,high,low,h=1e-4) :
   answer1=0
   answer2=0   
   iteration = (high-low)/h
   for i in range(int(iteration)) :
    answer1 += f(low+((i)*h))*h
    answer2 += f(high-((i)*h))*h
   return (answer1+answer2)/2

#Composite trapezoidal method
def integrate_trapezoidal(f,high,low,h=1e-4):
    iteration = (high-low)/h
    answer = low+high
    for i in range(1,int(iteration)-1):
        answer += 2*f(low+(i*h))
    answer = answer*h*0.5
    return answer

#Composite simpson's rule method
def integrate_simpson(f,high,low,h=1e-4):
    iteration = (high-low)/h
    answer = low
    for i in range(1,int(iteration)-1):
        if i%2==0:
            answer += 2*f(low+(i*h))
        else:
            answer += 4*f(low+(i*h))
    answer = answer*h*(1/3) 
    return answer       
      
#Derivative of function - asymmetric
def derivative1(f,x,h=1e-7) :
   answer = (f(x+h)-f(x))/h    
   return answer

#Derivative of function - symmetric   
def derivative2(f,x,h=1e-7) :
   answer = ((f(x+h)-f(x-h))/(2*h))
   return answer
   
#Here comes grapher, one stop solution to plot graph easily
def grapher(plots,types,legends,pos=0):
    for i in range(0,len(plots)):
        plt.plot(plots[i][0],plots[i][1],types[i])
    plt.legend(legends,loc=pos)
    plt.show()    
   
#Interpolations - Newtonian (Yet not a clean code)
def interpolation_newton(xData, yData, prec=0.01):
    m = len(xData)    
    a = yData.copy()
    for k in range(1,m):
        for i in range(k,m): 
         a[i] = (a[i]-a[k-1])/(xData[i]-xData[k-1])
    i = min(xData)-prec
    procY = []
    procX = []
    while i<=max(xData):
        i+=prec
        procX.append(i)
        p=a[m-1]
        for k in range(1,m):
         p = a[m-1-k] + (i-xData[m-1-k])*p
        procY.append(p)
    return [procX,procY]

#Interpolations - Lagrangian (Yet not a clean code)
def interpolation_lagrange(xData, yData, prec=0.01):
     m = len(xData)
     i = min(xData)-prec     
     procY = []
     procX = []
     while i<=max(xData):
      i+=prec
      procX.append(i)  
      ans=0
      for n in range(m):
       delta=1
       for j in range(m):
        if not n==j :
         delta*=(i-xData[j])/(xData[n]-xData[j])
       ans+=yData[n]*delta
      procY.append(ans) 
     return [procX,procY]

#Interpolation - Linear : (Bugs exists)
def interpolation_linear(xData, yData, prec=0.01):
     m = len(xData)
     i = min(xData)-prec     
     procY = []
     procX = []
     while i<=max(xData):
      i+=prec
      procX.append(i)  
      index = -1  
      for p in range(m-1):
       if xData[p]<=i and xData[p+1]>=i:
           index=p
           break
      a = (yData[index]*(xData[index+1]-i)/(xData[index+1]-xData[index]))
      b = (yData[index+1]*(i-xData[index])/(xData[index+1]-xData[index]))
      procY.append(a+b)
     return [procX,procY]         

#Supplementary methods for fits
def aAdd(x):
    out = 0    
    for i in range(0,len(x)):
        out += x[i]
    return out    

def mAdd(x,y):
    if len(x)==len(y):
        for i in range(0, len(x)):
            x[i] += y[i]
    else:
        return 0
    return x     
    
def mMult(x,y):
    if len(x)==len(y):
        for i in range(0, len(x)):
            x[i] *= y[i]
    else:
        return 0
    return x
#End of supplementary methods for fits

#This method does the linear fit for the given data 
def linearFit(x,y,err="~~"):
   (a0,a1)=(0,0)   
   if len(x)==len(y) :
       x1 = np.array(x)
       y1 = np.array(y)
       if err=="~~":
            w=np.array((x1/x1))
       else:
            w = np.array(err)
            w = (1/w)
            w = w*w
       a1 = (aAdd(x1*y1*w)*aAdd(w)-aAdd(x1*w)*aAdd(y1*w))/(aAdd(x1*x1*w)*aAdd(w)-aAdd(x1*w)**2)
       a0 = (aAdd(y1*w)-(a1*aAdd(x1*w)))/aAdd(w)     
   else:
       return 0    
   return [a0,a1]       

#This method doeas the exponential fit for the given data
def exponentialFit(x,y,err="~~"):
   (a0,a1)=(0,0)   
   if len(x)==len(y):
       x1 = np.array(x)
       y1 = np.array(y)
       z1 = np.log(y1)
       if err=="~~":
           w=y1
           w=w*w
       else:
           w = np.array(err)
           w = ((y1*y1)/(w))
           w=w*w
       a1 = (aAdd(x1*z1*w)*aAdd(w)-aAdd(x1*w)*aAdd(z1*w))/(aAdd(x1*x1*w)*aAdd(w)-aAdd(x1*w)**2)
       a0 = (aAdd(z1*w)-(a1*aAdd(x1*w)))/aAdd(w)        
   else:
       return 0    
   return [np.e**a0,a1] 
     
#Returns R value for the given fit     
def rValue(func,x,y,fitType="linear"):
    x1 = np.array(x)
    y1 = np.array(y)
    avY = aAdd(y1)/(len(x))
    if fitType=="linear":
        f=np.array(func[0]+(func[1]*x1))
    elif fitType=="exponential":
        f=np.array(np.log(func[0]+(func[1]*np.log(x1))))
    numer = aAdd((y1-f)**2)
    denom = aAdd((y1-avY)**2)
    return 1- numer/denom
    
#Returns X value for the given fit    
def xValue(func,x,y,e,fitType="linear"):
    x1 = np.array(x)
    y1 = np.array(y)
    e1 = np.array(e)
    if fitType=="linear":
        f=np.array(func[0]+(func[1]*x1))
    elif fitType=="exponential":
        f=np.array(func[0]*np.e**(func[1]*x1))
    numer = aAdd(((y1-f)/e1)**2)
    return numer

#Root finding - Brute Force method (BF) - Returns the interval where root lies
#if the int_mode=1
def root_BruteMethod(f,b,delta=1e-4,int_mode=0):
    a=b[0]    
    i=0
    while(a<b[1]):
        i+=1
        if(f(a)*f(a+delta)>0):
            a = a+delta
            #print(a)
        else:
            print(i)
            if(int_mode==1):    
                return [a,a+delta]   
            else:
                return (2*a+delta)*0.5
   
    
#Root finding - Bisection (B)
def root_Bisection(f,x1,x2,tol=1e-9,switch=0,int_mode=0):
    f1 = f(x1)
    t=0
    if f1==0: return x1
    f2 = f(x2)
    if f2==0: return x2
    if np.sign(f1)==np.sign(f2): return "~~"
    n = int(np.ceil(np.log(abs(x2 - x1)/tol)/np.log(2)))
    for i in range(n):
        t+=1
        x3=0.5*(x1+x2)
        f3 = f(x3)
        if switch==1 and abs(f3)>abs(f1) and abs(f3)>abs(f2):
            return "~~"
        if x3==0:
            return x3
        if np.sign(f2)!=np.sign(f3):
            x1=x3
            f1=f3
        else:
            x2=x3
            f2=f3
    print(t)            
    if int_mode==1:      
        return [x1,x2]
    else:
        return (x1+x2)*0.5
       
#Root finding - Newton Raphson method (NR)
def root_NewtonRaphson(f,df,a,b,tol=1e-9):
    x=a
    i=0
    while(x<=b):
        i+=1
        dx = -f(x)/df(x)
        x+=dx
        if abs(dx)<tol :
            print(i)
            return x
    print("Too many iterations") 
           
    
    
#Root finding - Combinations.. (NR <--> B) (Basically the same code as in book..
#as my earlier code was not good as this is given).
def root_NR_Bisection(f,df,a,b,iterations=30,tol=1e-9):
    fa = f(a)
    t=0
    if fa==0:
        return a
    fb = f(b)
    if fb==0:
        return b
    if np.sign(fa)==np.sign(fb):
        return "~~"
    x = (a+b)*0.5
    for i in range(iterations):
        t+=1
        fx=f(x)
        if fx==0:
            return x
        if np.sign(fa)!=np.sign(fx):
            b=x
        else:
            a=x
        dfx = df(x)
        try:
            dx=-fx/dfx 
        except ZeroDivisionError:
            dx = b-a
        x=x+dx
        if (b-x)*(a-x)<0:
            dx=0.5*(b-a)
            x=a+dx
        if abs(dx)<tol*max(abs(b),1.0):
            print(t)
            return x
    print("Too many iterations in Newton-Raphson")  


#**********************************************Welcome to PHY106 code segment>***********************************************

#There are some tools to supplement the code for solving linear equations
#The ones not using the dominant ones are still given above
def swapRows(A,a,b):
    x=A[a].copy()
    A[a]=A[b]
    A[b]=x
    return A
    
def scaled_matrix(A):
    s=[]
    for i in range(0,len(A)):
        s.append(max(abs(A[i])))
    return s         
    
def returnDominant(A,B):
    n = len(A)
    A = np.array(A)
    B = np.array(B)    
    s = scaled_matrix(A)    
    for k in range(0,n-1):
        p = int(np.argmax(abs(A[k:n,k]/s[k:n])))+k
        if p!=k:
            swapRows(A,k,p) 
            swapRows(B,k,p)
    return A,B  

#End of the code fragment here...

#These are the methods to solve the linear equations. Standard ones 
#------------------------------------------------------------------------------------------------------------
#Gauss method
def gauss_solution_n(a,b,expl_out=False):
    n=len(b)
    x=np.zeros(n)
    expl=""
    for k in range(0,n):
        for i in range(k+1,n):
          expl+="Using pivot as "+str(k)+" row."+" Transforming "+str(k+1)+" row on index "+str(i)+"\n"  
          if a[i][k]!=0.0:           
           lamb_da= a[i][k]/a[k][k]
           expl+=("\nLambda is "+str(lamb_da)+"\n")
           a[i][k+1:n] = np.array(a[i][k+1:n]) - lamb_da*np.array(a[k][k+1:n])
           b[i] = b[i] - lamb_da*b[k]
          expl+="The Matrix A : \n "+str(a)+"\nThe Matrix B : \n"+str(b)+"\n---\n"             
    for t in range(n-1,-1,-1):
        x[t]=(b[t]-np.dot(a[t][t+1:n],x[t+1:n]))/a[t][t]
    expl+="The solution matrix is :"+str(x)    
    if(expl_out):
        print(expl)
    return x
    
#LU_Decompositon method
def LU_decomposition_n(A,B):
    a=A.copy() 
    b=B.copy()  
    n=len(b)   
    
    for k in range(0,n-1):
        for i in range(k+1,n):
           if a[i][k]!=0.0:           
            lamb_da= float(a[i][k]/a[k][k])
            a[i][k+1:n] = np.array(a[i][k+1:n]) - lamb_da*np.array(a[k][k+1:n])
            a[i][k]=lamb_da
           
    for t in range(1,n):
        b[t]=(b[t]-np.dot(a[t][0:t],b[0:t]))
    b[n-1]=b[n-1]/a[n-1][n-1]    
    
    for q in range(n-2,-1,-1):
        b[q]=(b[q]-np.dot(a[q][q+1:n],b[q+1:n]))/a[q][q]  
    return b  
#--------------------------------------------------------------------------------------------------------------    


#Modified tools to solve the linear equations
#--------------------------------------------------------------------------------------------------------------
#Gauss method - Modified
def gauss_solution_m(a,b,expl_out=False):
    a,b = returnDominant(a,b)
    n=len(b)
    x=np.zeros(n)
    expl=""
    for k in range(0,n-1):
        for i in range(k+1,n):
          expl+="Using pivot as "+str(k)+" row."+" Transforming "+str(k+1)+" row on index "+str(i)+"\n"  
          if a[i][k]!=0.0:           
           lamb_da= a[i][k]/a[k][k]
           expl+=("\nLambda is "+str(lamb_da)+"\n")
           a[i][k+1:n] = np.array(a[i][k+1:n]) - lamb_da*np.array(a[k][k+1:n])
           b[i] = b[i] - lamb_da*b[k]
          expl+="The Matrix A : \n "+str(a)+"\nThe Matrix B : \n"+str(b)+"\n---\n"             
    for t in range(n-1,-1,-1):
        x[t]=(b[t]-np.dot(a[t][t+1:n],x[t+1:n]))/a[t][t]
    expl+="The solution matrix is :"+str(x)    
    if(expl_out):
        print(expl)
    return x

#LU Decompositon method
def LU_decomposition_m(A,B,pivot=True):
    a=A.copy() 
    b=B.copy()   
    if(pivot):
     a,b=returnDominant(a,b) 
    n=len(b)   
    for k in range(0,n-1):
        for i in range(k+1,n):
           if a[i][k]!=0.0:           
            lamb_da= float(a[i][k]/a[k][k])
            a[i][k+1:n] = np.array(a[i][k+1:n]) - lamb_da*np.array(a[k][k+1:n])
            a[i][k]=lamb_da
      
         
    for t in range(1,n):
        b[t]=(b[t]-np.dot(a[t][0:t],b[0:t]))
    b[n-1]=b[n-1]/a[n-1][n-1]    
     
    for q in range(n-2,-1,-1):
        b[q]=(b[q]-np.dot(a[q][q+1:n],b[q+1:n]))/a[q][q]  
    return b    
#--------------------------------------------------------------------------------------------------------------

#Generalized method - Better to refer them
def gauss_solution(A,B):
    if(type(B[0])!=list and type(B[0])!=np.ndarray):
        return gauss_solution_n(A,B) 
    else:
        Ans=[]
        for i in range(len(B)):
            Ans.append(LU_decomposition_m(A,B[i]))
        np.transpose(Ans)
        return Ans


def LU_decomposition(A,B,pivot=True):
    if(type(B[0])!=list and type(B[0])!=np.ndarray):
        return LU_decomposition_m(A,B,pivot) 
    else:
        Ans=[]
        for i in range(len(B)):
            Ans.append(LU_decomposition_m(A,B[i],pivot))
        np.transpose(Ans)
        return Ans

#To find inverse
def inverse(A): 
  n = len(A)  
  B=np.zeros((n,n))
  np.fill_diagonal(B,1)    
  return LU_decomposition(A,B,False)  

#End of linear equation solver code fragment

#A large component of eigenvalue component is here........................
#***Some Basic operations to be resolved***

#---operation of Matrix multiplied by vector---
def operation(A,v):
        X=[]
        for i in range (len(A)):
            X.append(np.dot(A[i],v))
        return X
#---Done---

#---norm of the vector---
def norm(v):
        return np.dot(v,v)**0.5
#---Done---

#Comparing two vectors to be similar, under the parameter 'eta'
def compare(v1,v2,eta=1e-2,neg=2):
   if(len(v1)!=len(v2)):
       return "~"
   v3=[]    
   for i in range (len(v1)):    
    v3.append(v1[i]-(-1)**neg*v2[i])
    
   if(np.dot(v3,v3)**0.5<eta):         
         return True
   return False
#--Done--

#**Basic operations are resolved, procceding to the power method

'''
Power function : Inputs - Matrix A, Guess v0, safetyClose sC
----------------------------------------------
Process -
0. eta is given, safetyClose (sC) is given 
1. Find z0 = A*v0
2. v1=v0/norm(v0)
3. Find z1 = A*v1
4. If compare(v0,v1,eta) returns True : eigenvalue = norm(v1), eigenvector = v1
5. If compare(v0,v1,eta) returns False : put v0=v1 and continue the loop
6. If loop exceeds the sC value, which by default is 50000, then loop crashes, 
   returning no eigenvalue 

'''
def power(A,v0,eta=1e-2,eigenVector=False,sC=50000):
    switch=True
    i=0
    neg=0
    while(switch):
        i+=1
        z0=operation(A,v0)
        v1=np.array(z0)/norm(z0)
        z1=operation(A,v1)       
        if(np.dot(z0,z1)<0):
          neg=1
        else:
          neg=2 
        #print(v0,v1,np.dot(z0,z1),neg) 
        if(compare(v0,v1,eta,neg)):
           if(eigenVector):
             return [norm(z1)*(-1)**neg,v1]           
           else:    
             return [norm(z1)*(-1)**neg]
        else:
           v0=v1
        if(i>sC):            
            return "~"
#--Done-- 

'''
Inverse power method : Inputs -  Matrix A, Guess v0, safetyClose sC
-----------------------------------------------------
Process - 
0. eta is given, safetyClose (sC) is given 
1. If A is invertible, great!!, and continue with power method by shift s
2. If not, return a crash

'''
def inversePower(A,v0,s=0,eta=1e-7,eigenVector=False,sC=50000):
        A = np.array(A)-np.eye(len(A))*s
        if(np.linalg.det(A)==0):
         return "~"
        a =power(inverse(A),v0,eta,eigenVector,sC)  
        if(a[0]=="~"):
            return("~")
        else:
            if(eigenVector):
             return([a[0]+s,np.array(a[1])])
            else:
             return ([a[0]+s])   
#--Done--

'''
Search for eigenvalues (eigenSearch) : Inputs - Matrix A, Guess v0, safetyClose sC
                                       shiftRange sR [Format : (min, max), step sT]
-------------------------------------------------------------------------------------
Process-
0.Everything is given
1.Check whether the inversePower works on the shift as progressed
2.List them, with eigenvectors also   
                                    
'''  
def eigenSearch(A,v0,sR,sT,eta=1e-7,eigenVector=False,sC=50):
    n = np.int((sR[1]-sR[0])/sT)
    eigenValues=[]
    eigenVectors=[]
    for i in range(n):
        a = inversePower(A,v0,sR[0]+i*sT,eta,eigenVector,sC)
        if(a[0]!="~"):
            eigenValues.append(a[0])
            if(eigenVector):
                eigenVectors.append(np.array(a[1]))                         
    return[eigenValues,eigenVectors]


'''
The jacobi method
----------------------------------------------------------------------------------------
Well, its funny. It is better to always use Jacobi
'''
def jacobi(Matrix,th):
    n=len(Matrix)
    A=Matrix.copy()
    P=np.eye(n)
    #Sub function to check all off diagonal values 
    #are lesser than threshold    
    def check():
        for i in range(n):
            for j in range(n):
                if(i!=j):
                    if(abs(A[i][j])>th):
                        return False
        return True            
    
    #The main process here                
    close=False
    while(not check()):
        for k in range(n):            
            #print("1",check())
            if(close):
                break
            for l in range(n):            
                if(k!=l):
                    if(abs(A[k][l])>th):
                        phi = -(A[k][k]-A[l][l])/(2*A[k][l])
                        if(A[k][k]!=A[l][l]):
                         t=np.sign(phi)/(abs(phi)+np.sqrt((phi**2)+1))
                        else:
                         t=1   
                        c=1/(np.sqrt(1+t**2))
                        s=t*c                                                               
                        tau=s/(1+c)    
                                                
                        #Transformation of A
                        A[k][k]=A[k][k]-t*A[k][l]
                        A[l][l]=A[l][l]+t*A[k][l]
                        A[k][l]=0
                        A[l][k]=0
                        for i in range(n):
                            if(i!=k and i!=l):
                                temp=A[k][i]
                                A[k][i]=A[k][i]-s*(A[l][i]+tau*A[k][i])
                                A[i][k]=temp
                                temp2=A[l][i]
                                A[l][i]=A[l][i]+s*(temp-tau*A[l][i])
                                A[i][l]=temp2
                                
                        #Transformation of P
                        for i in range(n):
                            temp=P[i][k]
                            P[i][k]=P[i][k]-s*(P[i][l]+tau*P[i][k])
                            P[i][l]=P[i][l]+s*(temp-tau*P[i][l])
                        
                        #Check whether off elements are fine, so that
                        #we can stop
                        if(check()):
                            #print("2")
                            close=True
                            break
    return [A,P]                

      

    
    
         
            
    
