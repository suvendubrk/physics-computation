#@ Suvendu Barik | B.Sc (Research) - Physics | Roll : 1510110412
#---Start of program---

'''
This is the pendulum wave simulation program made in vidle. It has two components
1. The simulation
2. The analyzed wave motion as seen

The first part would show you the actual simulation of how pendulum waves looks like.
The second part would allow us to study about the behaviour of pendulum waves by
demonstrating the wave motion in x-z plane.

The prerequisites to understand the simulation is to understand how pendulum works.

Please make sure that you have a good processor (Minimum i3 Core, 4th gen)
and atleast 4GB of RAM to run the simulation effectively.
'''

#Importing modules 
from visual import *
from visual.graph import *  
import numpy as np

#Scene windows related code
scene1=display(title="Pendulum Waves Simulation")
scene2 = display(title="Wave motion analyzer - (Top view of arrangement)",x=scene1.x+scene1.width)
scene2.userspin=false
scene2.autoscale=true

#Parameters defined (User changable)
#----------------------------------

g = 9.81
depression = np.deg2rad(15)
l = 10
step = 2

#These are the parameters where we would set for what time this
#simulation needs to run
dx=1e-2
x=0.0
xStop = 20 #Running simulation time

#Number of balls (balls = number given-1)
n=10

#Rate of running
r = 100

#Important, the inital conditions
#[program const ,computed theta ,velocity (in rad/s)]
y=np.array([1,0,0])

#----------------------------------

#Below are some of the basic tools needed for simulation..

#---They are the heart of the program---

#Runga Kutta method (customized)
def run_kut4(F,x,y,h,n):
    K0 = h*F(x,y,n)
    K1 = h*F(x+h/2,y+K0/2,n)
    K2 = h*F(x+h/2,y+K1/2,n)
    K3 = h*F(x+h,y+K2,n)
    return (K0+(K1+K2)*2+K3)/6

def RK4_diff(F,x,y,xStop,h,n):        
    X=[]
    Y=[]
    X.append(x)
    Y.append(y)
    while x<xStop:
        h = min(h,xStop-x)
        y = y + run_kut4(F,x,y,h,n)
        x=x+h
        X.append(x)
        Y.append(y)
    return np.array(X),np.array(Y)

#---End of section---

#---The process section---

#Length method - Used to find the length of pendulum w.r.t index n
def length(n):
    return l-(n-1)*np.tan(depression)*step

#Defining the S.H.M motion of the pendulum, for index n
def f(x,y,n):
    F= np.zeros(3)
    F[0]=y[1] #1st diff
    F[1]=-np.sin(y[0])*(length(n)/g)**0.5#2nd diff
    F[2]=y[0] #function
    return F

#Making the initial array C1 containing X and Y information
#Then using the theta array, we would create the position array C2
C1=[]
C2=[]

#Here compututaion for theta and position arrays are done
for i in range(1,n):
    X,Y=RK4_diff(f,x,y,xStop,dx,i)
    C1.append([X,Y])
    X = -length(1)*np.sin(C1[i-1][1][:,2])
    Y = -length(1)*np.cos(C1[i-1][1][:,2])
    Z = np.ones(len(X))*step*i
    C2.append([X,Y,Z])
print("Positions for each pendulum computed.")

#---End of process section---

#---Deployment of simulation---

#Assigning pendulums. Before that, we would make the pendulum array C3
#Also, we would make thread array C4
C3=[]
C4=[]
for i in range(1,n):
    C3.append(sphere(radius=0.1,display=scene1))
    C4.append(cylinder(pos=(0,length(i)-l,i*step),radius=0.02,display=scene1))

#Making a bar showing the handle of pendulum model
bar = cylinder(pos=(0,0,0), axis=(0,length(n)-l,n*step),radius=0.2,display=scene1)

#Making a reference line showing all balls in same level
reference = cylinder(pos=(0,-l,0),axis=(0,0,n*step),
                     radius=0.02,color=color.red,display=scene1)
print("Deployment of scene1 completed...")

#Making stuffs for scene 2
#We would make an object pos, for simulation of points
posArray=[]

#Length defined to set distance between points
stepLength = 5
siz=5
norm = (stepLength+siz)*n+10 
for i in range(1,n):
    posArray.append([step*i+stepLength*(i-1),0])    
point = points(pos=posArray,size=siz,color=color.yellow,display=scene2)

#Making axes
xAxis = curve(pos=[(-1,0),(norm,0)], radius=0.05,display=scene2)
yAxis = curve(pos=[(0,-l),(0,l)], radius=0.05,display=scene2)
scene2.center=(norm/2,0,0)

#Some text content
txt = label(pos=(norm/2,l+30,0),
            text="Along the horizontal line towards left, you will see bobs\n"+
            " of pendulum having progressively shorter lengths..."+
            "\nPress any key in pendulum window to continue",height=10)
print("Deployment of scene2 completed.")

#---End of deployment section---

#Data stimulant/recorder module (Yet to make for extension)
#Creation of pvalues pVal, and other essential arrays
#C2 already exists, so not much needed for it
#...
      
#Now, setting the positions and simulating the motion

for t in range(len(C2[0][0])):
   rate(r)
   pvalues=[]            
   for i in range(1,n):
       pvalues.append((C2[i-1][2][t]+stepLength*(i-1),
                       C2[i-1][0][t]))
       C3[i-1].pos = [C2[i-1][0][t],
                      C2[i-1][1][t],
                      C2[i-1][2][t]]       
       C4[i-1].axis= [C2[i-1][0][t],
                      C2[i-1][1][t]-length(i)+10,
                      0]
       
       #If uncommented, it will print the positions...
       #print(C3[i-1].pos)
       
   if(t==0):
    print("Press any key to continue in scene1")   
    ev = scene1.waitfor('keydown')
    txt.visible=false
    del txt
   point.pos = pvalues 

#---End of simulation loop---

#---End of program code---   
              
             
         
                
    
        
        

    
    








