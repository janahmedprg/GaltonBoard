import numpy as np
import math
import matplotlib.pyplot as plt
import random
from scipy import stats
import seaborn as sns
import csv


pi=np.pi

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def make_ngon(n):
    # Makes the heaxagonal bounds
    tabX=[1]
    tabY=[0]
    for i in range(1,n+1):
        tabX+=[np.cos(i*2*np.pi/n)]
        tabY+=[np.sin(i*2*np.pi/n)]
    # print(tabX)
    # print(tabY)
    return (tabX,tabY)

def getLines(tabX,tabY,n):
    # Gets the vertices and vectors of the hexagonal table for plotting
    lineEqs=[]
    for i in range(0,n):
        vX=tabX[i+1]-tabX[i]
        vY=tabY[i+1]-tabY[i]
        lineEqs.append([tabX[i],tabY[i],vX,vY])
    return lineEqs

def mat_mul(Rinv,T,R,vX,vY,vS):
    # Matrix multiplication function
    V=[[vS],[vX],[vY]]
    V= np.matmul(R,V)
    V = np.matmul(T,V)
    V = np.matmul(Rinv,V)
    return (float(V[1][0]),float(V[2][0]),float(V[0][0]))

def rotate(x):
    # Creates a rotational matrix
    return np.matrix([[1,0,0],[0,cos(x),sin(x)],[0,-sin(x),cos(x)]],dtype=float)

def reflect(pX,pY,vX,vY,vS):
    # No-slip reflection on the disperser
    Tr = [[-cos(eta*pi),sin(eta*pi),0],[sin(eta*pi),cos(eta*pi),0],[0,0,-1]]
    if pX==0 and pY>0:
        Rinv=rotate(0)
        R=rotate(0)
        (vX,vY,vS)=mat_mul(Rinv,Tr,R,vX,vY,vS)
    elif pX==0 and pY<0:
        Rinv=rotate(-np.pi)
        R=rotate(np.pi)
        (vX,vY,vS)=mat_mul(Rinv,Tr,R,vX,vY,vS)
    else:
        if (pX>0 and pY>=0):
            normalAng=np.arctan(pY/pX)
        elif(pX<0 and pY>=0):
            normalAng=np.pi+np.arctan(pY/pX)
        elif (pX<0 and pY<=0):
            normalAng=np.pi+np.arctan(pY/pX)
        else:
            normalAng=2*np.pi+np.arctan(pY/pX)
        Rinv=rotate(np.pi/2-normalAng)
        R=rotate(normalAng-np.pi/2)
        (vX,vY,vS)=mat_mul(Rinv,Tr,R,vX,vY,vS)
    return (vX,vY,vS)

def BilliardIte(pX,pY,vX,vY,vS,tabL,wall,r,isTorus,time):
    bestTime=1000
    bestTime1=1000

    if(not isTorus):
        # if we pass in a point on the disperser we call reflect
        vX,vY,vS=reflect(pX,pY,vX,vY,vS)
        # norm=(vX**2+vY**2+vS**2)**0.5
        # vX=vX/norm
        # vY=vY/norm
        # vS=vS/norm

    for i in range(0,len(tabL)):
        if i==wall:
            D = (tabL[i][3]*vX-vY*tabL[i][2])**2-4*(0.5*gravity*tabL[i][2])*(pX*tabL[i][3]-tabL[i][0]*tabL[i][3]-pY*tabL[i][2]+tabL[i][1]*tabL[i][2])
            if D<0:
                continue
            t1 = (-(tabL[i][3]*vX-vY*tabL[i][2])+D**0.5)/(2*(0.5*gravity*tabL[i][2]))
            t2 = (-(tabL[i][3]*vX-vY*tabL[i][2])-D**0.5)/(2*(0.5*gravity*tabL[i][2]))

            if t1>10**(-8) and t1>t2:
                ti = t1
            elif t2>10**(-8):
                ti=t2
            else:
                continue

            if ti<bestTime1:
                bestWall1 = i
                bestTime1 = ti
                bestxTravel1 = vX*ti
                bestyTravel1 = vY*ti - 0.5*(ti**2)*gravity
                bestPx1 = pX + vX*ti
                bestPy1 = pY + vY*ti - 0.5*(ti**2)*gravity
            continue

        D = (tabL[i][3]*vX-vY*tabL[i][2])**2-4*(0.5*gravity*tabL[i][2])*(pX*tabL[i][3]-tabL[i][0]*tabL[i][3]-pY*tabL[i][2]+tabL[i][1]*tabL[i][2])
        if D<0:
            continue
        t1 = (-(tabL[i][3]*vX-vY*tabL[i][2])+D**0.5)/(2*(0.5*gravity*tabL[i][2]))
        t2 = (-(tabL[i][3]*vX-vY*tabL[i][2])-D**0.5)/(2*(0.5*gravity*tabL[i][2]))

        for ti in [t1,t2]:
            if ti<0:
                continue
            if ti<bestTime1:
                bestWall1 = i
                bestTime1 = ti
                bestxTravel1= vX*ti
                bestyTravel1= vY*ti - 0.5*(ti**2)*gravity
                bestPx1 = pX + vX*ti
                bestPy1 = pY + vY*ti - 0.5*(ti**2)*gravity

    # This checks the collision with the disperser by solving for time t
    # We will have 4 different times because we have a parabolic trajectory
    coeff=[0.25*(gravity**2),-vY*gravity,-pY*gravity+vX**2+vY**2,2*pX*vX+2*pY*vY,pX**2+pY**2-r**2]
    intersect = np.roots(coeff)
    for ti in intersect:
        if np.iscomplex(ti) or ti<10**(-9):
            continue
        ti = float(ti.real)
        if ti < bestTime:
            bestTime = ti
            bestxTravel= vX*ti
            bestyTravel= vY*ti - 0.5*(ti**2)*gravity
            bestPx = pX + vX*ti
            bestPy = pY + vY*ti - 0.5*(ti**2)*gravity

    if bestTime1>bestTime and bestTime>10**(-9):
        time+=bestTime
        return (bestPx,bestPy,vX,vY-bestTime*gravity,vS,False,-1,time,bestxTravel,bestyTravel)
    else:
        time+=bestTime1
        return (bestPx1,bestPy1,vX,vY-bestTime1*gravity,vS,True,bestWall1,time,bestxTravel1,bestyTravel1)


def torus(pX,pY,wall):
    # Makes the torus effect. Instead of reflecting it starts on the opposite side of the hexagonal bounds.
    if(isTorus):
        if(wall==1 or wall==4):
            pY=-pY
            if wall==1:
                wall=4
            else:
                wall=1
        elif(wall==0 or wall==3):
            rDis=(pX**2+pY**2)**0.5
            if wall==0:
                ang=np.arctan(pY/pX)
                pX=rDis*np.cos(4/3*np.pi-ang)
                pY=rDis*np.sin(4/3*np.pi-ang)
                wall=3
            else:
                ang=np.arctan(pY/pX)
                pX=rDis*np.cos(1/3*np.pi-ang)
                pY=rDis*np.sin(1/3*np.pi-ang)
                wall=0
        elif(wall==2 or wall==5):
            rDis=(pX**2+pY**2)**0.5
            if wall==2:
                ang=np.arctan(pY/pX)
                pX=rDis*np.cos(10/6*np.pi-ang)
                pY=rDis*np.sin(10/6*np.pi-ang)
                wall=5
            else:
                ang=np.arctan(pY/pX)
                pX=rDis*np.cos(2/3*np.pi-ang)
                pY=rDis*np.sin(2/3*np.pi-ang)
                wall=2

        else:
            print("WARNING: VERTEX HIT")
    return (pX,pY,wall)

def getXYAng(r,epsilon,n):
    # Gets initial positions with angles.
    xyPos=[[],[],[],[]]
    for ii in range(n):
        sang = np.random.uniform(np.pi+np.pi/6,np.pi/6)
        x=r*cos(sang)
        y=r*sin(sang)
        if(y<0):
            y-=epsilon
        else:
            y+=epsilon
        sang-=np.pi/2
        j = np.random.uniform(sang,sang+np.pi)
        xyPos[0].append(x)
        xyPos[1].append(y)
        xyPos[2].append(j)
        xyPos[3].append(np.random.uniform(-0.5,0.5))
    return xyPos
################################################################################
################################ Interact ######################################
################################################################################
# r- radius of disperser; etaList- a list of eta's;
# timeCap- maximum allowed time before the program decides that the particle is trapped;
# Nparticles- number of particles; gravity- gravity;
# heightCap- the height the particle needs to reach to capture the distribution.
r=(1/1.125)/(0.5/(-0.5**2+1)**0.5)*0.5
sides=6
etaList= [np.arccos(1/3)/np.pi]
timeCap=500
Nparticles = 100000
gravity = 1
heightCap = -15
################################################################################
################################################################################
epsilon=0.0001
########################## TRAJECTORY MAP ######################################
xyang=getXYAng(r,epsilon,Nparticles)
sum = 0
fname='galton'+'_particles'+str(Nparticles)+'_grav'+str(round(gravity,3))+'_distv6'
trapped = 0
totalDist = []

for ETA in etaList:
    eta = ETA
    csvName = 'distData_eta'+str(round(eta,2))+'_timeCap'+str(timeCap)+'_particles'+str(Nparticles)+'_grav'+str(gravity)+'.csv'
    f = open(csvName, 'w', newline='')
    writer = csv.writer(f)
    writer.writerow(['xFrame','yFrame','pX','pY','vX','vY','vS','wall','isTorus','time','trapped'])
    xDist = []
    trapped = 0
    for px,py,startAng,spin in zip(xyang[0],xyang[1],xyang[2],xyang[3]):
        pX=px
        pY=py
        vX=np.cos(startAng)
        vY=np.sin(startAng)
        vS=spin
        norm=(vX**2+vY**2+vS**2)**0.5
        vX=vX/norm
        vY=vY/norm
        vS=vS/norm
        isTorus=True # Don't change unless we are starting on the disperses
        wall=-1
        time=0
        (tabX,tabY)=make_ngon(sides)
        tabLineEqs=getLines(tabX,tabY,sides)
        xFrame = pX
        yFrame = pY
        trap = True
        while time<timeCap:
            (pX,pY,vX,vY,vS,isTorus,wall,time,xTravel,yTravel)=BilliardIte(pX,pY,vX,vY,vS,tabLineEqs,wall,r,isTorus,time)
            xFrame += xTravel
            yFrame += yTravel
            if (yFrame<=heightCap):
                trap = False
                break
            if isTorus:
                (pX,pY,wall)=torus(pX,pY,wall)

        if trap:
            # print('trappe+d', trapped)
            trapped +=1
            writer.writerow([xFrame,yFrame,pX,pY,vX,vY,vS,wall,isTorus,time,True])
        else:
            D=vY**2-2*(yFrame-heightCap)*gravity
            if (D<0):
                print("D<0, this shouldn't happen check D eq and the Y finish line")
                continue

            trev1 = (-vY+D**0.5)/(-1*gravity)
            trev2 = (-vY-D**0.5)/(-1*gravity)
            if trev1<0 and trev2>0:
                trev = trev2
            elif trev2<0 and trev1>0:
                trev = trev1
            elif trev2>trev1:
                trev = trev1
            elif trev1> trev2:
                trev = trev2
            elif trev1<0 and trev2<0:
                print("NO TIME REVERSE")
                continue

            xDist.append(xFrame-vX*(trev-0.7))
            writer.writerow([xFrame,yFrame,pX,pY,vX,vY,vS,wall,isTorus,time,False])

    print(trapped)
    xDist.sort()
    totalDist.append(xDist)
    f.close()

sns.displot(totalDist, kind="kde", legend=False)
############################## Save of Show ###################################
plt.xlabel("x-coord")
plt.ylabel("proportrion")
plt.legend([round(e,2) for e in etaList])
# plt.show()
# plt.axis('off')
plt.savefig(fname+'.eps',transparent=True)
# plt.close('all')
