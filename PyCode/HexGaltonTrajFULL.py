import numpy as np
import math
import matplotlib.pyplot as plt
import random


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
    bestDistance=1000
    D=((2*pY*vY+2*pX*vX)**2-4*(vX**2+vY**2)*(pY**2+pX**2-r**2))
    if D>=0 and isTorus:
        # This if statement checks the collision with the disperser by solving for time t
        t1=(-2*pY*vY-2*pX*vX+D**0.5)/(2*(vX**2+vY**2))
        t2=(-2*pY*vY-2*pX*vX-D**0.5)/(2*(vX**2+vY**2))
        if t1<0 and t2<0:
            pass
        else:
            for i in [t1,t2]:
                if i<0:
                    continue
                newPx=pX+vX*i
                newPy=pY+vY*i
                newNorm=((pX-newPx)**2+(pY-newPy)**2)**0.5
                if newNorm<bestDistance:
                    bestDistance=newNorm
                    bestPx=newPx
                    bestPy=newPy
                    besttime=i
                    bestxTravel=vX*i
                    bestyTravel=vY*i
            time+=besttime
            return (bestPx,bestPy,vX,vY,vS,False,-1,time,bestxTravel,bestyTravel)

    elif(not isTorus):
        # if we pass in a point on the disperser we call reflect
        vX,vY,vS=reflect(pX,pY,vX,vY,vS)

    # If it misses the disperser and we are not passing in a point on the disperser
    # we check for collisions with the boundaries.
    for i in range(0,len(tabL)):
        if i==wall:
            continue
        if (tabL[i][2]*vY-tabL[i][3]*vX) == 0:
            print('WARNING DIVISION BY ZERO LINE 113')
            continue
        t=(tabL[i][1]*tabL[i][2]+tabL[i][3]*pX-tabL[i][3]*tabL[i][0]-tabL[i][2]*pY)/(tabL[i][2]*vY-tabL[i][3]*vX)
        if t<0:
            continue
        newPx=pX+vX*t
        newPy=pY+vY*t
        newNorm=((pX-newPx)**2+(pY-newPy)**2)**0.5
        if newNorm<bestDistance:
            bestDistance=newNorm
            bestPx=newPx
            bestPy=newPy
            bestWall=i
            besttime=t
            bestxTravel=vX*t
            bestyTravel=vY*t
    time+=besttime
    return (bestPx,bestPy,vX,vY,vS,True,bestWall,time,bestxTravel,bestyTravel)

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

def getXYAng(r,epsilon,n,m):
    # Gets initial positions with angles.
    xyPos=[[],[],[]]
    for sang in np.linspace(np.pi+np.pi/6,np.pi/6,n):
        x=r*cos(sang)
        y=r*sin(sang)
        if(y<0):
            y-=epsilon
        else:
            y+=epsilon
        sang-=np.pi/2
        for j in np.linspace(sang,sang+np.pi,m):
            xyPos[0].append(x)
            xyPos[1].append(y)
            xyPos[2].append(j)
    return xyPos
################################################################################
################################ Interact ######################################
################################################################################
r=(1/1.125)/(0.5/(-0.5**2+1)**0.5)*0.5
sides=6
startX=0.6
startY=0
spin=0
eta=0.9
N=32
timeCap=30000
etaRange=11
nXY=1
nAng=1
etaStart=0
etaEnd=1
################################################################################
################################################################################
epsilon=0.0001
########################## TRAJECTORY MAP ######################################
# xyang=getXYAng(r,epsilon,nXY,nAng)
xyang=[[0.777],[0.175],[2.4000]]

for px,py,startAng in zip(xyang[0],xyang[1],xyang[2]):
    fname='galton('+str(round(eta,2))+','+str(round(r,2))+'_x'+str(round(px,2))+'_y'+str(round(py,2))+'_ang'+str(round(startAng,5))+')_ite'+str(N)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    pX=px
    pY=py
    vX=np.cos(startAng)
    vY=np.sin(startAng)
    vS=spin
    norm=(vX**2+vY**2+vS**2)**0.5
    vX=vX/norm
    vY=vY/norm
    vS=vS/norm
    trajX=[]
    trajY=[]
    isTorus=True # Don't change unless we are starting on the disperses
    wall=-1
    time=0


    (tabX,tabY)=make_ngon(sides)
    tabLineEqs=getLines(tabX,tabY,sides)

    xFrame = pX
    yFrame = pY
    for i in range(0,N):
        trajX.append(xFrame)
        trajY.append(yFrame)
        (pX,pY,vX,vY,vS,isTorus,wall,time,xTravel,yTravel)=BilliardIte(pX,pY,vX,vY,vS,tabLineEqs,wall,r,isTorus,time)
        xFrame += xTravel
        yFrame += yTravel
        trajX.append(xFrame)
        trajY.append(yFrame)
        if isTorus:
            (pX,pY,wall)=torus(pX,pY,wall)

    # ax.plot(tabX, tabY,'k',linewidth=1)
    ax.plot(trajX, trajY,'k',linewidth=0.3)
    fullboard=[]
    for xCenter in np.linspace(-2*(1.5),2*(1.5),3):
        for yCenter in np.linspace(2*(-(1-(0.5**2))**0.5),2*((1-(0.5**2))**0.5),3):
            fullboard.append(plt.Circle((xCenter, yCenter), r, fill=False,color='k',linewidth=0.3))
    for xCenter in np.linspace(-3*1.5,3*1.5,4):
        for yCenter in np.linspace(3*(-(1-(0.5**2))**0.5),3*((1-(0.5**2))**0.5),4):
            fullboard.append(plt.Circle((xCenter, yCenter), r, fill=False,color='k',linewidth=0.3))
    for jj in range(0,3*3+4*4):
        ax.add_patch(fullboard[jj])
    ############################## Save of Show ###################################
    # plt.show()
    plt.axis('off')
    plt.savefig(fname+'.eps',transparent=True)
    plt.close('all')
