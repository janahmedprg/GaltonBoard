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
    if n == 4:
        return ([1, 1, -1, -1, 1],[-1, 1, 1, -1, -1])
    tabX=[1]
    tabY=[0]
    for i in range(1,n+1):
        tabX+=[np.cos(i*2*np.pi/n)]
        tabY+=[np.sin(i*2*np.pi/n)]
    # print(tabX)
    # print(tabY)
    print(tabX)
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
            # print('WARNING DIVISION BY ZERO LINE 113')
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

def box(pX,pY,wall):
    # Makes the torus effect. Instead of reflecting it starts on the opposite side of the hexagonal bounds.
    if(wall==0 or wall==2):
        pX=-pX
        if wall==0:
            wall=2
        else:
            wall=0
    elif(wall==1 or wall==3):
        pY=-pY
        if wall==1:
            wall=3
        else:
            wall=1
    else:
        print("WARNING: VERTEX HIT")
    return (pX,pY,wall)

def getXYAng(r,epsilon,n):
    # Gets initial positions with angles.
    xyPos=[[],[],[],[]]
    for sang in np.linspace(0,np.pi/2-0.05,n):
        x=r*cos(sang)
        y=r*sin(sang)
        if(y<0):
            y-=epsilon
        else:
            y+=epsilon
        angle = 0

        xyPos[0].append(x)
        xyPos[1].append(y)
        xyPos[2].append(angle)
        xyPos[3].append(sang)
    return xyPos
################################################################################
################################ Interact ######################################
################################################################################
r=0
sides=4
eta= np.arccos(1/3)/np.pi
N=5000
grid=100
perturbation = 0.01
################################################################################
################################################################################
epsilon=0.0001
########################## TRAJECTORY MAP ######################################


fname='sinayExp('+str(round(eta,3))+')_p'+str(round(perturbation,3))+'_grid' + str(grid)
fig, ax = plt.subplots()
ax.set_aspect('equal')

phiAx = np.linspace(0,np.pi/2-0.05,grid) #start and end needs to be same as in geXYAng() steps needs to be +1
rAx = np.linspace(0,0.99,grid) #start and end needs to be same as in r in linspace below steps needs to be +1
nColAx = np.zeros((grid,grid))
iCol = 0
iRow = 0

for r in rAx:
    xyang=getXYAng(r,epsilon,grid)
    for px,py,startAng,phiAng in zip(xyang[0],xyang[1],xyang[2],xyang[3]):
        pX=px
        pY=py
        vX=np.cos(startAng)
        vY=np.sin(startAng)
        vS=-vX*np.sin(phiAng)/(((1-np.cos(np.pi*eta))/(np.cos(np.pi*eta)+1))**0.5)+perturbation
        norm=(vX**2+vY**2+vS**2)**0.5
        vX=vX/norm
        vY=vY/norm
        vS=vS/norm
        isTorus=True # Don't change unless we are starting on the disperses
        wall=-1
        time=0
        nCol = -1

        (tabX,tabY)=make_ngon(sides)
        tabLineEqs=getLines(tabX,tabY,sides)

        for i in range(0,N):
            (pX,pY,vX,vY,vS,isTorus,wall,time,xTravel,yTravel)=BilliardIte(pX,pY,vX,vY,vS,tabLineEqs,wall,r,isTorus,time)
            if pY > r+0.1 or pY<-r-0.1:
                nCol = i
                break
            if isTorus:
                (pX,pY,wall)=box(pX,pY,wall)

        if nCol == -1:
            nCol = N
        nColAx[iRow][iCol] = nCol
        iCol += 1
        iCol = iCol%grid
    iRow += 1

ax.pcolormesh(phiAx, rAx, nColAx, shading='nearest', vmin=nColAx.min(), vmax=nColAx.max())
############################## Save of Show ###################################
plt.savefig(fname+'.eps')
plt.show()
plt.close('all')
