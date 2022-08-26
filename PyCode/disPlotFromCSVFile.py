import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import csv

############################ Interact ##########################################
particles = 100000
gravity = 1
etaList = [0,0.39]
timeCap = 500
heightCap = -15
freeFall = 1.4
################################################################################

totalDist = []
fname='galton'+'_particles'+str(particles)+'_grav'+str(round(gravity,2))+'_heightCap'+str(-heightCap)+'_timeCap'+str(timeCap)+'_dist'
for eta in etaList:
    xDist = []
    # for i in range(1000,1000+particles):
    csvName = r'D:\Users\janah\Desktop\Projects\GaltonBoard\PyCode\distData_eta' + str(eta) + '_timeCap' + str(timeCap)+'_particles' + str(particles) + '_grav' + str(gravity)+'.csv'
    # print(csvName)
    rf = open(csvName,'r')
    reader = csv.reader(rf)
    row = next(reader)
    count = 0
    for i in range(0,particles):
        row = next(reader)
        if row[10] == "False":
            count +=1
            xFrame = float(row[0])
            yFrame = float(row[1])
            vX = float(row[4])
            vY = float(row[5])
            # print(xFrame,yFrame,vX,vY)
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
            xDist.append(xFrame-vX*(trev-freeFall))
    rf.close()
    xDist.sort()
    totalDist.append(xDist)
    print(count,'/',particles)

sns.displot(totalDist, kind="kde", legend=False)
############################## Save of Show ###################################
plt.xlabel("x-coord")
plt.ylabel("proportrion")
plt.legend([round(e,2) for e in etaList])
plt.show()
# plt.axis('off')
# plt.savefig(fname+'.eps',transparent=True)
# plt.close('all')
