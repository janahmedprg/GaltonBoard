import csv

particles = 5
timeStart =1000
timeEnd = 8000
increment =1000
gravity = 1
eta = 0


for timeCap in range(timeStart,timeEnd+increment,increment):
    gdata = [[],[],[],[],[],[],[],[],[],[],[]]
    for i in range(0,particles):
        csvName = 'gdata_eta_' + str(eta) + 'time_' + str(timeCap)+'_particle' + str(i) + '_grav' + str(gravity)+'.csv'
        rf = open(csvName,'r')
        reader = csv.reader(rf)

        row = next(reader)
        gdata[0].append(float(row[0]))
        gdata[1].append(float(row[1]))
        gdata[2].append(float(row[2]))
        gdata[3].append(float(row[3]))
        gdata[4].append(float(row[4]))
        gdata[5].append(float(row[5]))
        gdata[6].append(float(row[6]))
        gdata[7].append(int(row[7]))
        gdata[8].append(bool(row[8]))
        gdata[9].append(float(row[9]))
        gdata[10].append(float(row[10]))
        rf.close()

    csvName = 'gdata_'+'eta_'+str(eta)+'time_'+str(timeCap)+'_particles'+str(particles)+'_grav'+str(gravity)+'.csv'
    f = open(csvName, 'w', newline='')
    writer = csv.writer(f)
    writer.writerow(['xFrame','yFrame','pX','pY','vX','vY','vS','wall','isTorus','time','particle'])
    for xframe,yframe,px,py,vx,vy,vs,wall,istorus,time,particle_idx in zip(gdata[0],gdata[1],gdata[2],gdata[3],gdata[4],gdata[5],gdata[6],gdata[7],gdata[8],gdata[9],gdata[10]):
        writer.writerow([xframe,yframe,px,py,vx,vy,vs,wall,istorus,time,particle_idx])
    f.close()
