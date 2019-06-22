import os
mycmd = "ls skeleton*.raw > rawlist.txt"
os.system(mycmd)
rf = open("rawlist.txt","r")
wf = open("rootfile.csv","w")
# For all the skeleton raw files
for line in rf:
    line = line.split(".")
    filename = line[0]+".raw"
    print("reading: "+filename)
    rawfile = open(filename,"r")
    points = []
    # Open the file and read the coordinates
    for point in rawfile:
        coords = point.split(' ')
        points.append((float(coords[0]),float(coords[1])))
    rawfile.close()
    # print(points)
    # Iterate through all the coordinates and
    for i in range(0,len(points)):
        # Look for the one that has occured 3 times (only root has degree 3)
        pcount=points.count(points[i])
        # print(i,"     ",pcount)
        if pcount == 3:
            towrite = filename+","+str(points[i][0])+","+str(points[i][1])
            break;
    wf.write(towrite)
    wf.write("\n")
rf.close()
wf.close()
