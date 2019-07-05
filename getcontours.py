# Created to get the bounding box data and segmentation data of each sketch to Created
# data for training

import os
mycmd = 'ls sketch*.png > imlist.txt'
os.system(mycmd)
rf = open("imlist.txt","r")
for name in rf:
    temp = name.split(".")
    name = temp[0]
    mycmd = 'mv '+name+'.png sketch.png'
    # print(mycmd)
    os.system(mycmd)
    mycmd = './trace2'
    # print(mycmd)
    os.system(mycmd)
    temp = name
    temp = temp.replace('sketch','')
    mycmd = 'cp triangle/contour.poly contours/contour'+temp+'.poly'
    # print(mycmd)
    os.system(mycmd)
    mycmd = 'mv sketch.png '+name+'.png';
    # print(mycmd)
    os.system(mycmd)
rf.close();
os.remove("imlist.txt");
