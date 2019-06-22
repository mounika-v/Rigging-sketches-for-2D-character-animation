import os
import re
import numpy as np
from scipy.linalg import solve
from shapely.geometry import Polygon
mycmd = 'ls sketch*.png > imlist.txt'
os.system(mycmd)
rf = open("imlist.txt","r")
# bfile = open("barylist.csv","w")
annid = 0;
num_keypoints=13;
annotations = []
images = []
for name in rf:
    temp = name.split('.')
    name = temp[0]
    print("Prepping ",name)
    temp = name
    temp = temp.replace('sketch','')

    #################################################################################################
    #for images attribute
    #------------------------------------------------------------------------------------------------
    license = 4
    file_name = name+".png"
    coco_url = "https://www.cse.iitb.ac.in/~mounika/dataset/images/"+file_name
    height = 512
    width = 512
    date_captured = "2019-05-01 00:00:00"
    id = temp

    image = {
        "license":license,
        "file_name":file_name,
        "coco_url":coco_url,
        "height":height,
        "width":width,
        "date_captured":date_captured,
        "id":id
    }

    images.append(image)
    # ----------------------------------------------------------------------------------------------
    #################################################################################################
    #################################################################################################
    # For annotations
    # ----------------------------------------------------------------------------------------------

    #Standard variables and initialization
    lx=512
    ly=512
    hx=0
    hy=0
    pixelc = 0
    vertexlist = []

    # Read segmentation from contour file reducecd to the (no. of vertices / 3)
    # ----------------------------------------------------------------------------------------------
    segmentation = []
    cf = open("contours/contour"+temp+".poly","r")

    for x in cf:
        coords = x.split(' ')
        if len(coords) == 3 and pixelc%3 == 0:
            coords[1] = float(coords[1])
            coords[2] = float(coords[2])
            segmentation.append(coords[1])
            segmentation.append(coords[2])
            vertexlist.append((coords[1],coords[2]))
            # For box dimensions
            if(lx > coords[1]):
                lx = coords[1]
            elif(hx < coords[1]):
                hx = coords[1]
            if(ly > coords[2]):
                ly = coords[2]
            elif(hy < coords[2]):
                hy = coords[2]
        pixelc+=1
    cf.close()
    segmentations = []
    imgpoly = Polygon(vertexlist)
    # print(imgpoly.area)
    segmentations.append(segmentation)

    iscrowd = 0

    # Creating key points list
    # ---------------------------------------------------------------------------------------------------------
    joints = []
    sortlist = []
    kpf = open("sem/2d/joints"+temp+".txt","r")
    dummyindex = 0
    for kp in kpf:
        kplist = kp.split(' ')
        joints.append((float(kplist[0]), float(kplist[1])))
        sortlist.append((float(kplist[1]),float(kplist[0]),dummyindex))
        dummyindex+=1
    sortlist.sort()
    kpf.close()

    neighborhood = []
    for i in range(0,15):
        neighborhood.append([])
    kpf = open("sem/1d/skeleton"+temp+".raw","r")
    evenflag = True
    firstpoint = (0,0)
    secondpoint = (0,0)
    for edl in kpf:
        edp = edl.split(' ')
        if(evenflag):
            firstpoint = (float(edp[0]),float(edp[1]))
            evenflag = False
        else:
            secondpoint = (float(edp[0]),float(edp[1]))
            jointsindexfirst = joints.index(firstpoint)
            jointsindexsecond = joints.index(secondpoint)
            neighborhood[jointsindexfirst].append(jointsindexsecond)
            neighborhood[jointsindexsecond].append(jointsindexfirst)
            evenflag = True
        # print("flag: ",evenflag,"  first point: ",firstpoint," secondpoint: ",secondpoint)

    # dfsqueue = []
    donequeue = []
    pointlabel = [""] * 15
    # dfsqueue.append(sortlist[0][2])
    # while len(dfsqueue) > 0:
    #     nodeindex = dfsqueue[0]
    #     neighs = []
    #     for neighbors in neighborhood[nodeindex]:
    #         if (len(donequeue)>0) and !(joints[neighbors] in donequeue):
    #             neighs.append(joints[neighbors])
    #     neighs.sort()
    def dfs(nodeindex):
        donequeue.append(nodeindex)
        neighs = []
        neighood = neighborhood[nodeindex]
        for neighbors in neighood:
            if(len(donequeue) > 0) and not(neighbors in donequeue):
                neighs.append((joints[neighbors][0],joints[neighbors][1]))
        neighs.sort()
        for neigh in neighs:
            dfs(joints.index(neigh))
    dfs(sortlist[0][2])
    # print(donequeue)

    labels = ["nose","neck","left_shoulder","left_elbow","left_wrist","root","left_hip","left_knee","left_ankle","right_hip","right_knee","right_ankle","right_shoulder","right_elbow","right_wrist"]
    originals = ["nose","left_eye","right_eye","left_ear","right_ear","left_shoulder","right_shoulder","left_elbow","right_elbow","left_wrist","right_wrist","left_hip","right_hip","left_knee","right_knee","left_ankle","right_ankle"]
    # ksd = 0
    # while ksd<=14:
    #     pointlabel[sortlist[ksd][2]] = labels[ksd]
    #     ksd+=1
    # # print(pointlabel)
    # ###########################################
    # Solving barycentric eq for Nodes
    #############################################
    # rootlabelindex = labels.index("root")
    # rootindex = donequeue[rootlabelindex]
    # rootx = joints[rootindex][0]
    # rooty = joints[rootindex][1]
    # necklabelindex = labels.index("neck")
    # neckindex = donequeue[necklabelindex]
    # neckx = joints[neckindex][0]
    # necky = joints[neckindex][1]
    # leftlabelindex = labels.index("left_hip")
    # leftindex = donequeue[leftlabelindex]
    # leftx = joints[leftindex][0]
    # lefty = joints[leftindex][1]
    # rightlabelindex = labels.index("right_hip")
    # rightindex = donequeue[rightlabelindex]
    # rightx = joints[rightindex][0]
    # righty = joints[rightindex][1]
    # A = np.array([[neckx,leftx,rightx],[necky,lefty,righty],[1,1,1]])
    # B = np.array([[rootx],[rooty],[1]])
    # x = solve(A,B)
    # writetobfile = file_name+","+str(x[0][0])+","+str(x[1][0])+","+str(x[2][0])+"\n"
    # bfile.write(writetobfile)
    ################################################
    keypoints = []
    ksd = 0
    while ksd <= 16:
        if(originals[ksd] in labels):
            target = labels.index(originals[ksd])
            jointindex = donequeue[target]
            keypoints.append(joints[jointindex][0])
            keypoints.append(joints[jointindex][1])
            keypoints.append(2)
        else:
            keypoints.append(0)
            keypoints.append(0)
            keypoints.append(0)
        ksd+=1
    image_id = temp
    category_id = 1
    bbox = [lx,ly,hx-lx,hy-ly]
    annotation = {
        "segmentation":segmentations,
        "num_keypoints":num_keypoints,
        "area":imgpoly.area,
        "iscrowd": iscrowd,
        "keypoints":keypoints,
        "image_id":image_id,
        "bbox":bbox,
        "category_id":category_id,
        "id":annid
    }
    annotations.append(annotation)
    annid+=1
    ##################################################################################################################
rf.close()
################# Writing data to json #############################
jfile = open("datahftext/headerdata.txt","r")
header = ""
for k in jfile:
    header += k
jfile.close()
jfile = open("datahftext/footerdata.txt","r")
footer = ""
for k in jfile:
    footer += k
jfile.close()
jfile = open("newdata.json","w+")
towrite = header
images = str(images).replace("\'","\"")
towrite += images
towrite += ",\"annotations\":"
annotations = str(annotations).replace("\'","\"")
towrite += annotations
towrite += footer
print(towrite)
jfile.write(towrite)
jfile.close()
# bfile.close() #File used for calculating barycentric coordinate values
# print(annotations)
# print(images)
