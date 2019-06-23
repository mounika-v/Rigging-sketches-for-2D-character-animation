import argparse

import cv2
import numpy as np
import torch

from models.with_mobilenet import PoseEstimationWithMobileNet
from modules.keypoints import extract_keypoints, group_keypoints, BODY_PARTS_KPT_IDS, BODY_PARTS_PAF_IDS
from modules.load_state import load_state
from val import normalize, pad_width

imagename = ""

class ImageReader(object):
    def __init__(self, file_names):
        self.file_names = file_names
        self.max_idx = len(file_names)

    def __iter__(self):
        self.idx = 0
        return self

    def __next__(self):
        if self.idx == self.max_idx:
            raise StopIteration
        img = cv2.imread(self.file_names[self.idx], cv2.IMREAD_COLOR)
        if img.size == 0:
            raise IOError('Image {} cannot be read'.format(self.file_names[self.idx]))
        self.idx = self.idx + 1
        return img


class VideoReader(object):
    def __init__(self, file_name):
        self.file_name = file_name
        try:  # OpenCV needs int to read from webcam
            self.file_name = int(file_name)
        except ValueError:
            pass

    def __iter__(self):
        self.cap = cv2.VideoCapture(self.file_name)
        if not self.cap.isOpened():
            raise IOError('Video {} cannot be opened'.format(self.file_name))
        return self

    def __next__(self):
        was_read, img = self.cap.read()
        if not was_read:
            raise StopIteration
        return img


def infer_fast(net, img, net_input_height_size, stride, upsample_ratio, cpu,
               pad_value=(0, 0, 0), img_mean=(128, 128, 128), img_scale=1/256):
    height, width, _ = img.shape
    scale = net_input_height_size / height

    scaled_img = cv2.resize(img, (0, 0), fx=scale, fy=scale, interpolation=cv2.INTER_CUBIC)
    scaled_img = normalize(scaled_img, img_mean, img_scale)
    min_dims = [net_input_height_size, max(scaled_img.shape[1], net_input_height_size)]
    padded_img, pad = pad_width(scaled_img, stride, pad_value, min_dims)

    tensor_img = torch.from_numpy(padded_img).permute(2, 0, 1).unsqueeze(0).float()
    if not cpu:
        tensor_img = tensor_img.cuda()

    stages_output = net(tensor_img)

    stage2_heatmaps = stages_output[-2]
    heatmaps = np.transpose(stage2_heatmaps.squeeze().cpu().data.numpy(), (1, 2, 0))
    heatmaps = cv2.resize(heatmaps, (0, 0), fx=upsample_ratio, fy=upsample_ratio, interpolation=cv2.INTER_CUBIC)

    stage2_pafs = stages_output[-1]
    pafs = np.transpose(stage2_pafs.squeeze().cpu().data.numpy(), (1, 2, 0))
    pafs = cv2.resize(pafs, (0, 0), fx=upsample_ratio, fy=upsample_ratio, interpolation=cv2.INTER_CUBIC)

    return heatmaps, pafs, scale, pad


def run_demo(net, image_provider, height_size, cpu):
    net = net.eval()
    if not cpu:
        net = net.cuda()

    stride = 8
    upsample_ratio = 4
    color = [255, 0, 0]
    colorcircle = [0, 0, 255]
    for img in image_provider:
        orig_img = img.copy()
        heatmaps, pafs, scale, pad = infer_fast(net, img, height_size, stride, upsample_ratio, cpu)

        total_keypoints_num = 0
        all_keypoints_by_type = []
        for kpt_idx in range(18):  # 19th for bg
            total_keypoints_num += extract_keypoints(heatmaps[:, :, kpt_idx], all_keypoints_by_type, total_keypoints_num)

        pose_entries, all_keypoints = group_keypoints(all_keypoints_by_type, pafs, demo=True)

        for kpt_id in range(all_keypoints.shape[0]):
            all_keypoints[kpt_id, 0] = (all_keypoints[kpt_id, 0] * stride / upsample_ratio - pad[1]) / scale
            all_keypoints[kpt_id, 1] = (all_keypoints[kpt_id, 1] * stride / upsample_ratio - pad[0]) / scale

        global imagename
        imagepath = imagename
        if imagepath != "":
            outputname = (imagepath[0].split('/')[-1])
            outputname = "output/"+(outputname.split('.')[0])+".skl"
        skeletonfile = open(outputname,"w")
        skeletonfile.write(imagepath[0]+"\n")
        # Write all the keypoints to output
        for kpt in all_keypoints:
            skeletonfile.write(str(kpt[0])+","+str(kpt[1])+",0.0\n")

        for n in range(len(pose_entries)):
            if len(pose_entries[n]) == 0:
                continue
            ######### calculate root point
            bneck = 0.406394346548266
            bleft = 0.295839233078765
            bright = 0.297766420372969
            neckx,necky = all_keypoints[int(pose_entries[n][1]), 0:2] #neck vertex
            leftx,lefty = all_keypoints[int(pose_entries[n][11]), 0:2]#Left leg vertex
            rightx,righty = all_keypoints[int(pose_entries[n][8]), 0:2]#right leg vertex
            rootx = bneck*neckx + bleft*leftx + bright*rightx
            rooty = bneck*necky + bleft*lefty + bright*righty
            # Add root point to keypoints
            skeletonfile.write(str(rootx)+","+str(rooty)+",0.0\n")
            ########### Plot root point
            cv2.circle(img,(int(rootx),int(rooty)),8,colorcircle,-1)
            ############ Draw lines
            writeedges = []
            for part_id in range(len(BODY_PARTS_PAF_IDS) - 2):
                kpt_a_id = BODY_PARTS_KPT_IDS[part_id][0]
                global_kpt_a_id = pose_entries[n][kpt_a_id]
                if global_kpt_a_id != -1:
                    x_a, y_a = all_keypoints[int(global_kpt_a_id), 0:2]
                    cv2.circle(img, (int(x_a), int(y_a)), 8, colorcircle, -1)
                kpt_b_id = BODY_PARTS_KPT_IDS[part_id][1]
                global_kpt_b_id = pose_entries[n][kpt_b_id]
                if global_kpt_b_id != -1:
                    x_b, y_b = all_keypoints[int(global_kpt_b_id), 0:2]
                    cv2.circle(img, (int(x_b), int(y_b)), 8, colorcircle, -1)
                keypointsize = len(all_keypoints)
                if global_kpt_a_id != -1 and global_kpt_b_id != -1:
                    if kpt_a_id == 1 and (kpt_b_id == 8 or kpt_b_id == 11):
                        if (global_kpt_a_id,keypointsize) not in writeedges:
                            writeedges.append((global_kpt_a_id,keypointsize))
                        if (keypointsize,global_kpt_b_id) not in writeedges:
                            writeedges.append((keypointsize,global_kpt_b_id))
                        # skeletonfile.write(str(global_kpt_a_id)+","+str(keypointsize)+"\n")
                        # skeletonfile.write(str(keypointsize)+","+str(global_kpt_b_id)+"\n")
                        cv2.line(img,(int(x_a), int(y_a)),(int(rootx),int(rooty)),color,3)
                        cv2.line(img,(int(rootx),int(rooty)),(int(x_b), int(y_b)),color,3)
                    else:
                        if (global_kpt_a_id,global_kpt_b_id) not in writeedges:
                            writeedges.append((global_kpt_a_id,global_kpt_b_id))
                        # skeletonfile.write(str(global_kpt_a_id)+","+str(global_kpt_b_id)+"\n")
                        cv2.line(img, (int(x_a), int(y_a)), (int(x_b), int(y_b)), color, 3)

            # Writing the edges to file
            for eachi in writeedges:
                skeletonfile.write(str(eachi[0])+","+str(eachi[1])+"\n")

        skeletonfile.close()

        if imagename != "":
            imagename = (imagename[0].split('/')[-1])
            imagename = "output/"+imagename
            # print(imagename)
            # cv2.imwrite(imagename, img)
        else:
            imagename = "output/tempop.png"
        cv2.imwrite(imagename,img)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Lightweight human pose estimation python demo.
                       This is just for quick results preview.
                       Please, consider c++ demo for the best performance.''')
    parser.add_argument('--checkpoint-path', type=str, required=True, help='path to the checkpoint')
    parser.add_argument('--height-size', type=int, default=256, help='network input layer height size')
    parser.add_argument('--video', type=str, default='', help='path to video file or camera id')
    parser.add_argument('--images', nargs='+', default='', help='path to input image(s)')
    parser.add_argument('--cpu', action='store_true', help='run network inference on cpu')
    args = parser.parse_args()

    if args.video == '' and args.images == '':
        raise ValueError('Either --video or --image has to be provided')

    net = PoseEstimationWithMobileNet()
    checkpoint = torch.load(args.checkpoint_path, map_location='cpu')
    load_state(net, checkpoint)

    frame_provider = ImageReader(args.images)
    if args.images != '':
        imagename = args.images
    if args.video != '':
        frame_provider = VideoReader(args.video)

    run_demo(net, frame_provider, args.height_size, args.cpu)
