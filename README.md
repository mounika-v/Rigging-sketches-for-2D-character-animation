# Parsing-sketches-for-2D-character-animation
Rigging is a crucial step in skeletal animation.  However, the current method doesn’tgive the ease of doing it automatically. It involves manual effort to perform rigging andskinning. We present a method to perform auto-rigging for 2D animation by parsing the2D model sketches of a character. The method we are suggesting will take an image ofthe character’s model sketch in T-pose as input and generates a skeleton and skinninginformation of the character for the animator to take it further. We adopted a 3D meshcontraction method for 2D mesh to achieve this goal.

# Steps to run mesh contraction method
Download libigl and triangle 

Extract them in the parent folder and install them as instructed in their respective homepages

```bash
g++ -o trace2 cimgdemo.cpp -lX11 -lpthread 
./trace2
cd triangle
triangle -p contour
triangle -rpa100 contour.1
showme contour.2.ele &
cd ..
g++ -o makeoff makeoff.cpp
./makeoff
sem/build ./example_bin
sem/1d ./draft1
sem/2d ./draft1
```
To generate the json file for training:
```bash
python3 getcontours.py
python3 genannotations.py
```
# Instructions for model.
go to model directory

Install all the requirements
```bash
pip -r requirements.txt
```
For training the model follow the instructions from here

For running the demo for individual image:
```bash
./rundemo.sh <path to image>
```
For running on the GUI tool:
```bash
python3 rigit.py
```
