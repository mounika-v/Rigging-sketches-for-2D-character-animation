# Parsing-sketches-for-2D-character-animation
Rigging is a crucial step in skeletal animation.  However, the current method doesn’tgive the ease of doing it automatically. It involves manual effort to perform rigging andskinning. We present a method to perform auto-rigging for 2D animation by parsing the2D model sketches of a character. The method we are suggesting will take an image ofthe character’s model sketch in T-pose as input and generates a skeleton and skinninginformation of the character for the animator to take it further. We adopted a 3D meshcontraction method for 2D mesh to achieve this goal. Using this data we fine-tuned an available pose estimation model to work for sketches as well. We are able to give a decent initialization which can be edited further.

Input:
![input](https://www.cse.iitb.ac.in/~mounika/finalmtp/ip/demo5.png) 

Output:
![output](https://www.cse.iitb.ac.in/~mounika/finalmtp/op/demo5.png)

### Steps to run mesh contraction method:
1. Download [libigl](https://github.com/libigl/libigl) and [triangle ](https://www.cs.cmu.edu/~quake/triangle.html)

2. Extract them in the parent folder and install them as instructed in their respective homepages

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
sem/build ./getskeleton_bin
sem/1d ./draft1
sem/2d ./draft1
```
3. To generate the json file for training:
```bash
python3 getcontours.py
python3 genannotations.py
```
4. Copy the json file to model directory

### Instructions for model:
1. Go to model directory

2. Install all the requirements
```bash
pip -r requirements.txt
```
3. For training the model follow the instructions from [here](https://github.com/Daniil-Osokin/lightweight-human-pose-estimation.pytorch)

4. For running the demo for individual image:
```bash
./rundemo.sh <path to image>
```
5. For running on the GUI tool:
```bash
python3 rigit.py

```
A complete report that was submitted towards the fulfilment of program can be found [here](https://www.cse.iitb.ac.in/~mounika/finalmtp/report.pdf) and links to bases papers on the [project home](https://www.cse.iitb.ac.in/~mounika/finalmtp)
