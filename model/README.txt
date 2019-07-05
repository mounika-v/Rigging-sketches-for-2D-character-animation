Contents of this folder
********************************************************
datasets
default_checkpoints
models
modules
new_checkpoints
output
_pycache_
scripts
sketches
demo.py
hi.png
LICENSE
newdata.json
prepared_train_annotation.json
README.txt
requirements.txt
rigit.py
rundemo.sh
train.py
val.py
*******************************************************
Most of the above contents are directly taken from the code made available by the publishers of base paper.
Please check git page of that project to get scripts for training and validating the project. 

new_checkpoints is the project where we have saved our new checkpoint after fine-tuning the existing model
The output of demo.py is saved to the output directory
sketches folder holds the samples sketches we have used to test the model
newdata.json is the new data we have used for finetuning the model. The data is similar to COCO Dataset.
train and val files are not changed.

Some changes are made to demo.py to accommodate our requirements. (Check the demo.py for help text of source code)
rigit.py shows a GUI tool that can be used to make modifications for the predicted skeleton.

rundemo.sh is to run the demo without GUI by passing path to any sketch.

