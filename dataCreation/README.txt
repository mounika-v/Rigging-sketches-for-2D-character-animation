Files in this directory
**********************************
CMakeLists.txt - make file
LICENSE
main.cpp - main file to run the skeleton contraction
README.txt
tutorial_shared_path.h - defines variables
1d/ - Files to determine the 1d skeleton will be saved here under pre-determined name (to be renamed as required)
2d/ - Files to determine the skeleton to mesh related saved here under pre-determined names (to be renamed as required)
build/ - Program is run from inside this directory
cmake/ - Project support files
offiles/ - .off files that are created from makeoff.cpp are saved here to be picked up by the project.
************************************************************************************************************************************

Check the libigl tutorial to see how to create a basic project.
To run the project:
-------------------------
Go to build directory
./getskeleton_bin
Skeleton files are generated in 1d directory. Change names as required
To view this go to 1d directory and run ./draft1
Skinning details are generated in 2d directory. Change names as required
To view this go to 2d directory and run ./draft1

In 1d directory, the skeleton can be edited and the respective skinning can be saved to 2d
Edit options:
I - Idle
M - Merge
L - Add leaf
V - Move
A - Add a node between 2 nodes
(To be done: merge A,L to be a single option but do both)
S - Option S will save the current skeleton and replace the existing files. So, running ./draft1 will load the saved skeleton

To enable any option, make sure the mode is Idle first.


