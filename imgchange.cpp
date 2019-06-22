#include<iostream>
#include<fstream>
#include "CImg.h"
using namespace std;
using namespace cimg_library;
int main()
{
CImg<unsigned char> isrc("sketch9.png");
CImg<unsigned char> dest = isrc.get_resize(512,512,1,3);
dest.save_png("finalimages/sketch9.png",0);
}