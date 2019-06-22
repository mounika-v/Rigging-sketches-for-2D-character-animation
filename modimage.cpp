#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

#include "CImg.h"
using namespace cimg_library;

int main(int argc,char *args[])
{
    string imgname;
    while(getline(cin, imgname))
    {
      cout<<"received: "<<imgname<<"   ";
      const char* dum = imgname.c_str();
      CImg<unsigned char> isrc(dum);
      // CImg<unsigned char> dest = isrc.get_resize(512,512,1,3);
      imgname = "finalimages/"+imgname;
      char const *dummy = imgname.data();
      cout<<"saving: "<<dummy<<endl;
      // dest.save(dummy);
    }
}
