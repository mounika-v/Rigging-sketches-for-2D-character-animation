import os
mycmd = 'ls sketch*.png > imlist.txt'
os.system(mycmd)
rf = open("imlist.txt","r")
for name in rf:
    temp = name.split(".")
    name = temp[0]
    print("reading "+name);
    wf = open("imgchange.cpp","w")
    wf.write("#include<iostream>\n")
    wf.write("#include<fstream>\n")
    wf.write("#include \"CImg.h\"\n")
    wf.write("using namespace std;\n")
    wf.write("using namespace cimg_library;\n")
    wf.write("int main()\n{\n")
    wf.write("CImg<unsigned char> isrc(\""+name+".png\");\n")
    wf.write("CImg<unsigned char> dest = isrc.get_resize(512,512,1,3);\n")
    wf.write("dest.save_png(\"finalimages/"+name+".png\",0);\n}")
    wf.close();
    # mycmd = 'cat imgchange.cpp';
    # os.system(mycmd)
    mycmd = 'g++ -o imgchange imgchange.cpp -lX11 -lpthread'
    os.system(mycmd)
    mycmd = './imgchange'
    os.system(mycmd)
    print("saved: "+name+"\n")
rf.close()
os.remove("imlist.txt")
