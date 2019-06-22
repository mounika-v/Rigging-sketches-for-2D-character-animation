#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

#include "CImg.h"
using namespace cimg_library;

bool entryexists(int curri,int currj,vector <int> contouri,vector <int> contourj);
void turn(int &dir, int &curri, int &currj, char turnto);
void fillneighbor(int dir,int &initi, int &initj, vector<int> &ni, vector<int> &nj);
void finddirection(int idiff,int jdiff,int &dir);

int main()
{
	CImg<unsigned char> isrc("sketch.png");//,512,512,1,3);

	int done,dir;
	vector <int> contouri,contourj;

	int width = isrc.width();
	int height = isrc.height();
	int bc = 0,c,r;
	bool breakflag = false;

	cout << width << "x" << height << endl;
	cout << isrc.depth() << " , " << isrc.spectrum()<<endl;

	CImg<unsigned char> src = isrc.get_resize(512,512,1,3);
	width = src.width();
	height = src.height();

	cout << width << "x" << height << endl;
	cout << src.depth() << " , " << src.spectrum()<<endl;


	src.channel(0);
	for (r = height-1 ; r>=0 ; r--)
	{
		for (c = 0; c < width; c++)
		{
			// cout<<r<<"  "<<c<<"   "<<(int)src(c,r,0,0)<<endl;
			if((int)src(c,r,0,0) < 175 )//&& (int)src(c,r,0,1) == 0 && (int)src(c,r,0,2) == 0)
			{
				breakflag = true;break;
			}
		}
		if(breakflag)
			break;
	}
	cout<<c<<" "<<r<<endl;
	// cout<<"black pixel count: "<<bc<<endl;

	/*---------------------------------- moores contour --------------------------------*/
	dir=2;
	done = 2;
	int initi=r,initj=c;
	vector<int> ni,nj;
	int idiff, jdiff,currc,currr;

	// done = 0;

	while(done>0)
	{
		if(!entryexists(initi,initj,contouri,contourj)) //If the contour pixel doesn't exist in the vector, push it.
		{
			contouri.push_back(initi);
			contourj.push_back(initj);
		}
		fillneighbor(dir,initi,initj,ni,nj); //Fill the neighbourhood vector in clockwise direction

		vector <int> :: iterator iiter,jiter;
		// cout<<"neighbor size: "<<ni.size();
		for(iiter = ni.begin(),jiter = nj.begin(); iiter != ni.end() ; iiter++,jiter++) //iterate through neighbourhood points
		{
			//cout<<*iiter<<" "<<*jiter<<endl;
			currc = *jiter; currr = *iiter;
			if((int)src(currc,currr,0,0) < 250) // && (int)src(currc,currr,0,1) == 0 && (int)src(currc,currr,0,2) == 0)
			{
				// cout<<"  -- met a black pixel"<<endl;
				initi = *iiter; initj = *jiter;
				break;
			}
		}
		if(initi == r && initj == c)
		{
			done--;
		}
		idiff = *iiter - *(iiter-1);
		jdiff = *jiter - *(jiter-1);
		finddirection(idiff,jdiff,dir);
	}
	cout<<"Contour size: "<<contouri.size()<<endl;

	/*------------------------------ Making the .poly file ---------------------------------*/

	vector <int> :: iterator ciiter,cjiter;
	ofstream outdata;
	outdata.open("triangle/contour.poly");
	if(!outdata)
	{
		cout<<"Error opening file"<<endl;return 0;
	}

	int csize = contouri.size();
	if(csize % 5 == 0 || (csize -1)%5 == 0)
		csize = (csize/5) + 1;
	else
		csize = (csize/5) + 2;

	//Start writing the file
	outdata<<csize<<" "<<2<<endl;
	int count=1,entryno=1;
	for(ciiter = contouri.begin(),cjiter=contourj.begin();ciiter!=contouri.end();ciiter++,cjiter++,count++)
	{
		if(count%5 == 1 || ciiter == contouri.end()-1)
		{
			outdata<<entryno<<" "<<*cjiter<<" "<<*ciiter<<endl;entryno++;
		}
	}
	//Adding segment data to poly file
	count = 1;
	outdata<<csize<<" "<<0<<endl;
	while(count <= csize)
	{
		outdata<<count<<" "<<((count+(csize-2))%csize)+1<<" "<<count<<endl;
		count++;
	}
	outdata<<0<<endl;
	outdata.close();

	//Silhouette image

	// CImg<unsigned char> silho(width,height,1,1,0);
	// const unsigned char black[1] = {0},white[3]={255,255,255};
	// silho.draw_fill(0,0,white);
	// for(ciiter = contouri.begin(),cjiter=contourj.begin();(ciiter+1)!=contouri.end();ciiter++,cjiter++)
	// {
	// 	cout<<*cjiter<<","<<*ciiter<<","<<*(cjiter+1)<<","<<*(ciiter+1)<<endl;
	// 	silho.draw_line(*cjiter,*ciiter,*(cjiter+1),*(ciiter+1),black);
	// }
	//
	// silho.draw_fill(c+1,r+1,black);
	// silho.save("silho.png");

	return 0;
}


void fillneighbor(int dir,int &initi, int &initj, vector<int> &ni, vector<int> &nj)
{
	int xd=0,yd=0;
	switch(dir)
	{
			case 0: xd = -1; yd = 0;break;
			case 1: xd = 0; yd = +1;break;
			case 2: xd = +1; yd = 0;break;
			case 3: xd = 0; yd = -1;break;
	}
	int x=initi+xd,y=initj+yd;
	for(int iter=-1;iter<8;)
	{
		while(((x-initi <= 1)&&(x-initi >= -1))&&((y-initj <= 1)&&(y-initj >= -1)))
		{
			iter++;
			if(iter>=8)
				break;
			if(ni.size() != 8)
			{
				ni.push_back(x);
				nj.push_back(y);
			}
			else
			{
				ni[iter]=x;nj[iter]=y;
			}
			x=x+xd;
			y=y+yd;
		}
		x-=xd;y-=yd;
		turn(dir,x,y,'r');
		switch(dir)
		{
			case 0: xd = -1; yd = 0;break;
			case 1: xd = 0; yd = +1;break;
			case 2: xd = +1; yd = 0;break;
			case 3: xd = 0; yd = -1;break;
		}
	}
}

void turn(int &dir, int &x, int &y, char turnto)
{
	dir = (dir+1+4)%4;
	switch(dir)
	{
		case 0:x--;break;
		case 1:y++;break;
		case 2:x++;break;
		case 3:y--;break;
	}
}

bool entryexists(int curri,int currj,vector <int> contouri,vector <int> contourj)
{
	vector <int> :: iterator i;
	for(i=contouri.begin();i != contouri.end(); i++)
	{
		if(*i == curri)
			if(currj == contourj.at(i-contouri.begin()))
				return true;
	}
	return false;
}

void finddirection(int idiff,int jdiff,int &dir)
{
	if(idiff == 0 && jdiff == 1)
		dir = 3;
	else if(idiff == 0 && jdiff == -1)
		dir = 1;
	else if(idiff == 1 && jdiff == 0)
		dir = 0;
	else if(idiff == -1 && jdiff == 0)
		dir = 2;
}

//6038
