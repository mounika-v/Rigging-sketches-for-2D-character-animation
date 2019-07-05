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

/*Implementing moores neighbourhood algorithm to extract the contour of the sketch*/
int main()
{
	/*To be implemented: To take any sketch as an input from arguments.
	The current implementation takes a sketch "sketch.png" always. so, before running this,
	rename the sketch to sketch.png"*/
	CImg<unsigned char> isrc("sketch.png");//,512,512,1,3);

	/* dir variable determines the direction of parsing after detecting a valid contour pixel */
	int done,dir;
	/*Contouri and contourj are the x,y coordinates of the contour*/
	vector <int> contouri,contourj;

	int width = isrc.width();
	int height = isrc.height();

	int bc = 0,c,r;
	bool breakflag = false;

	/* To resize the image for creating training data with images of constant size */
	CImg<unsigned char> src = isrc.get_resize(512,512,1,3);
	width = src.width();
	height = src.height();

	cout << width << "x" << height << endl;
	cout << src.depth() << " , " << src.spectrum()<<endl;

	/*Change it to a grey scale image*/
	src.channel(0);
	/*Detech the first valid pixel on the imagge scanning from bottom left*/
	for (r = height-1 ; r>=0 ; r--)
	{
		for (c = 0; c < width; c++)
		{
			if((int)src(c,r,0,0) < 175 )//&& (int)src(c,r,0,1) == 0 && (int)src(c,r,0,2) == 0)
			/*If the pixel value is less than 175 we consider it a valid pixel*/
			{
				breakflag = true;break;
			}
		}
		if(breakflag)
			break;
	}

	/*---------------------------------- moores contour --------------------------------*/
	dir=2;
	done = 2; /*Iteration counter for the no. of times the contour is parsed to make sure it is not stuck in a neighbourhood of single pixel */
	int initi=r,initj=c;
	vector<int> ni,nj;
	int idiff, jdiff,currc,currr;

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
			currc = *jiter; currr = *iiter; //Current column and row
			// done = 0;

			if((int)src(currc,currr,0,0) < 250) // && (int)src(currc,currr,0,1) == 0 && (int)src(currc,currr,0,2) == 0)
			//Any pixel that has grey value <250 because contours have very light areas.
			{
				initi = *iiter; initj = *jiter;
				break;
			}
		}
		if(initi == r && initj == c) //When the current pixel is same as the initial ones one iteration is done.
		{
			done--;
		}
		idiff = *iiter - *(iiter-1);
		jdiff = *jiter - *(jiter-1);
		/*Find the direction in which we are parsing currently to determine the direction in which we have to parse for neighbourhood*/
		finddirection(idiff,jdiff,dir);
	}
	cout<<"Contour size: "<<contouri.size()<<endl;

	/*------------------------------ Making the .poly file ---------------------------------*/
	/*This holds the countour of the image to parse it for triangulation*/

	vector <int> :: iterator ciiter,cjiter;
	ofstream outdata;
	outdata.open("triangle/contour.poly");
	if(!outdata)
	{
		cout<<"Error opening file"<<endl;return 0;
	}

	/*Traingulation needs a polygon with an edge big enough for triangulation. So, we can't use all the pixels on contour.
	We have taken pixels at a regular interval of 5(to keep a dense mesh). Can be changed as required.*/
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

	return 0;
}


void fillneighbor(int dir,int &initi, int &initj, vector<int> &ni, vector<int> &nj)
/*Find all the neighbourhood pixels of the current pixel*/
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
/*Change the direction to go around in clockwise direction*/
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
/*Check if the current valid pixel is already added to contour*/
/*To be done: Change the implementation to save ordered pairs instead of two different vectors
So, this could be avoided*/
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
/*To determine the initial direction of parsing based on
difference in x,y of current and previous contour pixel.
0 - right
1 - down
2 - left
3 - up
*/
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
