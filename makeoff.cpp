#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;

int main()
{
	FILE *nodef,*elef;
	nodef = fopen("triangle/contour.2.node","r");
	elef = fopen("triangle/contour.2.ele","r");
	
	ofstream outdata;
	outdata.open("sem/build/mesh1.off");
	outdata<<"OFF"<<endl;
	
	int n=0,ne=0,iter = 1;
	float x,y,z;

	fscanf(nodef,"%d %f %f %f\n",&n,&x,&y,&z);
	outdata<<n<<" ";
	fscanf(elef,"%d %f %f\n",&ne,&x,&y);
	outdata<<ne<<" "<<0<<endl;
	
	while(iter != n)
	{
		fscanf(nodef,"%d %f %f %f",&iter,&x,&y,&z);
		outdata<<x<<" "<<y<<" "<<z<<endl;
	}
	
	iter = 1;
	
	//fscanf(elef,"%d %f %f",&n,&x,&y);
	//printf(" n - %d  iter - %d\n",n,iter);
	//outdata<<n<<endl;
	
	while(iter != ne)
	{
		fscanf(elef,"%d %f %f %f",&iter,&x,&y,&z);
		outdata<<3<<" "<<x-1<<" "<<y-1<<" "<<z-1<<endl;
	}

	fclose(nodef);
	fclose(elef);
	outdata.close();
	return 0;
}
