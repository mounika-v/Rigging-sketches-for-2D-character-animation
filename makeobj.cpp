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
	outdata.open("sem/mesh1.obj");
	
	int n=0,iter = 1;
	float x,y,z;

	fscanf(nodef,"%d %f %f %f",&n,&x,&y,&z);
	//printf(" n - %d  iter - %d\n",n,iter);
	outdata<<n<<endl;
	while(iter != n)
	{
		fscanf(nodef,"%d %f %f %f",&iter,&x,&y,&z);
		outdata<<"v "<<x<<" "<<y<<" "<<z<<endl;
	}
	
	n=0;
	iter = 1;
	
	fscanf(elef,"%d %f %f",&n,&x,&y);
	//printf(" n - %d  iter - %d\n",n,iter);
	outdata<<n<<endl;
	
	while(iter != n)
	{
		fscanf(elef,"%d %f %f %f",&iter,&x,&y,&z);
		outdata<<"f "<<x<<" "<<y<<" "<<z<<endl;
	}

	fclose(nodef);
	fclose(elef);
	outdata.close();
	return 0;
}
