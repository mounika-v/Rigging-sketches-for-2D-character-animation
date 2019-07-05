/*
  Modified from An Introduction to OpenGL Programming,
  Ed Angel and Dave Shreiner, SIGGRAPH 2013
  Modified from 03_colorcube_rotate tutorial by Parag Chaudhuri, 2015
*/


#include "draft1.hpp"
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

GLuint shaderProgram;
GLuint vbo, vao;

glm::mat4 rotation_matrix;
glm::mat4 ortho_matrix;
glm::mat4 modelview_matrix;
GLuint uModelViewMatrix;

//For glu project and gluUnProject
GLdouble pos3D_x, pos3D_y, pos3D_z;
GLdouble model_view[16];
GLdouble projection[16];
GLint viewport[4];

//Skeleton data storage
vector<vector<int>> history,E,E2,stage2junct;
vector<int> stage2nodes;
vector<vector<double>> V,U;
double xpos1,ypos1,xpos2,ypos2;
int fetchedpoints = 0;
int scaleby = 200;
//-----------------------------------------------------------------

//6 faces, 2 triangles/face, 3 vertices/triangle
const int num_vertices = 5000;

int tri_idx=0;
bool fetched_first=false;
int srcnodeindex=0,destnodeindex=0;
glm::vec4 v_positions[num_vertices];
glm::vec4 v_colors[num_vertices];

glm::vec4 avail_colors[4] = {
    glm::vec4(1.0f,0.0f,0.0f,1.0f),
    glm::vec4(0.0f,1.0f,0.0f,1.0f),
    glm::vec4(0.0f,0.0f,1.0f,1.0f),
    glm::vec4(1.0f,1.0f,0.0f,1.0f)
};

int getvpositionsize()
{
    int count = 0;
    for(int i=0;i<num_vertices;i++)
    {
       if(v_positions[i].x != 0 || v_positions[i].y != 0)
        count++;
        else
        break;
    }
    return count;
}

void loadvertexdata()
{
    string inpline;
    ifstream infile("vertexdata.txt");
    getline(infile,inpline);
    istringstream iss1(inpline);
    int i;
    iss1 >> i;

    double x,y;
    while(i>0)
    {
        getline(infile,inpline);
        istringstream iss1(inpline);
        iss1 >> x >> y;
        vector<double> temp;
        temp.push_back(x);
        temp.push_back(y);
        V.push_back(temp);
        i--;
    }

    getline(infile,inpline);
    istringstream iss2(inpline);
    iss2 >> i;

    while(i>0)
    {
        getline(infile,inpline);
        istringstream iss2(inpline);
        iss2 >> x >> y;
        vector<double> temp;
        temp.push_back(x);
        temp.push_back(y);
        U.push_back(temp);
        i--;
    }

    getline(infile,inpline);
    istringstream iss3(inpline);
    iss3 >> i;

    int dumx,dumy;
    while(i>0)
    {
        getline(infile,inpline);
        istringstream iss3(inpline);
        iss3 >> dumx >> dumy;
        vector<int> temp;
        temp.push_back(dumx);
        temp.push_back(dumy);
        E.push_back(temp);
        i--;
    }

    getline(infile,inpline);
    istringstream iss4(inpline);
    iss4 >> i;
    while(i>0)
    {
        getline(infile,inpline);
        istringstream iss4(inpline);
        iss4 >> dumx >> dumy;
        vector<int> temp;
        temp.push_back(dumx);
        temp.push_back(dumy);
        E2.push_back(temp);
        i--;
    }
}

void loadSkeletondata()
{
    string inpline;
    ifstream infile("skeldata.txt");
    getline(infile,inpline);
    istringstream iss(inpline);
    int i;
    iss >> i;
    getline(infile, inpline);
    istringstream iss1(inpline);
    int x;
    while(i>0)
    {
        iss1 >> x;
        stage2nodes.push_back(x);
        vector<int> dummy;
        stage2junct.push_back(dummy);
        i--;
    }

    for(int i=0;i<E.size();i++)
    {
        if(!(E[i][0]==0 && E[i][1]==0))
        {
            int ind = find(stage2nodes.begin(),stage2nodes.end(),E[i][0]) - stage2nodes.begin();
            if(stage2junct.at(ind).size() == 0 || find(stage2junct.at(ind).begin(), stage2junct.at(ind).end(), E[i][1]) == stage2junct.at(ind).end())
            {
                stage2junct.at(ind).push_back(E[i][1]);
            }
            ind = find(stage2nodes.begin(),stage2nodes.end(),E[i][1]) - stage2nodes.begin();
            if(stage2junct.at(ind).size() == 0 || find(stage2junct.at(ind).begin(), stage2junct.at(ind).end(), E[i][0]) == stage2junct.at(ind).end())
            {
                stage2junct.at(ind).push_back(E[i][0]);
            }
        }
    }

    getline(infile, inpline);
    istringstream iss2(inpline);
    iss2 >> i;
    int nodex;
    while(i>0)
    {
        int hsize;
        getline(infile,inpline);
        istringstream iss3(inpline);
        iss3 >> hsize;
        getline(infile,inpline);
        istringstream iss4(inpline);
        vector<int> temp;
        while(hsize>0)
        {
            iss4 >> nodex;
            temp.push_back(nodex);
            hsize--;
        }
        history.push_back(temp);
        temp.size();
        temp.clear();
        i--;
    }
}

void loadSkeleton(void)
{
    FILE *fp;
    int i;
    double x,y,z;
    // fp = fopen("stdskel.raw","r");
    fp = fopen("skeleton.raw","r");
    // fp = fopen("stand.raw","r");
    if(fp == NULL)
    {
        printf("Unable to open file for reading\n");
        exit(1);
    }
    i=0;
    while(fscanf(fp,"%lf %lf %lf\n",&x,&y,&z)==3)
    {
        // v_positions[i] = glm::vec4(x/10,y/10,z/10,1.0f);
        v_positions[i]=glm::vec4(x/scaleby - 1.0, -y/scaleby + 1.0,z/scaleby,1.0f);
        //Scaleby because the window size could be smaller than the sketch skeleton size
        // v_positions[i]=glm::vec4(6*x/100, 6*y/100, 6*z/100, 1.0f);
        v_colors[i] = avail_colors[i%4];//glm::vec4(1.0f,1.0f,0.0f,1.0f);
        i++;
    }
    fclose(fp);
    // num_vertices = i;
}

void renderGL(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  rotation_matrix = glm::rotate(glm::mat4(1.0f), xrot, glm::vec3(1.0f,0.0f,0.0f));
  rotation_matrix = glm::rotate(rotation_matrix, yrot, glm::vec3(0.0f,1.0f,0.0f));
  rotation_matrix = glm::rotate(rotation_matrix, zrot, glm::vec3(0.0f,0.0f,1.0f));
  ortho_matrix = glm::ortho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);

  modelview_matrix = ortho_matrix * rotation_matrix;

  glUniformMatrix4fv(uModelViewMatrix, 1, GL_FALSE, glm::value_ptr(modelview_matrix));
  // Draw
  glDrawArrays(GL_LINES, 0, num_vertices);
  glDrawArrays(GL_POINTS, 0, num_vertices);

}

void changeObjects(){

  // //Copy the points into the current buffer
  // glBufferData (GL_ARRAY_BUFFER, sizeof (v_positions) + sizeof(v_colors), NULL, GL_DYNAMIC_DRAW);
  // glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(v_positions), v_positions );
  // glBufferSubData( GL_ARRAY_BUFFER, sizeof(v_positions), sizeof(v_colors), v_colors );


  //Ask GL for a Vertex Attribute Object (vao)
  glGenVertexArrays (1, &vao);
  //Set it as the current array to be used by binding it
  glBindVertexArray (vao);

  //Ask GL for a Vertex Buffer Object (vbo)
  glGenBuffers (1, &vbo);
  //Set it as the current buffer to be used by binding it
  glBindBuffer (GL_ARRAY_BUFFER, vbo);
  //Copy the points into the current buffer
  glBufferData (GL_ARRAY_BUFFER, sizeof (v_positions) + sizeof(v_colors), NULL, GL_DYNAMIC_DRAW);
  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(v_positions), v_positions );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(v_positions), sizeof(v_colors), v_colors );

  // Load shaders and use the resulting shader program
  std::string vertex_shader_file("03_vshader.glsl");
  std::string fragment_shader_file("03_fshader.glsl");

  std::vector<GLuint> shaderList;
  shaderList.push_back(csX75::LoadShaderGL(GL_VERTEX_SHADER, vertex_shader_file));
  shaderList.push_back(csX75::LoadShaderGL(GL_FRAGMENT_SHADER, fragment_shader_file));

  shaderProgram = csX75::CreateProgramGL(shaderList);
  glUseProgram( shaderProgram );

  // set up vertex arrays
  GLuint vPosition = glGetAttribLocation( shaderProgram, "vPosition" );
  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  GLuint vColor = glGetAttribLocation( shaderProgram, "vColor" );
  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(v_positions)) );

  uModelViewMatrix = glGetUniformLocation( shaderProgram, "uModelViewMatrix");

}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'M') // && mods == GLFW_MOD_SHIFT)
    {
        cout<<"Mouse clicked and mode is M and src picked? "<<fetched_first<<endl;
        if(!fetched_first)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            // gluUnProject(xpos1, ypos1, 0.01,model_view, projection, viewport,	&pos3D_x, &pos3D_y, &pos3D_z);
            //
            // xpos1 = 2 * pos3D_x;
            // ypos1 = 2 * pos3D_y;

            // cout<<"Clicked  and comparing with"<<endl;

            //Search for the clicked point in the existing nodes.
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);
                // cout<<i<<" : "<<xpos1<<" "<<ypos1<<" -- "<<chkx<<" "<<chky<<endl;

                if((xpos1 <= chkx+5 && xpos1 >= chkx - 5) && (ypos1 <= chky+5 && ypos1 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetched_first = !fetched_first; break;
                }
            }
        }
        else
        {
            glfwGetCursorPos(window, &xpos2, &ypos2);
            // gluUnProject(xpos1, ypos1, 0.01,model_view, projection, viewport,	&pos3D_x, &pos3D_y, &pos3D_z);
            //
            // xpos1 = 2 * pos3D_x;
            // ypos1 = 2 * pos3D_y;

            // cout<<"Clicked  and comparing with"<<endl;

            //Search for the clicked point in the existing nodes.
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);
                // cout<<i<<" : "<<xpos1<<" "<<ypos1<<" -- "<<chkx<<" "<<chky<<endl;

                if((xpos2 <= chkx+5 && xpos2 >= chkx - 5) && (ypos2 <= chky+5 && ypos2 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                {
                    destnodeindex = stage2nodes.at(i);
                    fetched_first = !fetched_first;
                    cout<<"fetched both the points. right click mouse to merge"<<endl;
                    cout<<fetched_first<<"   "<<srcnodeindex<<"    "<<destnodeindex<<endl;break;
                }
            }
        }
    }
    else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS && mode == 'M' && !fetched_first && srcnodeindex != destnodeindex)
    {
        //Merging the srcnode and destnodeindex

        //Correcting vpositions and vcolors matrix data.
        for(int i=0;i<2*stage2nodes.size();i+=2)
        {
            double vposx = v_positions[i].x;
            double vposy = v_positions[i].y;
            double vposx1 = v_positions[i+1].x;
            double vposy1 = v_positions[i+1].y;
            double srcvertx = U[srcnodeindex][0]/scaleby-1.0;
            double srcverty = 1.0 - (U[srcnodeindex][1]/scaleby);
            double destvertx = U[destnodeindex][0]/scaleby-1.0;
            double destverty = 1.0 - (U[destnodeindex][1]/scaleby);

            if((abs(vposx - srcvertx) < 0.01 && abs(vposy - srcverty) < 0.01) && !(abs(vposx1 - destvertx)<0.01 && abs(vposy1 - destverty)<0.01))
            {
                cout<<"Changing i"<<endl;
                v_positions[i].x = destvertx;
                v_positions[i].y = destverty;
            }
            else if((abs(vposx1 - srcvertx)<0.01 && abs(vposy1 - srcverty) < 0.01) && !(abs(vposx - destvertx) < 0.01 && abs(vposy - destverty) < 0.01 ))
            {
                cout<<"Changing i+1"<<endl;
                v_positions[i+1].x = destvertx;
                v_positions[i+1].y = destverty;
            }
            else if((abs(vposx - srcvertx) < 0.01 && abs(vposy - srcverty)<0.01 && abs(vposx1 - destvertx) < 0.01 && abs(vposy1 - destverty)<0.01) ||
                    (abs(vposx - destvertx)<0.01 && abs(vposy - destverty)<0.01 && abs(vposx1 - srcvertx) < 0.01 && abs(vposy1 - srcverty)<0.01))
            {
                cout<<"removing at "<<i<<endl;
                for(int ki=i; ki <= 2*stage2nodes.size()-1; ki+=2)
                {
                    v_positions[ki] = v_positions[ki+2];
                    v_positions[ki+1] = v_positions[ki+3];
                }
                i-=2;
            }
        }


        int positionindex = 0;
        //Removing the edges
        cout<<"removing edges"<<endl;
        for(int i=0;i<E.size();i++)
        {
            if(!(E[i][0] == 0 && E[i][1] == 0))
            {
                if(E[i][0] == srcnodeindex && E[i][1] != destnodeindex)
                {
                  E[i][0] = destnodeindex;
                }
                else if(E[i][1] == srcnodeindex && E[i][0] != destnodeindex)
                {
                  E[i][1] = destnodeindex;
                }
                else if((E[i][0] == srcnodeindex && E[i][1] == destnodeindex)|| (E[i][1] == srcnodeindex && E[i][0] == destnodeindex))
                {
                  E[i][1] = 0;  E[i][0] = 0;
                }
                positionindex++;
            }
        }

        //change history of destination node
        cout<<"changing vertex history"<<endl;
        for(int i=0;i<history.at(srcnodeindex).size();i++)
        {
            if(find(history.at(destnodeindex).begin(),history.at(destnodeindex).end(),history.at(srcnodeindex).at(i)) == history.at(destnodeindex).end())
            {
                history.at(destnodeindex).push_back(history.at(srcnodeindex).at(i));
            }
        }

        //Making vertex position as 0,0
        U[srcnodeindex][0] = 0.0; U.at(srcnodeindex).at(1) = 0.0;

        //Remove from nodeslist,junctlist;
        int targind = find(stage2nodes.begin(), stage2nodes.end(), destnodeindex) - stage2nodes.begin();
        stage2junct.at(targind).erase(find(stage2junct.at(targind).begin(), stage2junct.at(targind).end(),srcnodeindex));
        int indextoremove;
        indextoremove = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex)-stage2nodes.begin();
        stage2nodes.erase(stage2nodes.begin()+indextoremove);
        stage2junct.erase(stage2junct.begin()+indextoremove);

        cout<<"Merging done"<<num_vertices<<endl;
        srcnodeindex = 0; destnodeindex = 0;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'A')
    {
        cout<<"In add mode: fetched points: "<<fetchedpoints<<endl;
        if(fetchedpoints == 0)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);

                if((xpos1 <= chkx+5 && xpos1 >= chkx - 5) && (ypos1 <= chky+5 && ypos1 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetchedpoints = (fetchedpoints+1)%3;
                    cout<<"Fetched source"<<endl;break;
                }
            }
        }
        else if(fetchedpoints == 1)
        {
            glfwGetCursorPos(window, &xpos2, &ypos2);

            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);

                if((xpos2 <= chkx+5 && xpos2 >= chkx - 5) && (ypos2 <= chky+5 && ypos2 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                {
                    destnodeindex = stage2nodes.at(i);
                    fetchedpoints = (fetchedpoints+1)%3;
                    cout<<"fetched both the points. click anywhere on the edge to add new point"<<endl;
                }
            }
        }
        else if(fetchedpoints == 2)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            gluUnProject(xpos1,ypos1,0.0,model_view,projection,viewport,&pos3D_x,&pos3D_y,&pos3D_z);
            xpos1 = 2*pos3D_x;
            ypos1 = 2*pos3D_y;

            cout<<xpos1<<" "<<ypos1<<endl;
            double newnodex = (xpos1+1.0)*scaleby;
            double newnodey = (ypos1+1.0)*scaleby;

            double xdist1,xdist2,ydist1,ydist2;
            xdist1 = xpos1 - (U[srcnodeindex][0]/scaleby) + 1.0; xdist1 = xdist1 > 0? xdist1 : -1 * xdist1;
            ydist1 = ypos1 - (U[srcnodeindex][1]/scaleby) + 1.0; ydist1 = ydist1 > 0? ydist1 : -1 * ydist1;
            xdist2 = xpos1 - (U[destnodeindex][0]/scaleby) + 1.0; xdist2 = xdist2 > 0? xdist2 : -1 * xdist2;
            ydist2 = ypos1 - (U[destnodeindex][1]/scaleby) + 1.0; ydist2 = ydist2 > 0? ydist2 : -1 * ydist2;

            double slopeofedge = (U[destnodeindex][1] - U[srcnodeindex][1]) / (U[destnodeindex][0] - U[srcnodeindex][0]);

            if(xdist1 > ydist1 && xdist2 > ydist2)
            //y coordinate of new point to be changed
            {
               newnodey = (slopeofedge * (newnodex - U[srcnodeindex][0])) + U[srcnodeindex][1];
               ypos1 = (newnodey/scaleby) - 1.0;
            }
            else if(ydist1 > xdist1 && ydist2 > xdist2)
            //x coordinate to be changed
            {
               newnodex = ((newnodey - U[srcnodeindex][1]) / slopeofedge) + U[srcnodeindex][0];
               xpos1 = (newnodex/scaleby) - 1.0;
            }

            vector<double> temp;
            temp.push_back(newnodex); temp.push_back(newnodey);

            //Insert in original vertices and in modified vertices
            V.push_back(temp);
            U.push_back(temp);

            //index of new node to be used in edges and dneighborhood, history
            int newindex = V.size() - 1;

            //creating an entry in history.
            vector<int> dummy;
            dummy.push_back(newindex);
            history.push_back(dummy);
            dummy.clear();

            //Add new edges to edges list
            vector<int> temp1;
            temp1.push_back(srcnodeindex); temp1.push_back(newindex); E.push_back(temp1);
            temp1.clear();
            temp1.push_back(newindex); temp1.push_back(destnodeindex); E.push_back(temp1);

            cout<<"Added two new edges"<<endl;

            //Removing edge from E
            for(int i=0;i<E.size();i++)
            {
               if(!(E[i][0]==0 && E[i][1]))
               {
                   if((E[i][0] == srcnodeindex && E[i][1] == destnodeindex) || (E[i][0] == destnodeindex && E[i][1] == srcnodeindex))
                   {
                       E[i][0] = 0; E[i][1] = 0;
                       break;
                   }
               }
            }
            cout<<"removed an edge"<<endl;

            //Editing v_positions
            //Remove the original edge from v_positions and adding new edges
            int ki;
            for(int i=0;i<2*stage2nodes.size();i+=2)
            {
               double vposx = v_positions[i].x;
               double vposy = v_positions[i].y;
               double vposx1 = v_positions[i+1].x;
               double vposy1 = v_positions[i+1].y;
               double srcvertx = U[srcnodeindex][0]/scaleby-1.0;
               double srcverty = 1.0 - (U[srcnodeindex][1]/scaleby);
               double destvertx = U[destnodeindex][0]/scaleby-1.0;
               double destverty = 1.0 - (U[destnodeindex][1]/scaleby);

               if((abs(vposx - srcvertx) < 0.01 && abs(vposy - srcverty)<0.01 && abs(vposx1 - destvertx) < 0.01 && abs(vposy1 - destverty)<0.01) ||
                       (abs(vposx - destvertx)<0.01 && abs(vposy - destverty)<0.01 && abs(vposx1 - srcvertx) < 0.01 && abs(vposy1 - srcverty)<0.01))
               {
                   cout<<"removing at "<<i<<endl;
                   for(ki=i; ki <= 2*stage2nodes.size()-1; ki+=2)
                   {
                       v_positions[ki] = v_positions[ki+2];
                       v_positions[ki+1] = v_positions[ki+3];
                   }

                   ki-=4;

                   v_positions[ki] = glm::vec4(srcvertx, srcverty,0.0f,1.0f);
                   v_colors[ki] = avail_colors[ki%4]; ki++;
                   v_positions[ki] = glm::vec4(xpos1, -ypos1,0.0f,1.0f);
                   v_colors[ki] = avail_colors[ki%4]; ki++;
                   v_positions[ki] = glm::vec4(xpos1, -ypos1,0.0f,1.0f);
                   v_colors[ki] = avail_colors[ki%4]; ki++;
                   v_positions[ki] = glm::vec4(destvertx, destverty,0.0f,1.0f);
                   v_colors[ki] = avail_colors[ki%4]; ki++;
                   break;
               }
            }

            //Adjusting History of src and new node
            for(int i=0; i<history.at(srcnodeindex).size(); i++)
            {
                int nodeindexiter = history[srcnodeindex][i];
                //distance from node to src
                double dist1 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
                double dist2 = sqrt(pow((U[newindex][0] - V[nodeindexiter][0]),2)+pow((U[newindex][1] - V[nodeindexiter][1]),2));
                if(dist1 > dist2) // Node is farther from original node than it is from new node
                {
                    history.at(newindex).push_back(nodeindexiter);
                    history.at(srcnodeindex).erase(history.at(srcnodeindex).begin()+i);
                }
            }

            //Adjusting History of dest and new node
            for(int i=0; i<history.at(destnodeindex).size(); i++)
            {
                int nodeindexiter = history[destnodeindex][i];
                //distance from node to src
                double dist1 = sqrt(pow((U[destnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[destnodeindex][1] - V[nodeindexiter][1]),2));
                double dist2 = sqrt(pow((U[newindex][0] - V[nodeindexiter][0]),2)+pow((U[newindex][1] - V[nodeindexiter][1]),2));
                if(dist1 > dist2) // Node is farther from original node than it is from new node
                {
                    history.at(newindex).push_back(nodeindexiter);
                    history.at(destnodeindex).erase(history.at(destnodeindex).begin()+i);
                }
            }

            //Updating stage2nodes list and stage2junct (neighborhood) list`
            stage2nodes.push_back(V.size()-1);
            vector<int> tempj;
            tempj.push_back(srcnodeindex);
            tempj.push_back(destnodeindex);
            stage2junct.push_back(tempj);//neighbourhood of new node
            // //Neighborhood of source node
            // int targetind = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex) - stage2nodes.begin();
            // stage2junct.at(targetind).push_back(newindex);
            // //Neighborhood of destination node
            // targetind = find(stage2nodes.begin(),stage2nodes.end(),destnodeindex) - stage2nodes.begin();
            // stage2junct.at(targetind).push_back(newindex);


            // cout<<"new history"<<history.at(newindex).size()<<endl;
            fetchedpoints = (fetchedpoints+1)%3;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'L')
    {
        cout<<"In adding leaf mode: fetched points: "<<fetchedpoints<<endl;
        if(fetchedpoints == 0)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);

                if((xpos1 <= chkx+5 && xpos1 >= chkx - 5) && (ypos1 <= chky+5 && ypos1 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetchedpoints = (fetchedpoints+1)%2;
                    cout<<"Fetched source"<<endl;break;
                }
            }
        }
        else if(fetchedpoints == 1)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            gluUnProject(xpos1,ypos1,0.0,model_view,projection,viewport,&pos3D_x,&pos3D_y,&pos3D_z);
            xpos1 = 2*pos3D_x;
            ypos1 = 2*pos3D_y;

            cout<<xpos1<<" "<<ypos1<<endl;
            double newnodex = (xpos1+1.0)*scaleby;
            double newnodey = (ypos1+1.0)*scaleby;

            vector<double> temp;
            temp.push_back(newnodex); temp.push_back(newnodey);

            //Insert in original vertices and in modified vertices
            V.push_back(temp);
            U.push_back(temp);

            //index of new node to be used in edges and dneighborhood, history
            int newindex = V.size() - 1;

            //creating an entry in history.
            vector<int> dummy;
            dummy.push_back(newindex);
            history.push_back(dummy);
            dummy.clear();

            //Add new edges to edges list
            vector<int> temp1;
            temp1.push_back(srcnodeindex); temp1.push_back(newindex); E.push_back(temp1);
            // temp1.push_back(newindex); temp1.push_back(destnodeindex); E.push_back(temp1);

            //Editing v_positions
            int i;
            for(i=0;;i+=2)
            {
                double vposx = v_positions[i].x;
                double vposy = v_positions[i].y;
                double vposx1 = v_positions[i+1].x;
                double vposy1 = v_positions[i+1].y;

                if(abs(vposx-0.0) < 0.0001 && abs(vposy - 0.0) < 0.0001  && abs(vposx1 - 0.0) < 0.0001 && abs(vposy1 - 0.0) < 0.0001)
                {
                    break;
                }
            }
            double srcvertx = U[srcnodeindex][0]/scaleby-1.0;
            double srcverty = 1.0 - (U[srcnodeindex][1]/scaleby);
            v_positions[i] = glm::vec4(srcvertx,srcverty,0.0f,1.0f);
            v_colors[i] = avail_colors[i%4]; i++;
            v_positions[i] = glm::vec4(xpos1,-ypos1,0.0f,1.0f);
            v_colors[i] = avail_colors[i%4]; i++;

            //Adjusting History of src and new node
            for(int i=0; i<history.at(srcnodeindex).size(); i++)
            {
                int nodeindexiter = history[srcnodeindex][i];
                //distance from node to src
                double dist1 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
                double dist2 = sqrt(pow((U[newindex][0] - V[nodeindexiter][0]),2)+pow((U[newindex][1] - V[nodeindexiter][1]),2));
                if(dist1 > dist2) // Node is farther from original node than it is from new node
                {
                    history.at(newindex).push_back(nodeindexiter);
                    history.at(srcnodeindex).erase(history.at(srcnodeindex).begin()+i);
                }
            }

            //Updating stage2nodes list and stage2junct (neighborhood) list
            stage2nodes.push_back(V.size()-1);
            vector<int> tempj;
            tempj.push_back(srcnodeindex);

            stage2junct.push_back(tempj);//Neighborhood of new node
            //Neighborhood of source node
            int targetind = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex) - stage2nodes.begin();
            stage2junct.at(targetind).push_back(newindex);

            // cout<<"new history"<<history.at(newindex).size()<<endl;
            fetchedpoints = (fetchedpoints+1)%2;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'V' )
    {
        cout<<"fetched_first: "<<fetched_first<<endl;
        if(!fetched_first)
        {
            glfwGetCursorPos(window,&xpos1,&ypos1);
            // gluUnProject(xpos1,ypos1,0.0,model_view,projection,viewport,&pos3D_x,&pos3D_y,&pos3D_z);
            // xpos1 = 2*pos3D_x;
            // ypos1 = 2*pos3D_y;
            //
            // // cout<<xpos1<<" "<<ypos1<<endl;
            // double nodeposx = (xpos1+1.0)*scaleby;
            // double nodeposy = (ypos1+1.0)*scaleby;
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U[stage2nodes.at(i)][0]/scaleby - 1.0;
                curry = U[stage2nodes.at(i)][1]/scaleby - 1.0;

                //Project it to compare with clicked coordinate.
                double chkx,chky,chkz;
                gluProject(currx/2,curry/2,0.0,model_view,projection,viewport,&chkx,&chky,&chkz);

                if((xpos1 <= chkx+5 && xpos1 >= chkx - 5) && (ypos1 <= chky+5 && ypos1 >= chky - 5))
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetched_first = !fetched_first;
                    cout<<"Fetched point. Click where to move"<<endl;
                    cout<<fetched_first<<endl;
                }
            }
        }
        else
        {
            glfwGetCursorPos(window,&xpos1,&ypos1);
            gluUnProject(xpos1,ypos1,0.0,model_view,projection,viewport,&pos3D_x,&pos3D_y,&pos3D_z);
            xpos1 = 2*pos3D_x;
            ypos1 = 2*pos3D_y;

            cout<<xpos1<<" "<<ypos1<<endl;
            double nodeposx = (xpos1+1.0)*scaleby;
            double nodeposy = (ypos1+1.0)*scaleby;

            //Changing v_positions
            for(int i=0;i<2*stage2nodes.size();i+=2)
            {
               double vposx = v_positions[i].x;
               double vposy = v_positions[i].y;
               double vposx1 = v_positions[i+1].x;
               double vposy1 = v_positions[i+1].y;
               double srcvertx = (U[srcnodeindex][0]/scaleby)-1.0;
               double srcverty = 1.0 - (U[srcnodeindex][1]/scaleby);

               if(abs(vposx - srcvertx) < 0.01 && abs(vposy - srcverty) < 0.01)
               {
                  v_positions[i].x = xpos1;
                  v_positions[i].y = -ypos1;
               }
               else if(abs(vposx1 - srcvertx) < 0.01 && abs(vposy1 - srcverty) < 0.01)
               {
                  v_positions[i+1].x = xpos1;
                  v_positions[i+1].y = -ypos1;
               }
            }

            //Changing actual vertex positions
            U[srcnodeindex][0] = nodeposx;
            U[srcnodeindex][1] = nodeposy;
            // V[srcnodeindex][0] = nodeposx;
            // V[srcnodeindex][1] = nodeposy;

            int ind = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex)-stage2nodes.begin();
            for(int ji = 0;ji < stage2junct.at(ind).size(); ji++)
            {
                int neighbornode = stage2junct[ind][ji];
                for(int hi=0; hi<history.at(neighbornode).size(); hi++)
                {
                    int nodeindexiter = history[neighbornode][hi];
                    //distance from node to src
                    double dist1 = sqrt(pow((U[neighbornode][0] - V[nodeindexiter][0]),2)+pow((U[neighbornode][1] - V[nodeindexiter][1]),2));
                    double dist2 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
                    if(dist1 > dist2) // Node is farther from original node than it is from new node
                    {
                        history.at(srcnodeindex).push_back(nodeindexiter);
                        history.at(neighbornode).erase(history.at(neighbornode).begin()+hi);hi--;
                    }
                }
                for(int hi=0; hi<history.at(srcnodeindex).size(); hi++)
                {
                    int nodeindexiter = history[srcnodeindex][hi];
                    //distance from node to src
                    double dist1 = sqrt(pow((U[neighbornode][0] - V[nodeindexiter][0]),2)+pow((U[neighbornode][1] - V[nodeindexiter][1]),2));
                    double dist2 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
                    if(dist1 < dist2) // Node is farther from original node than it is from new node
                    {
                        history.at(neighbornode).push_back(nodeindexiter);
                        history.at(srcnodeindex).erase(history.at(srcnodeindex).begin()+hi);hi--;
                    }
                }
            }

            fetched_first = !fetched_first;
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'S' )
    {
          ofstream op;
          op.open("skeleton.raw");
          // int counterl=0;
          for(int i=0;i<E.size();i++)
          {
            if(!(E[i][0] == 0 && E[i][1] == 0))
            {
              // counterl++;
              op<<U[E[i][0]][0]<<" "<<U[E[i][0]][1]<<" "<<U[E[i][0]][2]<<endl;
              op<<U[E[i][1]][0]<<" "<<U[E[i][1]][1]<<" "<<U[E[i][1]][2]<<endl;
            }
          }
          op.close();
          op.open("vertexdata.txt");
          op<<V.size()<<endl;
          for(int kind = 0; kind < V.size(); kind++)
          {
              op<<V[kind][0]<<" "<<V[kind][1]<<endl;
          }
          op<<U.size()<<endl;
          for(int kind = 0; kind < U.size(); kind++)
          {
              op<<U[kind][0]<<" "<<U[kind][1]<<endl;
          }
          op<<E.size()<<endl;
          for(int kind = 0; kind < E.size(); kind++)
          {
              op<<E[kind][0]<<" "<<E[kind][1]<<endl;
          }
          op<<E2.size()<<endl;
          for(int kind=0;kind < E2.size(); kind++)
          {
              op<<E2[kind][0]<<" "<<E2[kind][1]<<endl;
          }
          op.close();

          op.open("skeldata.txt");
          op<<stage2nodes.size()<<endl;
          for(int indx=0;indx<stage2nodes.size();indx++)
          {
              op<<stage2nodes.at(indx)<<" ";
          }
          op<<endl;
          // for(int indx=0;indx<E.rows();indx++)
          // {
          // 		if(E(indx,0) != 0 && E(indx,1)!=0)
          // 		{
          // 				op<<E(indx,0)<<" "<<E(indx,1)<<endl;
          // 		}
          // }
          op<<history.size()<<endl;
          for(int indx=0;indx<history.size();indx++)
          {
              // int vertindex = stage2nodes.at(indx);
              op<<history.at(indx).size()<<endl;
              for(int kin = 0; kin<history.at(indx).size();kin++)
              {
                  op<<history.at(indx).at(kin)<<" ";
              }
              op<<endl;
          }
          op.close();
          cout<<"Saved content"<<endl;
    }
    changeObjects();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    renderGL();
}

void initBuffersGL(void)
{
  //cg::loadSkeleton();
  loadvertexdata();
  loadSkeletondata();
  loadSkeleton();


  //For glu project and glu unproject
  glGetDoublev(GL_MODELVIEW_MATRIX, model_view);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);

  //Ask GL for a Vertex Attribute Object (vao)
  glGenVertexArrays (1, &vao);
  //Set it as the current array to be used by binding it
  glBindVertexArray (vao);

  //Ask GL for a Vertex Buffer Object (vbo)
  glGenBuffers (1, &vbo);
  //Set it as the current buffer to be used by binding it
  glBindBuffer (GL_ARRAY_BUFFER, vbo);
  //Copy the points into the current buffer
  glBufferData (GL_ARRAY_BUFFER, sizeof (v_positions) + sizeof(v_colors), NULL, GL_DYNAMIC_DRAW);
  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(v_positions), v_positions );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(v_positions), sizeof(v_colors), v_colors );

  // Load shaders and use the resulting shader program
  std::string vertex_shader_file("03_vshader.glsl");
  std::string fragment_shader_file("03_fshader.glsl");

  std::vector<GLuint> shaderList;
  shaderList.push_back(csX75::LoadShaderGL(GL_VERTEX_SHADER, vertex_shader_file));
  shaderList.push_back(csX75::LoadShaderGL(GL_FRAGMENT_SHADER, fragment_shader_file));

  shaderProgram = csX75::CreateProgramGL(shaderList);
  glUseProgram( shaderProgram );

  // set up vertex arrays
  GLuint vPosition = glGetAttribLocation( shaderProgram, "vPosition" );
  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  GLuint vColor = glGetAttribLocation( shaderProgram, "vColor" );
  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(v_positions)) );

  uModelViewMatrix = glGetUniformLocation( shaderProgram, "uModelViewMatrix");
}


int main(int argc, char** argv)
{
  //! The pointer to the GLFW window
  GLFWwindow* window;

  //! Setting up the GLFW Error callback
  glfwSetErrorCallback(csX75::error_callback);

  //! Initialize GLFW
  if (!glfwInit())
    return -1;

  //We want OpenGL 4.0
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  //This is for MacOSX - can be omitted otherwise
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  //We don't want the old OpenGL
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  //! Create a windowed mode window and its OpenGL context
  window = glfwCreateWindow(900, 900, "Assignment 1", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    return -1;
  }

  //! Make the window's context current
  glfwMakeContextCurrent(window);

  //Initialize GLEW
  //Turn this on to get Shader based OpenGL
  glewExperimental = GL_TRUE;
  GLenum err = glewInit();
  if (GLEW_OK != err)
    {
      //Problem: glewInit failed, something is seriously wrong.
      std::cerr<<"GLEW Init Failed : %s"<<std::endl;
    }

  //Print and see what context got enabled
  std::cout<<"Vendor: "<<glGetString (GL_VENDOR)<<std::endl;
  std::cout<<"Renderer: "<<glGetString (GL_RENDERER)<<std::endl;
  std::cout<<"Version: "<<glGetString (GL_VERSION)<<std::endl;
  std::cout<<"GLSL Version: "<<glGetString (GL_SHADING_LANGUAGE_VERSION)<<std::endl;

  //Keyboard Callback
  glfwSetKeyCallback(window, csX75::key_callback);
  //Framebuffer resize callback
  glfwSetFramebufferSizeCallback(window, csX75::framebuffer_size_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

  //Initialize GL state
  csX75::initGL();
  initBuffersGL();
  cout<<"Current mode is idle. Press any option to start editing"<<endl;

  // Loop until the user closes the window
  while (glfwWindowShouldClose(window) == 0)
    {

      // Render here// int ind = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex)-stage2nodes.begin();
            // for(int ji = 0;ji < stage2junct.at(ind).size(); ji++)
            // {
            //     int neighbornode = stage2junct[ind][ji];
            //     for(int hi=0; hi<history.at(neighbornode).size(); hi++)
            //     {
            //         int nodeindexiter = history[neighbornode][hi];
            //         //distance from node to src
            //         double dist1 = sqrt(pow((U[neighbornode][0] - V[nodeindexiter][0]),2)+pow((U[neighbornode][1] - V[nodeindexiter][1]),2));
            //         double dist2 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
            //         if(dist1 > dist2) // Node is farther from original node than it is from new node
            //         {
            //             history.at(srcnodeindex).push_back(nodeindexiter);
            //             history.at(neighbornode).erase(history.at(srcnodeindex).begin()+hi);
            //         }
            //     }
            //     for(int hi=0; hi<history.at(srcnodeindex).size(); hi++)
            //     {
            //         int nodeindexiter = history[srcnodeindex][hi];
            //         //distance from node to src
            //         double dist1 = sqrt(pow((U[neighbornode][0] - V[nodeindexiter][0]),2)+pow((U[neighbornode][1] - V[nodeindexiter][1]),2));
            //         double dist2 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
            //         if(dist1 < dist2) // Node is farther from original node than it is from new node
            //         {
            //             history.at(neighbornode).push_back(nodeindexiter);
            //             history.at(srcnodeindex).erase(history.at(srcnodeindex).begin()+hi);
            //         }
            //     }
            // }
      renderGL();

      // Swap front and back buffers
      glfwSwapBuffers(window);

      // Poll for and process events
      glfwPollEvents();
    }

    ofstream op;

    float colorvs[6][3] =
    {
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0,
        1.0,1.0,0.0,
        1.0,0.0,1.0,
        0.0,1.0,1.0
    };
    vector<int> colorcodes(V.size(),0);
    int cc=0;
    for(int i=0; i<stage2nodes.size(); i++)
    {
        cout<<history.at(stage2nodes[i]).size()<<endl;
        for(int j=0; j< history.at(stage2nodes[i]).size();j++)
        {
              colorcodes[history[stage2nodes[i]][j]] = cc;
        }
        cc = (cc+1)%6;
    }

    op.open("../2d/skeleton.raw");
    for(int i=0; i<E2.size(); i++)
    {
        op<<V[E2[i][0]][0]<<" "<<V[E2[i][0]][1]<<" "<<V[E2[i][0]][2]<<endl;
        int q = colorcodes[E2[i][0]];
        op<<colorvs[q][0]<<" "<<colorvs[q][1]<<" "<<colorvs[q][2]<<endl;
        op<<V[E2[i][1]][0]<<" "<<V[E2[i][1]][1]<<" "<<V[E2[i][1]][2]<<endl;
        q = colorcodes[E2[i][1]];
        op<<colorvs[q][0]<<" "<<colorvs[q][1]<<" "<<colorvs[q][2]<<endl;
    }

    for(int i=0;i<E.size();i++)
    {
      if(!(E[i][0]== 0 && E[i][1]== 0))
      {
        op<<U[E[i][0]][0]<<" "<<U[E[i][0]][1]<<" "<<U[E[i][0]][2]<<endl;
        op<<1.0<<" "<<1.0<<" "<<1.0<<endl;
        op<<U[E[i][1]][0]<<" "<<U[E[i][1]][1]<<" "<<U[E[i][1]][2]<<endl;
        op<<1.0<<" "<<1.0<<" "<<1.0<<endl;
      }
    }
    op.close();

    op.open("../2d/vertices.txt");
    op<<V.size()<<endl;
    for(int kind = 0; kind < V.size(); kind++)
    {
        op<<V[kind][0]<<" "<<V[kind][1]<<" "<<0.0<<endl;
    }
    op.close();

    op.open("../2d/joints.txt");
    for(int kind = 0; kind < stage2nodes.size(); kind++)
    {
        op<<U[stage2nodes[kind]][0]<<" "<<U[stage2nodes[kind]][1]<<" "<<0.0<<endl;
    }
    op.close();

    op.open("../2d/mapping.txt");
    for(int kind=0;kind<stage2nodes.size();kind++)
    {
        int currind = stage2nodes[kind];
        for(int kk=0; kk < history[currind].size(); kk++)
        {
            op<<history[currind][kk]<<" ";
        }
        op<<endl;
    }
    op.close();


  glfwTerminate();
  return 0;
}

//-------------------------------------------------------------------------
