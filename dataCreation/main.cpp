#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>
#include <igl/outer_element.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/edge_lengths.h>
#include <igl/edge_flaps.h>
#include <fstream>
#include <algorithm>

using namespace std;

namespace mesh
{
	template <typename DerivedV, typename DerivedF, typename DeriveddblA>
	IGL_INLINE void doublearea(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedF> & F,
	Eigen::PlainObjectBase<DeriveddblA> & dblA)
	{
		const int dim = V.cols();
		const size_t m = F.rows();

		static int debugflag = 1;

		// Projected area helper
		const auto & proj_doublearea =
		[&V,&F](const int x, const int y, const int f)
		->typename DerivedV::Scalar
		{
			if( V(F(f,0),x)!=V(F(f,0),x) || V(F(f,2),x) != V(F(f,2),x) ||
			    V(F(f,1),x)!=V(F(f,1),x) || V(F(f,2),x) != V(F(f,2),x) ||
					V(F(f,0),y) != V(F(f,0),y) || V(F(f,2),y)!=V(F(f,2),y) ||
					V(F(f,1),y) != V(F(f,1),y) || V(F(f,2),y) != V(F(f,2),y))
			{
					// cout<<"there is a non-number returning zero"<<endl;
					return 0.0;
			}
			auto rx = V(F(f,0),x)-V(F(f,2),x);
			auto sx = V(F(f,1),x)-V(F(f,2),x);
			auto ry = V(F(f,0),y)-V(F(f,2),y);
			auto sy = V(F(f,1),y)-V(F(f,2),y);
			double area = rx*sy - ry*sx;
			// if(area<0)
			// 	return -1*area;
			// else
				return area;
		};

		dblA = DeriveddblA::Zero(m,1);
		for(size_t f = 0;f<m;f++)
		{
			for(int d = 0;d<2;d++)
			{
				const auto dblAd = proj_doublearea(d,(d+1)%2,f);
				//if(debugflag >12 && (dblAd > 0))
				//	cout<<dblAd<<"		";
				dblA(f) += dblAd*dblAd;

			}
		}
		dblA = dblA.array().sqrt().eval();debugflag++;
	}
}


int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
   MatrixXd OV;
   MatrixXi OF;
   igl::readOFF("mesh1.off",OV,OF);

	 VectorXi bnd,b(2,1);
   // igl::boundary_loop(OF,bnd);

	 cout << OF.rows()<<endl;
	 cout<<bnd<<endl<<endl;

   igl::opengl::glfw::Viewer viewer;

	 /*Some variable names are directly taken from the reference paper. Please refer to the paper*/
   MatrixXd V,U,tempU;
   MatrixXi F,tempF;
   SparseMatrix<double> L;
   VectorXd area,interarea;
   double areathreshold;
   int sl=2.0,vcount;
   MatrixXd wl,wh,whi;
   ArrayXf ringarea;
   int num_collapsed;
	 int stage = 1;
	 vector <MatrixXd> qmatrix;
	 vector <vector<int>> history;
	 VectorXi EMAP;
	 MatrixXi E,EF,EI;
	 MatrixXi E2;
	 double wa = 1.0, wb = 0.1;
	 VectorXd edgecosts;
	 MatrixXd fbedge;
	 MatrixXd homogeneousV;
	 MatrixXd fa;
	 vector<int> stage2nodes;
	 vector<vector<int>> stage2junct;
	 vector<vector<int>> vertextoedgemap;
	 vector<vector<int>> boundaryvertices;

   const auto & reset = [&]()
   {
     F=OF;
     V=OV;
     num_collapsed = 0;
     vcount = V.rows();
     U=V;

		 cout<<"Vertices: "<<V.rows()<<endl;
		 cout<<"Faces: "<<F.rows()<<endl;

     igl:cotmatrix(V,F,L);
     mesh::doublearea(V,F,area);
     area = area.array()/2;
     double area_avg = area.mean();

     cout<<"Area of mesh :"<<area.sum()<<endl;
     ringarea = ArrayXf::Zero(vcount);

     for(int i=0; i<F.rows(); i++)
     {
       ringarea(F(i,0)) += area(i);
       ringarea(F(i,1)) += area(i);
       ringarea(F(i,2)) += area(i);
     }

     cout<<"Ring area calculated"<<endl;
     int i;
     wl = MatrixXd::Zero(vcount, vcount);
     for(int i=0;i<vcount;i++)
     {
       wl(i,i) = 0.001 * sqrt(area_avg);
     }
     cout<<"wl"<<endl;

     wh = MatrixXd::Zero(vcount,vcount);
     for(i=0; i<vcount; i++)
     {
       wh(i,i) = 1.0;
     }
     cout<<"wh"<<endl;
     whi = wh;
     cout<<"Reset done"<<endl;
     viewer.data().clear();
     viewer.data().set_mesh(V,F);
     viewer.data().set_face_based(true);

   };

	 const auto &setstage2 = [&]()
	 {
		 		cout<<"In setstage2 function"<<endl;
		 		igl::edge_flaps(F,E,EMAP,EF,EI); //Creating an edge Matrix
				// cout<<EMAP<<endl<<endl<<endl;

				cout<<"EF: "<<EF(0,0)<<" "<<EF(0,1)<<endl;
				cout<<"EI: "<<EI(0,0)<<" "<<EI(0,1)<<endl;
				E2 = E;

				//History matrix to maintain the list of vertices collapsed.
				for(int i=0;i<U.rows(); i++)
				{
					vector<int> temp;
					temp.push_back(i);
					history.push_back(temp);
				}


				//Creating a Q matrix for all the vertices.
				MatrixXd dummy = MatrixXd::Zero(4,4);
				for (int i = 1; i <= vcount; i++) //Initializing them all with zero matrix.
					qmatrix.push_back(dummy);

				//For each edge, we create Kij and add it to the respective Vi qmatrix
				for(int i=0; i<E.rows(); i++)
				{
						for(int j=0; j<2; j++)
						{
							if(!(E(i,0) == 0 && E(i,1) == 0))
							{
									MatrixXd kij = MatrixXd::Zero(3,4);
									Vector3d edgij;
									int k=(j+1)%2;
									edgij << (U(E(i,k),0)-U(E(i,j),0)), (U(E(i,k),1)-U(E(i,j),1)), (U(E(i,k),2)-U(E(i,j),2));
									//VectorXd edgij(3);
									edgij = edgij.normalized();
									MatrixXd a = MatrixXd::Zero(3,3);
									a(0,1) = -edgij(2); a(1,0) = -a(0,1);
									a(0,2) = edgij(1); a(2,0) = -a(0,2);
									a(1,2) = -edgij(0); a(2,1) = -a(1,2);
									Vector3d vi;
									vi << U(E(i,j),0), U(E(i,j),1), U(E(i,j),2);//= U(E(i,j)).row();
									Vector3d b;
									b = -(edgij.cross(vi));
									kij<<a,b;
									MatrixXd ktk = kij.transpose() * kij;
									qmatrix.at(E(i,j)) += ktk;
							}
						}
				}

				cout<<"Done creating Qmatrix"<<endl;
		};

		const auto &fafbupdate = [&]()
		{

				//Creating fa for each vertex. For each edge we do fa(i) + fa(j)
				homogeneousV = U.rowwise().homogeneous();
				fa = MatrixXd::Zero(E.rows(),2);
				for(int i =0; i< E.rows(); i++)
				{
					VectorXd pt(4);
					pt << (homogeneousV(E(i,0),0)+homogeneousV(E(i,1),0))/2, (homogeneousV(E(i,0),1)+homogeneousV(E(i,1),1))/2, (homogeneousV(E(i,0),2)+homogeneousV(E(i,1),2))/2, (homogeneousV(E(i,0),3)+homogeneousV(E(i,1),3))/2;
					fa(i,0) = pt.transpose() * qmatrix.at(E(i,0)) * pt;
					fa(i,1) = pt.transpose() * qmatrix.at(E(i,1)) * pt;
				}

				//Creating fb. We are taking two cols of fb because the direction of collapse matters for fb.
				//calculating norms of all the edges.
				VectorXd normvec = VectorXd::Zero(E.rows());
				for(int i=0;i<E.rows();i++)
				{
					if(!(E(i,1)== 0 && E(i,0)== 0))
					{
							double dist = sqrt( pow(( U(E(i,1),0) - U(E(i,0),0) ),2) +  pow(( U(E(i,1),1) - U(E(i,0),1) ),2) + pow(( U(E(i,1),2) - U(E(i,0),2) ),2) );
							normvec(i) = dist;
					}
				}

				// cout<<"Creating FB for vertices"<<endl;
				fbedge = MatrixXd::Zero(E.rows(),2);

				//creating fb entry as product of currentedge and sum of all the other edges.
				for(int i=0; i<E.rows(); i++)
				{
					for(int j=0;j<2;j++)
					{
							double magv = 0;
							//cout<<"Calculating magv"<<endl;
							for(int k=0; k<E.rows(); k++)
							{
								if( (E(k,0)==E(i,j) || E(k,1)==E(i,j)) && (k != i))
									magv += normvec(k);
							}
							//cout<<"magv calculated"<<endl;
							//cout<<i<<" , "<<j<<endl;
							double dummy = normvec(i) * magv;
							fbedge(i,j) = dummy;
					}
				}
				// cout<<"Done creating fb"<<endl;
				viewer.core.is_animating = true;
	 };

	 const auto &removeDupFaces = [&]()
	 //Removing duplicating faces
	 {
		 	vector<vector<int>> faceslist;
		 	for(int i=0; i<F.rows(); i++)
			{
					vector<int> tempvec;
					tempvec.push_back(F(i,0));
					tempvec.push_back(F(i,1));
					tempvec.push_back(F(i,2));
					sort(tempvec.begin(),tempvec.end());
					faceslist.push_back(tempvec);
					tempvec.clear();
			}
			int iter_pos;
			for(int i=0; i<F.rows(); i++)
			{
					if(!(F(i,0) == 0 && F(i,1) == 0 && F(i,2) == 0))
					{
							std::vector<vector<int>>::iterator iter = faceslist.begin()+i+1;
							while(iter != faceslist.end())
							{
									vector<vector<int>>::iterator dummy_iter;
									dummy_iter = find(iter,faceslist.end(),faceslist.at(i));
									if(dummy_iter != faceslist.end())
									{
											iter_pos = dummy_iter - faceslist.begin();
											vector<int> zerovec = {0,0,0};
											faceslist.at(iter_pos) = zerovec;
											cout<<"Making "<<F(iter_pos,0)<<" "<<F(iter_pos,1)<<" "<<F(iter_pos,2)<<"  zero for duplication"<<endl;
											F(iter_pos,0) = 0; F(iter_pos,1) = 0; F(iter_pos,2) = 0;
											iter = dummy_iter +1;
									}
										else
											iter = dummy_iter;
							}
					}
			}
			// cout<<"done"<<endl;
	 };

	 const auto &removeDupEdges = [&]()
	 /*Removing duplicate edges*/
	 {
		 	vector<vector<int>> edgeslist;
		 	for(int i=0; i<E.rows(); i++)
			{
					vector<int> tempvec;
					tempvec.push_back(E(i,0));
					tempvec.push_back(E(i,1));
					sort(tempvec.begin(),tempvec.end());
					edgeslist.push_back(tempvec);
					tempvec.clear();
			}
			int iter_pos;
			for(int i=0; i<E.rows(); i++)
			{
					if(!(E(i,0) == 0 && E(i,1) == 0))
					{
							std::vector<vector<int>>::iterator iter = edgeslist.begin()+i+1;
							while(iter != edgeslist.end())
							{
									vector<vector<int>>::iterator dummy_iter;
									dummy_iter = find(iter,edgeslist.end(),edgeslist.at(i));
									if(dummy_iter != edgeslist.end())
									{
											iter_pos = dummy_iter - edgeslist.begin();
											vector<int> zerovec = {0,0};
											edgeslist.at(iter_pos) = zerovec;
											// cout<<"Making "<<E(iter_pos,0)<<" "<<E(iter_pos,1)<<"  zero for duplication"<<endl;
											E(iter_pos,0) = 0; E(iter_pos,1) = 0;
											iter = dummy_iter +1;
									}
										else
											iter = dummy_iter;
							}
					}
			}
			// cout<<"done"<<endl;
	 };

	 const auto &isedgeExists = [&](int v1,int v2) -> bool
	 {
		 	for(int iter=0; iter<vertextoedgemap.at(v1).size(); iter++)
			{
					int edgeindex = vertextoedgemap.at(v1).at(iter);
					if(find(vertextoedgemap.at(v2).begin(), vertextoedgemap.at(v2).end(), edgeindex) != vertextoedgemap.at(v2).end())
							return true;
			}
	 };

	 const auto &setstage3 = [&]()
	 {
		 for(int iter=0;iter<V.rows();iter++)
		 {
			 	vector<int> temp;
				// temp.push_back(0);
				vertextoedgemap.push_back(temp);
		 }
		 //For each vertex we are saving all the edges that are sharing the vertex.
		 for(int iter=0; iter<E2.rows(); iter++)
		 {
			 		vertextoedgemap.at(E2(iter,0)).push_back(iter);
					vertextoedgemap.at(E2(iter,1)).push_back(iter);
		 }

		 for(int iter=0;iter<E.rows();iter++)
		 {
			 if(!(E(iter,0)==0&&E(iter,1)==0))
			 {
				 cout<<"Edge: "<<E(iter,0)<<" "<<E(iter,1)<<"\nVertices: "<<U(E(iter,0),0)<<" "<<U(E(iter,0),1)<<"  ---  "<<U(E(iter,1),0)<<" "<<U(E(iter,1),1)<<endl;
				 //cout<<"Non zero face is: "<<F(iter,0)<<","<<F(iter,1)<<","<<F(iter,2)<<endl;
				 //zeroflag = true; break;
				 if(stage2nodes.size() == 0 || find(stage2nodes.begin(), stage2nodes.end(), E(iter,0)) == stage2nodes.end())
				 {
							 stage2nodes.push_back(E(iter,0));
							 vector<int> temp;
							 temp.push_back(E(iter,1));
							 stage2junct.push_back(temp);
				 }
				 else
				 {
							 int temp = find(stage2nodes.begin(), stage2nodes.end(),E(iter,0))-stage2nodes.begin();
							 stage2junct.at(temp).push_back(E(iter,1));
				 }
				 if(stage2nodes.size() == 0 || find(stage2nodes.begin(), stage2nodes.end(), E(iter,1)) == stage2nodes.end())
				 {
							 stage2nodes.push_back(E(iter,1));
							 vector<int> temp;
							 temp.push_back(E(iter,0));
							 stage2junct.push_back(temp);
				 }
				 else
				 {
							 int temp = find(stage2nodes.begin(), stage2nodes.end(), E(iter,1)) - stage2nodes.begin();
							 stage2junct.at(temp).push_back(E(iter,0));
				 }
			 }
		 }

		 for(int iter=0;iter<stage2nodes.size(); iter++)
		 {
			 	int nodeindex = stage2nodes.at(iter);
				vector<int> tempboundary;
			 	for(int jiter=0;jiter<history.at(nodeindex).size(); jiter++)
				{
						int curvertex = history.at(nodeindex).at(jiter);
						for(int kiter = 0; kiter < vertextoedgemap.at(curvertex).size(); kiter++)
						{
								int curedge = vertextoedgemap.at(curvertex).at(kiter);
								int endpoint1 = E2(curedge,0);
								int endpoint2 = E2(curedge,1);
								if(!(find(history.at(nodeindex).begin(), history.at(nodeindex).end(), endpoint1) != history.at(nodeindex).end() &&
									 find(history.at(nodeindex).begin(), history.at(nodeindex).end(), endpoint2) != history.at(nodeindex).end()))
								{
										if( find(history.at(nodeindex).begin(), history.at(nodeindex).end(), endpoint1) != history.at(nodeindex).end() )
										{
												if(find(tempboundary.begin(), tempboundary.end(), endpoint1) == tempboundary.end())
														tempboundary.push_back(endpoint1);
										}
										else
										{
												if(find(tempboundary.begin(), tempboundary.end(), endpoint2) == tempboundary.end())
														tempboundary.push_back(endpoint2);
										}
								}
						}
				}
				boundaryvertices.push_back(tempboundary);
		 }
		 for(int iter=0; iter<boundaryvertices.size(); iter++)
		 {
			 		cout<<stage2nodes.at(iter)<<"  -  ";
			 		for(int jiter=0; jiter<boundaryvertices.at(iter).size(); jiter++)
					{
								cout<<boundaryvertices.at(iter).at(jiter)<<"  ";
					}
					cout<<endl;
		 }

	 };

   const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool
   {
      if(viewer.core.is_animating)
      {
					if(stage == 1)
					{
							bool something_collapsed = false;

        			MatrixXd a(vcount,vcount);
        			a=wl*L;
        			MatrixXd lhs(2*vcount, vcount);
        			lhs << a,wh;

        			MatrixXd b(vcount, vcount);
        			b = wh * U;
        			ArrayXXd zro = ArrayXXd :: Zero(vcount,3);
        			MatrixXd rhs(2*vcount,3);
        			rhs <<zro,b;

        			tempU = U;
        			tempF = F;

        			ColPivHouseholderQR<MatrixXd> solver(lhs);
        			U = solver.solve(rhs);
        			something_collapsed = true;
        			cout<<"Solving done"<<endl;

        			mesh::doublearea(U,F,interarea);
        			interarea = interarea.array()/2;

        			areathreshold = interarea.sum()/area.sum();
        			cout<<"Areathreshold: "<<areathreshold<<endl;

        			if(areathreshold < 0.0001)
        			{
          			U = tempU;
          			F= tempF;
								stage = 2;
								setstage2();
								fafbupdate();
								// viewer.core.is_animating = false;
        			}

        			igl::cotmatrix(U,F,L);
        			wl = wl*sl;
        			ArrayXf interringarea = ArrayXf::Zero(vcount);
        			for(int i=0;i<F.rows(); i++)
        			{
          			interringarea(F(i,0)) += interarea(i);
          			interringarea(F(i,1)) += interarea(i);
          			interringarea(F(i,2)) += interarea(i);
        			}

        			for(int i=0; i<vcount; i++)
        			{
          			wh(i,i) = whi(i,i) *sqrt(ringarea(i)/interringarea(i));
        			}

        			if(something_collapsed)
        			{
								// cout<<U<<endl<<endl;
          			viewer.data().clear();
          			viewer.data().set_mesh(U,F);
          			viewer.data().set_face_based(true);
          			// viewer.core.is_animating = false;
        			}
					}
					else if(stage == 2)
					{
							int rowi=0,collapse_i=0,collapse_j=0;
							double mincost=0;
							edgecosts = VectorXd::Zero(E.rows());
							for(int i=0; i<E.rows(); i++)
							{
								if(!(E(i,0)== 0 && E(i,1)== 0))
								{
									//collapse collapse_i to collapse_j
									//Min of fbedge(i,0) holds fb of (0->1) collapse and (i,1) holds (1->0) collapse.
									//So we chose the smallest to be ci and cj and as collapse_i and collapse_j
									int ci = fbedge(i,0) < fbedge(i,1) ? 0 : 1;
									int cj = fbedge(i,0) < fbedge(i,1) ? 1 : 0;
									double totalcost = wa * (fa(i,0)+fa(i,1)) + wb * fbedge(i,ci);

									if(EF(i,1)==-1)
										totalcost *= 2;

									if(i == 0 || mincost == 0)
									{
										mincost = totalcost;
										collapse_i=ci;collapse_j=cj;rowi=i;
									}
									else if(totalcost < mincost)
									{
										mincost = totalcost;
										rowi = i;collapse_i = ci;collapse_j=cj;
									}
									edgecosts(i) = totalcost;
									// cout<<totalcost<<"    ";
								}
							}
							// cout<<endl<<"Mincost: "<<mincost<<endl<<endl;


						//Collapsing the edge
						// cout<<"Collapsing Edge: "<<rowi<<" - ("<<E(rowi,0)<<","<<E(rowi,1)<<") from "<<E(rowi,collapse_i)<<" to "<<E(rowi,collapse_j)<<endl;
						// cout<<"   Faces: ";
						//Faces with the edge is made zero and the face with only source is changed to have the destination

						//Calculating mid-point for half edge collapse. x1+x2 / 2 and y1+y2/2
						double midx = (U(E(rowi,collapse_i),0) + U(E(rowi,collapse_j),0))/2;
						double midy = (U(E(rowi,collapse_i),1) + U(E(rowi,collapse_j),1))/2;
						double midz = (U(E(rowi,collapse_i),2) + U(E(rowi,collapse_j),2))/2;

						//Update the vertex at E(rowi,collapse_j) position to mid point of collapse_i and collapse_j.
						U(E(rowi,collapse_j),0) = midx;
						U(E(rowi,collapse_j),1) = midy;
						U(E(rowi,collapse_j),2) = midz;

						for(int i=0;i<F.rows();i++)
						{
							if(!(F(i,0)== 0 && F(i,1)== 0 && F(i,2)== 0))
							{
								if((F(i,0)==E(rowi,0)||F(i,1)==E(rowi,0)||F(i,2)==E(rowi,0))&&(F(i,0)==E(rowi,1)||F(i,1)==E(rowi,1)||F(i,2)==E(rowi,1)))
								{
									// cout<<i<<" ";
									// cout<<"**********    "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
									F(i,0) = 0;F(i,1)= 0;F(i,2) = 0;
								}
								else
								{
										if(F(i,0) == E(rowi,collapse_i))// && F(i,1) != E(rowi,collapse_j) && F(i,2) != E(rowi,collapse_j))
										{
												F(i,0) = E(rowi,collapse_j);
												// cout<<"$$$$$$$$$$   "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
										}
										else if(F(i,1) == E(rowi,collapse_i))// && F(i,0) != E(rowi,collapse_j) && F(i,2) != E())
										{
												F(i,1) = E(rowi,collapse_j);
												// cout<<"$$$$$$$$$   "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
										}
										else if(F(i,2) == E(rowi,collapse_i))
										{
											  F(i,2) = E(rowi,collapse_j);
												// cout<<"$$$$$$$$$   "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
										}
								}
							}
						}

						//Have to remove duplicates from faces.
						// removeDupFaces();
						// There will not be any duplicate faces.

						// Looking up if any of collapse_i entries of history already exists in history at collapse_j. If not we push them all in history.
						for(int iter=0;iter<history.at(E(rowi,collapse_i)).size();iter++)
						{
								if(find(history.at(E(rowi,collapse_j)).begin(),history.at(E(rowi,collapse_j)).end(),history.at(E(rowi,collapse_i)).at(iter)) == history.at(E(rowi,collapse_j)).end())
								{
										history.at(E(rowi,collapse_j)).push_back(history.at(E(rowi,collapse_i)).at(iter));
								}
						}

							// cout<<"    Modified: ";
							//make the edge zero and change existing edges that has i to j

							// cout<<"Edge to be REMOVED: "<<E(rowi,0)<<","<<E(rowi,1)<<"   "<<E(rowi,collapse_i)<<" to "<<E(rowi,collapse_j)<<endl;

							for(int i=0; i<E.rows(); i++)
							{
								if(i != rowi && !(E(i,0) == 0 && E(i,1)== 0))
								{
										if(E(i,0) == E(rowi,collapse_i) && E(i,1) != E(rowi,collapse_j))
										{
												E(i,0) = E(rowi,collapse_j);
												// cout<<"^^^^^^^^^^^^   "<<E(i,0)<<" "<<E(i,1)<<endl;
										}
										else if(E(i,1) == E(rowi,collapse_i) && E(i,0) != E(rowi, collapse_j))
										{
												E(i,1) = E(rowi, collapse_j);
												// cout<<"^^^^^^^^^^^   "<<E(i,0)<<" "<<E(i,1)<<endl;
										}
								}
							}


							//Removing duplicate edges
							removeDupEdges();

							// cout<<"Back in the predraw"<<endl;

							//Update the qmatrix and the fa of destination.
							// VectorXd pt(4);
							int ind = E(rowi,collapse_j);
							// homogeneousV = U.rowwise().homogeneous();
							// pt << homogeneousV(ind,0), homogeneousV(ind,1), homogeneousV(ind,2), homogeneousV(ind,3);
							qmatrix.at(ind)+= qmatrix.at(E(rowi,collapse_i));
							// fa(ind) = pt.transpose() * qmatrix.at(ind) * pt;

							fafbupdate();

							//zeroflag = false; //When all faces are zero flag remains false and while loop breaks. If a face is present, the flag becomes true.
							int nonzerofaces = 0;
							for(int iter=0;iter<F.rows();iter++)
							{
								if(!(F(iter,0)== 0 && F(iter,1)== 0 && F(iter,2)== 0))
								{
									nonzerofaces++;
									//cout<<"Non zero face is: "<<F(iter,0)<<","<<F(iter,1)<<","<<F(iter,2)<<endl;
									//zeroflag = true; break;
								}
							}

							//Removing the edges
							E(rowi,0) = 0; E(rowi,1) = 0;

							// viewer.core.is_animating=false;

							viewer.data().clear();
							viewer.data().set_mesh(U,F);
							// viewer.data().set_face_based(true);
							// viewer.core.is_animating = false;

							int nonzeroedges = 0;
							for(int iter=0;iter<E.rows();iter++)
							{
								if(!(E(iter,0)== 0 && E(iter,1)== 0))
								{
									nonzeroedges++;
									//cout<<"Non zero face is: "<<F(iter,0)<<","<<F(iter,1)<<","<<F(iter,2)<<endl;
									//zeroflag = true; break;
								}
							}



							cout<<"faces: "<<F.rows()<<"   rem: "<<nonzerofaces;
							cout<<"     edges: "<<E.rows()<<"   rem: "<<nonzeroedges<<endl;

							if(nonzerofaces == 0)
							{
								stage = 3;
								setstage3();

								// MatrixXi E1;
								// igl::edge_flaps(F,E1,EMAP,EF,EI);
								//
								// cout<<F<<endl<<endl;
								// cout<<E<<endl<<endl;
								// cout<<E1<<endl<<endl;
								// cout<<EF<<endl<<endl;
								// cout<<EI<<endl<<endl;

								// viewer.core.is_animating = false;
								viewer.data().clear();
								viewer.data().set_mesh(U,F);
								// viewer.data().set_face_based(true);
									// zeroflag = true;
								 // cout<<"Number of faces: "<<nonzerofaces<<endl;//<<endl<<endl;

							}
					}
					else if(stage == 3)
					{
						cout<<"In stage 3: "<<endl;
						cout<<"History: "<<history.size()<<endl;
						cout<<"Vertices: "<<V.rows()<<endl;
						cout<<"size: "<<stage2nodes.size()<<endl;
						VectorXd distances = VectorXd::Zero(V.rows());

						for(int i=0;i<stage2nodes.size();i++)  // For each skeleton node
						{
								if(stage2junct.at(i).size() >= 2) // If it is a non-terminal node
								{
										double alllooplen = 0.0;
										Vector3d junctdisp = {0,0,0};
										for(int nodeiter = 0; nodeiter < stage2junct.at(i).size(); nodeiter++) // For each boundary
										{
												vector<int> currentboundary;
												double lensum = 0;
												Vector3d dispavg={0,0,0};
												int adjnode = stage2junct.at(i).at(nodeiter); //Current adjacent node for calculating boundary dj

												for(int biter=0; biter < boundaryvertices.at(i).size(); biter++) // For each boundary vertex of the current skeleton node,
												{
														int stageindex = find(stage2nodes.begin(), stage2nodes.end(), adjnode) - stage2nodes.begin();
														for(int niter=0; niter < boundaryvertices.at(stageindex).size(); niter++)	// For each boundary vertex of the adjancent node,
														{
																//If an edge exists, then the boundary vertex of skeleton node in the boundary under consideration.
																if(isedgeExists(boundaryvertices.at(i).at(biter), boundaryvertices.at(stageindex).at(niter)) &&
																		(currentboundary.size() == 0 || find(currentboundary.begin(), currentboundary.end(), boundaryvertices.at(i).at(biter)) == currentboundary.end()))
																{
																		currentboundary.push_back(boundaryvertices.at(i).at(biter));
																}
														}
												}

												for(int viter = 0; viter < currentboundary.size(); viter++) // For each vertex in the current boundary
												{
														int vindex = currentboundary.at(viter);
														Vector3d vi = {V(vindex,0),V(vindex,1),V(vindex,2)}, ui={U(vindex,0),U(vindex,1),U(vindex,2)}; // Vi' and vi
														double totaledgelen = 0.0;
														//caclulating L(i,j)
														int ecount = 0;
														for(int iviter=0; iviter < currentboundary.size(); iviter++) // For each of remaining vertices in current boundary
														{
																if(iviter != viter && isedgeExists(currentboundary.at(viter), currentboundary.at(iviter))) // If edge exists, add it to l(i,j)
																{
																		int uindex = currentboundary.at(iviter);
																		double len = sqrt(pow(V(vindex,0)-V(uindex,0),2)+pow(V(vindex,1)-V(uindex,1),2)+pow(V(vindex,2)-V(uindex,2),2));
																		totaledgelen += len;
																		ecount++;
																}
														}
														if(ecount==1)
															totaledgelen*=2;
														dispavg+= totaledgelen * (ui - vi);
														lensum += totaledgelen;
												}
												//We have sigma(l(i,j) * (vi' - vi)), sigma(l(i,j)) for the current boundary.
												// For testing only node with 2 boundaries being considered for now. Since avg of 2 nodes has to be subtracted. Dividing by 2 here it self.
												dispavg = dispavg / lensum;
												if(stage2junct.at(i).size() == 2)
												{
														dispavg = dispavg / 2;
														U(stage2nodes.at(i),0) -= dispavg(0);
						 							 	U(stage2nodes.at(i),1) -= dispavg(1);
						 							 	U(stage2nodes.at(i),2) -= dispavg(2);
												}
												else if(stage2junct.at(i).size() > 2)
												{
														double looplength = 0.0;
														for(int viter = 0; viter<currentboundary.size(); viter++) //For each vertex in the current boundary
														{
																int vindex = currentboundary.at(viter);
																for(int iviter = viter+1; iviter < currentboundary.size()-1; iviter++) //FOr all the subsequent vertices in the boundary
																{
																		if(isedgeExists(currentboundary.at(viter),currentboundary.at(iviter)))
																		{
																				int uindex = currentboundary.at(iviter);
																				double len=sqrt(pow(V(vindex,0)-V(uindex,0),2)+pow(V(vindex,1)-V(uindex,1),2) + pow(V(vindex,2)-V(uindex,2),2));
																				looplength += len;
																		}
																}
														}
														dispavg *= looplength;
														junctdisp += dispavg;
														alllooplen += looplength;
												}
										}
										if(stage2junct.at(i).size() > 2)
										{
												junctdisp = junctdisp / alllooplen;
												U(stage2nodes.at(i),0) += junctdisp(0);
												U(stage2nodes.at(i),1) += junctdisp(1);
												U(stage2nodes.at(i),2) += junctdisp(2);
										}
								}

								else
								{ //If terminal node.

							 			Vector3d dispavg={0,0,0};
							 			for(int j=0;j<history.at(stage2nodes.at(i)).size();j++)
							 			{
								 				int vindex = history.at(stage2nodes.at(i)).at(j);
												// cout<<"Stagenodes: "<<stage2nodes.at(i)<<"   History entry: "<<history.at(stage2nodes.at(i)).at(0)<<"   vindex: "<<vindex<<endl;
												// cout<<V(vindex,0)<<" "<<V(vindex,1)<<" "<<V(vindex,2)<<" ----  "<<U(vindex,0)<<" "<<U(vindex,1)<<" "<<U(vindex,2)<<endl;
												Vector3d vi = {V(vindex,0),V(vindex,1),V(vindex,2)}, ui={U(vindex,0),U(vindex,1),U(vindex,2)};
												// cout<<vi-ui<<endl;
												dispavg += vi-ui;
							 			}
							 			// cout<<dispavg.transpose()<<endl;
							 			// cout<<U(stage2nodes.at(i),0)<<" "<<U(stage2nodes.at(i),1)<<" "<<U(stage2nodes.at(i),2)<<endl;
							 			// cout<<dispavg.transpose()<<"   ";
							 			dispavg = dispavg/history.at(stage2nodes.at(i)).size();
							 			// cout<<dispavg.transpose()<<endl;
							 			U(stage2nodes.at(i),0) -= dispavg(0);
							 			U(stage2nodes.at(i),1) -= dispavg(1);
							 			U(stage2nodes.at(i),2) -= dispavg(2);
							 			// cout<<dispavg.transpose() <<endl;
							 			// cout<<"After: "<<U(stage2nodes.at(i),0)<<" "<<U(stage2nodes.at(i),1)<<" "<<U(stage2nodes.at(i),2)<<endl;
							 			// U(stage2nodes.at(i)) = U(stage2nodes.at(i)) - dispavg;
								}
					}



						// for(int vindex=0; vindex<V.rows(); vindex++)
						// {
						// 		double dist = sqrt(pow(V(vindex,0)-U(vindex,0),2) + pow(V(vindex,1)-U(vindex,1),2) + pow(V(vindex,2)-U(vindex,2),2));
						// 		distances(vindex) = dist;
						// }
						//
						//Distances from node to Vertices
						for(int ind=0;ind<stage2nodes.size();ind++)
						{
								for(int vindex=0; vindex < history.at(stage2nodes.at(ind)).size(); vindex++)
								{
										int nodei = stage2nodes.at(ind);
										int nodej = history.at(nodei).at(vindex);
										double dist = sqrt(pow((V(nodej,0)-U(nodei,0)),2) + pow((V(nodej,1) - U(nodei,1)),2) + pow((V(nodej,2) - U(nodei,2)),2));
										distances(nodej) = dist;
								}
						}

						// cout<<stage2junct.size()<<endl;
						// for(int i=0; i<stage2junct.size();i++)
						// {
						// 		for(int j=0; j<stage2junct.at(i).size();j++)
						// 		{
						// 				cout<<stage2junct.at(i).at(j)<<" ";
						// 		}
						// 		cout<<endl;
						// }
						// cout<<stage2nodes.size()<<endl;

						cout<<"Junction nodes: "<<endl;
						for(int i=0; i<stage2nodes.size(); i++) //Iterate through all the nodes in 1d skeleton
						{
								bool fwdflag = true;
								// cout<<"Node: "<<stage2nodes.at(i)<<endl;
								while(fwdflag)
								{
								if(stage2junct.at(i).size() > 2) // If a node i is junction node
								{
										double stdDev = 0.0,distsum=0,distavg=0;
										for(int j=0; j<history.at(stage2nodes.at(i)).size();j++) //Sum of distances from nodes that are merged to node i to the node i
										{
												double disp=0.0;
												int ind = history.at(stage2nodes.at(i)).at(j);
												distsum += distances(history.at(stage2nodes.at(i)).at(j));
										}
										distavg = distsum / history.at(stage2nodes.at(i)).size(); //Average of distances
										double deviationsqr = 0.0;
										for(int j=0; j<history.at(stage2nodes.at(i)).size(); j++)
										{
												deviationsqr += pow((distavg - distances(history.at(stage2nodes.at(i)).at(j))),2);
										}
										deviationsqr /= history.at(stage2nodes.at(i)).size();
										stdDev = sqrt(deviationsqr); //Standard deviation of node distances of node i
										int nodeindex = -1;
										double minstddev=0;
										for(int j=0; j<stage2junct.at(i).size(); j++) //Iterate through all the nodes adjacent to node i
										{
												double distavg2 = distsum,stdDevM=0;
												deviationsqr = 0.0;
												VectorXd distances2 = VectorXd::Zero(history.at(stage2junct.at(i).at(j)).size());
												for(int k=0; k<history.at(stage2junct.at(i).at(j)).size(); k++)
												{
														int nodei = stage2nodes.at(i);
														int nodej = history.at(stage2junct.at(i).at(j)).at(k);
														double dist = sqrt(pow((V(nodej,0)-U(nodei,0)),2) + pow((V(nodej,1) - U(nodei,1)),2) + pow((V(nodej,2) - U(nodei,2)),2));
														distances2(k) = dist;
														// distances(nodej) = dist;
														distavg2 += dist;
														//distances(history.at(stage2junct.at(i).at(j)).at(k));
												}
												distavg2 /= (history.at(stage2nodes.at(i)).size() + history.at(stage2junct.at(i).at(j)).size());
												double deviationsqr = 0.0;
												for(int k=0; k<history.at(stage2nodes.at(i)).size(); k++)
												{
														deviationsqr += pow((distavg2 - distances(history.at(stage2nodes.at(i)).at(k))),2);
												}
												for(int k=0; k<history.at(stage2junct.at(i).at(j)).size(); k++)
												{
														deviationsqr += pow((distavg2 - distances2(k)),2);
												}
												deviationsqr /= (history.at(stage2nodes.at(i)).size() + history.at(stage2junct.at(i).at(j)).size());
												stdDevM = sqrt(deviationsqr);
												cout << stage2junct.at(i).at(j)<<"      "<<stdDev << "     " <<stdDevM << "     "<< 0.96 * stdDev << endl;

												if(stdDevM < 0.96 * stdDev)
												{
														if(nodeindex == -1 || minstddev > stdDevM)
														{
																nodeindex = j; minstddev = stdDevM;
														}
												}
										}
										cout<<"Merge: node index "<<nodeindex<<"   with std dev: "<<minstddev<<endl;
										if(nodeindex >= 0)
										{
												cout<<"Merging nodes"<<endl;
												int nodevertex = stage2junct.at(i).at(nodeindex);

												for(int kk=0; kk < E.rows(); kk++)
												{
														if(!(E(kk,0)==0 && E(kk,1)==0))
														{
																if(E(kk,0) == nodevertex && E(kk,1) != stage2nodes.at(i))
																{
																		E(kk,0) = stage2nodes.at(i);
																}
																else if(E(kk,1) == nodevertex && E(kk,0) != stage2nodes.at(i))
																{
																		E(kk,1) = stage2nodes.at(i);
																}
																else if((E(kk,0) == nodevertex && E(kk,1) == stage2nodes.at(i)) || (E(kk,1) == nodevertex && E(kk,0) == stage2nodes.at(i)))
																{
																		E(kk,0) = 0; E(kk,1) = 0;
																}
														}
												}

												//Merge the vertices collapsed in History.
												for(int kk=0; kk<history.at(nodevertex).size(); kk++)
												{
														if(find(history.at(stage2nodes.at(i)).begin(), history.at(stage2nodes.at(i)).end(),history.at(nodevertex).at(kk)) == history.at(stage2nodes.at(i)).end())
														{
																history.at(stage2nodes.at(i)).push_back(history.at(nodevertex).at(kk));
														}
												}

												int mergedindex = find(stage2nodes.begin(), stage2nodes.end(), nodevertex)-stage2nodes.begin();
												for(int kk=0; kk<stage2junct.at(mergedindex).size(); kk++)
												{
														if(stage2junct.at(mergedindex).at(kk) != stage2nodes.at(i) &&
															find(stage2junct.at(i).begin(), stage2junct.at(i).end(), stage2junct.at(mergedindex).at(kk)) == stage2junct.at(i).end())
														{
																stage2junct.at(i).push_back(stage2junct.at(mergedindex).at(kk));
														}
												}
												stage2junct.at(i).erase(stage2junct.at(i).begin()+nodeindex);
												stage2junct.erase(stage2junct.begin()+mergedindex);
												stage2nodes.erase(stage2nodes.begin()+mergedindex);
												// distances.erase();

												for(int ind=0;ind<stage2nodes.size();ind++)
												{
														for(int vindex=0; vindex < history.at(stage2nodes.at(ind)).size(); vindex++)
														{
																int nodei = stage2nodes.at(ind);
																int nodej = history.at(nodei).at(vindex);
																double dist = sqrt(pow((V(nodej,0)-U(nodei,0)),2) + pow((V(nodej,1) - U(nodei,1)),2) + pow((V(nodej,2) - U(nodei,2)),2));
																distances(nodej) = dist;
														}
												}
										}
										else
										{
												fwdflag = false;
										}
									}
									else
									{
										fwdflag = false;
									}
								}
								// cout<<stage2nodes.at(i);
								// cout<<" "<<stage2junct.at(i).size()<<endl;
						}

						int minindexv = 0;
						int maxindexv = 0;
						int minindexh = 0;
						int maxindexh = 0;
						for(int ind=0; ind< stage2nodes.size(); ind++)
						{
								if(U(stage2nodes.at(minindexv),1) > U(stage2nodes.at(ind),1))
								{
										minindexv = ind;
								}
								else if(U(stage2nodes.at(maxindexv),1) < U(stage2nodes.at(ind),1))
								{
										maxindexv = ind;
								}
								if(U(stage2nodes.at(minindexh),0) > U(stage2nodes.at(ind),0))
								{
										minindexh = ind;
								}
								else if(U(stage2nodes.at(maxindexh),0) < U(stage2nodes.at(ind),0))
								{
										maxindexh = ind;
								}
						}

						cout<<"Minh: "<<minindexh<<endl;
						cout<<"Maxh: "<<maxindexh <<endl;
						cout<<"Minv: "<<minindexv <<endl;
						cout<<"Maxv: "<<maxindexv <<endl;
						cout<<"Degree minh and maxh: "<<stage2junct.at(minindexh).size()<<" "<<stage2junct.at(maxindexh).size()<<endl;

						//remove unnecessary leaf nodes by removing them if their neighbor's degree is >2
						// for(int ind=0; ind < stage2junct.size(); ind++)
						// {
						// 		if(stage2junct.at(ind).size()==1 )
						// 		{
						// 				cout<<"one degree"<<endl;
						// 				int neighbornode = stage2junct.at(ind).at(0);
						// 				int neighborindex = find(stage2nodes.begin(),stage2nodes.end(),neighbornode)-stage2nodes.begin();
						// 				cout<<stage2junct.size()<<" n "<<neighborindex<<endl;
						//
						// 				if(stage2junct.at(neighborindex).size() > 2) //neighbor is a junction node.
						// 				{
						// 						//Update edges
						// 						int nodevertex = stage2nodes.at(ind);
						// 						for(int kk=0; kk < E.rows(); kk++)
						// 						{
						// 								if(!(E(kk,0)==0 && E(kk,1)==0))
						// 								{
						// 										if(E(kk,0) == nodevertex && E(kk,1) != neighbornode)
						// 										{
						// 												E(kk,0) = neighbornode;
						// 										}
						// 										else if(E(kk,1) == nodevertex && E(kk,0) != neighbornode)
						// 										{
						// 												E(kk,1) = neighbornode;
						// 										}
						// 										else if((E(kk,0) == nodevertex && E(kk,1) == neighbornode) || (E(kk,1) == nodevertex && E(kk,0) == neighbornode))
						// 										{
						// 												E(kk,0) = 0; E(kk,1) = 0;
						// 										}
						// 								}
						// 						}
						// 						for(int hind=0;hind<history.at(stage2nodes.at(ind)).size();hind++)
						// 						{
						// 								if(find(history.at(neighbornode).begin(), history.at(neighbornode).end(), history.at(stage2nodes.at(ind)).at(hind)) == history.at(neighbornode).end())
						// 								{
						// 									history.at(neighbornode).push_back(history.at(stage2nodes.at(ind)).at(hind));
						// 								}
						//
						// 						}
						//
						// 						//We need not update the neighborhood of newly merged index because, the neighborhood of the source node is only 1
						//
						// 						cout<<"Testing testing: "<<ind<<endl;
						// 						cout<<"Testing testing: "<<stage2nodes.at(3)<<endl;
						// 						int indexinneigh = find(stage2junct.at(neighborindex).begin(),stage2junct.at(neighborindex).end(),stage2nodes.at(ind))-stage2junct.at(neighborindex).begin();
						// 						stage2junct.at(neighborindex).erase(stage2junct.at(neighborindex).begin()+indexinneigh);
						//
						// 						stage2junct.erase(stage2junct.begin()+ind);
						// 						stage2nodes.erase(stage2nodes.begin()+ind);
						// 						ind--;
						// 						cout<<stage2junct.size()<<"  "<<ind+1<<endl;
						// 				}
						// 		}
						// }


						ofstream op;
						op.open("../1d/skeleton.raw");
						for(int i=0;i<E.rows();i++)
						{
							if(!(E(i,0)== 0 && E(i,1)== 0))
							{
								op<<U(E(i,0),0)<<" "<<U(E(i,0),1)<<" "<<U(E(i,0),2)<<endl;
								op<<U(E(i,1),0)<<" "<<U(E(i,1),1)<<" "<<U(E(i,1),2)<<endl;
							}
						}
						op.close();
						float colorvs[6][3] =
						{
								1.0,0.0,0.0,
								0.0,1.0,0.0,
								0.0,0.0,1.0,
								1.0,1.0,0.0,
								1.0,0.0,1.0,
								0.0,1.0,1.0
						};
						VectorXi colorcodes = VectorXi::Zero(V.rows());
						int cc=0;
						for(int i=0; i<stage2nodes.size(); i++)
						{
								for(int j=0; j< history.at(stage2nodes.at(i)).size();j++)
								{
											colorcodes(history.at(stage2nodes.at(i)).at(j)) = cc;
								}
								cc = (cc+1)%6;
						}

						op.open("../2d/skeleton.raw");
						for(int i=0; i<E2.rows(); i++)
						{
								op<<V(E2(i,0),0)<<" "<<V(E2(i,0),1)<<" "<<V(E2(i,0),2)<<endl;
								int q = colorcodes(E2(i,0));
								op<<colorvs[q][0]<<" "<<colorvs[q][1]<<" "<<colorvs[q][2]<<endl;
								op<<V(E2(i,1),0)<<" "<<V(E2(i,1),1)<<" "<<V(E2(i,1),2)<<endl;
								q = colorcodes(E2(i,1));
								op<<colorvs[q][0]<<" "<<colorvs[q][1]<<" "<<colorvs[q][2]<<endl;
						}

						for(int i=0;i<E.rows();i++)
						{
							if(!(E(i,0)== 0 && E(i,1)== 0))
							{
								op<<U(E(i,0),0)<<" "<<U(E(i,0),1)<<" "<<U(E(i,0),2)<<endl;
								op<<1.0<<" "<<1.0<<" "<<1.0<<endl;
								op<<U(E(i,1),0)<<" "<<U(E(i,1),1)<<" "<<U(E(i,1),2)<<endl;
								op<<1.0<<" "<<1.0<<" "<<1.0<<endl;
							}
						}

						op.close();

						op.open("../1d/vertexdata.txt");
						op<<V.rows()<<endl;
						for(int kind = 0; kind < V.rows(); kind++)
						{
								op<<V(kind,0)<<" "<<V(kind,1)<<endl;
						}
						op<<U.rows()<<endl;
						for(int kind = 0; kind < U.rows(); kind++)
						{
								op<<U(kind,0)<<" "<<U(kind,1)<<endl;
						}
						op<<E.rows()<<endl;
						for(int kind = 0; kind < E.rows(); kind++)
						{
								op<<E(kind,0)<<" "<<E(kind,1)<<endl;
						}
						op<<E2.rows()<<endl;
						for(int kind=0;kind < E2.rows(); kind++)
						{
								op<<E2(kind,0)<<" "<<E2(kind,1)<<endl;
						}
						op.close();

						op.open("../1d/skeldata.txt");
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
						viewer.core.is_animating = false;
					}
      }
      return false;

   };

   const auto &key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod) -> bool
   {
      switch(key)
      {
        case ' ':
          viewer.core.is_animating ^= 1;
          break;
        default: return false;
      }
      return true;
   };

   reset();
   viewer.core.background_color.setConstant(1);
   viewer.core.is_animating = true;
   viewer.callback_key_down = key_down;
   viewer.callback_pre_draw = pre_draw;
   return viewer.launch();

}
