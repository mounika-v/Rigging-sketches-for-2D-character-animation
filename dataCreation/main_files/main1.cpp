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


#include <igl/boundary_facets.h>
#include <igl/unique.h>


using namespace std;
using namespace Eigen;

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
		// Compute edge lengths
		//Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 3> l;
		//igl::edge_lengths(V,F,l);

		//igl::doublearea(l,0.0,dblA);

		// Projected area helper
		const auto & proj_doublearea =
		[&V,&F](const int x, const int y, const int f)
		->typename DerivedV::Scalar
		{
			auto rx = V(F(f,0),x)-V(F(f,2),x);
			auto sx = V(F(f,1),x)-V(F(f,2),x);
			auto ry = V(F(f,0),y)-V(F(f,2),y);
			auto sy = V(F(f,1),y)-V(F(f,2),y);
			double area = rx*sy - ry*sx;
			// if(isnan(area))
			// 	return 0.0;
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


MatrixXd V,U,tempU;
MatrixXi F,tempF;
SparseMatrix<double> L;
VectorXd area,interarea;
igl::opengl::glfw::Viewer viewer;
int sl = 2.0,vcount;
double areathreshold;


//testing edge collapse
VectorXi EMAP;
MatrixXi E,EF,EI;
// MatrixXd V,OV;
// MatrixXi F,OF;
// typedef std::set<std::pair<double,int> > PriorityQueue;
// PriorityQueue Q;
// std::vector<PriorityQueue::iterator > Qit;
// // If an edge were collapsed, we'd collapse it to these points:
// MatrixXd C;
// int num_collapsed;


int main(int argc, char * argv[])
{
					// Load a mesh in OFF format
					igl::readOFF("mesh1.off", V, F);

					vcount = V.rows();
					U= V;

					//L matrix
					igl::cotmatrix(V,F,L);

					//find average area for Wl
					igl::doublearea(V,F,area);
					area = area.array() / 2;
					double area_avg   = area.mean();

					cout<<"Total area of original mesh: "<<area.sum()<<endl;
					cout<<"V:: "<<V.rows()<<" "<<V.cols()<<endl;
					cout<<"F:: "<<F.rows()<<" "<<F.cols()<<endl;

					//RowVectorXd ringarea(vcount);
					ArrayXf ringarea = ArrayXf::Zero(vcount);

					//Define a one_ring area matrix for V
					//We need this to udpate the Wh for each iteration
					//For each face we are adding its area to all 3 vertices. So, all the vertices will have total area of faces sharing it.
					for(int i = 0; i < F.rows(); i++)
					{
						ringarea(F(i,0)) += area(i);
						ringarea(F(i,1)) += area(i);
						ringarea(F(i,2)) += area(i);
					}


					//Create initial wl
					int i;
					MatrixXd wl=MatrixXd::Zero(vcount, vcount);
					for(i=0;i<vcount;i++)
					{
						wl(i,i) = 0.001 * sqrt(area_avg);
					}

					//Initial wh
					MatrixXd wh=MatrixXd::Zero(vcount, vcount);
					for(i=0;i<vcount;i++)
					{
						wh(i,i) = 1.0;
					}

//-----------------------------------------

MatrixXi E1;
igl::boundary_facets(F,E1);
// Find boundary vertices
VectorXi b,IA,IC;
igl::unique(E1,b,IA,IC);

cout<<endl;cout<<"Size of E1: "<<E1.rows()<<"x"<<E1.cols();
cout<<endl;cout<<"Size of b: "<<b.rows()<<"x"<<b.cols();
cout<<endl;cout<<"Size of IA: "<<IA.rows()<<"x"<<IA.cols();
cout<<endl;cout<<"Size of IC: "<<IC.rows()<<"x"<<IC.cols();

//-----------------------------------------


					MatrixXd whi = wh;
					while(true)
					{
						MatrixXd a(vcount,vcount);
						a = wl * L;
						MatrixXd lhs(2*vcount,vcount);
						lhs<< a, wh;
						MatrixXd b(vcount,vcount);
						b = wh * U;
						ArrayXXd zro = ArrayXXd::Zero(vcount,3);
						MatrixXd rhs(2*vcount,3) ;
						rhs << zro, b;

						//Preserving the mesh incase the area threshold become 0 (work around)
						tempU = U;
						tempF = F;

						//Solving for V(t+1)
						ColPivHouseholderQR<MatrixXd> solver(lhs);
						U= solver.solve(rhs);


						//Compare surface area of mesh with original
						mesh::doublearea(U,F,interarea);

						interarea = interarea.array() / 2;

						//cout<<"Intermediate area: "<<interarea.sum()<<endl;

						areathreshold=interarea.sum()/area.sum();
						cout<<"Areathreshold: "<<areathreshold<<endl;
						//cout<<"inter size: "<<interarea.size()<<endl;

						if(areathreshold < 0.0001)
						{

							U = tempU;
							F = tempF;
							break;
						}
						igl::cotmatrix(U,F,L);
						//if(areathreshold < 0.3)
						//	cout<<L<<endl;

						//Update Wl for next iteration
						//Sl is suggested to be taken 2.0
						wl = wl * sl;

						//Calculate new one ring area
						ArrayXf interringarea = ArrayXf::Zero(vcount);
						for(int i = 0; i < F.rows(); i++)
						{
							interringarea(F(i,0)) += interarea(i);
							interringarea(F(i,1)) += interarea(i);
							interringarea(F(i,2)) += interarea(i);
						}

						//Update wh for next iteration
						for(i = 0;i < vcount; i++)
						{
							wh(i,i) = whi(i,i) * sqrt(ringarea(i)/interringarea(i));
						}
					}

//------------------------------------------------------------------------------------------------------
					//Removing all the zero area faces.
					//Compare surface area of mesh with original
					/*mesh::doublearea(U,F,interarea);
					interarea = interarea.array() / 2;


					for(int i=0; i<interarea.size(); i++)
					{
							if(interarea(i) < 0.0001)
							{
									double dist01 = sqrt( pow(( U(F(i,1),0) - U(F(i,0),0) ),2) +  pow(( U(F(i,1),1) - U(F(i,0),1) ),2) + pow(( U(F(i,1),2) - U(F(i,0),2) ),2));
									double dist02 = sqrt( pow(( U(F(i,2),0) - U(F(i,0),0) ),2) +  pow(( U(F(i,2),1) - U(F(i,0),1) ),2) + pow(( U(F(i,2),2) - U(F(i,0),2) ),2));
									double dist12 = sqrt( pow(( U(F(i,1),0) - U(F(i,2),0) ),2) +  pow(( U(F(i,1),1) - U(F(i,2),1) ),2) + pow(( U(F(i,1),2) - U(F(i,2),2) ),2));
									int s,t;
									if(dist01 < dist02 && dist01 < dist12)
									{
										s=0;t=1;
									}
									else if(dist02 <dist01 && dist02 < dist12)
									{
										s=0; t=2;
									}
									else
									{
										s=1; t=2;
									}

									for(int k=0; k<E.rows(); k++)
									{
										if((E(k,1)==F(i,s) && E(k,0)==F(i,t)) || (E(k,1)==F(i,t) && E(k,0) == F(i,s)))
										{
											E(k,0) = 0; E(k,1) = 0;
										}
										else if(E(k,1) == F(i,s))
										{
											E(k,1) = F(i,t);
										}
										else if(E(k,0) == F(i,s))
										{
											E(k,0) = F(i,t);
										}
									}
									int e1 = F(i,s),e2 = F(i,t);
									for(int k = 0; k<F.rows(); k++)
									{
											if((F(k,0)==e1 || F(k,1) == e1 || F(k,2)==e1) && (F(k,0) == e2 || F(k,2) == e2 || F(k,1)==e2))
											{
												F(k,0) = 0; F(k,1)=0; F(k,2) = 0;
											}
											else if(F(k,0) == e1 && F(k,1) != e2 && F(k,2) != e2)
											{
												F(k,0) = e2;
											}
											else if(F(k,1) == e1 && F(k,0) != e2 && F(k,2) != e2)
											{
												F(k,1) = e2;
											}
											else if(F(k,2) == e1 && F(k,1) != e2 && F(k,0) != e2)
											{
												F(k,2) = e2;
											}
									}

							}
					}

//-------------------------------------------------------------------------------------------------------
					/*Second step:*/

				vector <MatrixXd> qmatrix;
				vector <vector<int>> history;

					igl::edge_flaps(F,E,EMAP,EF,EI); //Creating an edge Matrix
					cout<<"Edges: "<<E.rows()<<"  Faces: "<<F.rows()<<endl;
					cout<<"EMAP: "<<EMAP.rows()<<","<<EMAP.cols()<<endl;
					cout<<"EF: "<<EF.rows()<<","<<EF.cols()<<endl;
					cout<<"EI: "<<EI.rows()<<","<<EI.cols()<<endl;

					// for(int i=0;i<EF.rows();i++)
					// {
					// 	cout<<EF(i,0)<<"   "<<EF(i,1)<<endl;
					// }

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

					//Creating fa for each vertex. For each edge we do fa(i) + fa(j)
					MatrixXd normV = U.rowwise().homogeneous();
					VectorXd fa = VectorXd::Zero(vcount);
					for(i =0; i< vcount; i++)
					{
						VectorXd pt(4);
						pt << normV(i,0), normV(i,1), normV(i,2), normV(i,3);
						fa(i) = pt.transpose() * qmatrix.at(i) * pt;
					}


					// bool zeroflag = true;
					// while(zeroflag)
					// {

					//Creating fb. We are taking two cols of fb because the direction of collapse matters for fb.
					//calculating norms of all the edges.
					VectorXd normvec = VectorXd::Zero(E.rows());
					for(int i=0;i<E.rows();i++)
					{
						if(!(E(i,1)==0 && E(i,0)==0))
						{
								double dist = sqrt( pow(( U(E(i,1),0) - U(E(i,0),0) ),2) +  pow(( U(E(i,1),1) - U(E(i,0),1) ),2) + pow(( U(E(i,1),2) - U(E(i,0),2) ),2) );
						 		normvec(i) = dist;
						}
					}
					//creating fb entry as product of currentedge and sum of all the other edges.
					MatrixXd fbedge(E.rows(),2);
					for(int i=0; i<E.rows(); i++)
					{
						for(int j=0;j<2;j++)
						{
								double magv = 0;
								for(int k=0; k<E.rows(); k++)
								{
									if( (E(k,0)==E(i,j) || E(k,1)==E(i,j)) && (k != i))
										magv += normvec(k);
								}
								fbedge(i,j) = normvec(i) * magv;
						}
					}

					double wa = 1.0, wb = 0.1;
int iterations = 0;

					//For each iteration we will isolate the least cost edge. coli is source. colj is destination.
					bool zeroflag = true;
					VectorXd edgecosts(E.rows());
					while(zeroflag)
					{

						int rowi=0,coli=-1,colj=-1;
						double mincost=0;
						for(int i=0; i<E.rows(); i++)
						{
							if(!(E(i,0)==0 && E(i,1)==0))
							{
								int ci = fbedge(i,0) < fbedge(i,1) ? 0 : 1;
								int cj = fbedge(i,0) < fbedge(i,1) ? 1 : 0;
								double totalcost = wa * (fa(E(i,0))+fa(E(i,1))) + wb * fbedge(i,ci);
								if(i == 0 || mincost == 0)
								{
									mincost = totalcost;
									coli=ci;colj=cj;rowi=i;
								}
								else if(totalcost < mincost)
								{
									mincost = totalcost;
									rowi = i;coli = ci;colj=cj;
								}
								edgecosts(i) = totalcost;
							}
						}


						int mincount = 0;
						for(int riter=0; riter<E.rows(); riter++)
						{
							  if(edgecosts(riter) == mincost)
									mincount++;
						}
						if(mincount>1)
						cout<<"**** Entries with min count: "<<mincount<<endl;

						// //Collapsing the edge
						// cout<<"Collapsing Edge: "<<rowi<<" - ("<<E(rowi,0)<<","<<E(rowi,1)<<")";
						// cout<<"   Faces: ";
						//Faces with the edge is made zero and the face with only source is changed to have the destination


						for(i=0;i<F.rows();i++)
						{
							if(!(F(i,0)==0 && F(i,1)==0 && F(i,2)==0))
							{
								if((F(i,0)==E(rowi,0)||F(i,1)==E(rowi,0)||F(i,2)==E(rowi,0))&&(F(i,0)==E(rowi,1)||F(i,1)==E(rowi,1)||F(i,2)==E(rowi,1)))
								{
									// cout<<i<<" ";
									F(i,0) = 0;F(i,1)=0;F(i,2) = 0;
								}
								else if((F(i,0)==E(rowi,coli)||F(i,1)==E(rowi,coli)||F(i,2)==E(rowi,coli))&&(F(i,0)!=E(rowi,colj)&&F(i,1)!=E(rowi,colj)&&F(i,2)!=E(rowi,colj)))
								{
									if(F(i,0) == E(rowi,coli))
									{
											// bool dupfound=false;
											// for(int riter=0; riter<F.rows(); riter++)
											// {
											// 		if( (E(rowi,colj) == F(riter,0) || E(rowi,colj) == F(riter, 1) || E(rowi,colj) == F(riter,2)) &&
											// 				(F(i,1) == F(riter,0) || F(i,1) == F(riter,1) || F(i,1) == F(riter,2)) &&
											// 				(F(i,2) == F(riter,0) || F(i,2) == F(riter,1) || F(i,2) == F(riter,2)) )
											// 		{
											// 			cout<<"Dup face found"<<F(riter,0)<<","<<F(riter,1)<<","<<F(riter,2)<<endl;
											// 			dupfound = true; break;
											// 		}
											// }
											// if(dupfound)
											// {
											// 	F(i,0) = 0; F(i,1) = 0; F(i,2) = 0;
											// }
											// else
												F(i,0) = E(rowi,colj);
									}
									else if(F(i,1) == E(rowi,coli))
									{
											F(i,1) = E(rowi,colj);
									}
									else if(F(i,2) == E(rowi,coli))
									{
											F(i,2) = E(rowi,colj);
									}
									// cout<<"  -"<<i<<"  ";
								}
							}
						}

						// cout<<"Merging "<<E(rowi,coli)<<" to "<<E(rowi,colj)<<" : "<<history.at(E(rowi,coli)).size()<<" to "<<history.at(E(rowi,colj)).size()<<endl;
						//Looking up if any of coli entries of history already exists in history at colj. If not we push them all in history.
						for(int iter=0;iter<history.at(E(rowi,coli)).size();iter++)
						{
								if(find(history.at(E(rowi,colj)).begin(),history.at(E(rowi,colj)).end(),history.at(E(rowi,coli)).at(iter)) == history.at(E(rowi,colj)).end())
								{
										history.at(E(rowi,colj)).push_back(history.at(E(rowi,coli)).at(iter));
								}
						}

							// cout<<"    Modified: ";
							//make the edge zero and change existing edges that has i to j

							//cout<<"Edge to be REMOVED: "<<E(rowi,0)<<","<<E(rowi,1)<<"   "<<E(rowi,coli)<<" to "<<E(rowi,colj)<<endl;
							for(i=0;i<E.rows();i++)
							{
								if(i!=rowi && !(E(i,0)==0 && E(i,1)==0))
								{
									//cout<<"Edge working on is: "<<E(i,0)<<","<<E(i,1);
									if(E(i,0)==E(rowi,coli) && E(i,1)!=E(rowi,colj))
									{

										//cout<<"  to "<<E(rowi,colj)<<","<<E(i,1);
										// //If duplicate do not add the edge. testing
										// bool existsflag = false;
										// for(int riter = 0; riter< E.rows(); riter++)
										// {
										// 		if((E(riter,0) == E(rowi,colj) && E(riter,1) == E(i,1)) || (E(riter,1) == E(rowi,colj) && E(riter,0) == E(i,1)))
										// 		{	//cout<<"    Edge found as duplicate: "<<E(riter,0)<<","<<E(riter,1);
										// 			existsflag = true; break; }
										// }
										// if(existsflag)
										// {
										// 		E(i,0) = 0; E(i,1) = 0;
										// }
										// else
												E(i,0) = E(rowi,colj);
									}
									else if(E(i,1)==E(rowi,coli) && E(i,0) != E(rowi,colj))
									{

										//cout<<"   to "<<E(i,0)<<","<<E(rowi,colj);
										//If duplicate do not add the edge. testing
										// bool existsflag = false;
										// for(int riter = 0; riter< E.rows(); riter++)
										// {
										// 		if((E(riter,0) == E(rowi,colj) && E(riter,1) == E(i,0)) || (E(riter,1) == E(rowi,colj) && E(riter,0) == E(i,0)))
										// 		{	//cout<<"    Edge found as duplicate: "<<E(riter,0)<<","<<E(riter,1);
										// 			existsflag = true; break;
										// 		}
										// }
										// if(existsflag)
										// {
										// 		E(i,0) = 0; E(i,1) = 0;
										// }
										// else
												E(i,1) = E(rowi,colj);
									}
									// cout<<endl;
								}
							}
							//cout<<endl<<endl;
							E(rowi,0) = 0; E(rowi,1) = 0;
							// cout<<endl<<endl<<endl;

							//Update the qmatrix and the fa of destination.
							VectorXd pt(4);
							int ind = E(rowi,colj);
							pt << normV(ind,0), normV(ind,1), normV(ind,2), normV(ind,3);
							qmatrix.at(ind)+= qmatrix.at(E(rowi,coli));
							fa(ind) = pt.transpose() * qmatrix.at(ind) * pt;

							zeroflag = false; //When all faces are zero flag remains false and while loop breaks. If a face is present, the flag becomes true.
							int nonzerofaces = 0;
							for(int iter=0;iter<F.rows();iter++)
							{
								if(!(F(iter,0)==0&&F(iter,1)==0&&F(iter,2)==0))
								{
									nonzerofaces++;
									//cout<<"Non zero face is: "<<F(iter,0)<<","<<F(iter,1)<<","<<F(iter,2)<<endl;
									//zeroflag = true; break;
								}
							}
							if(nonzerofaces > 0)
							{
									zeroflag = true;
								 // cout<<"Number of faces: "<<nonzerofaces<<endl;//<<endl<<endl;
							}


							// int nonzeroedges = 0;
							// for(int riter = 0; riter < E.rows(); riter++)
							// {
							// 		if(!(E(riter,0) == 0 && E(riter,1) == 0))
							// 			nonzeroedges++;
							// }
							// cout<<"Non zero edges: "<< nonzeroedges <<endl;


iterations++;
					}
					cout<<endl;

cout<<"No. of iterations: "<<iterations<<endl;


					int count=0;
					vector<int> remain;
					for(int i=0;i<E.rows();i++)
					{
						if(!(E(i,0)==0 && E(i,1)==0))
						{
							if(find(remain.begin(),remain.end(),E(i,0)) == remain.end())
							{
								remain.push_back(E(i,0));
							}
							else if(find(remain.begin(),remain.end(),E(i,1)) == remain.end())
							{
								remain.push_back(E(i,1));
							}
							// count++;
							// cout<<count<<"     "<<history.at(E(i,0)).size()<<" , "<<history.at(E(i,1)).size()<<endl;
						}
					}

					/*Third step

					for(int i=0; i<remain.size(); i++)
					{
							float xdist = 0,ydist=0;
							float dist = 0;
							vector<int> indlist = history.at(remain.at(i));
							for(int j=1; j<indlist.size(); j++)
							{
								dist += sqrt(pow((U(indlist.at(j),0) - V(indlist.at(j),0)),2) + pow((U(indlist.at(j),1) - V(indlist.at(j),1)),2));
								xdist += (U(indlist.at(j),0) - V(indlist.at(j),0));
								ydist += (U(indlist.at(j),1) - V(indlist.at(j),1));
							}

							// option1
							xdist /= (indlist.size()-1);ydist /= (indlist.size()-1);
							U(remain.at(i),0) -= xdist; U(remain.at(i),1) -= ydist;

							// option2
							//dist /= (length of sum of two adjacent edges) (?)
					} */



					ofstream op;
					op.open("../1d/skeleton.raw");
					for(int i=0;i<E.rows();i++)
					{
						if(!(E(i,0)==0 && E(i,1)==0))
						{
							op<<U(E(i,0),0)<<" "<<U(E(i,0),1)<<" "<<U(E(i,0),2)<<endl;
							op<<U(E(i,1),0)<<" "<<U(E(i,1),1)<<" "<<U(E(i,1),2)<<endl;
						}
					}
					op.close();


		// Plot the mesh (Error: The new mesh has a different number of vertices/faces. Please clear the mesh before plotting.) doubt
		viewer.data().clear();
		viewer.data().set_mesh(U, F);
		return viewer.launch();
}
