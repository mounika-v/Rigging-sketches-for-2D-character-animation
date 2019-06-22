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

   igl::opengl::glfw::Viewer viewer;

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
	 double wa = 1.0, wb = 0.1;
	 VectorXd edgecosts;
	 MatrixXd fbedge;
	 MatrixXd homogeneousV;
	 MatrixXd fa;

   const auto & reset = [&]()
   {
     F=OF;
     V=OV;
     num_collapsed = 0;
     vcount = V.rows();
     U=V;

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

				//History matrix to maintain the list of vertices collapsed.
				// for(int i=0;i<U.rows(); i++)
				// {
				// 	vector<int> temp;
				// 	temp.push_back(i);
				// 	history.push_back(temp);
				// }


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
									VectorXd vi;
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


				// cout<<"Creating FA for vertices"<<endl;

				// bool zeroflag = true;
				// while(zeroflag)
				// {

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
	 {
		 	// cout<<"Removing duplicating faces ^^^^^^^^^^"<<endl;
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
	 {
		 	// cout<<"Removing duplicating edges ^^^^^^^^^^"<<endl;
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
          			viewer.core.is_animating = false;

								setstage2();
								fafbupdate();
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
          			viewer.core.is_animating = false;
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
						cout<<"Collapsing Edge: "<<rowi<<" - ("<<E(rowi,0)<<","<<E(rowi,1)<<") from "<<E(rowi,collapse_i)<<" to "<<E(rowi,collapse_j)<<endl;
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
										if(F(i,0) == E(rowi,collapse_i))// && F(i,1) != E(rowi,colj) && F(i,2) != E(rowi,colj))
										{
												F(i,0) = E(rowi,collapse_j);
												// cout<<"$$$$$$$$$$   "<<F(i,0)<<" "<<F(i,1)<<" "<<F(i,2)<<endl;
										}
										else if(F(i,1) == E(rowi,collapse_i))// && F(i,0) != E(rowi,colj) && F(i,2) != E())
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

						// cout<<"Merging "<<E(rowi,coli)<<" to "<<E(rowi,colj)<<" : "<<history.at(E(rowi,coli)).size()<<" to "<<history.at(E(rowi,colj)).size()<<endl;
						//Looking up if any of coli entries of history already exists in history at colj. If not we push them all in history.
						// for(int iter=0;iter<history.at(E(rowi,coli)).size();iter++)
						// {
						// 		if(find(history.at(E(rowi,colj)).begin(),history.at(E(rowi,colj)).end(),history.at(E(rowi,coli)).at(iter)) == history.at(E(rowi,colj)).end())
						// 		{
						// 				history.at(E(rowi,colj)).push_back(history.at(E(rowi,coli)).at(iter));
						// 		}
						// }

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
							viewer.data().set_face_based(true);
							viewer.core.is_animating = false;

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



							cout<<"Total: "<<F.rows()<<"      Non zero faces: "<<nonzerofaces<<endl;
							cout<<"Total: "<<E.rows()<<"      Non zero edges: "<<nonzeroedges<<endl;

							if(nonzerofaces == 0)
							{
								stage = 3;
								for(int iter=0;iter<E.rows();iter++)
								{
									if(!(E(iter,0)==0&&E(iter,1)==0))
									{
										cout<<"Edge: "<<E(iter,0)<<" "<<E(iter,1)<<endl;
										//cout<<"Non zero face is: "<<F(iter,0)<<","<<F(iter,1)<<","<<F(iter,2)<<endl;
										//zeroflag = true; break;
									}
								}

								// viewer.core.is_animating = false;
								// viewer.data().clear();
								// viewer.data().set_mesh(U,F);
								// viewer.data().set_face_based(true);
									// zeroflag = true;
								 // cout<<"Number of faces: "<<nonzerofaces<<endl;//<<endl<<endl;

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
							}
							/*else
							{
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
										//viewer.core.is_animating^=1;
							}*/
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
