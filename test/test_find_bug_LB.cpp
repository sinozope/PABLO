/*
 * test_find_LB_bug.cpp
 *
 *  Created on: 7 sept. 2015
 *      Author: marco
 */

#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if NOMPI==0
        MPI::Init(argc, argv);

        {
#endif
                int iter = 0;
                int dim = 2;
                int ga = 2;

                /**<Instantation of a 2D para_tree object.*/
                Class_Para_Tree<2> pablo;
                pablo.setBalance((uint32_t)0,false);

                /**<Refine globally four level and write the para_tree.*/
                for (iter=0; iter<ga; iter++){
                        pablo.adaptGlobalRefine();
                }

#if NOMPI==0
                /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
                pablo.loadBalance();
#endif

                /**<Define a center point and a radius.*/
                double xc, yc;
                xc = yc = 0.75;
                double radius = 0.25;

                /**<Define vectors of data.*/
                uint32_t nocts = pablo.getNumOctants();
                uint32_t nghosts = pablo.getNumGhosts();
                vector<double> weight;
                vector<double> weight2;
                vector<double> weight_ghosts;


                /**<Set marker inside circle*/
                for (int i=0; i<nocts; i++){
                        vector<vector<double> > nodes = pablo.getNodes(i);
                        for (int j=0; j<global2D.nnodes; j++){
                                double x = nodes[j][0];
                                double y = nodes[j][1];
                                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                                        pablo.setMarker(i,1);
                                }
                        }
                }

                pablo.adapt();
                if(pablo.rank==3 || pablo.rank==4)
                        pablo.setMarker((uint32_t)0,1);
                pablo.adapt();
                if(pablo.rank==3 || pablo.rank==4)
                        pablo.setMarker((uint32_t)0,1);
                pablo.adapt();
                if(pablo.rank==3)
                        pablo.setMarker((uint32_t)0,1);
                pablo.adapt();
                if(pablo.rank==3)
                        pablo.setMarker((uint32_t)0,1);
                pablo.adapt();
                nocts = pablo.getNumOctants();
                cout << pablo.rank << " nocts " << nocts << endl;
                weight.resize(nocts, 1.0);
                weight2.resize(nocts, 1.0);
                for (int i=0; i<nocts; i++){
                        vector<vector<double> > nodes = pablo.getNodes(i);
                        for (int j=0; j<global2D.nnodes; j++){
                                double x = nodes[j][0];
                                double y = nodes[j][1];
                                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                                    weight[i] = 2;
                                    weight2[i] = 2;
                                }
                        }
                }
                double sum = 0.0;
                for (int i=0; i<nocts; i++){
                	sum += weight[i];
                }
                cout << pablo.rank << " nocts " << nocts << endl;
                cout << pablo.rank << " sum weight " << sum << endl;


#if NOMPI==0
                User_Data_LB<vector<double> > data_lb(weight, weight_ghosts);
				pablo.loadBalance(data_lb, &weight2);
				weight2 = weight;
                nocts = pablo.getNumOctants();
                for (int i=0; i<nocts; i++){
                        vector<vector<double> > nodes = pablo.getNodes(i);
                        for (int j=0; j<global2D.nnodes; j++){
                                double x = nodes[j][0];
                                double y = nodes[j][1];
                                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                                    weight[i] = 5;
                                    weight2[i] = 5;
                                }
                        }
                }
				pablo.loadBalance(data_lb, &weight2);
//				pablo.loadBalance(data_lb);
#endif

                nocts = pablo.getNumOctants();
                sum = 0.0;
                for (int i=0; i<nocts; i++){
                	sum += weight[i];
                }
                cout << pablo.rank << " nocts " << nocts << endl;
                cout << pablo.rank << " size weight " << weight.size() << endl;
                cout << pablo.rank << " sum weight " << sum << endl;

                /**<Update the connectivity and write the para_tree.*/
                iter = 0;
                pablo.updateConnectivity();
                pablo.writeTest("PabloLB_iter"+to_string(static_cast<unsigned long long>(iter)), weight);

#if NOMPI==0
        }

        MPI::Finalize();
#endif
}




