#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "User_Data_Comm.hpp"
#include "User_Data_LB.hpp"

using namespace std;

// =================================================================================== //

class bubble{
public:
	double c[3];
	double r;


};


int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{
		int iter = 0;

		/**<Instantation of a 3D para_tree object.*/
		Class_Para_Tree<3> pabloBB;

		/**<Set 2:1 balance for the octree.*/
		int idx = 0;
		pabloBB.setBalance(idx,false);

		/**<Refine globally four level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pabloBB.adaptGlobalRefine();
		}

		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pabloBB.loadBalance();

		/**<Define a set of bubbles.*/
		cout << time(NULL) << endl;
		srand(time(NULL));
//		srand(1418143772);


		//8 proc deadlock 1418143772  at iter 51 proc

		int nb = 100;
		vector<bubble> BB;
		vector<bubble> BB0;
		vector<double> DZ;
		vector<double> OM;
		vector<double> AA;

		for (int i=0; i<nb; i++){
			double randc[3];
			randc[0] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[1] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[2] = (double) (rand()) /  RAND_MAX - 0.5;
			double randr = (0.05 * (double) (rand()) / RAND_MAX + 0.04);
			//double randr = 0.005 * (double) (rand()) / RAND_MAX + 0.002;
			double dz = 0.005 + 0.05 * (double) (rand()) / RAND_MAX;
			double omega = 0.5 * (double) (rand()) / RAND_MAX;
			double aa = 0.15 * (double) (rand()) / RAND_MAX;
			bubble bb;
			bb.c[0] = randc[0];
			bb.c[1] = randc[1];
			bb.c[2] = randc[2];
			bb.r = randr;
			BB.push_back(bb);
			BB0.push_back(bb);
			DZ.push_back(dz);
			OM.push_back(omega);
			AA.push_back(aa);
		}
		double t0 = 0;
		double t = t0;
		double Dt = 0.5;
//		double Dy = 0.01;
//		double omega = 2.0*3.14/0.001;

		/**<Define vectors of data.*/
		uint32_t nocts = pabloBB.getNumOctants();
		uint32_t nghosts = pabloBB.getNumGhosts();
		vector<double> oct_data(nocts, 99), ghost_data(nghosts, 0.0);

		/**<Adapt itend times with data injection on new octants.*/
		int itstart = 1;
		int itend = 200;


		int nrefperiter = 3;

		for (iter=itstart; iter<itend; iter++){
			if(pabloBB.rank==0) cout << "iter " << iter << endl;
			t += Dt;

			for (int i=0; i<nb; i++){
				BB[i].c[0] = BB0[i].c[0];// + AA[i]*cos(OM[i]*t);
				BB[i].c[1] = BB0[i].c[1];// + AA[i]*cos(OM[i]*t);
				BB[i].c[2] = BB[i].c[2]+ Dt*DZ[i];
			}

			for (int iref=0; iref<nrefperiter; iref++){

				for (int i=0; i<nocts; i++){
					bool inside = false;
					vector<vector<double> > nodes = pabloBB.getNodes(i);
					vector<double> center = pabloBB.getCenter(i);
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double zc = BB[ib].c[2];
						double radius = BB[ib].r;
						//oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
						oct_data[i] = 99;
						for (int j=0; j<global3D.nnodes; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							double z = nodes[j][2];
							if ((pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) <= 1.15*pow(radius,2.0) &&
									pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) >= 0.85*pow(radius,2.0)) ||
									(pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0) <= 1.15*pow(radius,2.0) &&
											pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0) >= 0.85*pow(radius,2.0))){
								if (pabloBB.getLevel(i) < 8){
									/**<Set to refine inside the sphere.*/
									pabloBB.setMarker(i,1);
								}
								else{
									pabloBB.setMarker(i,0);
								}
								inside = true;
							}
						}
						ib++;
					}
					if (pabloBB.getLevel(i) > 5 && !inside){
						/**<Set to coarse if the octant has a level higher than 5.*/
						pabloBB.setMarker(i,-1);
					}
				}

				/**<Adapt the octree.*/
				//vector<uint32_t> mapidx;
				//bool adapt = pabloBB.adapt(mapidx);
				bool adapt = pabloBB.adapt();

//				/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
//				pabloBB.loadBalance();

				nocts = pabloBB.getNumOctants();
				nghosts = pabloBB.getNumGhosts();
				oct_data.resize(nocts, 99);

			}

			/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
			pabloBB.loadBalance();
			nocts = pabloBB.getNumOctants();
			nghosts = pabloBB.getNumGhosts();
			oct_data.resize(nocts, 99);

			/**<Assign to the new octant the data after an adaption.*/
//			for (int i=0; i<nocts; i++){
//				vector<vector<double> > nodes = pabloBB.getNodes(i);
//				vector<double> center = pabloBB.getCenter(i);
//				oct_data[i] = 99;
//				int ib = 0;
//				while (ib<nb){
//					double xc = BB[ib].c[0];
//					double yc = BB[ib].c[1];
//					double zc = BB[ib].c[2];
//					double radius = BB[ib].r;
//					if (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0) <= pow(radius,2.0)){
//						if (oct_data[i] == 99) oct_data[i]= -(pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0));
//
//					}
//				}
//				ib++;
//			}


			/**<Update the connectivity and write the para_tree.*/
			pabloBB.updateConnectivity();
			pabloBB.writeTest("PabloBubble_iter"+to_string(iter), oct_data);
			//pabloBB.write("PabloBubble_iter"+to_string(iter));
		}
	}

	MPI::Finalize();

}

