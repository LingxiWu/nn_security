#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

// conductances should be dumped from main.cpp not test.cpp. 

double WireResistance(int x, int y=0){ // x col, y row

	double wireResistanceRow = 0;
	double wireResistanceCol = 0; // array.h
	double resistanceAccess = 15000.0;
	double arrayRowSize = 400.0;

	double Rho = 0.0000000273;
	double wireWidth = 100.0;
	double AR = 2.30;
	double unitLengthWireResistance =  Rho / ( wireWidth*1e-9 * wireWidth*1e-9 * AR );
	
	//double wireLength = wireWidth * 1e-9 * 2;	// 2F
	//wireResistanceRow = unitLengthWireResistance * wireLength;
	//wireResistanceCol = unitLengthWireResistance * wireLength;

	printf("unitLengthWireResistance: %.20f \n", unitLengthWireResistance);

	double numRow = 400;
	double numCol = 100;
	double lengthRow = 0.0000128000;
	double lengthCol = 0.0000512000;
	

	wireResistanceRow = lengthRow / numCol * unitLengthWireResistance;
	wireResistanceCol = lengthCol / numRow * unitLengthWireResistance;

	//printf("wireresistancerow: %.20f, wireResistanceCol: %.20f \n", wireResistanceRow, wireResistanceCol);
	double totalWireResistance = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol + resistanceAccess;

	return totalWireResistance;
}

double ReadCell(int x, int y, double conductance){ // x: col y: row

	double totalWireResistance = WireResistance(x, y);
	//printf("totalWireResistance: %.20f \n", totalWireResistance);
	double cellCurrent = 0;
	double readVoltage = 0.5;
	char readNoise = 0;

	if(readNoise){

	} else {
		cellCurrent = readVoltage / (1.0/conductance + totalWireResistance);
	}

	return cellCurrent;
}

double DeltaEnergy(double Isum){ 

	double delta_e;
	double readVoltage = 0.5;
	double readPulseWidth = 0.000000005;
	double wireCapRow = 0.00000000000000256;
	double activation_e = wireCapRow * readVoltage * readVoltage;
	delta_e = Isum * readVoltage * readPulseWidth;
	delta_e += activation_e;

	return delta_e;
}

int main()
{
	// conductance --> isum --> energy
	// if energy != recorded, try another conductance

	double readVoltage = 0;
	double totalWireResistance = 0;

	double min_cond = 100e-9;
	double max_cond = 5e-6;
	double conductance_levels = 63;

	// create a list of conductances for an ideal device
	vector<double> conductance_vec;
	for (int i=0;i<=conductance_levels;i++){
		double conductance = min_cond + 1/conductance_levels * (max_cond - min_cond) * i;
		conductance_vec.push_back(conductance);
	}

	/* attack starts here */

	int num_synapses = 6;
	double accumulated_energy = 0; // accumulate energy caused by num_synapse weights

	FILE *fp_delta = fopen("MLP_NeuroSim_V3.0/delta_energy.txt", "r");
	int i = 1;
	double delta = 0;
	while (fscanf(fp_delta, "%lf", &delta) != EOF){
		accumulated_energy += delta;
		i += 1;
		if(i>num_synapses*2){
			break;
		}
	}
	printf("accumulated_energy: %.20f \n", accumulated_energy);
	cout << " " << endl;

	// print recorded conductance and isum
	FILE *fp_cond = fopen("MLP_NeuroSim_V3.0/cond_1.txt", "r");
	double cond = 0;
	double ene = 0;
	int j = 0;
	list<double> victim;
	while (fscanf(fp_cond, "%lf", &cond) != EOF){
		victim.push_back(cond);
		
		ene += DeltaEnergy(ReadCell(j,0,cond));
		j += 1;
		if(j>num_synapses-1){
			victim.sort();
			break;
		}
	}
	list<double>::iterator it;
	for (it=victim.begin();it!=victim.end();it++){
		printf("victim conductance: %.20f \n", *it);
	}

cout << " --- " << endl;



/*
	vector<double> stolen_conductance;
	double epsilon = 0.00000000000000000001;
	// first synapse
	for(int c0=0;c0<=conductance_levels;c0++){
		double cond_0 = conductance_vec[c0];
		double isum_0 = ReadCell(0,0,cond_0);
		double energy_0 = DeltaEnergy(isum_0);
		// second synapse
		for(int c1=0;c1<=conductance_levels;c1++){
			double cond_1 = conductance_vec[c1];
			double isum_1 = ReadCell(1,0,cond_1);
			double energy_1 = DeltaEnergy(isum_1);

			double total_energy = energy_0 + energy_1;
			if(fabs(total_energy - accumulated_energy) < epsilon){
				printf("cond_0: %.20f, cond_1: %.20f \n", cond_0, cond_1);
			}

		}
	}
*/

/*
	vector<double> stolen_conductance;
	double epsilon = 0.00000000000000000001;
	// first synapse
	for(int c0=0;c0<=conductance_levels;c0++){
		double cond_0 = conductance_vec[c0];
		double isum_0 = ReadCell(0,0,cond_0);
		double energy_0 = DeltaEnergy(isum_0);
		// second synapse
		for(int c1=0;c1<=conductance_levels;c1++){
			double cond_1 = conductance_vec[c1];
			double isum_1 = ReadCell(1,0,cond_1);
			double energy_1 = DeltaEnergy(isum_1);
			// third synapse
			for(int c2=0;c2<=conductance_levels;c2++){
				double cond_2 = conductance_vec[c2];
				double isum_2 = ReadCell(2,0,cond_2);
				double energy_2 = DeltaEnergy(isum_2);

				double total_energy = energy_0 + energy_1 + energy_2;
				if(fabs(total_energy - accumulated_energy) < epsilon){

					if(find(stolen_conductance.begin(), stolen_conductance.end(), cond_0) == stolen_conductance.end()
						&& find(stolen_conductance.begin(), stolen_conductance.end(), cond_1) == stolen_conductance.end()
						&& find(stolen_conductance.begin(), stolen_conductance.end(), cond_2) == stolen_conductance.end()){

							stolen_conductance.push_back(cond_0);
							stolen_conductance.push_back(cond_1);
							stolen_conductance.push_back(cond_2);
							printf("cond_0: %.20f, cond_1: %.20f, cond_2: %.20f \n", cond_0, cond_1, cond_2);

					}
				}
			}
		}
	}
*/

/*

	//vector<double> stolen_conductance;
	vector<list<double>> stolen_conductance;
	double epsilon = 0.00000000000000000001;
	// first synapse
	for(int c0=0;c0<=conductance_levels;c0++){
		double cond_0 = conductance_vec[c0];
		double isum_0 = ReadCell(0,0,cond_0);
		double energy_0 = DeltaEnergy(isum_0);
		// second synapse
		for(int c1=0;c1<=conductance_levels;c1++){
			double cond_1 = conductance_vec[c1];
			double isum_1 = ReadCell(1,0,cond_1);
			double energy_1 = DeltaEnergy(isum_1);
			// third synapse
			for(int c2=0;c2<=conductance_levels;c2++){
				double cond_2 = conductance_vec[c2];
				double isum_2 = ReadCell(2,0,cond_2);
				double energy_2 = DeltaEnergy(isum_2);
				// fourth synapse
				for(int c3=0;c3<=conductance_levels;c3++){
					double cond_3 = conductance_vec[c3];
					double isum_3 = ReadCell(3,0,cond_3);
					double energy_3 = DeltaEnergy(isum_3);

					double total_energy = energy_0 + energy_1 + energy_2 + energy_3;
					
					// evaluate
					if(fabs(total_energy - accumulated_energy) < epsilon){
						//printf("cond_0: %.20f, cond_1: %.20f, cond_2: %.20f, cond_3: %.20f \n", cond_0, cond_1, cond_2, cond_3);
						list<double> s;
						s.push_back(cond_0);s.push_back(cond_1);s.push_back(cond_2);s.push_back(cond_3);
						s.sort();
						if(stolen_conductance.size() == 0){
							stolen_conductance.push_back(s);
						}
						else{
							char added = 0;
							for(int i=0;i<stolen_conductance.size();i++){
								if(stolen_conductance[i] == s){
									added = 1;
								} 
							}
							if(!added){
								stolen_conductance.push_back(s);
							}
						}
					}
				}
			}
		}
	}


	for(int i=0;i<stolen_conductance.size();i++){
		list<double>::iterator itr;
		for(itr=stolen_conductance[i].begin();itr!=stolen_conductance[i].end();itr++){
			cout << *itr << endl;
		}
		cout << " " << endl;
	}

*/

/*
	vector<list<double>> stolen_conductance;
	double epsilon = 0.00000000000000000001;
	// first synapse
	for(int c0=0;c0<=conductance_levels;c0++){
		double cond_0 = conductance_vec[c0];
		double isum_0 = ReadCell(0,0,cond_0);
		double energy_0 = DeltaEnergy(isum_0);
		// second synapse
		for(int c1=0;c1<=conductance_levels;c1++){
			double cond_1 = conductance_vec[c1];
			double isum_1 = ReadCell(1,0,cond_1);
			double energy_1 = DeltaEnergy(isum_1);
			// third synapse
			for(int c2=0;c2<=conductance_levels;c2++){
				double cond_2 = conductance_vec[c2];
				double isum_2 = ReadCell(2,0,cond_2);
				double energy_2 = DeltaEnergy(isum_2);
				// fourth synapse
				for(int c3=0;c3<=conductance_levels;c3++){
					double cond_3 = conductance_vec[c3];
					double isum_3 = ReadCell(3,0,cond_3);
					double energy_3 = DeltaEnergy(isum_3);
					// fifth synapse
					for(int c4=0;c4<=conductance_levels;c4++){
						double cond_4 = conductance_vec[c4];
						double isum_4 = ReadCell(4,0,cond_4);
						double energy_4 = DeltaEnergy(isum_4);
					
						double total_energy = energy_0 + energy_1 + energy_2 + energy_3 + energy_4;

						if(fabs(total_energy - accumulated_energy) < epsilon){
							//printf("%.20f | %.20f | %.20f | %.20f | %.20f \n", cond_0, cond_1, cond_2, cond_3, cond_4);
							list<double> s;
							s.push_back(cond_0);s.push_back(cond_1);s.push_back(cond_2);s.push_back(cond_3);s.push_back(cond_4);
							s.sort();
							if(stolen_conductance.size() == 0){
								stolen_conductance.push_back(s);
							}
							else{
								char added = 0;
								for(int i=0;i<stolen_conductance.size();i++){
									if(stolen_conductance[i] == s){
										added = 1;
									} 
								}
								if(!added){
									stolen_conductance.push_back(s);
								}
							}
				 		}
					}
				}
			}
		}
	}

*/
	

	vector<list<double>> stolen_conductance;
	double epsilon = 0.00000000000000000001;
	// first synapse
	for(int c0=0;c0<=conductance_levels;c0++){
		double cond_0 = conductance_vec[c0];
		double isum_0 = ReadCell(0,0,cond_0);
		double energy_0 = DeltaEnergy(isum_0);
		// second synapse
		for(int c1=0;c1<=conductance_levels;c1++){
			double cond_1 = conductance_vec[c1];
			double isum_1 = ReadCell(1,0,cond_1);
			double energy_1 = DeltaEnergy(isum_1);
			// third synapse
			for(int c2=0;c2<=conductance_levels;c2++){
				double cond_2 = conductance_vec[c2];
				double isum_2 = ReadCell(2,0,cond_2);
				double energy_2 = DeltaEnergy(isum_2);
				// fourth synapse
				for(int c3=0;c3<=conductance_levels;c3++){
					double cond_3 = conductance_vec[c3];
					double isum_3 = ReadCell(3,0,cond_3);
					double energy_3 = DeltaEnergy(isum_3);
					// fifth synapse
					for(int c4=0;c4<=conductance_levels;c4++){
						double cond_4 = conductance_vec[c4];
						double isum_4 = ReadCell(4,0,cond_4);
						double energy_4 = DeltaEnergy(isum_4);
						// sixth synapse
						for(int c5=0;c5<=conductance_levels;c5++){
							double cond_5 = conductance_vec[c5];
							double isum_5 = ReadCell(5,0,cond_5);
							double energy_5 = DeltaEnergy(isum_5);
					
							double total_energy = energy_0 + energy_1 + energy_2 + energy_3 + energy_4 + energy_5;

							if(fabs(total_energy - accumulated_energy) < epsilon){
								//printf("%.20f | %.20f | %.20f | %.20f | %.20f \n", cond_0, cond_1, cond_2, cond_3, cond_4);
								list<double> s;
								s.push_back(cond_0);s.push_back(cond_1);s.push_back(cond_2);s.push_back(cond_3);s.push_back(cond_4);s.push_back(cond_5);
								s.sort();
								if(stolen_conductance.size() == 0){
									stolen_conductance.push_back(s);
								}
								else{
									char added = 0;
									for(int i=0;i<stolen_conductance.size();i++){
										if(stolen_conductance[i] == s){
											added = 1;
										} 
									}
									if(!added){
										stolen_conductance.push_back(s);
									}
								}
					 		}
					 	}
					}
				}
			}
		}
	}


	cout << "found: " << to_string(stolen_conductance.size()) << " candidates" << endl;
	/*
	for(int i=0;i<stolen_conductance.size();i++){
		list<double>::iterator itr;
		for(itr=stolen_conductance[i].begin();itr!=stolen_conductance[i].end();itr++){
			cout << *itr << endl;
		}
		cout << " " << endl;
	}
	*/

	return 0;

}
