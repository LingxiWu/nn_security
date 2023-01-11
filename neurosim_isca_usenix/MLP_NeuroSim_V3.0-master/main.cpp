/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
*   
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer. 
*   
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
*   
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "Cell.h"
#include "Array.h"
#include "formula.h"
#include "NeuroSim.h"
#include "Param.h"
#include "IO.h"
#include "Train.h"
#include "Test.h"
#include "Mapping.h"
#include "Definition.h"
#include "omp.h"
 
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

	//printf("unitLengthWireResistance: %.20f \n", unitLengthWireResistance);

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

double ReadCell_lingxi(int x, int y, double conductance){ // x: col y: row

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


int baseline_attack(){

	int num_thermo_levels = 256; // ADC resolution steps
	// double COND_LEVEL = num_thermo_levels; // original code
	double COND_LEVEL = 32;

	/* cell conductance range */
	// ideal device
	// double MAX_COND = 5e-6;
	// double MIN_COND = 100e-9; 
	// Ag:a-Si
	// double MAX_COND = 3.8462e-8;
	// double MIN_COND = 3.0769e-9; 
	// FeFET --- done
	// double MAX_COND = 1.788e-6;		
	// double MIN_COND = 3.973e-8;
	// PCM
	// double MAX_COND = 2.123e-4;		
	// double MIN_COND = 1.072e-5;
	// TaOx/HfOx --- done
	double MAX_COND = 1e-5;		
	double MIN_COND = 1e-6;
	// 3T1C (PCM)
	// double MAX_COND = 2e-5;
	// double MIN_COND = 4e-6;


	/* ADC range */
	// Parallel architecture with resistive crosspoint array for dictionary learning acceleration
	// long double min_cur = 0;
	// long double max_cur = 10e-6;
	// ideal device
	// long double min_cur = 0; 
	// long double max_cur = 8e-6;
	// Ag:a-Si
	long double min_cur = 0; 
	long double max_cur = 7.96e-6;
	// PCM - test
	// long double min_cur = 0;
	// long double max_cur = 0.00010615;


	// array dimensions
	double ROW_SIZE = 400;
	double COL_SIZE = 100;

	/*
		dump trained cond of layer one into a file. row by row
	*/
	ofstream cond_file_1;
	if(remove("trained_cond_1.csv")==0)
 		cout << "previous trained_cond_1.csv removed ..." << endl;
 	cond_file_1.open("trained_cond_1.csv", ios::app);

 	for(int row=0;row<param->nInput;row++){
 		for(int col=0;col<param->nHide;col++){
 			cond_file_1 << setprecision(20) << static_cast<eNVM*>(arrayIH->cell[col][row])->conductance << ",";
 		}
 		cond_file_1 << endl;
 	}
 	cond_file_1.close();


	/*
		generate current for each cell, save it to a file
			   col_0	col_1    col_2	col_3    ...
		row_0
		row_1
		row_2
		...
	*/
	string cur_file_name = "cur_sweep.csv";
	if(remove(cur_file_name.c_str())==0)
		cout << "previous " << cur_file_name << " removed ..." << endl;

	// output Isum of each column to a file
	ofstream cur_file_out;
	cur_file_out.open(cur_file_name, ios::app);

	for(int row=0;row<ROW_SIZE;row++){
		for(int col=0;col<COL_SIZE;col++){
			cur_file_out << setprecision(10) << ReadCell_lingxi(col, row, static_cast<eNVM*>(arrayIH->cell[col][row])->conductance) << ",";
		}
		cur_file_out << endl;
	}

	cur_file_out.close();

	/*
		scan individual current, and determine: 
		(1). if they can be stolen by single row activation
		(2). if so, which one of the ADC level they fell into
	*/
	vector<double> thermo_cur_levels(num_thermo_levels);
	
	long double cur_increment = (max_cur - min_cur) / num_thermo_levels;

	for(int i=1;i<=num_thermo_levels;i++){
		thermo_cur_levels[i-1] = cur_increment * i + min_cur;
		// printf("thermo_cur_levels[%d]: %.4e\n",i-1, thermo_cur_levels[i-1]);
	}

	// input cur sweep file
	ifstream cur_file_in(cur_file_name);
	string cur_line;
	vector<vector<double>> cur_vector(ROW_SIZE, vector<double>(COL_SIZE)); // row-major store current at each cell location
	
	if(cur_file_in.is_open()){
		string delimiter = ",";
		int row = 0;

		while(getline(cur_file_in,cur_line)){ // cond_1.csv stores cond row-major way. each line is a row of curs
			// fill up cur_vector row by row
			vector<string> tokens{}; // row of curs
			size_t pos = 0;
			while((pos = cur_line.find(delimiter)) != string::npos){
				tokens.push_back(cur_line.substr(0, pos));
				cur_line.erase(0, pos + delimiter.length());
			}
			for(int col=0;col<COL_SIZE;col++){
				cur_vector[row][col] = stod(tokens[col]);
			}
			row += 1;
		}
	}

	cur_file_in.close();

	/*
		thermo level histogram 
	*/
	cout << "total num weights: " << (int)(ROW_SIZE * COL_SIZE) << endl;
	vector<int> above_vec(num_thermo_levels);

	for(int row=0;row<ROW_SIZE;row++){
		for(int col=0;col<COL_SIZE;col++){
			for(int i=(num_thermo_levels-1);i>=0;i--){
				if(cur_vector[row][col] >= thermo_cur_levels[i]){
					above_vec[i]+=1;
					break;
				}
			}

		}
	}
	// histogram printout
	int total_above = 0;
	for(int i=0;i<num_thermo_levels;i++){
		total_above += above_vec[i];
		cout << "above level " << i << "(" << thermo_cur_levels[i] << " uA)" << ": " << above_vec[i] << endl;
	}

	cout << "total detectable using single ADC side channel: " << total_above << endl;

	/*
		baseline way, recover only recoverable. if below min thermo level, replace with smallest cond
		note: ADC range needs to be calibrated (thermo_cur_levels)
		step-1: pre-calculate the current steps generated by conductances: cond_cur_vec[COND_LEVEL]: [..., 0.5(v) / (1/cond_i), ...]
		step-1: for each cell, identify the corresponding thermo current cur_thermo
		step-2: from cond_cur_vec, find the a current value that is closest to cur_thermo, get the index, and retrieve the conductance accordingly
	*/

	// step-1 fill cond_cur_vector
	vector<double> cond_vector(COND_LEVEL);

	// build a list of all possible conductances
	double cond_increment = (MAX_COND - MIN_COND) / (COND_LEVEL-1);
	for(int i=0;i<COND_LEVEL;i++){
		// printf("conductances[%d]: %.4e\n", i, cond_increment * i + MIN_COND);
		// printf("%.4e\n",cond_increment * i + MIN_COND);
		cond_vector[i] = cond_increment * i + MIN_COND;
	}
	cout << " "  << endl;
	// populate cond_cur_vector
	double readVoltage = 0.5;
	vector<double> cond_cur_vec(COND_LEVEL);
	for(int i=0;i<COND_LEVEL;i++){
		cond_cur_vec[i] = readVoltage / (1/cond_vector[i]);
		printf("cond_cur_vec[%d]: %.4e\n",i,cond_cur_vec[i]);
	}

	vector<vector<double>> stolen_cond_vec(ROW_SIZE, vector<double>(COL_SIZE)); // row-major
	cout << "\nBaseline attack ... " << endl;


	/////////////////// attack procedure
	double clone_accuracy = 0;
	for(int iteration = -3; iteration<3; iteration++){
		// step-2: find closest thermo_cur at each cell location, "cur_thermo" is thermo_cur_levels[i]
		for(int row=0;row<ROW_SIZE;row++){
			for(int col=0;col<COL_SIZE;col++){
			
				bool found_thermo_cur = false;
				for(int i=(num_thermo_levels-1);i>=0;i--){ // find the min ADC res step
					if(cur_vector[row][col] >= thermo_cur_levels[i]){ // check if detectable by ADC at that location
					
						if(row==399 && col > 95){ // DEBUG
							// printf("cur_vector[%d][%d]: %.4e uA, thermo_cur_levels[%d]: %.4e uA\n",row, col, cur_vector[row][col], i, thermo_cur_levels[i]);
						}
					
						// step-3: compare thermo_cur_levels[i] with cond_cur_vec, find the closest one, record the index in cond_cur_vec
						found_thermo_cur = true;

						for(int j=(COND_LEVEL-1);j>=0;j--){ // find a suitable cond_level
							if(thermo_cur_levels[i] >= cond_cur_vec[j]){
								// stolen_cond_vec[row][col] = cond_vector[j];
								// ad-hoc: consider calibration 
								stolen_cond_vec[row][col] = cond_vector[j+iteration];
								if(row==399 && col > 95){ // DEBUG
									// printf("found cond_cur: %.4e \n", i, j, cond_vector[j]);
								}
								break;
							}
						}
						break;
					}
				}
				// substitute with min cond
				if(found_thermo_cur == false){
					stolen_cond_vec[row][col] = cond_vector[0];
				}
			}
		}

		/*
			Validate accuracy of the stolen conductances
		*/
		for(int row=0;row<param->nInput;row++){
			for(int col=0;col<param->nHide;col++){
				static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = stolen_cond_vec[row][col];
			}
		}
		
		Validate();
		if(clone_accuracy < (double)correct/param->numMnistTestImages*100){
			clone_accuracy = (double)correct/param->numMnistTestImages*100;
		}
		printf("baseline side-channel attack, iteration %d Accuracy: %.2f%\n", iteration, (double)correct/param->numMnistTestImages*100);

	}
	printf("final clone accuracy: %.2f%\n",clone_accuracy);


	/*
		store baseline stolen conds row by row
	*/ 
	string baseline_stolen_file_name = "baseline_stolen_cond_1.csv";
	if(remove(baseline_stolen_file_name.c_str())==0)
		cout << "previous " << baseline_stolen_file_name << " removed ..." << endl;
	
	ofstream baseline_stolen_file;
	baseline_stolen_file.open(baseline_stolen_file_name, ios::app);

	for(int row=0;row<ROW_SIZE;row++){
		for(int col=0;col<COL_SIZE;col++){
			baseline_stolen_file << setprecision(20) << stolen_cond_vec[row][col] << ",";
		}
		baseline_stolen_file << endl;
	}
	
	baseline_stolen_file.close();

}



int main() {
	gen.seed(0);
	
	/* Load in MNIST data */
	ReadTrainingDataFromFile("patch60000_train.txt", "label60000_train.txt");
	ReadTestingDataFromFile("patch10000_test.txt", "label10000_test.txt");

	/* Initialization of synaptic array from input to hidden layer */
	// arrayIH->Initialization<IdealDevice>();
	arrayIH->Initialization<RealDevice>(); 
	// arrayIH->Initialization<MeasuredDevice>();
	//arrayIH->Initialization<SRAM>(param->numWeightBit);
	//arrayIH->Initialization<DigitalNVM>(param->numWeightBit,true);
	// arrayIH->Initialization<HybridCell>(); // the 3T1C+2PCM cell
	//arrayIH->Initialization<_2T1F>();

	
	/* Initialization of synaptic array from hidden to output layer */
	// arrayHO->Initialization<IdealDevice>();
	arrayHO->Initialization<RealDevice>();
	// arrayHO->Initialization<MeasuredDevice>();
	//arrayHO->Initialization<SRAM>(param->numWeightBit);
	//arrayHO->Initialization<DigitalNVM>(param->numWeightBit,true);
	// arrayHO->Initialization<HybridCell>(); // the 3T1C+2PCM cell
	//arrayHO->Initialization<_2T1F>();

    omp_set_num_threads(16);
	/* Initialization of NeuroSim synaptic cores */
	param->relaxArrayCellWidth = 0;
	NeuroSimSubArrayInitialize(subArrayIH, arrayIH, inputParameterIH, techIH, cellIH);
	param->relaxArrayCellWidth = 1;
	NeuroSimSubArrayInitialize(subArrayHO, arrayHO, inputParameterHO, techHO, cellHO);
	/* Calculate synaptic core area */
	NeuroSimSubArrayArea(subArrayIH);
	NeuroSimSubArrayArea(subArrayHO);
	
	/* Calculate synaptic core standby leakage power */
	NeuroSimSubArrayLeakagePower(subArrayIH);
	NeuroSimSubArrayLeakagePower(subArrayHO);
	
	/* Initialize the neuron peripheries */
	NeuroSimNeuronInitialize(subArrayIH, inputParameterIH, techIH, cellIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
	NeuroSimNeuronInitialize(subArrayHO, inputParameterHO, techHO, cellHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);

	/* Calculate the area and standby leakage power of neuron peripheries below subArrayIH */
	double heightNeuronIH, widthNeuronIH;
	NeuroSimNeuronArea(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH, &heightNeuronIH, &widthNeuronIH);
	double leakageNeuronIH = NeuroSimNeuronLeakagePower(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
	/* Calculate the area and standby leakage power of neuron peripheries below subArrayHO */
	double heightNeuronHO, widthNeuronHO;
	NeuroSimNeuronArea(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO, &heightNeuronHO, &widthNeuronHO);
	double leakageNeuronHO = NeuroSimNeuronLeakagePower(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);
	
	/* Print the area of synaptic core and neuron peripheries */
	double totalSubArrayArea = subArrayIH->usedArea + subArrayHO->usedArea;
	double totalNeuronAreaIH = adderIH.area + muxIH.area + muxDecoderIH.area + dffIH.area + subtractorIH.area;

	cout << " " << endl;
	printf("1st layer neuron peripheries breakdown --> adderIH.area: %.4e, muxIH.area: %.4e, mucDecoderIH.area: %.4e, dffIH.area: %.4e, subtractorIH.area: %.4e\n", adderIH.area, muxIH.area, muxDecoderIH.area, dffIH.area, subtractorIH.area);
	cout << " " << endl;

	cout << " " << endl;
	printf("2nd layer neuron peripheries breakdown --> adderIH.area: %.4e, muxIH.area: %.4e, mucDecoderIH.area: %.4e, dffIH.area: %.4e, subtractorIH.area: %.4e\n", adderHO.area, muxHO.area, muxDecoderHO.area, dffHO.area, subtractorHO.area);
	cout << " " << endl;

	subArrayIH->CalculateLatency(NULL);
	subArrayHO->CalculateLatency(NULL);

	printf("subarrayIH readEnergy: %.4e, latency: %.4e\n",NeuroSimSubArrayReadEnergy(subArrayIH),subArrayIH->readLatency);
	printf("subArrayHO readEnergy: %.4e, latency: %.4e\n",NeuroSimSubArrayReadEnergy(subArrayHO),subArrayHO->readLatency);

	// latency report
	subArrayIH->wlDecoder.CalculateLatency(1e20, subArrayIH->wlDecoderOutput.capNorInput, NULL, 1, 1);	// Don't care write
	printf("subarrayIH latency: wl: %.4e\n", subArrayIH->wlDecoder.readLatency);
	subArrayHO->wlDecoder.CalculateLatency(1e20, subArrayHO->wlDecoderOutput.capNorInput, NULL, 1, 1);	// Don't care write
	printf("subarrayHO latency: wl: %.4e\n", subArrayHO->wlDecoder.readLatency);

	subArrayIH->wlDecoderOutput.CalculateLatency(subArrayIH->wlDecoder.rampOutput, subArrayIH->capRow2, subArrayIH->resRow, 1, 1);
	printf("subarrayIH latency: wl_OUTPUT: %.4e\n", subArrayIH->wlDecoderOutput.readLatency);
	subArrayHO->wlDecoderOutput.CalculateLatency(subArrayHO->wlDecoder.rampOutput, subArrayHO->capRow2, subArrayHO->resRow, 1, 1);// Don't care write
	printf("subarrayHO latency: wl_OUTPUT: %.4e\n", subArrayHO->wlDecoderOutput.readLatency);

	subArrayIH->blSwitchMatrix.CalculateLatency(1e20, subArrayIH->capRow1, subArrayIH->resRow, subArrayIH->numReadPulse, 1); 
	printf("subArrayIH latency: blSwitchMatrix: %.4e\n", subArrayIH->blSwitchMatrix.readLatency);
	subArrayHO->blSwitchMatrix.CalculateLatency(1e20, subArrayHO->capRow1, subArrayHO->resRow, subArrayHO->numReadPulse, 1); 
	printf("subArrayHO latency: blSwitchMatrix: %.4e\n", subArrayHO->blSwitchMatrix.readLatency);

	subArrayIH->slSwitchMatrix.CalculateLatency(1e20, subArrayIH->capCol, subArrayIH->resCol, 1, 1);
	printf("subArrayIH latency: SlSwitchMatrix: %.4e\n", subArrayIH->slSwitchMatrix.readLatency);
	subArrayHO->slSwitchMatrix.CalculateLatency(1e20, subArrayHO->capCol, subArrayHO->resCol, 1, 1);
	printf("subArrayHO latency: SlSwitchMatrix: %.4e\n", subArrayHO->slSwitchMatrix.readLatency);

	printf("NeuroSimSubArrayReadLatency(subArrayIH): %.4e\n",NeuroSimSubArrayReadLatency(subArrayIH));
	printf("NeuroSimSubArrayReadLatency(subArrayHO): %.4e\n",NeuroSimSubArrayReadLatency(subArrayHO));

	subArrayIH->subtractor.CalculateLatency(1e20, 0, subArrayIH->numReadPulse);
	printf("subarrayIH.subtractor.latency: %.4e\n", subArrayIH->subtractor.readLatency);
	subArrayHO->subtractor.CalculateLatency(1e20, 0, subArrayHO->numReadPulse);
	printf("subarrayHO.subtractor.latency: %.4e\n", subArrayHO->subtractor.readLatency);

	// subArrayIH->mux.CalculateLatency(0, 0, 1);
	// printf("subarrayIH.mux.latency: %.4e\n", subArrayIH->mux.readLatency);
	// subArrayHO->mux.CalculateLatency(0, 0, 1);
	// printf("subarrayHO.mux.latency: %.4e\n", subArrayHO->mux.readLatency);

	adderIH.CalculateLatency(1e20, subArrayIH->mux.capTgDrain, 1);
	subArrayIH->mux.CalculateLatency(subArrayIH->adder.rampOutput, subArrayIH->dff.capTgDrain, 1);
	subArrayIH->muxDecoder.CalculateLatency(1e20, subArrayIH->mux.capTgGateN * subArrayIH->adder.numAdder, 
	subArrayIH->mux.capTgGateP * subArrayIH->adder.numAdder, 1, 1);	// Don't care write
	subArrayIH->subtractor.CalculateLatency(1e20, subArrayIH->mux.capTgDrain, 1);
	dffIH.CalculateLatency(1e20, 1);
	printf("neuron peripheral IH latency --> adderIH: %.4e, mux: %.4e, muxDeco: %.4e, subtractor: %.4e, dffIH: %.4e\n", 
		adderIH.readLatency, subArrayIH->mux.readLatency, subArrayIH->muxDecoder.readLatency, 
		subArrayIH->subtractor.readLatency, dffIH.readLatency);

	adderHO.CalculateLatency(1e20, subArrayHO->mux.capTgDrain, 1);
	subArrayHO->mux.CalculateLatency(subArrayHO->adder.rampOutput, subArrayHO->dff.capTgDrain, 1);
	subArrayHO->muxDecoder.CalculateLatency(1e20, subArrayHO->mux.capTgGateN * subArrayHO->adder.numAdder, 
	subArrayHO->mux.capTgGateP * subArrayHO->adder.numAdder, 1, 1);	// Don't care write
	subArrayHO->subtractor.CalculateLatency(1e20, subArrayHO->mux.capTgDrain, 1);
	dffHO.CalculateLatency(1e20, 1);
	printf("neuron peripheral HO latency --> adder: %.4e, mux: %.4e, muxDeco: %.4e, subtractor: %.4e, dffHO: %.4e\n", 
		adderHO.readLatency, subArrayHO->mux.readLatency, subArrayHO->muxDecoder.readLatency, 
		subArrayHO->subtractor.readLatency, dffHO.readLatency);

	// subArrayIH->wlDecoderOutput.CalculateLatency(subArrayIH->wlDecoder.rampOutput, subArrayIH->capRow2, subArrayIH->resRow, 1, 1);	// Don't care write
	// subArrayIH->blSwitchMatrix.CalculateLatency(1e20, subArrayIH->capRow1, subArrayIH->resRow, subArrayIH->numReadPulse, 1); 

cout << " " << endl;
	double totalNeuronAreaHO = adderHO.area + muxHO.area + muxDecoderHO.area + dffHO.area + subtractorHO.area;
	printf("Total SubArray (synaptic core) area=%.4e m^2\n", totalSubArrayArea);
	printf("Total Neuron (neuron peripheries) area=%.4e m^2\n", totalNeuronAreaIH + totalNeuronAreaHO);
	printf("Total area=%.4e m^2\n", totalSubArrayArea + totalNeuronAreaIH + totalNeuronAreaHO);

	/* Print the standby leakage power of synaptic core and neuron peripheries */
	printf("Leakage power of subArrayIH is : %.4e W\n", subArrayIH->leakage);
	printf("Leakage power of subArrayHO is : %.4e W\n", subArrayHO->leakage);
	printf("Leakage power of NeuronIH is : %.4e W\n", leakageNeuronIH);
	printf("Leakage power of NeuronHO is : %.4e W\n", leakageNeuronHO);
	printf("Total leakage power of subArray is : %.4e W\n", subArrayIH->leakage + subArrayHO->leakage);
	printf("Total leakage power of Neuron is : %.4e W\n", leakageNeuronIH + leakageNeuronHO);
	
	/* Initialize weights and map weights to conductances for hardware implementation */
	WeightInitialize();
	if (param->useHardwareInTraining)
    	WeightToConductance();
	srand(0);	// Pseudorandom number seed
	

	ofstream mywriteoutfile;
	mywriteoutfile.open("output.csv");                                      
////////

	// cout << "test previous trained cond accuracy" << endl;
	// // validate on cond_file_1
	// ifstream cond_1_file("../cond_1.csv");
	// string cur_line;
	// vector<vector<double>> cond_vec(param->nInput, vector<double>(param->nHide)); // row-major

	// // read out baseline stolen conductances to a vector
	// if(cond_1_file.is_open()){
	// 	string delimiter = ",";
	// 	int row = 0;

	// 	while(getline(cond_1_file,cur_line)){ // cond stored row-major way. each line is a row of curs

	// 		vector<string> tokens{}; // row of conductances
	// 		size_t pos = 0;
	// 		while((pos = cur_line.find(delimiter)) != string::npos){
	// 			tokens.push_back(cur_line.substr(0, pos));
	// 			cur_line.erase(0, pos + delimiter.length());
	// 		}
	// 		for(int col=0;col<param->nHide;col++){
	// 			cond_vec[row][col] = stod(tokens[col]);
	// 		}
	// 		row += 1;
	// 	}
	// }

	// // cond_1_file.close();
	// for(int row=0;row<param->nInput;row++){
	// 	for(int col=0;col<param->nHide;col++){
	// 		static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = cond_vec[row][col];
	// 	}
	// }
	// Validate();
	// printf(" Accuracy of prev conds: %.2f%\n", (double)correct/param->numMnistTestImages*100);

////////

	// for (int i=1; i<=param->totalNumEpochs/param->interNumEpochs; i++){
	for (int i=1; i<=5; i++){
		Train(param->numTrainImagesPerEpoch, param->interNumEpochs,param->optimization_type);
		if (!param->useHardwareInTraining && param->useHardwareInTestingFF) { WeightToConductance(); }
		Validate();
        if (HybridCell *temp = dynamic_cast<HybridCell*>(arrayIH->cell[0][0]))
            WeightTransfer();
        else if(_2T1F *temp = dynamic_cast<_2T1F*>(arrayIH->cell[0][0]))
            WeightTransfer_2T1F();
                
		mywriteoutfile << i*param->interNumEpochs << ", " << (double)correct/param->numMnistTestImages*100 << endl;
		
		printf("Accuracy at %d epochs is : %.2f%\n", i*param->interNumEpochs, (double)correct/param->numMnistTestImages*100);
		/* Here the performance metrics of subArray also includes that of neuron peripheries (see Train.cpp and Test.cpp) */
		// printf("\tRead latency=%.4e s\n", subArrayIH->readLatency + subArrayHO->readLatency);
		// printf("\tWrite latency=%.4e s\n", subArrayIH->writeLatency + subArrayHO->writeLatency);
		// printf("\tRead energy=%.4e J\n", arrayIH->readEnergy + subArrayIH->readDynamicEnergy + arrayHO->readEnergy + subArrayHO->readDynamicEnergy);
		// printf("\tWrite energy=%.4e J\n", arrayIH->writeEnergy + subArrayIH->writeDynamicEnergy + arrayHO->writeEnergy + subArrayHO->writeDynamicEnergy);
		// if(HybridCell* temp = dynamic_cast<HybridCell*>(arrayIH->cell[0][0])){
  //           printf("\tTransfer latency=%.4e s\n", subArrayIH->transferLatency + subArrayHO->transferLatency);
  //           printf("\tTransfer latency=%.4e s\n", subArrayIH->transferLatency);	
  //           printf("\tTransfer energy=%.4e J\n", arrayIH->transferEnergy + subArrayIH->transferDynamicEnergy + arrayHO->transferEnergy + subArrayHO->transferDynamicEnergy);
  //       }
  //       else if(_2T1F* temp = dynamic_cast<_2T1F*>(arrayIH->cell[0][0])){
  //           printf("\tTransfer latency=%.4e s\n", subArrayIH->transferLatency);	
  //           printf("\tTransfer energy=%.4e J\n", arrayIH->transferEnergy + subArrayIH->transferDynamicEnergy + arrayHO->transferEnergy + subArrayHO->transferDynamicEnergy);
  //       }
        // printf("\tThe total weight update = %.4e\n", totalWeightUpdate);
        // printf("\tThe total pulse number = %.4e\n", totalNumPulse);
	}
	
	baseline_attack();
	
	//////////////////////////////////////////////////////////////////////////////////////////
	// cout << "Training done. Testing conductance scramble ..." << endl;

	// /* get layer one conductance, for every CHUNK in the col, shuffle cond  */
	// vector<vector<double>> cond_array(param->nHide, vector<double>(param->nInput)); // 100 x 400
	// int CHUNK_SIZE = 9;
	// vector<double> chunk(CHUNK_SIZE);
	// for(int col=0;col<param->nHide;col++){ // 100
	// 	int c = 0;
	// 	for(int row=0;row<param->nInput;row++){ // 400 (row: 0 ~ 399)
	// 		/* gather a chunk of cond along a col, shuffle chunk */
	// 		if(c<CHUNK_SIZE){
	// 			chunk[c] = static_cast<eNVM*>(arrayIH->cell[col][row])->conductance;
	// 		} else if(c==CHUNK_SIZE) { // obtained a chunk, ready to shuffle

	// 			/* one chunk ready to be scrambled */
	// 			random_shuffle(chunk.begin(), chunk.end());

	// 			/* put the chunk into cond_array */
	// 			int row_index = row - CHUNK_SIZE;
	// 			for(int i=0;i<CHUNK_SIZE;i++){
	// 				cond_array[col][row_index+i] = chunk[i];
	// 			}

	// 			// start to gather a new chunk
	// 			c = 0;
	// 			chunk[c] = static_cast<eNVM*>(arrayIH->cell[col][row])->conductance;
	// 		}
	// 		c += 1;
	// 	}
	// 	/* deal with the last chunk 380 ~ 399*/
	// 	random_shuffle(chunk.begin(), chunk.end());
	// 	int row_index = param->nInput - CHUNK_SIZE;
	// 	/* put last chunk into cond_array*/
	// 	for(int i=0;i<CHUNK_SIZE;i++){
	// 		cond_array[col][row_index+i] = chunk[i];
	// 	}
	// }

	// /* remap the conductances*/
	// for(int col=0;col<param->nHide;col++){
	// 	for(int row=0;row<param->nInput;row++){
	// 		static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = cond_array[col][row];
	// 	}
	// }

	// /* re-test accurracy */  
	// Validate();
	// printf("chunk_size: %d, Accuracy: %.2f%\n", CHUNK_SIZE, (double)correct/param->numMnistTestImages*100);
	//////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////
	// print out old conductances into file
	// ofstream og_cond_col1;
	// og_cond_col1.open("og_cond_col1.txt", ios::app);
	// ofstream shuf_cond_col1;
	// shuf_cond_col1.open("shuf_cond_col1.txt", ios::app);
	// for(int i=0;i<400;i++){
	// 	og_cond_col1 << setprecision(10) << static_cast<eNVM*>(arrayIH->cell[0][i])->conductance << endl;
	// }
	// for(int j=0;j<400;j++){
	// 	shuf_cond_col1 << setprecision(10) << cond_array[0][j] << endl;
	// }
	// og_cond_col1.close();
	// shuf_cond_col1.close();
	/////////////////////////////////////////////


	/////////////////////////////////////////////
	/*
		dump trained cond of layer one into a file. row by row
	*/
	// ofstream cond_file_1;
	// if(remove("../cond_1.csv")==0)
 // 		cout << "previous cond_1.csv removed ..." << endl;
 // 	cond_file_1.open("../cond_1.csv", ios::app);

 // 	for(int row=0;row<param->nInput;row++){
 // 		for(int col=0;col<param->nHide;col++){
 // 			cond_file_1 << setprecision(20) << static_cast<eNVM*>(arrayIH->cell[col][row])->conductance << ",";
 // 		}
 // 		cond_file_1 << endl;
 // 	}
 // 	cond_file_1.close();
	/////////////////////////////////////////////

	/////////////////////////////////////////////
	/*  === test ===
		to simulate baseline attack, replace weight below thermo_level[0]: 0.125 uA with min_cond
	*/
	// for(int row=0;row<param->nInput;row++){
	// 	for(int col=0;col<param->nHide;col++){
	// 		if(static_cast<eNVM*>(arrayIH->cell[col][row])->conductance < 1.78e-07){
	// 			static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = 1e-7;
	// 		} else {
	// 			if(row == 0)
	// 				TO-DO
	// 				// printf("org: %.6e, round: %.6e \n", static_cast<eNVM*>(arrayIH->cell[col][row])->conductance, setprecision(3)static_cast<eNVM*>(arrayIH->cell[col][row])->conductance)));
	// 			static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = static_cast<eNVM*>(arrayIH->cell[col][row])->conductance;
	// 		}
	// 	}
	// }
	/////////////////////////////////////////////

 	/////////////////////////////////////////////
 	/*
		test baseline stolen cond accuracy
 	*/
	// ifstream baseline_stolen_cond_file("../baseline_stolen.csv");
	// string cur_line;
	// vector<vector<double>> baseline_cond_vec(param->nInput, vector<double>(param->nHide)); // row-major

	// // read out baseline stolen conductances to a vector
	// if(baseline_stolen_cond_file.is_open()){
	// 	string delimiter = ",";
	// 	int row = 0;

	// 	while(getline(baseline_stolen_cond_file,cur_line)){ // cond stored row-major way. each line is a row of curs

	// 		vector<string> tokens{}; // row of conductances
	// 		size_t pos = 0;
	// 		while((pos = cur_line.find(delimiter)) != string::npos){
	// 			tokens.push_back(cur_line.substr(0, pos));
	// 			cur_line.erase(0, pos + delimiter.length());
	// 		}
	// 		for(int col=0;col<param->nHide;col++){
	// 			baseline_cond_vec[row][col] = stod(tokens[col]);
	// 		}
	// 		row += 1;
	// 	}
	// }

	// baseline_stolen_cond_file.close();

	// // remap stolen conductances to 
	// for(int row=0;row<param->nInput;row++){
	// 	for(int col=0;col<param->nHide;col++){
	// 		static_cast<eNVM*>(arrayIH->cell[col][row])->conductance = baseline_cond_vec[row][col];
	// 	}
	// }

	// check accuracy
	// Validate();
	// printf("baseline side-channel attack (test), Accuracy: %.2f%\n", (double)correct/param->numMnistTestImages*100);

 	/////////////////////////////////////////////
	printf("\n");
	return 0;
}


