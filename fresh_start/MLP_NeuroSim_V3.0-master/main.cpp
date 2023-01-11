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
#include <stdio.h>
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

using namespace std;

int main() {
	gen.seed(0);
	
	/* Load in MNIST data */
	ReadTrainingDataFromFile("patch60000_train.txt", "label60000_train.txt");
	ReadTestingDataFromFile("patch10000_test.txt", "label10000_test.txt");

	/* Initialization of synaptic array from input to hidden layer */
	arrayIH->Initialization<IdealDevice>();
	// arrayIH->Initialization<RealDevice>();
	//arrayIH->Initialization<MeasuredDevice>();
	//arrayIH->Initialization<SRAM>(param->numWeightBit);
	// arrayIH->Initialization<DigitalNVM>(param->numWeightBit,true); // true: consider refColumn

	
	/* Initialization of synaptic array from hidden to output layer */
	arrayHO->Initialization<IdealDevice>();
	// arrayHO->Initialization<RealDevice>();
	//arrayHO->Initialization<MeasuredDevice>();
	//arrayHO->Initialization<SRAM>(param->numWeightBit);
	// arrayHO->Initialization<DigitalNVM>(param->numWeightBit,true);


	/* Initialization of NeuroSim synaptic cores */
	// synaptic cores --> subArray, which contains cells*** and peripherals
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

	// int numAdder = (int)ceil(((double)subArrayIH->numCol / subArrayIH->numCellPerSynapse) / subArrayIH->numColMuxed);
	// printf("numCol: %d, numColMuxed: %d, numCellPerSynapse: %d, numAdder: %d\n", subArrayIH->numCol,subArrayIH->numColMuxed,subArrayIH->numCellPerSynapse,numAdder);
	
	/* Initialize weights and map weights to conductances for hardware implementation */
	WeightInitialize();
	param->useHardwareInTraining = false;
	if (param->useHardwareInTraining) { WeightToConductance(); }

	// srand(0);	// Pseudorandom number seed
	
	ofstream mywriteoutfile;
	mywriteoutfile.open("my_log.csv");                                                                                                            
	
	if(remove("neuron_col.txt")==0)
		cout << "previous neuron_col.txt removed ..." << endl;

	// for (int i=1; i<=param->totalNumEpochs/param->interNumEpochs; i++) {
	for (int i=1; i<=1; i++) {
        //cout << "Training Epoch : " << i << endl;
		Train(param->numTrainImagesPerEpoch, param->interNumEpochs,param->optimization_type);
		if (!param->useHardwareInTraining && param->useHardwareInTestingFF) { WeightToConductance(); }
		// wipe out the zero file before each validation

		Validate();
		mywriteoutfile << i*param->interNumEpochs << ", " << (double)correct/param->numMnistTestImages*100 << endl;
		
		printf("Accuracy at %d epochs is : %.2f%\n", i*param->interNumEpochs, (double)correct/param->numMnistTestImages*100);
		// printf("\tRead latency=%.4e s\n", subArrayIH->readLatency + subArrayHO->readLatency);
		// printf("\tWrite latency=%.4e s\n", subArrayIH->writeLatency + subArrayHO->writeLatency);
		// printf("\tRead energy=%.4e J\n", arrayIH->readEnergy + subArrayIH->readDynamicEnergy + arrayHO->readEnergy + subArrayHO->readDynamicEnergy);
		// printf("\tWrite energy=%.4e J\n", arrayIH->writeEnergy + subArrayIH->writeDynamicEnergy + arrayHO->writeEnergy + subArrayHO->writeDynamicEnergy);
	}

	cout << "training done ... " << endl;

	// save weights, clicks, and binary weights to files
	if(remove("weight_1.txt")==0)
        cout << "previous weight_1.txt removed ..." << endl;
    if(remove("clicks_1.txt")==0)
        cout << "previous clicks_1.txt removed ..." << endl;
    if(remove("bits_1.txt")==0)
        cout << "previous bits_1.txt removed ..." << endl;

	ofstream weight_file_1;
	weight_file_1.open("weight_1.txt", ios::app);
	// layer one
	for(int j=0;j<param->nInput;j++){ // 400
		for(int k=0;k<param->nHide;k++){ // 100
			weight_file_1 << weight1[k][j] << endl;
		}
	}

	// convert weight to clicks (0 ~ 64)
	ofstream clicks_file_1;
	clicks_file_1.open("clicks_1.txt", ios::app);	
	int numCellPerSynapse = 6;	
	int numLevel = pow(2, numCellPerSynapse);

	for(int j=0;j<param->nInput;j++){ // 400
		for(int k=0;k<param->nHide;k++){ // 100
			int targetWeightDigits = (int)((weight1[k][j] + 1)/2 * (numLevel-1)); 
			clicks_file_1 << targetWeightDigits << endl;
		}
	}

	// convert clicks to binary string
	ofstream bits_file_1;
	bits_file_1.open("bits_1.txt", ios::app);
	for(int j=0;j<param->nInput;j++){ // 400
		for(int k=0;k<param->nHide;k++){ // 100
			int targetWeightDigits = (int)((weight1[k][j] + 1)/2 * (numLevel-1)); 
			bits_file_1 << targetWeightDigits << ",";
			for (int n=0; n<numCellPerSynapse; n++){
				int bit = ((targetWeightDigits >> n) & 1);
				bits_file_1 << bit;
			}
			bits_file_1 << endl;
		}
	}


	// // save bits to file
	// int bits[600][400]; // bits[col][row]
	// for (int x=0; x<param->nHide; x++) { // hidden neuron numbers 100

	// 	for (int y=0; y<param->nInput; y++) { // 400 rows
	// 		int numCellPerSynapse = 6;
	// 		for (int n=0; n<numCellPerSynapse; n++) { 
	// 			// for each neuron, iterate through each bit
	// 			int colIndex = (x+1) * numCellPerSynapse - (n+1);
	// 			double readVoltage = 0.5;
	// 			double totalWireResistance;
				
	// 			double wireResistanceRow = 0;
	// 			double wireResistanceCol = 0; // array.h
	// 			double resistanceAccess = 15000.0;
	// 			double arrayRowSize = 400.0;

	// 			double Rho = 0.0000000273;
	// 			double wireWidth = 100.0;
	// 			double AR = 2.30;
	// 			double unitLengthWireResistance =  Rho / ( wireWidth*1e-9 * wireWidth*1e-9 * AR );

	// 			double numRow = 400;
	// 			double numCol = 100;
	// 			double lengthRow = 0.0000128000;
	// 			double lengthCol = 0.0000512000;
				
	// 			wireResistanceRow = lengthRow / numCol * unitLengthWireResistance;
	// 			wireResistanceCol = lengthCol / numRow * unitLengthWireResistance;

	// 			totalWireResistance = (x + 1) * wireResistanceRow + (arrayRowSize - y) * wireResistanceCol + resistanceAccess;

	// 			double cellCurrent;
	// 			cellCurrent = readVoltage / (1/static_cast<eNVM*>(arrayIH->cell[colIndex][y])->conductance + totalWireResistance);
	// 			int bit;
	// 			if (cellCurrent >= static_cast<DigitalNVM*>(arrayIH->cell[colIndex][y])->refCurrent) 
 //                	bit = 1;
 //                else 
	// 				bit = 0;
	// 			bits[colIndex][y] = bit;
	// 		}
	// 	}
	// }

	// if(remove( "bitsFile.txt" ) == 0 )
 //    	printf("bitsFile.txt removed\n");
 //    if(remove( "bitsRatio.txt" ) == 0 )
 //    	printf("bitsRatio.txt removed\n");

	// ofstream bitsFile;
	// ofstream bitsRatio;
 //  	bitsFile.open ("bitsFile.txt",ios::app);
 //  	bitsRatio.open ("bitsRatio.txt",ios::app);
  	
	// for (int row =0;row<400;row++){
	// 	string row_bits = "";
	// 	int zero_counter = 0;
	// 	int one_counter = 0;
	// 	for (int col=0;col<600;col++){
	// 		row_bits += to_string(bits[col][row]);
	// 		row_bits += ",";
	// 		if(bits[col][row] == 0){
	// 			zero_counter += 1;
	// 		} else {
	// 			one_counter += 1;
	// 		}
	// 	}
	// 	bitsFile << row_bits << endl;
	// 	bitsRatio << "zero counter: " << to_string(zero_counter) << " one counter: " << to_string(one_counter) << endl;
	// }

	// bitsFile.close();
	// bitsRatio.close();




	printf("\n");
	return 0;
}



