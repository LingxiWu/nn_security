#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <list>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

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

void current_gen(vector<int> input_bits, vector<vector<double>> conductances, double ROW_SIZE, double COL_SIZE, string file_name){
	if(remove(file_name.c_str())==0)
		cout << "previous " << file_name << " removed ..." << endl;
	
	// output Isum of each column to a file
	ofstream isum_file;
	isum_file.open(file_name, ios::app);
	
	// output current from each column to the file
	for(int col=0;col<COL_SIZE;col++){ // col 
		double Isum = 0;
		for(int row=0;row<ROW_SIZE;row++){ // row
			if(input_bits[row] == 1){
				Isum += ReadCell(col,row,conductances[col][row]);
			}
		}
		isum_file << setprecision(20) << Isum << endl; 
	}

	isum_file.close();
	cout << "New " << file_name << " generated ..." << endl; 
}



int main()
{
	// 
	double READ_VOLTAGE = 0;

	// cell properties

	// ideal device
	double MIN_COND = 100e-9; 
	double MAX_COND = 5e-6;
	double COND_LEVEL = 63;

	// real device
	// double MIN_COND = 3.8462e-8;
	// double MAX_COND = 3.0769e-9;
	// double COND_LEVEL = ;
	
	double cond_click = (MAX_COND - MIN_COND)/COND_LEVEL;

	// array dimensions
	double ROW_SIZE = 400;
	double COL_SIZE = 100;

	// memory array of conductances
	vector<vector<double> > 
	conductances(COL_SIZE, vector<double>(ROW_SIZE));

	// input vector indicates which rows are open
	vector<int> input_bits(ROW_SIZE);
	for(int i=0;i<ROW_SIZE;i++){
		input_bits[i] = 0; // default inactive
	}

	/*
		assign random conductances
	*/
	/////////////////////////////////////////
	// srand(time(0));
	// for(int col=0; col<COL_SIZE; col++){ // 100
	// 	for(int row=0; row<ROW_SIZE; row++){ // 400
	// 		conductances[col][row] = MIN_COND + (rand() % (int)COND_LEVEL) * cond_click;
	// 	}
	// }
	/////////////////////////////////////////

	/*
		load conductances with file IO
	*/
	/////////////////////////////////////////
	string line;
	ifstream cond_file_1("cond_1.csv");
	int row = 0;
	int col = 0;
	int total = 0;
	if(cond_file_1.is_open()){

		while(getline(cond_file_1,line)){
			// cout << col << ", " << row << endl;
			conductances[col][row] = stod(line);
			row += 1;
			if(row > 399){
				row = 0;
				col += 1;
			}
		}
	}
	cond_file_1.close();
	/////////////////////////////////////////

	// generate current
	for(int i=20; i<= 400; i+=20){ 
		
		string file_name = "isum_";
		file_name += to_string(i);
		file_name += ".csv";

		for(int j=0;j<i;j++){
			input_bits[j] = 1;
		}

		current_gen(input_bits, conductances, ROW_SIZE, COL_SIZE, file_name);

	}

	// // first 200 rows are open
	// for(int i=0;i<200;i++){
	// 	input_bits[i] = 1;
	// }
	// current_gen(input_bits, conductances, ROW_SIZE, COL_SIZE, "isum_0_200.csv");

	// for(int i=200;i<400;i+=20){

	// }

	
	return 0;

}
