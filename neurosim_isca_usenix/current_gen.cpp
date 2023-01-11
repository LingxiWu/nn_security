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

void current_gen_chunk(vector<int> input_bits, vector<vector<double>> conductances, double ROW_SIZE, double COL_SIZE, string file_name){
	
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

	// ideal device
	double MIN_COND = 100e-9; 
	double MAX_COND = 5e-6;
	double COND_LEVEL = 64;

	// real device
	// double MIN_COND = 3.8462e-8;
	// double MAX_COND = 3.0769e-9;
	// double COND_LEVEL = ;
	
	double cond_click = (MAX_COND - MIN_COND)/COND_LEVEL;

	// array dimensions
	double ROW_SIZE = 400;
	double COL_SIZE = 100;

	// array of conductances
	vector<vector<double> > 
	conductances(ROW_SIZE, vector<double>(COL_SIZE));

	// input vector indicates which rows are open
	// vector<int> input_bits(ROW_SIZE);
	// for(int i=0;i<ROW_SIZE;i++){
	// 	input_bits[i] = 0; // default inactive
	// }


	/////////////////////////////////////////
	/*
		assign random conductances
	*/
	// srand(time(0));
	// for(int col=0; col<COL_SIZE; col++){ // 100
	// 	for(int row=0; row<ROW_SIZE; row++){ // 400
	// 		conductances[col][row] = MIN_COND + (rand() % (int)COND_LEVEL) * cond_click;
	// 	}
	// }
	/////////////////////////////////////////


	/////////////////////////////////////////
	/*
		load conductances with file IO
	*/
	// string line;
	// ifstream cond_file_1("cond_1.csv");
	// int row = 0;
	// int col = 0;
	// int total = 0;
	// if(cond_file_1.is_open()){

	// 	while(getline(cond_file_1,line)){ // cond_1.csv stores cond row-major way, scan file row-by-row
	// 		conductances[col][row] = stod(line);
	// 		col += 1;
	// 		if(col > 99){
	// 			col = 0;
	// 			row += 1;
	// 		}
	// 	}
	// }
	// cond_file_1.close();

	ifstream cond_file("cond_1.csv");
	string cond_line;
	
	if(cond_file.is_open()){
		string delimiter = ",";
		int row = 0;

		while(getline(cond_file,cond_line)){ // cond_1.csv stores cond row-major way. each line is a row of curs
			// fill up conductances row by row
			vector<string> tokens{}; // row of curs
			size_t pos = 0;
			while((pos = cond_line.find(delimiter)) != string::npos){
				tokens.push_back(cond_line.substr(0, pos));
				cond_line.erase(0, pos + delimiter.length());
			}
			for(int col=0;col<COL_SIZE;col++){
				conductances[row][col] = stod(tokens[col]);
			}
			row += 1;
		}
	}

	cond_file.close();


	/////////////////////////////////////////



	/////////////////////////////////////////
	/*
		generate current, sweep the accumulative current of CHUNK rows
	
	*/
	// CHUNK = 20
	// for(int i=CHUNK; i<= 400; i+=CHUNK){ 
		
	// 	string file_name = "isum_";
	// 	file_name += to_string(i);
	// 	file_name += ".csv";

	// 	for(int j=0;j<i;j++){
	// 		input_bits[j] = 1;
	// 	}

	// 	current_gen_chunk(input_bits, conductances, ROW_SIZE, COL_SIZE, file_name);
	// }
	/////////////////////////////////////////

	/////////////////////////////////////////////
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
			cur_file_out << setprecision(10) << ReadCell(col, row, conductances[row][col]) << ",";
		}
		cur_file_out << endl;
	}

	cur_file_out.close();
	/////////////////////////////////////////////

    /////////////////////////////////////////////
	/*
		scan individual current, and determine: 
		(1). if they can be stolen by single row activation
		(2). if so, which one of the 64 levels they fell into
		(3). can we map them to their conductnce
	*/
	char num_thermo_levels = 64;
	vector<double> thermo_cur_levels(num_thermo_levels);
	double min_cur = 0; // micro ampere
	double max_cur = 8e-6;
	double cur_increment = (max_cur - min_cur) / num_thermo_levels;
	for(int i=1;i<=num_thermo_levels;i++){
		thermo_cur_levels[i-1] = cur_increment * i + min_cur;
		printf("thermo_cur_levels[%d]: %.4e\n",i-1, thermo_cur_levels[i-1]);
	}

	// input cur sweep file
	ifstream cur_file_in(cur_file_name);
	string cur_line;
	vector<vector<double>> cur_vector(ROW_SIZE, vector<double>(COL_SIZE)); // row-major
	
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
	/////////////////////////////////////////////

	/////////////////////////////////////////////
	/*
		try to reverse engineer cond based on thermo code
		(1). Each thermo code level correspond to a current range
		(2). for each current range, "guess" the corresponding conductance --> coarse-grain, lower precision weight
		(3). optimization: for consecutive chunks that are made of several "under powered" 
	*/
	double min_cond = 100e-9;	    // ideal Minimum cell conductance (S)
	double max_cond = 5e-6;		// ideal Maximum cell conductance (S)
	double cond_level = 64;
	vector<double> cond_vector(cond_level);

	// fill cond_vector  
	double cond_increment = (max_cond - min_cond) / (cond_level-1);
	for(int i=0;i<cond_level;i++){
		printf("conductances[%d]: %.4e\n", i, cond_increment * i + min_cond);
		// printf("%.4e\n",cond_increment * i + min_cond);
		cond_vector[i] = cond_increment * i + min_cond;
	}


	vector<vector<double>> stolen_cond_vec(ROW_SIZE, vector<double>(COL_SIZE)); // row-major
	
	/////////////////////////////////////////////
	// baseline way, recover only recoverable. if below min thermo level, replace with smallest cond
	// note: ADC range needs to be calibrated (thermo_cur_levels)
	// step-1: pre-calculate the current steps generated by conductances: cond_cur_vec[cond_level]: [..., 0.5(v) / (1/cond_i), ...]
	// step-1: for each cell, identify the corresponding thermo current cur_thermo
	// step-2: from cond_cur_vec, find the a current value that is closest to cur_thermo, get the index, and retrieve the conductance accordingly
	cout << "\nBaseline attack ... " << endl;

	// step-1: 
	double readVoltage = 0.5;
	vector<double> cond_cur_vec(cond_level);
	for(int i=0;i<cond_level;i++){
		cond_cur_vec[i] = readVoltage / (1/cond_vector[i]);
	}


	// step-2: find appropriate thermo_cur, "cur_thermo" is thermo_cur_levels[i]
	for(int row=0;row<ROW_SIZE;row++){
		for(int col=0;col<COL_SIZE;col++){
			
			bool found_thermo_cur = false;
			for(int i=(num_thermo_levels-1);i>=0;i--){
				if(cur_vector[row][col] >= thermo_cur_levels[i]){
					if(row==399 && col > 95){ // DEBUG
						printf("cur_vector[%d][%d]: %.4e uA, thermo_cur_levels[%d]: %.4e uA\n",row, col, cur_vector[row][col], i, thermo_cur_levels[i]);
					}
					// step-3: find corresponding (closest) current from cond_cur_vec
					found_thermo_cur = true;
					for(int j=(cond_level-1);j>=0;j--){
						if(thermo_cur_levels[i] >= cond_cur_vec[j]){
							// stolen_cond_vec[row][col] = cond_vector[j];
							// ad-hoc: consider calibration 
							stolen_cond_vec[row][col] = cond_vector[j+3];
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

	// store baseline stolen conds row by row
	string baseline_stolen_file_name = "baseline_stolen.csv";
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
	/////////////////////////////////////////////


	




	/////////////////////////////////////////////




	
	return 0;

}
