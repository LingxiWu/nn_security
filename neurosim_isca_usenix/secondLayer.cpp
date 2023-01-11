	// for the first 14 col of conductances, see if its possible to find rows to activate 
	// such that 1 col is 1 and 13 cols are zero
	// quick and dirty, gather large conductance in each cols until cross the current threshold
	

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

int main(){

	// array dimensions
	double ROW_SIZE = 400;
	double COL_SIZE = 100;

	// array of conductances
	vector<vector<double> > 
	conductances(ROW_SIZE, vector<double>(COL_SIZE));

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

	// printf("cond[%d][%d]: %.4e\n", 0,0,conductances[0][0]);
	// printf("cond[%d][%d]: %.4e\n", 0,1,conductances[0][1]);
	// printf("cond[%d][%d]: %.4e\n", 1,0,conductances[1][0]);
	// printf("cond[%d][%d]: %.4e\n", 1,1,conductances[1][1]);
	// printf("cond[%d][%d]: %.4e\n", 399,0,conductances[399][0]);
	// printf("cond[%d][%d]: %.4e\n", 0,99,conductances[0][99]);
	// printf("cond[%d][%d]: %.4e\n", 399,99,conductances[399][99]);

	


	return 0;
}