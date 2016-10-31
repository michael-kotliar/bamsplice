//
// Created by kot4or on 8/8/16.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>


#ifndef TEST_1_RPKM_CALCULATION_H
#define TEST_1_RPKM_CALCULATION_H
#endif //TEST_1_RPKM_CALCULATION_H

using namespace std;

void print_weight_array(const vector<vector<double> > & weight_array, const string & title = "");
void transform_to_density ( vector <vector <double> > & weight_array);
double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row);
double get_average_by_row (const vector <vector <double> > & weight_array, const int & row);
double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row = 1);
double get_average_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row = 1);