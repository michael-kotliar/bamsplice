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
double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row, int & count);
double get_average_by_row (const vector <vector <double> > & weight_array, const int & row);
double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int & count, int start_row = 1);
double get_average_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row = 1);
vector <double> get_average_density_by_all_isoforms (const vector <vector <double> > & weight_array, int start_row = 1);
vector <double> get_sum_density_by_all_intervals (const vector <vector <double> > & weight_array, int start_row = 1);
void update_isoforms_density_to_average_for_isoform (vector <vector <double> > & weight_array);
void print_array (const vector <double> & intput_array, const string & title, string tab_spacer = "   ");
void adjust_isoforms_density_by_coef (vector <vector <double> > & weight_array, const vector <double> & original_densities_by_isoforms,  const vector <double> & new_densities_by_isoforms);
void run_cycle (vector <vector <double> > & weight_array);


