//
// Created by kot4or on 8/8/16.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "interval_map.h"
#include <boost/thread/mutex.hpp>

#ifndef TEST_1_RPKM_CALCULATION_H
#define TEST_1_RPKM_CALCULATION_H
#endif //TEST_1_RPKM_CALCULATION_H


static boost::mutex iso_var_map_mutex;

using namespace std;
void print_isoform_by_name (const vector<vector<double> > & data_array,
                            std::map <string, std::map <string, Isoform> > iso_var_map,
                            string chr,
                            string isoform_name,
                            std::ostream & out = cout);
void print_weight_array(const vector<vector<double> > & weight_array, const string & title = "");
void transform_to_density ( vector <vector <double> > & weight_array);
double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row, int & count);
double get_average_by_row (const vector <vector <double> > & weight_array, const int & row);
double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int & count);
double get_average_by_column (const vector <vector <double> > & weight_array, const int & column);
vector <double> get_average_density_by_all_isoforms (const vector <vector <double> > & weight_array);
vector <double> get_sum_density_by_all_intervals (const vector <vector <double> > & weight_array);
void update_isoforms_density_to_average_for_isoform (vector <vector <double> > & weight_array);
void print_array (const vector <double> & intput_array, const string & title, string tab_spacer = "   ");
void adjust_isoforms_density_by_coef (vector <vector <double> > & weight_array, const vector <double> & original_densities_by_isoforms,  const vector <double> & new_densities_by_isoforms);
int run_cycle (vector <vector <double> > & weight_array);
double sum_all (const vector <vector <double> > dens_matrix);
void subtract_matrix (vector <vector <double> > & first, const vector <vector <double> > & second);
void calculate_totReads_density (const vector<vector<double> > & weight_array, std::map <string, Isoform> & iso_map,
                                 const std::map <string, int> & correspondence_map, int cycles = 0, string bin_id = ".");
void calculate_rpkm (std::map <string, std::map <string, Isoform> > & iso_var_map, const long & aligned);
