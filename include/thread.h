//
// Created by kot4or on 11/26/16.
//

#ifndef GEEP_THREAD_H
#define GEEP_THREAD_H

#endif //GEEP_THREAD_H

#include <iostream>
#include <fstream>
#include <cstdio>

#include <iomanip>

#include <set>
#include <map>

#include <string.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "rpkm_calculation.h"
#include <boost/thread.hpp>
#include "test.h"

using namespace std;

void process (   vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > chrom_vector,
                 std::map <string, pair <int, int> > chromosome_info_map,
                 std::map <string, std::map <string, Isoform> > & iso_var_map,
                 string bam_full_path_name,
                 int thread_number,
                 string test_results_path,
                 int min_length,
                 bool keep_unique);

void filter_weight_array (  vector<vector<double> > & weight_array,
                            const interval_map<long, MapElement> & gtf_records_splitted,
                            const std::map <string, int> & correspondence_map,
                            double min_weight,
                            double min_length );