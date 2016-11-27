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

using namespace std;

void process (   vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > chrom_vector,
                 std::map <string, pair <int, int> > chromosome_info_map,
                 std::map <string, std::map <string, Isoform> > & iso_var_map,
                 string bam_full_path_name,
                 int thread_number);