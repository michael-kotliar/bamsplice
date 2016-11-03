//
// Created by kot4or on 8/8/16.
//

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <forward_list>
#include <boost/icl/interval_map.hpp>
#include <list>

#include "api/BamReader.h"
#include "api/BamMultiReader.h"

#ifndef TEST_1_BAM_READER_H
#define TEST_1_BAM_READER_H
#endif //TEST_1_BAM_READER_H


using namespace std;
using namespace BamTools;

class BamRecord;
typedef boost::shared_ptr<BamRecord> BamRecordPtr;

struct BamGeneralInfo {
    long total;
    long not_aligned;
    BamGeneralInfo ():
            total (0),
            not_aligned (0)
    {
    }
};


// CLASSes

// to save current read
class BamRecord {
public:
    long start_pose; // start position of the read
    long end_pose; // stop position of the read
    string read_id; // text identificator of the read
    short slices; // if the read is spliced - set here total amount of parts, that form this read, default = 1
    // ADD HERE OTHER INMPORTANT FIELDS

    // Constructor with parameters
    BamRecord (long start, long end, string read, short parts_number)
            : start_pose (start)
            , end_pose (end)
            , read_id (read)
            , slices (parts_number)
    {}

    // Empty constructor
    BamRecord ()
            : start_pose (0)
            , end_pose (0)
            , read_id (" ")
            , slices (1)
    {}

    void print ();
};

// FUNCTIONs

list <BamRecordPtr> split_to_single_reads (const BamAlignment & current_alignment);

bool get_bam_record (BamReader & bam_reader, BamRecordPtr & bam_record, BamGeneralInfo & bam_general_info, bool freeze = false);

void print_ref_info (const std::map <string, pair <int, int> > & info_map);

std::map <string, pair <int, int> > get_chromosome_map_info (const BamReader & reader);

bool make_index (BamReader & bam_reader);
bool flag_check (const BamAlignment & al, BamGeneralInfo & bam_general_info);


