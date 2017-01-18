//
// Created by kot4or on 8/8/16.
//

#include "bam_reader.h"
#include <iomanip>
#include <iostream>
#include <fstream>



#ifndef TEST_1_ANNOTATION_READER_H
#define TEST_1_ANNOTATION_READER_H
#endif //TEST_1_ANNOTATION_READER_H

using namespace std;
using namespace string_tools;


enum cds_stat {
    none,
    unk,
    incmpl,
    cmpl
};


class GffRecord;
typedef boost::shared_ptr<GffRecord> GffRecordPtr;

class GffRecord {
public:
    long start_pose; // start position of annotation
    long end_pose; // stop position of annotation
    string exon_id; // text-identificator of current annotation. Looks like we will use it only for debug
    string isoform_id; // set the name of the isoform to which current annotation belongs
    vector < BamRecordPtr > bam_records; // array of pointers to all of the reads, which belongs to this exon. Need this for debug, than we can delete this field
    int reads_count; // total number of reads, which belongs to this exon
    bool strand; // true for +
    bool start_exon;
    bool stop_exon;
    GffRecordPtr previous_gff; // ptr to the previous annotation in the same isoform. In NULL - first annotation in current isoform

    // CONSTRUCTOR WITH PARAMETERS
    GffRecord (long start, long end, string exon, string isoform, GffRecordPtr pre_gff, bool strnd, bool start_ex = false, bool stop_ex = false)
            : start_pose (start)
            , end_pose (end)
            , exon_id (exon)
            , isoform_id (isoform)
            , previous_gff (pre_gff)
            , reads_count (0)
            , strand (strnd)
            , start_exon (start_ex)
            , stop_exon (stop_ex)
    {}

    // EMPTY CONSTRUCTOR
    GffRecord ()
            : start_pose (0)
            , end_pose (0)
            , exon_id ("")
            , isoform_id ("")
            , previous_gff (NULL)
            , reads_count (0)
            , strand (false)
            , start_exon (false)
            , stop_exon (false)
    {}

};

class Isoform {
public:
    int bin; // 0
    string name; // 1
    string chrom; // 2
    bool strand;  // true - +, false  - - // 3
    long tx_start; // 4
    long tx_end; // 5
    long cds_start; // 6
    long cds_end; // 7
    int exon_count; // 8
    int length; // length of isofom = sum of all exons' lengths
    int total_reads; // total number of reads in isoform
    double density; // total_reads/length
    double rpkm; // rpkm
    int cycles; // number of iterations to balance the table
    string bin_id; // id to be able to find out which isoform have been processed in the same matrix
    int index; // index that will be used in weigth array to save reads there
    set <long> exon_starts;
    set <long> exon_ends;
//    set <long> exon_frames; // pointers to exon frames, not necessary to be sorted // 15

    double score; //11
    string name2; // 12
    cds_stat cds_start_stat; // 13
    cds_stat cds_end_stat; // 14

    void print ();

    // Constructor
    Isoform (string line, bool gtf = false);

    // Empty constructor
    Isoform ();

    Isoform& operator+=(const Isoform& other_iso);
};

bool str_to_cds_stat(const string &value, cds_stat &result);

void print_iso_var_map (const std::map <string, std::map <string, Isoform> > & iso_var_map);
void print_iso_var_map_to_file (const std::map <string, std::map <string, Isoform> > & iso_var_map, const string path);

// global_annotation_map_ptr : key - chromosome name, value - multimap of annotations, sorted by not-unique key - start pose of annotation
// NOTE : forward list of annotations should be sorted by start pose with rule a<b
bool load_annotation (const string & full_path_name,
                      std::map <string, multimap <long, GffRecordPtr> > & global_annotation_map_ptr,
                      std::map <string, std::map <string, Isoform> > & iso_var_map);

bool is_duplicate (const Isoform & original_isoform, const Isoform & new_isoform, bool gtf);

map <string, string> split_attributes (string line);