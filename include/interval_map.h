//
// Created by kot4or on 8/8/16.
//

#include "annotation_reader.h"

#ifndef TEST_1_INTERVAL_MAP_H
#define TEST_1_INTERVAL_MAP_H
#endif //TEST_1_INTERVAL_MAP_H


using namespace boost::icl;
static boost::thread_specific_ptr< bool > allow_skip_rest;


// Class for interval_map only
// if we are not planning to add any other fields, we can change it
// into the simple forward list of GffRecordPtr
class MapElement {
public:
    forward_list < GffRecordPtr > gtf_records; // by default unsorted list of pointers GffRecordPtr to annotations

    // when we add items due to operator += we don't care about the order
    inline MapElement& operator+=(const MapElement& other_element) {
        for (auto it = other_element.gtf_records.begin(); it != other_element.gtf_records.end(); ++it){
            this->gtf_records.push_front(*it);
        }
        return *this;
    }

    inline bool operator==(const MapElement& other_element) const
    {
        return this->gtf_records == other_element.gtf_records;
    }

};


// Class to save not only pointer to annotation, but its start and stop iterators from interval map
class GffAndStartStopIt {
public:
    GffRecordPtr annotation;
    interval_map<long, MapElement>::iterator start_it;
    interval_map<long, MapElement>::iterator stop_it;
    GffAndStartStopIt (GffRecordPtr new_annotation, interval_map<long, MapElement>::iterator new_start_it, interval_map<long, MapElement>::iterator new_stop_it)
            : annotation(new_annotation)
            , start_it (new_start_it)
            , stop_it (new_stop_it)
    {
    }
    bool operator<(const GffAndStartStopIt& other) const
    {
        return annotation < other.annotation;  //assume that you compare the record based on a
    }
};


// FUNCTIONs

// Updates iterator current_gtf_records_splitted_it for the interval map segment which includes the starting position of the current_bam_read
// If fails - returns false
// If we need to skip bam_record, we set freeze = false to get the new current_bam_record
// If we need to skip interval map segment we set freeze to true to avoid changing current_bam_record, when iterating over the loop,
// that get new value from bam records array
// If skip interval map segment - increment iterator current_gtf_records_splitted_it
bool find_start_segment_annotation (BamRecordPtr current_bam_record,
                                    BamRecord previous_bam_record,
                                    interval_map<long, MapElement>::iterator & current_gtf_records_splitted_it,
                                    bool & freeze);

// Updates iterator temp_gtf_records_splitted_it for interval map segment which includes stop position of current_bam_record
// If fail - returns false
bool find_stop_segment_annotation (BamRecordPtr current_bam_record,
                                   interval_map<long,MapElement>::iterator & temp_gtf_records_splitted_it,
                                   interval_map<long,MapElement>::iterator max_segment_annotation,
                                   bool & freeze);

// FOR DEBUG ONLY
void print_segment_annotation (const string & title, interval_map<long, MapElement>::iterator current_gtf_records_splitted_it);

// Returns set of annotation pointers which are made as a result of intersection between two annotation pointers sets from two interval map segments
set<GffRecordPtr> get_intersection (interval_map<long, MapElement>::iterator input_1, interval_map<long, MapElement>::iterator input_2);


// Check if current read is a part of the spliced read. If false - returns true and doesn't check any additional conditions
// If true, check if current annotation fits current read with following conditions:
// 1. If read is the first part of big spliced read and its end position is equal to end position of current annotation - return true;
// 2. If read is the last part of big spliced read and its start position is equal to start position of current annotation - return true;
// 3. If read is the middle part of big spliced read and both its start and end position are equal to start and end position of current annotation correspondingly - return true
// In all other cases return false
// NOTE temp_set has only one record. We can change it for GffRecordPtr and it send to the function just temp_set.begin()
bool fit_spliced_read_condition(const long & current_slice, const long & slice_number, const set <GffAndStartStopIt> & temp_set, BamRecordPtr current_bam_record);


// if input_set forms linked list, for each of the element, except the first one, I can find the previous one.
// If not - return false;
bool form_line (const set<GffAndStartStopIt> & complete_input_set);