//
// Created by kot4or on 8/8/16.
//

#include "interval_map.h"


bool find_start_segment_annotation (BamRecordPtr current_bam_record, BamRecord previous_bam_record, interval_map<long, MapElement>::iterator & current_gtf_records_splitted_it, bool & freeze){
    static bool allow_skip_rest;
    assert (current_bam_record.use_count() > 0);
    if (current_bam_record->read_id == previous_bam_record.read_id and allow_skip_rest){
        freeze = false;
        return false;
    }

    if (current_bam_record->start_pose < current_gtf_records_splitted_it->first.lower()) {
        cout << "   Skip read " << current_bam_record->read_id << " [" <<
             current_bam_record->start_pose << "," <<
             current_bam_record->end_pose << "]" << endl;
        cout << "     " << current_bam_record->start_pose << " < " << current_gtf_records_splitted_it->first.lower() << endl;
        freeze = false; // Set freeze to false to change current_bam_record
        allow_skip_rest = true;
        return false;
    }
    allow_skip_rest = false;
    if (current_bam_record->start_pose >= current_gtf_records_splitted_it->first.upper()) {
//        cout << current_bam_record->start_pose << " > " << current_gtf_records_splitted_it->first.upper() << endl;
        cout << "   Skip segment annotation : " << "[" <<  current_gtf_records_splitted_it->first.lower() << ","
             << current_gtf_records_splitted_it->first.upper() << "]" << endl;
        current_gtf_records_splitted_it++;
        freeze = true; // Set freeze to true to prevent changing current_bam_record
        return false;
    }
    return true;



//    if ( current_bam_record->start_pose >= current_gtf_records_splitted_it->first.lower() &&
//         current_bam_record->start_pose <= current_gtf_records_splitted_it->first.upper() ){
//        allow_skip_rest = false;
//        return true;
//    } else if ( current_bam_record->start_pose < current_gtf_records_splitted_it->first.lower() ){
//        cout << "   Skip read " << current_bam_record->read_id << " [" <<
//             current_bam_record->start_pose << "," <<
//             current_bam_record->end_pose << "]" << endl;
//        cout << "     " << current_bam_record->start_pose << " < " << current_gtf_records_splitted_it->first.lower() << endl;
//        freeze = false; // Set freeze to false to change current_bam_record
//        allow_skip_rest = true;
//        if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
//            cerr << "   Skip read " << current_bam_record->read_id << " [" <<
//                 current_bam_record->start_pose << "," <<
//                 current_bam_record->end_pose << "]" << endl;
//            cerr << "     " << current_bam_record->start_pose << " < " << current_gtf_records_splitted_it->first.lower() << endl;
//        }
//        return false;
//    } else if ( current_bam_record->start_pose > current_gtf_records_splitted_it->first.upper() ){
//        //        cout << current_bam_record->start_pose << " > " << current_gtf_records_splitted_it->first.upper() << endl;
//        cout << "   Skip segment annotation : " << "[" <<  current_gtf_records_splitted_it->first.lower() << ","
//             << current_gtf_records_splitted_it->first.upper() << "]" << endl;
//        if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
//            cerr << "   Skip segment annotation : " << "[" <<  current_gtf_records_splitted_it->first.lower() << ","
//                 << current_gtf_records_splitted_it->first.upper() << "]" << endl;
//        }
//        current_gtf_records_splitted_it++;
//        freeze = true; // Set freeze to true to prevent changing current_bam_record
//        allow_skip_rest = false;
//        return false;
//    }

}

bool find_stop_segment_annotation (BamRecordPtr current_bam_record,
                                   interval_map<long,MapElement>::iterator & temp_gtf_records_splitted_it,
                                   interval_map<long,MapElement>::iterator max_segment_annotation,
                                   bool & freeze){
    assert (current_bam_record.use_count() > 0);
    while (true){
        if ( current_bam_record->end_pose >= temp_gtf_records_splitted_it->first.lower() &&
             current_bam_record->end_pose <= temp_gtf_records_splitted_it->first.upper() ){
            return true;
        } else if (temp_gtf_records_splitted_it == max_segment_annotation){
            freeze = false;
            return false;
        } else {
            temp_gtf_records_splitted_it++;
        }
        if (current_bam_record->end_pose < temp_gtf_records_splitted_it->first.lower()){
            freeze = false;
            return false;
        }
    }
//    while (current_bam_record->end_pose < temp_gtf_records_splitted_it->first.lower() or
//           current_bam_record->end_pose > temp_gtf_records_splitted_it->first.upper()){
//        // check if we reached the end of the gtf_records array. If true, change go to the next bam_record
//        if (temp_gtf_records_splitted_it == max_segment_annotation){
//            freeze = false;
//            return false;
//        }
//        temp_gtf_records_splitted_it++;
//    }
//    return true;
}

void print_segment_annotation (const string & title, interval_map<long, MapElement>::iterator current_gtf_records_splitted_it){
    cout << title << " " << "[" << current_gtf_records_splitted_it->first.lower() << "," << current_gtf_records_splitted_it->first.upper() << "] :";
    for (auto start_segment_annotation_it = current_gtf_records_splitted_it->second.gtf_records.begin();
         start_segment_annotation_it != current_gtf_records_splitted_it->second.gtf_records.end(); ++start_segment_annotation_it){
        assert (start_segment_annotation_it->use_count() > 0);
        cout << " " << start_segment_annotation_it->get()->exon_id << " (" << start_segment_annotation_it->get()->isoform_id << ")";
    }
    cout << endl;
}

// Returns set of annotation pointers which are made as a result of intersection between two annotation pointers sets from two interval map segments
set<GffRecordPtr> get_intersection (interval_map<long, MapElement>::iterator input_1, interval_map<long, MapElement>::iterator input_2){
    set<GffRecordPtr> intersection;
    input_1->second.gtf_records.sort();
    input_2->second.gtf_records.sort();
    set_intersection(input_1->second.gtf_records.begin(),input_1->second.gtf_records.end(),
                     input_2->second.gtf_records.begin(),input_2->second.gtf_records.end(),
                     std::inserter(intersection, intersection.begin()));
    cout << "   Intersection : ";
    for (auto intersection_segment_annotation_it = intersection.begin();
         intersection_segment_annotation_it != intersection.end(); ++intersection_segment_annotation_it){
        assert (intersection_segment_annotation_it->use_count() > 0);
        cout << " " << intersection_segment_annotation_it->get()->exon_id << " ("<< intersection_segment_annotation_it->get()->isoform_id << "), ";
    }
    cout << endl;
    return intersection;
}

bool fit_spliced_read_condition(const long & current_slice, const long & slice_number, const set <GffAndStartStopIt> & temp_set, BamRecordPtr current_bam_record){
    assert (temp_set.begin()->annotation.use_count() > 0);
    if (slice_number == 1) return true; // added just in case to make sure that we are not trying ot check unspliced read
    if (current_slice == 1 and temp_set.begin()->annotation->end_pose != current_bam_record->end_pose ) {
        return false;
    } else
    if (current_slice == slice_number and temp_set.begin()->annotation->start_pose != current_bam_record->start_pose){
        return false;
    } else
    if (current_slice > 1 and
        current_slice < slice_number and
        temp_set.begin()->annotation->end_pose != current_bam_record->end_pose and
        temp_set.begin()->annotation->start_pose != current_bam_record->start_pose){
        return false;
    }
    return true;
}

// if input_set forms linked list, for each of the element, except the first one, I can find the previous one.
// If not - return false;
bool form_line (const set<GffAndStartStopIt> & complete_input_set){
    cout << "Check if set of annotation forms linked list" << endl;

    set<GffRecordPtr> input_set;
    for (auto it = complete_input_set.begin(); it != complete_input_set.end(); ++it){
        input_set.insert(it->annotation);
    }

    if (input_set.size() == 1) return true;
    int counter = 0;
    for (auto it = input_set.begin(); it != input_set.end(); ++it){
        assert (it->use_count() > 0);
        if (it->get()->previous_gff.use_count() > 0){
            auto it_check = input_set.find(it->get()->previous_gff);
            if (it_check == input_set.end()) {
                cout << "for element " << it->get()->exon_id << " from " << it->get()->isoform_id << " cannot find any previous" << endl;
                counter++;
            } else{
                assert (it_check->use_count() > 0);
                cout << "for element " << it->get()->exon_id << " from " << it->get()->isoform_id
                     << " found previous " << it_check->get()->exon_id << " from " << it_check->get()->isoform_id <<  endl;
            }
        } else {
            cout << "Found first element of isoform: " << it->get()->exon_id << " from " << it->get()->isoform_id << endl;
            counter++;
        }
    }
    if (counter > 1){
        cout << "counter > 1. Impossible to have two or more unlinked elements" << endl;
        return false;
    }
    return true;
}