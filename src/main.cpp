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

//#include "interval_map.h"
#include "rpkm_calculation.h"

#include "test.h"


using namespace std;
using namespace boost::icl;
using namespace BamTools;

bool test_mode = false;


int main(int argc, char **argv) {
    // Read the paths from arguments
    if (argc < 3){
        cerr << "Set <full path to bam-file> <full path to tab-delimited file>" << endl;
        return 0;
    }

    // if put --test instead of path to the log file
    if ( argc > 3 && string(argv[3]) == "--test" ){
        test_mode = true;
    }

    string test_results_path = "";
    if (argc > 4 && test_mode ){
        test_results_path = string (argv[4]);
    }


    if (argc > 3 && (!test_mode) ){
        string log_filename = string(argv[3]);
        cerr << "Log file " << log_filename << endl;
        freopen(log_filename.c_str(), "a", stdout); // TODO Check what happens when filename is not correct
        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        cout << (now->tm_year + 1900) << '-'
             << (now->tm_mon + 1) << '-'
             << now->tm_mday << "   "
             << now->tm_hour << ":" << now->tm_min
             << endl << endl;
    }

    string results_path = "";
    if (argc > 4 && (!test_mode) ){
        results_path = string (argv[4]);
    }

    // Set paths to bam and annotation files
    string bam_full_path_name = string(argv[1]);
    string annotation_full_path_name = string(argv[2]);

    // read from BAM file
    BamReader bam_reader;
    if (not bam_reader.Open(bam_full_path_name)) {
        cout << "Couldn't open file " << bam_full_path_name << endl;
        return 0;
    } else cout << "Open " << bam_reader.GetFilename() << endl;

    cout << endl << endl;

    // key - chromosome name, value - <RefId, Length> for corresponding chromosome from the BAM file
    // TODO chromosome_info_map - saves correspondence between chromosome name and RefId from the BamReader object
    std::map <string, pair <int, int> > chromosome_info_map = get_chromosome_map_info (bam_reader);

    print_ref_info (chromosome_info_map); // Only for DEBUG

    // Check if current bam file is indexed (and that index data is loaded into program)
    if (not make_index(bam_reader)){
        return 0;
    }

    // read from tab delimited file
    // map < chromosome_key, multimap < exon_start_pose, exon_pointer> >
    // exon_start_pose - is needed for sorting mutlimap by the start pose of exon; we cannot do if we add only pointers
    // TODO global_annotation_map_ptr - map to store all annotation data from tab delimited file
    std::map <string, multimap <long, GffRecordPtr> > global_annotation_map_ptr;

    // map to save <chromosome name, <isoform name, correspondent Isoform object> >
    std::map <string, std::map <string, Isoform> > iso_var_map;
    if (not load_annotation (annotation_full_path_name, global_annotation_map_ptr, iso_var_map)){
        return 0;
    }
    cout << endl;
    print_iso_var_map (iso_var_map); // for DEBUG only
    cout << endl;

    // FOR DEBUG USE ONLY
    //    cout << "ANNOTATIONS" << endl;
    //    for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
    //        cout << "Chromosome: " << chrom_it->first << endl;
    //        for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
    //            assert (start_it->second.use_count() > 0);
    //            cout << "  " << start_it->second->exon_id  << " " << start_it->second->isoform_id << " - [";
    //            cout << start_it->second->start_pose << "," << start_it->second->end_pose << "]" << endl;
    //            for (auto bam_record_it = start_it->second->bam_records.begin(); bam_record_it != start_it->second->bam_records.end(); ++bam_record_it){
    //                assert (bam_record_it->use_count() > 0);
    //                cout << bam_record_it->get()->read_id << " ";
    //            }
    //        }
    //    }

    cerr << "Gathering info about bam file" << endl;
    BamGeneralInfo bam_general_info;
    // Need to run through the whole bam file, because when we align reads according to exons, we can skip some of its parts
    get_bam_info (bam_reader, bam_general_info);

    for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it) {

        cout << endl;
        string chrom = chrom_it->first;
        cout << "Current chromosome from annotation file: " << chrom << endl;
        auto map_it = chromosome_info_map.find(chrom);
        if (map_it == chromosome_info_map.end()){
            cout << "Cannot locate RefId for " << chrom << endl;
            cout << "Skip the whole chromosome from annotation file" << endl;
            continue;
        }
        int ref_id = map_it->second.first;
        int length = map_it->second.second;
        cout << "Found corresponding RefId:  " << ref_id << endl;
        if (not bam_reader.HasIndex()){
            cout << "ERROR: Current bam file isn't indexed" << endl;
            cout << "EXIT. Find a bug in a code. You suppose to index BAM file right after opening" << endl;
            return 0;
        }
        cout << "Current BAM file is indexed" << endl;
        cout << "Trying to set region limited by current chromosome: " << chrom << endl;


        // TODO Could couse a problem of mapped reads count, because we skip regions
        if (not bam_reader.SetRegion(ref_id, 0, ref_id, length)){
            cout << "Cannot set region. Exit" << endl;
            cout << bam_reader.GetErrorString(); // added just in case
            return 0;
        }
        cout << "Region from BAM file is succesfully set for chromosome " << chrom << " with RefId = " << ref_id << endl << endl;

        // Making an interval map on the base of the annotation.
        // By default for each current_map_element we add one the corresponding annotation pointer
        interval_map<long, MapElement> gtf_records_splitted;
        for (auto it = chrom_it->second.begin(); it != chrom_it->second.end(); ++it) {
            assert (it->second.use_count() > 0);
            MapElement current_map_element;
            current_map_element.gtf_records.push_front(it->second);
            assert (it->second->start_pose < it->second->end_pose);
            gtf_records_splitted.add( make_pair(interval<long>::right_open(it->second->start_pose, it->second->end_pose), current_map_element) );
        }

        // create an empty matrix: column - one interval from interval map, row - isoforms, initialize it with -1
        vector <vector <double> > weight_array (iso_var_map[chrom].size()+1, vector <double> (gtf_records_splitted.iterative_size(), 0));

        // Set the length of intervals into the first line of weight_array
        // put min_weight value instead of the 0 in all of the intervals where we have exons
        double min_weight = 1.0e-29; // TODO put it as argument

        int temp_n = 0;
        for (auto temp_it = gtf_records_splitted.begin(); temp_it != gtf_records_splitted.end(); ++temp_it) {
            double length = temp_it->first.upper() - temp_it->first.lower();
            assert (length > 0);
            weight_array[0][temp_n] = length;
            for (auto  gtf_it = temp_it->second.gtf_records.begin(); gtf_it != temp_it->second.gtf_records.end(); ++gtf_it){
                GffRecordPtr temp_gtf_ptr = *gtf_it;
                weight_array [ iso_var_map[chrom][temp_gtf_ptr->isoform_id].index ] [temp_n] = min_weight;
            }
            temp_n++;
        }

        // FOR DEBUG ONLY
//        print_weight_array (weight_array);



        // FOR DEBUG USE ONLY
        //        cerr << "GENERATE INTERVAL MAP ONLY FOR NM_001198798" << endl;
        //        for (auto temp_it = gtf_records_splitted.begin(); temp_it != gtf_records_splitted.end(); ++temp_it) {
        //            bool print_it = false;
        //            std::stringstream list_of_exons;
        //            for (auto it = temp_it->second.gtf_records.begin(); it != temp_it->second.gtf_records.end(); ++it) {
        //                assert (it->use_count() > 0);
        //                list_of_exons << "      ex:" << it->get()->exon_id << " (" << it->get()->isoform_id << ")";
        //                if (it->get()->isoform_id  == "NM_001198798"){
        //                    print_it = true;
        //                }
        //            }
        //            if (print_it){
        //                cerr << "[" << temp_it->first.lower() << "," << temp_it->first.upper() << "]   ";
        //                cerr << list_of_exons.str();
        //                cerr << endl;
        //            }
        //        }
        //        cerr << endl;


        // The main iterator to iterate over interval map segment array
        interval_map<long, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();
        // Additinal iterator to backup current_gtf_records_splitted_it iterator when we are processing spliced read]
        // If the read is standard the backup iterator doesn't influence on the result
        interval_map<long, MapElement>::iterator backup_current_gtf_records_splitted_it = gtf_records_splitted.begin();

        BamRecordPtr current_bam_record; // pointer to save current bam record, which we are trying to align
        bool freeze = false; // if true - calling the get_bam_record function return's the same read as it it did it before
        int slice_number = 1; // if read is spliced slice_number > 1
        int current_slice = 0; // defines the position of current read as a part of a big spliced read
        std::map<string, set<GffAndStartStopIt> > iso_map; // map to arrange exons according to the isoform key
        BamRecord previous_bam_record; // temporal bam record to detect the moment when next bam record isn't a part of big scpliced read

        cerr << "Processing reads" << endl;
        int reads_tem_count = 0;
        static int lower_index = 0;
        while (get_bam_record(bam_reader, current_bam_record, freeze)) { // Check if I can get new record from BAM file
            //            cerr << "*" << endl;
            //            if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
            //                cerr << "Get read A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549" << endl;
            //                cerr << " [" << current_bam_record->start_pose << "," <<
            //                                current_bam_record->end_pose << "] - " <<
            //                                current_bam_record->slices << " parts" << endl;
            //            }
            reads_tem_count++;
            if (reads_tem_count % 1000 == 0){
                cerr << "*";
            }
            if (reads_tem_count > 6000){
                cerr << "\r";
                reads_tem_count = 0;
            }
            // Check if gtf records array is already empty. Break the while loop
            if (current_gtf_records_splitted_it == gtf_records_splitted.end()){
                cerr << "reached the end of interval map 1" << endl;
                break;
            }

            assert (current_bam_record.use_count() > 0);

            // FOR DEBUG USE ONLY
            //            cout << endl << "Process read " << current_bam_record->read_id << " [" <<
            //            current_bam_record->start_pose << "," <<
            //            current_bam_record->end_pose << "] - " <<
            //            current_bam_record->slices << " parts" << endl;

            // get the number of parts from current bam record
            slice_number = current_bam_record->slices;

            // if the next read isn't a part of the same big spliced read (in this case we check if read_id is equal)
            // we set current_slice equal to 0 and clear iso_map, because we don't want to have annotation which are relevant for the previous read
            // make a backup for interval map iterator.
            // We need to make a backup because the next read could be a part of spliced read
            // and in that case we will change current_gtf_records_splitted_it which we will restore from
            // its backup version when we rich the last part in spliced read
            // TODO maybe it should be the other way to check it
            if (previous_bam_record.read_id != current_bam_record->read_id) {
                current_slice = 0;
                iso_map.clear();
                backup_current_gtf_records_splitted_it = current_gtf_records_splitted_it; // update the backup for gff iterator
                cout << "iso_map is cleared because of the new read" << endl;
            }

            // increment current_slice (position inside a spliced read)
            if (not freeze){
                current_slice++;
            }

//            if (backup_current_gtf_records_splitted_it->first.lower() != lower_index){
//                cerr << "backup_current_gtf_records_splitted_it = " << backup_current_gtf_records_splitted_it->first.lower() << endl;
//                assert (lower_index <= backup_current_gtf_records_splitted_it->first.lower());
//                lower_index = backup_current_gtf_records_splitted_it->first.lower();
//            }

            // Find start segment annotation
            // If false call get_bam_record again. freeze is true if we need to skip the read and false if we want to skip annotation
            if (not find_start_segment_annotation(current_bam_record, previous_bam_record, current_gtf_records_splitted_it, freeze)){
                previous_bam_record = *current_bam_record; // update previous_bam_record with current value
                //                if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
                //                    cerr << "For read A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549 didn't find any start interval" << endl;
                //                }
                if (current_gtf_records_splitted_it == gtf_records_splitted.end()){
                    cerr << "reached the end of interval map 2" << endl;
                    break;
                }
                if (current_slice == current_bam_record->slices && current_bam_record->slices > 1 && (not freeze)){ // TODO Double check it
                    current_gtf_records_splitted_it = backup_current_gtf_records_splitted_it;
                }
                continue;
            }
            previous_bam_record = *current_bam_record; // update previous_bam_record with current value

            // FOR DEBUG USE ONLY
            //            print_segment_annotation("   start segment annotation", current_gtf_records_splitted_it);

            // Find stop segment annotation
            // If false call get_bam_record with freeze = false ===> skip current bam record and get the next one
            if (current_gtf_records_splitted_it == gtf_records_splitted.end()){
                cerr << "reached the end of interval map 3" << endl;
                break;
            }
            interval_map<long, MapElement>::iterator temp_gtf_records_splitted_it = current_gtf_records_splitted_it;

            if (not find_stop_segment_annotation(current_bam_record, temp_gtf_records_splitted_it, gtf_records_splitted.end(), freeze)){
                //                if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
                //                    cerr << "For read A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549 didn't find any stop interval" << endl;
                //                };
                continue;
            }

            // FOR DEBUG USE ONLY
            //            print_segment_annotation("   stop segment annotation", temp_gtf_records_splitted_it);



            // find intersection of two sets
            // we receive set of annotation which is similar between start and end interval map segments
            // in other words we receive pointers to the annotations which includes current bam read
            set<GffRecordPtr> gff_intersection = get_intersection(current_gtf_records_splitted_it, temp_gtf_records_splitted_it);
            // FOR DEBUG USE ONLY
            //            if (current_bam_record->read_id == "A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549"){
            //                cerr << "Intersection set for A9DF6E9B-909F-CBB5-A534-DDFB5895BE2A.fastq.12378549" << endl;
            //                for (auto intersect_iter = gff_intersection.begin(); intersect_iter != gff_intersection.end(); ++intersect_iter){
            //                    cerr <<  intersect_iter->get()->isoform_id << " - "<< intersect_iter->get()->exon_id << endl;
            //                }
            //            };


            // TODO maybe we'll need it when have to filter by strand
            // Filter gff_intersection for only those GffRecordPtr, who has the same strand as current_bam_record
//            for (auto intersect_iter = gff_intersection.begin(); intersect_iter != gff_intersection.end(); ) {
//                if (intersect_iter->get()->strand != current_bam_record->strand){
//                    intersect_iter = gff_intersection.erase (intersect_iter);
//                    cerr << "erased GffRecord because of the strand" << endl;
//                } else {
//                    ++intersect_iter;
//                }
//            }

            // iterate over the intersection set
            for (auto gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it) {
                set<GffAndStartStopIt> temp_set; // we save only one pointer to this set, but we need to use set, because we want to add it into iso_map
                GffRecordPtr temp_ptr (*gff_it);
                temp_set.insert(GffAndStartStopIt (temp_ptr, current_gtf_records_splitted_it, temp_gtf_records_splitted_it));
                // check if slice_number == 1 which is equal that read isn't spliced or
                // current annotation fits the condition which are required for specific part of the spliced read (start, middle or end part of spliced read)
                if (slice_number == 1 or fit_spliced_read_condition(current_slice, slice_number, temp_set, current_bam_record)) {
                    // Add new element into map
                    pair<std::map<string, set<GffAndStartStopIt> >::iterator, bool> ret;
                    pair<string, set<GffAndStartStopIt> > input_pair(gff_it->get()->isoform_id, temp_set);
                    ret = iso_map.insert(input_pair);
                    // If already exist with the same isoform key - add annotation to the corresponding set of the map
                    if (ret.second == false) {
                        cout << "Isoform " << input_pair.first << " already exists" << endl;
                        cout << "Updating the set with a value " << gff_it->get()->exon_id << endl;
                        ret.first->second.insert(temp_set.begin(), temp_set.end());
                    }
                }
            }

            // If we reached the end of the big spliced read (equal to current_slice >= slice_number)
            if (current_slice >= slice_number) {
                // iterating over the isoforms from iso_map to get divider
                double global_koef = slice_number; // defines the weight koef of each slice of the spliced read
                double vertical_koef = 0; // defines the weight koef if read is aligned on few isoforms
                for (auto map_iterator = iso_map.begin(); map_iterator != iso_map.end(); map_iterator++) {
                    if (map_iterator->second.size() == slice_number){
                        if (not form_line (map_iterator->second)) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                    vertical_koef++;
                }

                // iterating over the isoforms from iso_map to fill exons' data
                for (auto map_iterator = iso_map.begin(); map_iterator != iso_map.end(); map_iterator++) {
                    // if current isoform has the same number of exons as number of parts from the spliced read
                    // we caught the correct spliced read
                    // if not - go to next isoform
                    if (map_iterator->second.size() == slice_number) {
                        cout << "Found possible location of spliced read" << endl;
                        if (not form_line (map_iterator->second)) {
                            cout << "Current set of annotations doesn't form linked list" << endl;
                            continue; // go to next isoform if annotation in current set don't form an "unbreakable" line
                        }

                        // iterating iver annotation of one specific isoform from iso_map and add to each of them
                        // pointer to a current bam record
                        for (auto gff_it = map_iterator->second.begin(); gff_it != map_iterator->second.end(); ++gff_it) {
                            BamRecordPtr bam_record_to_put_in_array = current_bam_record; // put it just in case. should be the same as just push original pointer to array, 'cos push copies
                            gff_it->annotation->bam_records.push_back(bam_record_to_put_in_array);
                            gff_it->annotation->reads_count++;
                            cout <<"[***]" << endl;
                            // FOR DEBUG USE ONLY
                            //                            if (gff_it->annotation->isoform_id == "NM_001198798"){
                            //                                cerr << "to NM_001198798 added: " << bam_record_to_put_in_array->read_id << endl;
                            //                            }
                            assert (gff_it->annotation.use_count() > 0);
                            assert (bam_record_to_put_in_array.use_count() > 0);
//                            cout << "   into annotation " << gff_it->annotation->exon_id << " from " << gff_it->annotation->isoform_id << " added read " << bam_record_to_put_in_array->read_id << endl;
                            // updated weight array
                            string isoform = gff_it->annotation->isoform_id;
                            int j = iso_var_map[chrom][isoform].index;
                            cout << "j = " << j << endl;
                            long start_i = distance (gtf_records_splitted.begin(), gff_it->start_it);
                            long stop_i = distance (gtf_records_splitted.begin(), gff_it->stop_it);
                            int horizontal_koef = distance (gff_it->start_it, gff_it->stop_it) + 1; // defines the weight koef if read occupies few intervals of some exon
//                            cout << "start_i = " << start_i << endl;
//                            cout << "stop_i = " << stop_i << endl;
//                            cout << "length = " << horizontal_koef << endl;
                            // FOR DEBUG ONLY
                            if ((global_koef * vertical_koef * (double)horizontal_koef) == 0) {
                                cerr << "Something went wrong. Find a bug" << endl;
                                throw ("Error: dividing by zero");
                            }
                            double weight = 1 / (global_koef * vertical_koef * (double)horizontal_koef);
//                            cout << "global_koef = " << global_koef << endl;
//                            cout << "vertical_koef = " << vertical_koef << endl;
//                            cout << "horizontal_koef = " << horizontal_koef << endl;
//                            cout << "weight = " << weight << endl;
                            for (int i = start_i; i <= stop_i; i++) {
                                weight_array[j][i] += weight;
                            }
                        }
                    }
                }
                // at the end of the iteration over the isoform map, we need to revert value of current_gtf_records_splitted_it
                // from its backup version
                current_gtf_records_splitted_it = backup_current_gtf_records_splitted_it; // get iterator from the backup
            }
            // set freeze to false to get new read from the bam file when calling function get_bam_record
            freeze = false;
        }

        if ( test_mode ) print_weight_array_test (weight_array, "Weight array", test_results_path);


        // Original weight array
        //        print_weight_array(weight_array, "Original weight array");

//        cerr << "Original weight array" << endl;
        //        print_isoform_by_name (weight_array, iso_var_map, "chr10", "NM_001198798", cout);

        // Transormed to density
        transform_to_density (weight_array);
//        print_weight_array(weight_array, "Original density array");

//        cerr << "Original density array" << endl;
//        print_isoform_by_name (weight_array, iso_var_map, "chr10", "NM_001198798", cout);

        cerr << "Started to run cycles" << endl;
        int cycles = run_cycle (weight_array);
        cerr << "Finished to run cycles : " << cycles << endl;
//        print_weight_array(weight_array, "Final density array");
        cout << endl;
        calculate_totReads_density (weight_array, iso_var_map[chrom]);
        cout << endl << endl;


    }

    cout << "mapped_reads_counter: " << bam_general_info.aligned <<endl;
    cout << "total_reads_counter: " << bam_general_info.total << endl;
    cerr << endl;
    cerr << "Calculate rpkm" << endl;

    calculate_rpkm (iso_var_map, bam_general_info.aligned);
    print_iso_var_map (iso_var_map);

    if (results_path.length() > 0){
        cerr << "Exporting results" << endl;
        print_iso_var_map_to_file (iso_var_map, results_path);
    }


    // FOR DEBUG USE ONLY
//    cout << endl << "RESULTS" << endl;
//    for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
//        cout << "Chrom: " << chrom_it->first << endl;
//        for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
//            cout << "exon " << (*start_it->second).exon_id  << " from " << (*start_it->second).isoform_id << " - [";
//            cout << (*start_it->second).start_pose << "," << (*start_it->second).end_pose << "]" << " : ";
//            for (auto bam_record_it = (*start_it->second).bam_records.begin(); bam_record_it != (*start_it->second).bam_records.end(); ++bam_record_it){
//                cout << (*bam_record_it)->read_id << " ";
//            }
//            cout << endl;
//        }
//    }

    return 0;
}
