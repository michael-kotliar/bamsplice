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
//#include "rpkm_calculation.h"

#include "test.h"
#include "thread.h"



using namespace std;
using namespace boost::icl;
using namespace BamTools;

bool test_mode = false;


int main(int argc, char **argv) {
    // Read the paths from arguments
    int threads_number = 1;

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
        cerr << "Couldn't open file " << bam_full_path_name << endl;
        return 0;
    } else cerr << "Open " << bam_reader.GetFilename() << endl;

    // key - chromosome name, value - <RefId, Length> for corresponding chromosome from the BAM file
    // chromosome_info_map - saves correspondence between chromosome name and RefId from the BamReader object
    std::map <string, pair <int, int> > chromosome_info_map = get_chromosome_map_info (bam_reader);

    //    print_ref_info (chromosome_info_map); // Only for DEBUG

    // Check if current bam file is indexed (and that index data is loaded into program)
    if (not make_index(bam_reader)){
        return 0;
    }

    cerr << "Gathering info about bam file" << endl;
    BamGeneralInfo bam_general_info;
    // Need to run through the whole bam file, because when we align reads according to exons, we can skip some of its parts
    get_bam_info (bam_reader, bam_general_info);
    bam_reader.Close();


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
//    cout << endl;
//    print_iso_var_map (iso_var_map); // for DEBUG only
//    cout << endl;

//     FOR DEBUG USE ONLY
//        cout << "ANNOTATIONS" << endl;
//        for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
//            cout << "Chromosome: " << chrom_it->first << endl;
//            for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
//                assert (start_it->second.use_count() > 0);
//                cout << " exon: [" << start_it->second->exon_id << "]"
//                     << " isoform: [" << start_it->second->isoform_id << "]"
//                     << " start/stop: [" << start_it->second->start_pose << "," << start_it->second->end_pose << "]"
//                     << " strand: [" << start_it->second->strand << "]"
//                     << " readscount: [" << start_it->second->reads_count << "]";
//                if (start_it->second->previous_gff.use_count() > 0){
//                    cout << " previous gtf: [" << start_it->second->previous_gff.get()->exon_id <<", " << start_it->second->previous_gff.get()->isoform_id << "]";
//                } else {
//                    cout << " previous gtf: [" << "null" << "]";
//                }
//                cout << endl;
//                cout << " Reads: [";
//                for (auto bam_record_it = start_it->second->bam_records.begin(); bam_record_it != start_it->second->bam_records.end(); ++bam_record_it){
//                    assert (bam_record_it->use_count() > 0);
//                    cout << bam_record_it->get()->read_id << " ";
//                }
//                cout << "]" << endl;
//            }
//        }
//        for (auto ext_it = iso_var_map.begin(); ext_it != iso_var_map.end(); ++ext_it){
//            cout << "Chromosome: " << ext_it->first << endl;
//            for (auto int_it = ext_it->second.begin(); int_it != ext_it->second.end(); ++int_it){
//                cout << setw(10) << "  isoform: " << setw(15) << int_it->first << endl;
//                int_it->second.print();
//            }
//        }





    // TODO run thread
    std::map <string, multimap <long, GffRecordPtr> >::iterator start_it = global_annotation_map_ptr.begin();
    std::map <string, multimap <long, GffRecordPtr> >::iterator stop_it = global_annotation_map_ptr.end();
    boost::thread_group process_threads;
    for (int t = 0; t < threads_number; t++){
        cerr << "Adding new thread: " << t << endl;
        process_threads.add_thread(new boost::thread(process, start_it, stop_it, chromosome_info_map, boost::ref(iso_var_map), bam_full_path_name));
    }
    process_threads.join_all();








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
