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


bool verify_params (cxxopts::Options params){
    // check mandatorary parameters
    if (!params.count("bam") || !params.count("annotation") || !params.count("output")){
        cerr << "You should set --bam, --annotation, --output parameters" << endl;
        return false;
    }

    try { // It will never fail, because we already checked it on the parameters parsing step
        cerr << "Parameters: " << endl;
        cerr << "  bam: " << params["bam"].as<std::string>() << endl;
        cerr << "  annotation: " << params["annotation"].as<std::string>() << endl;
        cerr << "  log: " << params["log"].as<std::string>() << endl;
        cerr << "  output: " << params["output"].as<std::string>() << endl;
        cerr << "  threshold: " << params["threshold"].as<double>() << endl;
        cerr << "  cutoff: " << params["cutoff"].as<double>() << endl;
        cerr << "  minIntLen: " << params["minIntLen"].as<int>() << endl;
        cerr << "  minReadLen: " << params["minReadLen"].as<int>() << endl;
        cerr << "  exclude: " << params["exclude"].as<std::string>() << endl;
        cerr << "  threads: " << params["threads"].as<int>() << endl;
        cerr << "  keepUnique: " << params["keepUnique"].as<bool>() << endl;
        cerr << "  dutp: " << params["dutp"].as<bool>() << endl;
    } catch (...){
        cerr << "Some of the parameters have a mistake" << endl;
        return false;
    }
    if (params["threshold"].as<double>() < 0){
        cerr << "threshold cannot be less than 0" << endl;
        return false;
    }
    if (params["cutoff"].as<double>() < 0){
        cerr << "cutoff cannot be less than 0" << endl;
        return false;
    }
    if (params["minIntLen"].as<int>() < 0){
        cerr << "minIntLen cannot be less than 0" << endl;
        return false;
    }
    if (params["minReadLen"].as<int>() < 0){
        cerr << "minReadLen cannot be less than 0" << endl;
        return false;
    }
    if (params["threads"].as<int>() < 1){
        cerr << "threads cannot be less than 1" << endl;
        return false;
    }
    return true;
}



int main(int argc, char **argv) {
    cxxopts::Options params("GEEP", "Gene Expression Evaluation Processor");

    params.add_options()
            ("b,bam", "path to the BAM file", cxxopts::value<std::string>(),
             "Set the path to the BAM file")
            ("g,annotation", "path to the GTF or TAB-delimited annotation file", cxxopts::value<std::string>(),
             "Set the path to the GTF or TAB-delimited file")
            ("l,log", "path to the LOG file", cxxopts::value<std::string>()->default_value("/dev/null"),
             "Set the path to save LOG file")
            ("o,output", "prefix for the output files", cxxopts::value<std::string>(),
             "Set the prefix to save output files")
            ("t,threshold", "RPKM threshold", cxxopts::value<double>()->default_value("0"),
             "Set rpkm cutoff threshold, below which everything will be changed to value set with --cutoff")
            ("c,cutoff", "RPKM cutoff value", cxxopts::value<double>()->default_value("0"),
             "Set rpkm cutoff value to be used for all rpkms below --threshold")
            ("x,exclude", "Ignore chromosome list", cxxopts::value<std::string>()->default_value(""),
             "Coma separated list of chromosomes to be ignored")
            ("i,minIntLen", "minimal interval length", cxxopts::value<int>()->default_value("0"),
             "Set the minimal interval length. All shorter intervals will be discarded")
            ("r,minReadLen", "minimal read length", cxxopts::value<int>()->default_value("0"),
             "Set the minimal read length. All parts of spliced reads that intersect with exon in less than minReadLen nucleotides will be discarded")
            ("p,threads", "number of threads", cxxopts::value<int>()->default_value("1"),
             "Set the number of threads")
            ("u,keepUnique", "flag to fix unique reads for the specific isoromf interval", cxxopts::value<bool>(),
             "Set this flag if you want prevent distributing the isoform unique reads among other isoforms")
            ("d,dutp", "set dutp enabled", cxxopts::value<bool>(),
             "Set this dutp flag if strand specific analysys should be made");

    try {
        params.parse(argc, argv);
    } catch (...){
        cerr << "Parameter parsing error" << endl;
        return 0;
    }

    if (!verify_params (params)){
        return 0;
    }

    // Get the log file path or default /dev/null
    string log_filename = params["log"].as<std::string>();
    freopen(log_filename.c_str(), "a", stdout); // TODO Check what happens when filename is not correct
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         << now->tm_mday << "   "
         << now->tm_hour << ":" << now->tm_min
         << endl << endl;

    // read from BAM file
    BamReader bam_reader;
    if (not bam_reader.Open(params["bam"].as<std::string>())) {
        cerr << "Couldn't open file " << params["bam"].as<std::string>() << endl;
        return 0;
    } else cerr << "Open " << bam_reader.GetFilename() << endl;

    // key - chromosome name, value - <RefId, Length> for corresponding chromosome from the BAM file
    // chromosome_info_map - saves correspondence between chromosome name and RefId from the BamReader object
    std::map <string, pair <int, int> > chromosome_info_map = get_chromosome_map_info (bam_reader);

//        print_ref_info (chromosome_info_map); // Only for DEBUG

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
    if (not load_annotation (params["annotation"].as<std::string>(), params["exclude"].as<std::string>(), global_annotation_map_ptr, iso_var_map)){
        return 0;
    }

//    cout << endl;
//    print_iso_var_map (iso_var_map); // for DEBUG only
//    cout << endl;

//     FOR DEBUG USE ONLY
//        cout << endl << "ANNOTATIONS" << endl;
//        for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
//            cout << "Chromosome: " << chrom_it->first << endl;
//            for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
//                assert (start_it->second.use_count() > 0);
//                cout << " isoform: [" << start_it->second->isoform_id << "]"
//                     << " exon: [" << start_it->second->exon_id << "]"
//                     << " start/stop: [" << start_it->second->start_pose << "," << start_it->second->end_pose << "]"
//                     << " strand: [" << start_it->second->strand << "]"
//                     << " readscount: [" << start_it->second->reads_count << "]"
//                     << " start_exon: [" << start_it->second->start_exon << "]"
//                     << " stop_exon: [" << start_it->second->stop_exon << "]" << endl;
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

    // Get intersection of chrom names in bam and annotaion file, save iterators from global_annotation_map_ptr
    vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > intersection_array;
    for (auto bam_it = chromosome_info_map.begin(); bam_it != chromosome_info_map.end(); ++bam_it){
        for (auto gtf_it = global_annotation_map_ptr.begin(); gtf_it != global_annotation_map_ptr.end(); ++gtf_it){
            if (bam_it->first == gtf_it->first){
                intersection_array.push_back(gtf_it);
            }
        }
    }

    int in_each_thread = (int) floor ((double)intersection_array.size() / fmin(params["threads"].as<int>(), intersection_array.size()));
    cerr << "on each thread: " << in_each_thread << endl;

    boost::thread_group process_threads;

    for (int t = 0; t < fmin (params["threads"].as<int>(), intersection_array.size()); t++){
//        cerr << "Adding new thread" << endl;
//
//        vector < std::map <string, multimap <long, GffRecordPtr> >::iterator >::iterator start_subvector = intersection_array.begin() + t * in_each_thread;
//        vector < std::map <string, multimap <long, GffRecordPtr> >::iterator >::iterator stop_subvector = start_subvector;
//        while ( stop_subvector != intersection_array.begin() + (t+1)*in_each_thread && stop_subvector != intersection_array.end()){
//            ++stop_subvector;
//        }
//        if (t == fmin (threads_number, intersection_array.size()) - 1 ){
//            stop_subvector = intersection_array.end();
//        }
//        vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > chrom_vector;
//        chrom_vector = vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > (start_subvector, stop_subvector);
//        process_threads.add_thread(new boost::thread(process, chrom_vector, chromosome_info_map, boost::ref(iso_var_map), bam_full_path_name, t));

            cerr << "Adding new thread: " << t <<  endl;
            int start_subvector = t * in_each_thread;
            int stop_subvector = start_subvector;
            while ( stop_subvector != (t+1)*in_each_thread && stop_subvector != intersection_array.size()){
                ++stop_subvector;
            }
            if (t == fmin (params["threads"].as<int>(), intersection_array.size()) - 1 ){
                stop_subvector = intersection_array.size();
            }
            vector < std::map <string, multimap <long, GffRecordPtr> >::iterator > chrom_vector;
            for (int j = start_subvector; j < stop_subvector; j++){
                chrom_vector.push_back(intersection_array[j]);
            }
            for (int j = 0; j < chrom_vector.size(); j++){
                cerr << start_subvector + j << ". " << chrom_vector[j]->first << endl;
            }
            process_threads.add_thread(new boost::thread(process,
                                                         chrom_vector,
                                                         chromosome_info_map,
                                                         boost::ref(iso_var_map),
                                                         t,
                                                         params) );
    }
    process_threads.join_all();


//            vector < int > chromosome_info_map_1;
//            chromosome_info_map_1.push_back(1);
//            chromosome_info_map_1.push_back(2);
//            chromosome_info_map_1.push_back(3);
//            chromosome_info_map_1.push_back(5);
//            chromosome_info_map_1.push_back(6);
//            vector < int > global_annotation_map_ptr_1;
//            global_annotation_map_ptr_1.push_back(1);
//            global_annotation_map_ptr_1.push_back(2);
//            global_annotation_map_ptr_1.push_back(3);
//            vector < int > intersection_array_1;
//            for (auto bam_it = chromosome_info_map_1.begin(); bam_it != chromosome_info_map_1.end(); ++bam_it){
//                for (auto gtf_it = global_annotation_map_ptr_1.begin(); gtf_it != global_annotation_map_ptr_1.end(); ++gtf_it){
//                    if (*bam_it == *gtf_it){
//                        intersection_array_1.push_back(*gtf_it);
//                    }
//                }
//            }
//
//            for (int i = 0; i < intersection_array_1.size(); i++)
//                cerr << i << ". " << intersection_array_1[i] << endl;
//
//
//            threads_number = 8;
//            cerr << "threads_number: " << threads_number << endl;

//            int in_each_thread = (int) floor ((double)intersection_array_1.size() / (double)fmin(threads_number, intersection_array_1.size()));
//            cerr << "in_each_thread: " << in_each_thread << endl;
//
//
//            for (int t = 0; t < fmin (threads_number, intersection_array_1.size()); t++){
//                vector < int > chrom_vector_1;
//                vector < int >::iterator start = intersection_array_1.begin() + t * in_each_thread;
//                vector < int >::iterator stop = start;
//                while ( stop != intersection_array_1.begin() + (t+1)*in_each_thread && stop != intersection_array_1.end())
//                    ++stop;
//                if (t == fmin (threads_number, intersection_array_1.size())-1 ){
//                    cerr << "last thread" << endl;
//                    stop = intersection_array_1.end();
//                }
//
//                chrom_vector_1 = vector < int > (start, stop);
//                cerr << "chrom_vector.size(): " << chrom_vector_1.size() << endl;
//                    for (int k = 0; k < chrom_vector_1.size(); k++)
//                        cerr << "   " << k << ". " << chrom_vector_1[k] << endl;
//            }


    cerr << "Calculate rpkm" << endl;

    calculate_rpkm (iso_var_map, bam_general_info.aligned);
    adjust_threshold (iso_var_map, params["threshold"].as<double>(), params["cutoff"].as<double>(), bam_general_info.aligned);
    print_iso_var_map (iso_var_map);


    cerr << "Exporting results" << endl;
    print_iso_var_map_to_file (iso_var_map, params["output"].as<std::string>()+"raw.txt");

    export_isoform_group (iso_var_map, params["output"].as<std::string>()+"isoforms.csv");
    export_gene_group (iso_var_map, params["output"].as<std::string>()+"genes.csv");


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
