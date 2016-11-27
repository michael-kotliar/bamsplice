//
// Created by kot4or on 8/8/16.
//

#include "bam_reader.h"

// CLASS BamRecord

using namespace string_tools;


void BamRecord::print (){
    cout << "Read_id: " << read_id << endl;
    cout << "[start, end]: [" << start_pose << ", " << end_pose << "]" << endl;
    cout << "number of parts: " << slices << endl;
};


list <BamRecordPtr> split_to_single_reads (const BamAlignment & current_alignment){
    // Parse CIGAR
    // If splice read - add all of them into array
    // NOTE we need to put it in that list in a right order. use push_back
    list <BamRecordPtr> single_read_array;
    bool strand = ! current_alignment.IsReverseStrand();
    vector<CigarOp> cigar_data = current_alignment.CigarData;
    long start_pose = current_alignment.Position;
    string read_id = current_alignment.Name;
    int slices = 1;
    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'N' || cigar_data[i].Type == 'D'){
            slices++;
        }
    }
    int shift = 0;
    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'M'
//            or
//            cigar_data[i].Type == 'I' or
//            cigar_data[i].Type == 'S' or
//            cigar_data[i].Type == '=' or
//            cigar_data[i].Type == 'X'
             ){
            shift += cigar_data[i].Length;
        } else if (cigar_data[i].Type == 'N' || cigar_data[i].Type == 'D') {
            BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
            single_read_array.push_back(single_read);
            start_pose += shift;
            if (cigar_data[i].Type == 'N'){
                start_pose += cigar_data[i].Length;
            }
            shift = 0;
        }
    }
    BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
    single_read_array.push_back(single_read);
    return single_read_array;
}


bool flag_check (const BamAlignment & al, BamGeneralInfo & bam_general_info){
    // TODO add other checks for pair-end reads
    bam_general_info.total++;
    if(al.IsMapped()) {
        if(al.IsPaired() && (!al.IsProperPair())) {
            bam_general_info.not_aligned++;
            return false;
        }
        if(al.IsPaired() && al.IsProperPair() && al.IsMateMapped() ) {
            // Looks like to exclude chimeric reads
            if( ( (al.Position<al.MatePosition) && al.IsReverseStrand() ) || ( (al.MatePosition < al.Position) && al.IsMateReverseStrand() )) {
                bam_general_info.not_aligned++;
                return false;
            }
        }
    } else {
        bam_general_info.not_aligned++;
        return false;
    }
    return true;
}

bool flag_check (const BamAlignment & al){
    BamGeneralInfo bam_general_info;
    return flag_check (al, bam_general_info);
}


// Gets the new read fro the BAM file through BamReader object
bool get_bam_record (BamReader & bam_reader, BamRecordPtr & bam_record, bool freeze){

    static boost::thread_specific_ptr< list <BamRecordPtr> > saved_reads_tls;
    if( ! saved_reads_tls.get() ) {
        saved_reads_tls.reset( new list <BamRecordPtr> );
    }

    if (freeze and bam_record){
        return true;
    }
    // try to get from the previously stored list (works in case of spliced reads)
    if (not saved_reads_tls->empty()){
        bam_record = saved_reads_tls->front();
        saved_reads_tls->pop_front();
        return true;
    }
    BamAlignment current_alignment;
    if (bam_reader.GetNextAlignment(current_alignment)){
//        cout << "DEBUG: " << current_alignment.Position << " " << current_alignment.Length << endl;
        if (not flag_check (current_alignment)) {
            bam_record.reset();
            return false;
        }
        saved_reads_tls.reset (new list <BamRecordPtr> (split_to_single_reads (current_alignment)) );
        bam_record = saved_reads_tls->front();
        saved_reads_tls->pop_front();
        return true;
    } else {
        bam_record.reset();
        return false;
    }

//    static list <BamRecordPtr> saved_reads; // save all of single reads, which we got from the spliced read
//    if (freeze and bam_record){
//        return true;
//    }
//    // try to get from the previously stored list (works in case of spliced reads)
//    if (not saved_reads.empty()){
//        bam_record = saved_reads.front();
//        saved_reads.pop_front();
//        return true;
//    }
//    BamAlignment current_alignment;
//    if (bam_reader.GetNextAlignment(current_alignment)){
////        cout << "DEBUG: " << current_alignment.Position << " " << current_alignment.Length << endl;
//        if (not flag_check (current_alignment)) {
//            bam_record.reset();
//            return false;
//        }
//        saved_reads = split_to_single_reads (current_alignment);
//        bam_record = saved_reads.front();
//        saved_reads.pop_front();
//        return true;
//    } else {
//        bam_record.reset();
//        return false;
//    }
}


void print_ref_info (const std::map <string, pair <int, int> > & info_map){
    cout << "PRINT REFERENCE DATA FORM BAM FILE" << endl;
    for (auto it = info_map.begin(); it != info_map.end(); ++it){
        cout << "  Chromosome: " << it->first << endl;
        cout << "  RefId: " << it->second.first << endl;
        cout << "  Length: " << it->second.second << endl;
        cout << endl;
    }
}


std::map <string, pair <int, int> > get_chromosome_map_info (const BamReader & reader){
    std::map <string, pair <int, int> > output_map;
    RefVector ref_data_vector = reader.GetReferenceData();
    for (int i = 0; i < ref_data_vector.size(); i++ ){
        RefData current_ref_data = ref_data_vector[i];
        string chrom_name = current_ref_data.RefName;
        int ref_id = reader.GetReferenceID(chrom_name);
        if (ref_id == -1){
            throw ("BUG: it shouldn't be like this");
        }
        int length = current_ref_data.RefLength;

        pair <int, int> internal_pair (ref_id, length);
        pair <string, pair <int, int> > external_pair (chrom_name, internal_pair);
        output_map.insert (external_pair);
    }
    return output_map;
}

// Check if current bam file is indexed (and that index data is loaded into program)
bool make_index (BamReader & bam_reader){
    if (not bam_reader.HasIndex()){
        cerr << "Current BAM file isn't indexed" << endl;
        cerr << "Trying to find index files in the same directory" << endl;
        // Trying to load index data from the filesystem
        // If BamIndex::STANDARD we are looking for BAI file
        if (bam_reader.LocateIndex(BamIndex::STANDARD)){
            cerr << "Located and loaded index file from disk" << endl;
        } else {
            cerr << "Couldn't locate index file" << endl;
            cerr << "Trying to create the new one" << endl;
            // Trying to create index data ourself
            // BamIndex::STANDARD - we are trying to create BAI file
            if (not bam_reader.CreateIndex(BamIndex::STANDARD)){
                cerr << "Cannot create index for current bam file. Exiting" << endl;
                return false;
            };
            cerr << "Index file for current BAM file is succesfully created" << endl;
        }
    }
    return true;
}




bool load_from_file (const string & full_filename, BamGeneralInfo & bam_general_info){
    BamGeneralInfo new_info;
    ifstream input_stream (full_filename);
    if (!input_stream) {
        cerr << "Cannot open file " << full_filename << endl;
        return false;
    }
    string line;
    while (getline(input_stream, line)) {
        if (include_key(line, "Number of input reads")){
            if (not str_to_long(new_info.total, split_line(line)[1]) ){
                return false;
            }
        }
        if (include_key(line, "Uniquely mapped reads number")){
            if (not str_to_long(new_info.aligned, split_line(line)[1]) ){
                return false;
            }
        }
    }
    if (new_info.total == 0 || new_info.aligned == 0){
        return false;
    }
    new_info.not_aligned = new_info.total - new_info.aligned;
    bam_general_info = new_info;
    return true;
}


void get_bam_info(BamReader & bam_reader, BamGeneralInfo & bam_general_info){
    // check if we can get it from the STAR output Log.final.out in the same folder.
    // If not - iterate over the file
    string log_filename_sufix = "Log.final.out";   // TODO put it as parameter
    string bam_filename (bam_reader.GetFilename(), bam_reader.GetFilename().find_last_of('/')+1);
    string bam_wo_ext (bam_filename, 0, bam_filename.find_last_of('.')+1);
    string path (bam_reader.GetFilename(), 0,  bam_reader.GetFilename().find_last_of('/')+1);
    string full_log_filename = path + bam_wo_ext + log_filename_sufix;
    if (not load_from_file (full_log_filename, bam_general_info)){
        cerr << "Couldn't find any file with mapping statistics. Calculating ..." << endl;
        BamAlignment al;
        while (bam_reader.GetNextAlignment(al)){
            flag_check (al, bam_general_info);
        }
        bam_general_info.aligned = bam_general_info.total - bam_general_info.not_aligned;
        bam_reader.Rewind();
    }
    cerr << "BAM alignment statistics:" << endl;
    cerr << "   Total: " << bam_general_info.total << endl;
    cerr << "   Aligned: " << bam_general_info.aligned << endl;
    cerr << "   Not aligned: " << bam_general_info.not_aligned << endl;
}