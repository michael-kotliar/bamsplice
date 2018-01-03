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

//static boost::thread_specific_ptr< list <BamRecordPtr> > saved_reads_tls;
//list <BamRecordPtr> split_to_single_reads (const BamAlignment & current_alignment, int min_read_segment_length){
//    list <BamRecordPtr> single_read_array;
//    bool strand = ! current_alignment.IsReverseStrand();
//    vector<CigarOp> cigar_data = current_alignment.CigarData;
//    long start_pose = current_alignment.Position;
//    string read_id = current_alignment.Name;
//    int slices = 1;
//    for (int i = 0; i < cigar_data.size(); i++){
//        if (cigar_data[i].Type == 'N' || cigar_data[i].Type == 'D'){
//            slices++;
//        }
//    }
//    int shift = 0;
//    for (int i = 0; i < cigar_data.size(); i++){
//        if (cigar_data[i].Type == 'M'){
//            shift += cigar_data[i].Length;
//        } else if (cigar_data[i].Type == 'N' || cigar_data[i].Type == 'D') {
//            BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
//            single_read_array.push_back(single_read);
//            start_pose += shift;
//            if (cigar_data[i].Type == 'N'){
//                start_pose += cigar_data[i].Length;
//            }
//            shift = 0;
//        }
//    }
//    BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
//    single_read_array.push_back(single_read);
//    return single_read_array;
//}




//
list <BamRecordPtr> split_to_single_reads (const BamAlignment & current_alignment, int min_read_segment_length){

    list <BamRecordPtr> single_read_array;

    bool strand = ! current_alignment.IsReverseStrand();
    vector<CigarOp> cigar_data = current_alignment.CigarData;
    long start_pose = current_alignment.Position;
    string read_id = current_alignment.Name;

    int slices = 1;
    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'N'){ // Spliced reads are separated only by N
            slices++;
        }
    }

    int shift = 0;
    int original_slices = slices;

    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'M' or cigar_data[i].Type == 'D'){
            assert (cigar_data[i].Length > 0);
            shift += cigar_data[i].Length;
        } else if (cigar_data[i].Type == 'N') {
            if (shift > min_read_segment_length){
                BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
                single_read_array.push_back(single_read);
            } else {
                slices--;
            }
            start_pose += shift;
            start_pose += cigar_data[i].Length;
            shift = 0;
        }
    }

    // in a case when is wasn't spliced read
    BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices, strand));
    single_read_array.push_back(single_read);

    // updated number of reads if some of the segments were excluded
    if (slices < original_slices){
        for (auto iter = single_read_array.begin(); iter != single_read_array.end(); ++iter){
            iter->get()->slices = slices;
        }
    }

    return single_read_array;
}

bool flag_check (BamAlignment & al, bool dUTP){
    if (!al.IsMapped() || !al.IsPrimaryAlignment() || al.IsDuplicate()){
        return false;
    }
    if (al.IsPaired() && (!al.IsProperPair() || !al.IsMateMapped()) ){
        return false;
    }
  // TODO Think if we really should do like this?
  // TODO check if it switch it back when is called twice
    if (dUTP && al.IsPaired() && al.IsSecondMate()){
        if(al.IsReverseStrand()) {
            al.SetIsReverseStrand(false);
        } else {
            al.SetIsReverseStrand(true);
        }
    }
    return true;
}


void put_bam_record_back (BamRecordPtr bam_record){
    if( ! saved_reads_tls.get() ) {
        saved_reads_tls.reset( new list <BamRecordPtr> );
    }
    saved_reads_tls->push_back(bam_record);
    cout << "Pushed back read: " << bam_record->read_id << endl;
}

void reset_saved_reads (){
    cout << "reset_saved_reads" << endl;
    saved_reads_tls.reset( new list <BamRecordPtr> );
}

// Gets the new read from the BAM file through BamReader object
bool get_bam_record (BamReader & bam_reader, BamRecordPtr & bam_record, int min_read_segment_length, bool dUTP, bool freeze){

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
    while (bam_reader.GetNextAlignment(current_alignment)){
        if (not flag_check (current_alignment, dUTP)) {
            continue;
        }
        saved_reads_tls.reset (new list <BamRecordPtr> (split_to_single_reads (current_alignment, min_read_segment_length)) );
        bam_record = saved_reads_tls->front();
        saved_reads_tls->pop_front();
        return true;
    }

    // reached the end of file
    bam_record.reset();
    return false;

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


std::map <string, pair <int, int> > get_chromosome_map_info (const BamReader & reader, const vector<string> & exclude_chr){
    std::map <string, pair <int, int> > output_map;
    RefVector ref_data_vector = reader.GetReferenceData();
    for (int i = 0; i < ref_data_vector.size(); i++ ){
        RefData current_ref_data = ref_data_vector[i];
        string chrom_name = current_ref_data.RefName;
        if (std::find(exclude_chr.begin(), exclude_chr.end(), boost::to_lower_copy(chrom_name)) != exclude_chr.end()){
            cerr << "Skip excluded chromosome [" << chrom_name << "]" << endl;
            continue;
        }

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



//
//bool load_from_file (const string & full_filename, BamGeneralInfo & bam_general_info){
//    BamGeneralInfo new_info;
//    ifstream input_stream (full_filename);
//    if (!input_stream) {
//        cerr << "Cannot open file " << full_filename << endl;
//        return false;
//    }
//    string line;
//    while (getline(input_stream, line)) {
//        if (include_key(line, "Number of input reads")){
//            if (not str_to_long(new_info.total, split_line(line)[1]) ){
//                return false;
//            }
//        }
//        if (include_key(line, "Uniquely mapped reads number")){
//            if (not str_to_long(new_info.aligned, split_line(line)[1]) ){
//                return false;
//            }
//        }
//    }
//    if (new_info.total == 0 || new_info.aligned == 0){
//        return false;
//    }
//    new_info.not_aligned = new_info.total - new_info.aligned;
//    bam_general_info = new_info;
//    return true;
//}


void get_bam_info(BamReader & bam_reader, BamGeneralInfo & bam_general_info, const std::map <string, pair <int, int> > & chromosome_info_map){
    // Refactore chromosome_info_map to allow easy search by RefIf
    std::map <int, string> allowed_chr_ids; // < RefID, chromosome name >
    for (auto it = chromosome_info_map.begin(); it != chromosome_info_map.end(); ++it){
        allowed_chr_ids[it->second.first] = it->first;
    }

    cerr << "Calculating mapping statistics" << endl;
    BamAlignment al;
    while (bam_reader.GetNextAlignment(al)){
        bam_general_info.total++;
        if ( allowed_chr_ids.count(al.RefID) != 1){
            bam_general_info.excluded++;
            continue;
        }
        if (not flag_check (al, false)){
            bam_general_info.not_aligned++;
        };
    }
    bam_general_info.aligned = bam_general_info.total - bam_general_info.not_aligned - bam_general_info.excluded;
    bam_reader.Rewind();

    cerr << "BAM alignment statistics:" << endl;
    cerr << "   Total: " << bam_general_info.total << endl;
    cerr << "   Excluded: " << bam_general_info.excluded << endl;
    cerr << "   Aligned: " << bam_general_info.aligned << endl;
    cerr << "   Not aligned: " << bam_general_info.not_aligned << endl;
}