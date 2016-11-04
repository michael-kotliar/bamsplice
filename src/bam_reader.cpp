//
// Created by kot4or on 8/8/16.
//

#include "bam_reader.h"

// CLASS BamRecord

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
    vector<CigarOp> cigar_data = current_alignment.CigarData;
    long start_pose = current_alignment.Position;
    string read_id = current_alignment.Name;
    int slices = 1;
    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'N'){
            slices++;
        }
    }
    int shift = 0;
    for (int i = 0; i < cigar_data.size(); i++){
        if (cigar_data[i].Type == 'M' or
            cigar_data[i].Type == 'I' or
            cigar_data[i].Type == 'S' or
            cigar_data[i].Type == '=' or
            cigar_data[i].Type == 'X'){
            shift += cigar_data[i].Length;
        } else {
            BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices));
            single_read_array.push_back(single_read);
            start_pose += shift;
            if (cigar_data[i].Type == 'N'){
                start_pose += cigar_data[i].Length;
            }
            shift = 0;
        }
    }
    BamRecordPtr single_read (new BamRecord (start_pose, start_pose + shift, read_id, slices));
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
    static list <BamRecordPtr> saved_reads; // save all of single reads, which we got from the spliced read
    if (freeze and bam_record){
        return true;
    }
    // try to get from the previously stored list (works in case of spliced reads)
    if (not saved_reads.empty()){
        bam_record = saved_reads.front();
        saved_reads.pop_front();
        return true;
    }
    BamAlignment current_alignment;
    if (bam_reader.GetNextAlignment(current_alignment)){
        if (not flag_check (current_alignment)) {
            bam_record.reset();
            return false;
        }
        saved_reads = split_to_single_reads (current_alignment);
        bam_record = saved_reads.front();
        saved_reads.pop_front();
        return true;
    } else {
        bam_record.reset();
        return false;
    }
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
        cout << "Current BAM file isn't indexed" << endl;
        cout << "Trying to find index files in the same directory" << endl;
        // Trying to load index data from the filesystem
        // If BamIndex::STANDARD we are looking for BAI file
        if (bam_reader.LocateIndex(BamIndex::STANDARD)){
            cout << "Located and loaded index file from disk" << endl;
        } else {
            cout << "Couldn't locate index file" << endl;
            cout << "Trying to create the new one" << endl;
            // Trying to create index data ourself
            // BamIndex::STANDARD - we are trying to create BAI file
            if (not bam_reader.CreateIndex(BamIndex::STANDARD)){
                cout << "Cannot create index for current bam file. Exiting" << endl;
                return false;
            };
            cout << "Index file for current BAM file is succesfully created" << endl;
        }
    }
    return true;
}

void get_bam_info(BamReader & bam_reader, BamGeneralInfo & bam_general_info){
    BamAlignment al;
    while (bam_reader.GetNextAlignment(al)){
        flag_check (al, bam_general_info);
    }
    // return read pointer to the beginnig of the file
    bam_reader.Rewind();
}