#include <iostream>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>
#include <map>
#include <forward_list>
#include <boost/icl/interval_map.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "include/api/BamReader.h"
#include "include/api/BamMultiReader.h"

using namespace std;
using namespace boost::icl;
using namespace BamTools;

class BamRecord;
class GffRecord;
class MapElement;

typedef boost::shared_ptr<BamRecord> BamRecordPtr;
typedef boost::shared_ptr<GffRecord> GffRecordPtr;

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
            , read_id ("")
            , slices (1)
    {}

    void print (){
        cout << "Read_id: " << read_id << endl;
        cout << "[start, end]: [" << start_pose << ", " << end_pose << "]" << endl;
        cout << "number of parts: " << slices << endl;
    }
};

class GffRecord {
public:
    long start_pose; // start position of annotation
    long end_pose; // stop position of annotation
    string exon_id; // text-identificator of current annotation. Looks like we will use it only for debug
    string isoform_id; // set the name of the isofom=rm to which current annotation belongs
    vector < BamRecordPtr > bam_records; // array of pointers to all of the reads, which belongs to this annotation. Need this for debug, than we can delete this field
    long reads_count; // total number of reads, which belongs to this annotation
    GffRecordPtr previous_gff; // ptr to the previous annotation in the same isoform. In NULL - first annotation in current isiform

    // CONSTRUCTOR WITH PARAMETERS
    GffRecord (long start, long end, string exon, string isoform, GffRecordPtr pre_gff)
            : start_pose (start)
            , end_pose (end)
            , exon_id (exon)
            , isoform_id (isoform)
            , previous_gff (pre_gff)
            , reads_count (0)
    {}

    // EMPTY CONSTRUCTOR
    GffRecord ()
            : start_pose (0)
            , end_pose (0)
            , exon_id ("")
            , isoform_id ("")
            , previous_gff (NULL)
            , reads_count (0)
    {}

};

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



list <BamRecordPtr> saved_reads; // save all of single reads, which we got from the spliced read
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

// This function should be upfated to the similar one when we load the real bam/sam file
bool get_bam_record (BamReader & bam_reader, BamRecordPtr & bam_record, bool freeze = false){
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
        saved_reads = split_to_single_reads (current_alignment);
        bam_record = saved_reads.front();
        saved_reads.pop_front();
        return true;
    } else {
        bam_record.reset();
        return false;
    }
}

// Updates iterator current_gtf_records_splitted_it for the interval map segment which includes the starting position of the current_bam_read
// If fails - returns false
// If we need to skip bam_record, we set freeze = false to get the new current_bam_record
// If we need to skip interval map segment we set freeze to true to avoid changing current_bam_record, when iterating over the loop,
// that get new value from bam records array
// If skip interval map segment - increment iterator current_gtf_records_splitted_it
bool allow_skip_rest = false;
bool find_start_segment_annotation (BamRecordPtr current_bam_record, BamRecord previous_bam_record, interval_map<long, MapElement>::iterator & current_gtf_records_splitted_it, bool & freeze){
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
}

// Updates iterator temp_gtf_records_splitted_it for interval map segment which includes stop position of current_bam_record
// If fail - returns false
bool find_stop_segment_annotation (BamRecordPtr current_bam_record,
                                   interval_map<long,MapElement>::iterator & temp_gtf_records_splitted_it,
                                   interval_map<long,MapElement>::iterator max_segment_annotation,
                                   bool & freeze){
    while (current_bam_record->end_pose < temp_gtf_records_splitted_it->first.lower() or
           current_bam_record->end_pose > temp_gtf_records_splitted_it->first.upper()){
        // check if we reached the end of the gtf_records array. If true, change go to the next bam_record
        if (temp_gtf_records_splitted_it == max_segment_annotation){
            freeze = false;
            return false;
        }
        temp_gtf_records_splitted_it++;
    }
    return true;
}

// FOR DEBUG ONLY
void print_segment_annotation (const string & title, interval_map<long, MapElement>::iterator current_gtf_records_splitted_it){
    cout << title << " " << "[" << current_gtf_records_splitted_it->first.lower() << "," << current_gtf_records_splitted_it->first.upper() << "] :";
    for (auto start_segment_annotation_it = current_gtf_records_splitted_it->second.gtf_records.begin();
         start_segment_annotation_it != current_gtf_records_splitted_it->second.gtf_records.end(); ++start_segment_annotation_it){
        cout << " " << (*start_segment_annotation_it)->exon_id << " (" << (*start_segment_annotation_it)->isoform_id << ")";
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
                    cout << " " << (*intersection_segment_annotation_it)->exon_id << " ("<< (*intersection_segment_annotation_it)->isoform_id << "), ";
                }
                cout << endl;
    return intersection;
}


// Check if current read is a part of the spliced read. If false - returns true and doesn't check any additional conditions
// If true, check if current annotation fits current read with following conditions:
// 1. If read is the first part of big spliced read and its end position is equal to end position of current annotation - return true;
// 2. If read is the last part of big spliced read and its start position is equal to start position of current annotation - return true;
// 3. If read is the middle part of big spliced read and both its start and end position are equal to start and end position of current annotation correspondingly - return true
// In all other cases return false
// NOTE temp_set has only one record. We can change it for GffRecordPtr and it send to the function just temp_set.begin()
bool fit_spliced_read_condition(const long & current_slice, const long & slice_number, const set <GffRecordPtr> & temp_set, BamRecordPtr current_bam_record){
    if (slice_number == 1) return true; // added just in case to make sure that we are not trying ot check unspliced read
    if (current_slice == 1 and (*temp_set.begin())->end_pose != current_bam_record->end_pose ) {
        return false;
    } else
        if (current_slice == slice_number and (*temp_set.begin())->start_pose != current_bam_record->start_pose){
            return false;
        } else
            if (current_slice > 1 and
                current_slice < slice_number and
                (*temp_set.begin())->end_pose != current_bam_record->end_pose and
                (*temp_set.begin())->start_pose != current_bam_record->start_pose){
                return false;
            }
    return true;
}





namespace string_tools {

    void remove_spaces(vector<string> &array) {
        vector<string> array_filtered;
        for (int i = 0; i < array.size(); i++)
            if (array[i].size() > 0 and array[i] != "\t") {
                array_filtered.push_back(array[i]);
            }
        array = array_filtered;
    }

    vector<string> split_line(const string &line, const string & key = "\t") {
        vector<string> line_splitted;
        boost::split(line_splitted, line, boost::is_any_of(key));
        remove_spaces(line_splitted);
        return line_splitted;
    }

    bool include_key(const string &line, const string &key_sequence) {
        vector<string> key_splitted = split_line(key_sequence);
        for (int i = 0; i < key_splitted.size(); i++) {
            if (line.find(key_splitted[i]) != std::string::npos) {
                return true;
            }
        }
        return false;
    }

}


enum cds_stat {
    none,
    unk,
    incmpl,
    cmpl
};

bool str_to_cds_stat(const string &value, cds_stat &result){
    if (value == "none") {
        result = none;
    } else
    if (value == "unk"){
        result = unk;
    } else
    if (value == "incmpl"){
        result = incmpl;
    } else
    if (value == "cmpl"){
        result = cmpl;
    } else {
        cout << "Cannot evalueate " << value << endl;
        return false;
    }
    return true;
}

bool str_to_long_ptr(boost::shared_ptr<long> &ptr, const string &value){
    long temp;
    try {
        temp = boost::lexical_cast<long>(value);
    }
    catch(...){
        cout << "Bad lexical cast of " << value << " as long" << endl;
        return false;
    }
    ptr.reset (new long (temp));
    return true;
}

bool str_to_long(long &var, const string &value){
    long temp;
    try {
        temp = boost::lexical_cast<long>(value);
    }
    catch(...){
        cout << "Bad lexical cast of " << value << " as long" << endl;
        return false;
    }
    var = temp;
    return true;
}

bool str_to_int_ptr(boost::shared_ptr<int> &ptr, const string &value){
    int temp;
    try {
        temp = boost::lexical_cast<int>(value);
    }
    catch(...){
        cout << "Bad lexical cast of " << value << " as int" << endl;
        return false;
    }
    ptr.reset (new int (temp));
    return true;
}

bool str_to_int(int &var, const string &value){
    int temp;
    try {
        temp = boost::lexical_cast<int>(value);
    }
    catch(...){
        cout << "Bad lexical cast of " << value << " as int" << endl;
        return false;
    }
    var = temp;
    return true;
}

bool str_array_to_long_ptr_array(const vector<string> &input, vector<boost::shared_ptr<long> > &output){
    for (int i = 0; i < input.size(); i++){
        boost::shared_ptr <long> temp;
        if (not str_to_long_ptr(temp, input[i])){
            return false;
        }
        output.push_back(temp);
    }
    return true;
}


bool str_array_to_long_array(const vector<string> &input, vector<long> &output){
    for (int i = 0; i < input.size(); i++){
        long temp;
        try {
            temp = boost::lexical_cast<long>(input[i]);
        }
        catch(...){
            cout << "Bad lexical cast of " << input[i] << " as long" << endl;
            return false;
        }
        output.push_back(temp);
    }
    return true;
}

void print_vector (const vector <string> & in, string title = ""){
    cout << title << endl;
    for (int i = 0; i < in.size(); i++){
        cout << i <<") " << in[i] << endl;
    }
}

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

    vector <long> exon_starts; // not necessary to be sorted // 9
    vector <long> exon_ends; // saves pointers to ends of the exons. not necessary to be sorted // 10
    vector <long> exon_frames; // pointers to exon frames, not necessary to be sorted // 15

    int score; //11
    string name2; // 12
    cds_stat cds_start_stat; // 13
    cds_stat cds_end_stat; // 14

    void print (){
        cout << "isoform: " << name << endl;
        cout << "chrom: " << chrom << endl;
        cout << "strand: " << strand << endl;
        cout << "[tx_start, tx_end]: [" << tx_start << ", " << tx_end << "]" << endl;
        cout << "[cds_start, cds_end]: [" << cds_start << ", " << cds_end << "]" << endl;
        cout << "exon_count: " << exon_count << endl;
        cout << "exon_starts.size(): " << exon_starts.size() << endl;
        cout << "exon_ends.size(): " << exon_ends.size() << endl;
        cout << "exon_frames.size(): " << exon_frames.size() << endl;
        cout << "Exons:" << endl;
        for (int i = 0; i < exon_count; i++){
            cout << "  " << i << ") " << "[" << exon_starts[i] << ", "<< exon_ends[i] << "] - " << exon_frames[i] << endl;
        }
        cout << "score: " << score << endl;
        cout << "name2: " << name2 << endl;
        cout << "cds_start_stat: " << cds_start_stat << endl;
        cout << "cds_end_stat: " << cds_end_stat << endl;
    }

    Isoform (string line){
        vector<string> line_splitted = string_tools::split_line(line);

        // BIN
        if (not str_to_int(bin, line_splitted[0])){
            throw ("Isoform constructor error");
        }

        // NAME
        name = line_splitted[1];

        // CHROM
        chrom = line_splitted[2];

        // STRAND
        strand = (line_splitted[3] == "+") ? true : false;

        // TX_START
        if (not str_to_long(tx_start, line_splitted[4])){
            throw ("Isoform constructor error");
        }

        // TX_END
        if (not str_to_long(tx_end, line_splitted[5])){
            throw ("Isoform constructor error");
        }

        // CDS_START
        if (not str_to_long(cds_start, line_splitted[6])){
            throw ("Isoform constructor error");
        }

        //CDS_END
        if (not str_to_long(cds_end, line_splitted[7])){
            throw ("Isoform constructor error");
        }

        // EXON_COUNT
        if (not str_to_int(exon_count, line_splitted[8])){
            throw ("Isoform constructor error");
        }

        //SCORE
        if (not str_to_int(score, line_splitted[11])){
            throw ("Isoform constructor error");
        }

        // NAME2
        name2 = line_splitted[12];

        // CDS_START_STAT
        if (not str_to_cds_stat(line_splitted[13], cds_start_stat)){
            throw ("Cannot set the value for cds_start_stat");
        }

        // CDS_END_STAT
        if (not str_to_cds_stat(line_splitted[14], cds_end_stat)){
            throw ("Cannot set the value for cds_end_stat");
        }

        vector<string> exon_starts_str = string_tools::split_line(line_splitted[9], ",");
        vector<string> exon_ends_str = string_tools::split_line(line_splitted[10], ",");
        vector<string> exon_frames_str = string_tools::split_line(line_splitted[15], ",");

//        print_vector (exon_starts_str, "exon_starts_str");
//        print_vector(exon_ends_str, "exon_ends_str");
//        print_vector(exon_frames_str, "exon_frames_str");


        if (not str_array_to_long_array(exon_starts_str, exon_starts)){
            throw ("Isoform class constructor fail");
        };

        if (not str_array_to_long_array(exon_ends_str, exon_ends)){
            throw ("Isoform class constructor fail");
        };

        if (not str_array_to_long_array(exon_frames_str, exon_frames)){
            throw ("Isoform class constructor fail");
        };

    }
};

void print_iso_var_map (const std::map <string, std::map <string, int> > & iso_var_map){
    for (auto ext_it = iso_var_map.begin(); ext_it != iso_var_map.end(); ++ext_it){
        cout << "Chromosome: " << ext_it->first << endl;
        for (auto int_it = ext_it->second.begin(); int_it != ext_it->second.end(); ++int_it){
            cout << "  isoform: " << int_it->first << " index: " << int_it->second << endl;
        }
    }
}


// global_annotation_map_ptr : key - chromosome name, value - multimap of annotations, sorted by not-unique key - start pose of annotation
// NOTE : forward list of annotations should be sorted by start pose with rule a<b
bool load_annotation (const string & full_path_name, std::map <string, multimap <long, GffRecordPtr> > & global_annotation_map_ptr, std::map <string, std::map <string, int> > & iso_var_map){
    ifstream input_stream (full_path_name);
    if (!input_stream) {
        cout << "Cannot open file " << full_path_name << endl;
        return false;
    }
    string line;
    while (getline(input_stream, line)) {
        if (string_tools::include_key(line, "name")) {
            continue;
        }

        Isoform current_isoform(line);

        pair <string, int> internal_pair_for_iso_var_map (current_isoform.name, iso_var_map[current_isoform.chrom].size());
        std::map <string, int> internal_iso_var_map;
        internal_iso_var_map.insert(internal_pair_for_iso_var_map);
        pair <std::map <string, std::map <string, int> >::iterator, bool> res;
        pair <string, std::map <string, int> > external_pair_for_iso_var_map (current_isoform.chrom, internal_iso_var_map);
        res = iso_var_map.insert (external_pair_for_iso_var_map);
        if (res.second == false){
            res.first->second.insert(internal_pair_for_iso_var_map);
        }

        GffRecordPtr previous_annotation;
        previous_annotation.reset();

        for (int i = 0; i < current_isoform.exon_count; i++){
            stringstream ss;
            ss << i+1;
            string exon_id = ss.str();

            GffRecordPtr current_gff (new GffRecord (current_isoform.exon_starts[i], current_isoform.exon_ends[i], exon_id, current_isoform.name, previous_annotation) );
            previous_annotation = current_gff;
            pair <long, GffRecordPtr> internal_pair (current_isoform.exon_starts[i], current_gff);
            multimap <long, GffRecordPtr> internal_multimap;
            internal_multimap.insert (internal_pair);
            pair <std::map <string, multimap <long, GffRecordPtr> >::iterator, bool> ret;
            pair <string, multimap <long, GffRecordPtr> > external_pair (current_isoform.chrom, internal_multimap);
            ret = global_annotation_map_ptr.insert (external_pair);
            if (ret.second == false) {
                ret.first->second.insert (internal_pair);
            }
        }
    }


    return true;
}


void print_ref_info (const std::map <string, pair <int, int> > & info_map){
    cout << "PRINT REFERENCE DATA FORM BAM FILE" << endl;
    for (auto it = info_map.begin(); it != info_map.end(); ++it){
        cout << "Chromosome: " << it->first << endl;
        cout << "RefId: " << it->second.first << endl;
        cout << "Length: " << it->second.second << endl;
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

// if input_set forms linked list, for each of the element, except the first one, I can find the previous one.
// If not - return false;
bool form_line (set<GffRecordPtr> input_set){
    cout << "Check if set of annotation forms linked list" << endl;
    if (input_set.size() == 1) return true;
    int counter = 0;
    for (auto it = input_set.begin(); it != input_set.end(); ++it){
        if ((*it)->previous_gff){
            auto it_check = input_set.find((*it)->previous_gff);
            if (it_check == input_set.end()) {
                cout << "for element " << (*it)->exon_id << " from " << (*it)->isoform_id << " cannot find any previous" << endl;
                counter++;
            } else{
                cout << "for element " << (*it)->exon_id << " from " << (*it)->isoform_id << " found previous " << (*it_check)->exon_id << " from " << (*it_check)->isoform_id <<  endl;
            }
        } else {
            cout << "Found first element of isoform: " << (*it)->exon_id << " from " << (*it)->isoform_id << endl;
            counter++;
        }
    }
    if (counter > 1){
        cout << "counter > 1. Impossible to have two or more unlinked elements" << endl;
        return false;
    }
    return true;
}

int main() {

    // read from BAM/SAM file
    string bam_full_path_name = "/Users/kot4or/ClionProjects/test_1/bam_ex_1.bam";
//    string bam_full_path_name = "/Users/kot4or/ClionProjects/samtools_primer/tutorial/alignments/sim_reads_aligned.bam";
    BamReader bam_reader;
    if (not bam_reader.Open(bam_full_path_name)) {
        cout << "Couldn't open file " << bam_full_path_name << endl;
        return 0;
    } else cout << "Open " << bam_reader.GetFilename() << endl;

    cout << endl << endl;

    // key - chromosome name, value - <RefId, Length> for corresponding chromosome from the BAM file
    std::map <string, pair <int, int> > chromosome_info_map = get_chromosome_map_info (bam_reader);

    print_ref_info (chromosome_info_map); // Only for DEBUG

    if (not bam_reader.HasIndex()){
        cout << "Current BAM file isn't indexed" << endl;
        cout << "Trying to find index files in the same directory" << endl;
        if (bam_reader.LocateIndex(BamIndex::STANDARD)){
            cout << "Located and loaded index file from disk" << endl;
        } else {
            cout << "Couldn't locate index file" << endl;
            cout << "Trying to create the new one" << endl;
            if (not bam_reader.CreateIndex(BamIndex::STANDARD)){
                cout << "Cannot create index for current bam file. Exiting" << endl;
                return 0;
            };
            cout << "Index file for current BAM file is succesfully created" << endl;
        }
    }

    // read from tab delimited file
    string annotation_full_path_name = "/Users/kot4or/ClionProjects/test_1/tab_del_ex_1";
    std::map <string, multimap <long, GffRecordPtr> > global_annotation_map_ptr;

    // map to save <chromosome name, <isoform name, correspondent index in array> >
    std::map <string, std::map <string, int> > iso_var_map;
    if (not load_annotation (annotation_full_path_name, global_annotation_map_ptr, iso_var_map)){
        return 0;
    }
    cout << endl;
    print_iso_var_map (iso_var_map);
    cout << endl;

    // FOR DEBUG USE ONLY
    cout << "ANNOTATIONS" << endl;
    for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
        for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
            cout << (*start_it->second).exon_id  << " " << (*start_it->second).isoform_id << " - [";
            cout << (*start_it->second).start_pose << "," << (*start_it->second).end_pose << "]" << endl;
            for (auto bam_record_it = (*start_it->second).bam_records.begin(); bam_record_it != (*start_it->second).bam_records.end(); ++bam_record_it){
                cout << (*bam_record_it)->read_id << " ";
            }
        }
    }


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
            cout << "EXIT. Find a bug in a code. You suppose to index BAM file after opening" << endl;
            return 0;
        }
        cout << "Current BAM file is indexed" << endl;
        cout << "Trying to set region limited by current chromosome: " << chrom << endl;

        if (not bam_reader.SetRegion(ref_id, 0, ref_id, length)){
            cout << "Cannot set region. Exit" << endl;
                cout << bam_reader.GetErrorString(); // added just in case
            return 0;
        }
        cout << "Region from BAM file is succesfully set for chromosome " << chrom << " with RefId = " << ref_id << endl << endl;

        // Making an interval map on the base of the annotation.
        // By default for each current_map_element we add the corresponding annotation pointer
        interval_map<long, MapElement> gtf_records_splitted;
        for (auto it = chrom_it->second.begin(); it != chrom_it->second.end(); ++it) {
            GffRecordPtr temp_ptr = it->second;
            MapElement current_map_element;
            current_map_element.gtf_records.push_front(temp_ptr);
            gtf_records_splitted.add(make_pair(interval<long>::closed((*it->second).start_pose, (*it->second).end_pose),
                                               current_map_element)
            );
        }

        // create an empty matrix: column - one interval from interval map, row - isoforms
        vector <vector <double> > weight_array (gtf_records_splitted.iterative_size(), vector <double> (iso_var_map[chrom].size(), 0));
        // FOR DEBUG ONLY
        cout << "WEIGHT ARRAY" << endl;
        for (int i = 0; i < weight_array.size(); i++) {
            cout <<  i <<") ";
            for (int j = 0; j < weight_array[i].size(); j++) {
                cout << weight_array[i][j] << " ";
            }
            cout << endl;
        }




        // FOR DEBUG USE ONLY
        cout << "GENERATE INTERVAL MAP" << endl;
        for (auto temp_it = gtf_records_splitted.begin(); temp_it != gtf_records_splitted.end(); ++temp_it) {
            cout << "[" << temp_it->first.lower() << "," << temp_it->first.upper() << "]" << " : ";
            for (auto it = temp_it->second.gtf_records.begin(); it != temp_it->second.gtf_records.end(); ++it) {
                cout << (*it)->exon_id << " (" << (*it)->isoform_id << "), ";
            }
            cout << endl;
        }
        cout << endl;

        // The main iterator to iterate over interval map segment array
        interval_map<long, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();
        // Additinal iterator to backup current_gtf_records_splitted_it iterator when we are processing spliced read]
        // If the read is standard the backup iterator doesn't influence on the result
        interval_map<long, MapElement>::iterator backup_current_gtf_records_splitted_it = gtf_records_splitted.begin();

        BamRecordPtr current_bam_record; // pointer to save current bam record, which we are trying to align
        bool freeze = false; // if true - calling the get_bam_record function return's the same read as it it did it before
        long slice_number = 1; // if read is spliced slice_number > 1
        long current_slice = 1; // defines the position of cuurent read as a part of a big spliced read
        std::map<string, set<GffRecordPtr> > iso_map; // map to arrange annotation according to the isoform key
        BamRecord previous_bam_record; // temporal pointer to bam record to detect the moment when next bam record isn't a part of big scpliced read

        while (get_bam_record(bam_reader, current_bam_record, freeze)) { // Check if I can get new record from BAM file
            // Check if gtf records array is already empty. Break the while loop
            if (current_gtf_records_splitted_it == gtf_records_splitted.end()) break;

            // FOR DEBUG USE ONLY
            cout << endl << "Process read " << current_bam_record->read_id << " [" <<
            current_bam_record->start_pose << "," <<
            current_bam_record->end_pose << "] - " <<
            current_bam_record->slices << " parts" << endl;

            // get the number of parts from current bam record
            slice_number = current_bam_record->slices;

            // if the next read isn't a part of the same big spliced read (in this case we check if read_id is equal)
            // we set current_slice equal to 1 and clear iso_map, because we don't want to have annotation which are relevant for the previous read
            // make a backup for interval map iterator.
            // We need to make a backup because the next read could be a part of spliced read
            // and in that case we will change current_gtf_records_splitted_it which we will restore from
            // its backup version when we rich the last part in spliced read
            // TODO maybe it should be the other way to check it
            if (previous_bam_record.read_id != current_bam_record->read_id) {
                current_slice = 1;
                iso_map.clear();
                backup_current_gtf_records_splitted_it = current_gtf_records_splitted_it; // update the buckup for gff iterator
                cout << "iso_map is cleared" << endl;
            }



            // Find start segment annotation
            // If false call get_bam_record again. freeze is true if we need to skip the read and false if we want to skip annotation
            if (not find_start_segment_annotation(current_bam_record, previous_bam_record, current_gtf_records_splitted_it, freeze)){
                previous_bam_record = *current_bam_record; // update previous_bam_record with current value
                continue;
            }
            previous_bam_record = *current_bam_record; // update previous_bam_record with current value


            // FOR DEBUG USE ONLY
            print_segment_annotation("   start segment annotation", current_gtf_records_splitted_it);

            // Find stop segment annotation
            // If false call get_bam_record with freeze = false ===> skip current bam record and get the next one
            interval_map<long, MapElement>::iterator temp_gtf_records_splitted_it = current_gtf_records_splitted_it;
            if (not find_stop_segment_annotation(current_bam_record, temp_gtf_records_splitted_it, gtf_records_splitted.end(), freeze))
                continue;

            // FOR DEBUG USE ONLY
            print_segment_annotation("   stop segment annotation", temp_gtf_records_splitted_it);



            // find intersection of two sets
            // we receive set of annotation which is similar between start and end interval map segments
            // in other words we receive pointers to the annotations which includes current bam read
            set<GffRecordPtr> gff_intersection = get_intersection(current_gtf_records_splitted_it,
                                                                  temp_gtf_records_splitted_it);

            // iteratore over the intersection set
            for (auto gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it) {
                set<GffRecordPtr> temp_set; // we save only one pointer to this set, but we need to use set, because we want to add it into iso_map
                temp_set.insert((*gff_it));
                // check if slice_number == 1 which is equal that read isn't spliced or
                // current annotation fits the condition which are required for specific part of the spliced read (start, middle or end part of spliced read)
                if (slice_number == 1 or fit_spliced_read_condition(current_slice, slice_number, temp_set, current_bam_record)) {
                    // Add new element into map
                    pair<std::map<string, set<GffRecordPtr> >::iterator, bool> ret;
                    pair<string, set<GffRecordPtr> > input_pair((*gff_it)->isoform_id, temp_set);
                    ret = iso_map.insert(input_pair);
                    // If already exist with the same isoform key - add annotation to the corresponding set of the map
                    if (ret.second == false) {
                        cout << "Isoform " << input_pair.first << " already exists" << endl;
                        cout << "Updating the set with a value " << (*gff_it)->exon_id << endl;
                        ret.first->second.insert(temp_set.begin(), temp_set.end());
                    }
                }
            }
            // increment current_slice (position inside a spliced read)
            current_slice++;

            // If we reached the end of the big spliced read (equal to current_slice > slice_number)
            // In case of unspliced read slice_number = 1 and current_slice will be 2 after incrementing in the previous line
            if (current_slice > slice_number) {
                // iterating over the isoforms from iso_map
                for (auto map_iterator = iso_map.begin(); map_iterator != iso_map.end(); map_iterator++) {
                    // if current isoform has the same number of annotation as number of parts from the spliced read
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
                        for (auto gff_it = map_iterator->second.begin();
                             gff_it != map_iterator->second.end(); ++gff_it) {
                            (*gff_it)->bam_records.push_back(current_bam_record);
                            (*gff_it)->reads_count++;
                            cout << "   into annotation " << (*gff_it)->exon_id << " from " << (*gff_it)->isoform_id << " added read " <<
                            current_bam_record->read_id << endl;
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

    }

    // FOR DEBUG USE ONLY

    cout << endl << "RESULTS" << endl;
    for (auto chrom_it = global_annotation_map_ptr.begin(); chrom_it != global_annotation_map_ptr.end(); ++chrom_it){
        cout << "Chrom: " << chrom_it->first << endl;
        for (auto start_it = chrom_it->second.begin(); start_it!=chrom_it->second.end(); ++start_it){
            cout << "exon " << (*start_it->second).exon_id  << " from " << (*start_it->second).isoform_id << " - [";
            cout << (*start_it->second).start_pose << "," << (*start_it->second).end_pose << "]" << " : ";
            for (auto bam_record_it = (*start_it->second).bam_records.begin(); bam_record_it != (*start_it->second).bam_records.end(); ++bam_record_it){
                cout << (*bam_record_it)->read_id << " ";
            }
            cout << endl;
        }
    }

    return 0;
}
