//
// Created by kot4or on 8/8/16.
//

#include "annotation_reader.h"



void Isoform::print (){
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

Isoform::Isoform (string line):
         bin (0)
        ,name ("")
        ,chrom("")
        ,strand(true)
        ,tx_start(0)
        ,tx_end(0)
        ,cds_start(0)
        ,cds_end(0)
        ,exon_count(0)
        ,length (0)
        ,total_reads(0)
        ,density(0)
        ,rpkm(0)
        ,index(0)
        ,score(0)
        ,name2("")
        ,cds_start_stat (cds_stat::none)
        ,cds_end_stat (cds_stat::none)
{
    vector<string> line_splitted = split_line(line);

    // BIN
    if (not str_to_long(bin, line_splitted[0])){
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
    if (not str_to_long(exon_count, line_splitted[8])){
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

    vector<string> exon_starts_str = split_line(line_splitted[9], ",");
    vector<string> exon_ends_str = split_line(line_splitted[10], ",");
    vector<string> exon_frames_str = split_line(line_splitted[15], ",");

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

    // calculating the length of isiform as sum of all exons' lengths
    for (int i = 0; i < exon_count; i++){
        length += (long)(exon_ends[i] - exon_starts[i]);
    }

}

Isoform::Isoform ():
        bin (0)
        ,name ("")
        ,chrom("")
        ,strand(true)
        ,tx_start(0)
        ,tx_end(0)
        ,cds_start(0)
        ,cds_end(0)
        ,exon_count(0)
        ,length (0)
        ,total_reads(0)
        ,density(0)
        ,rpkm(0)
        ,index(0)
        ,score(0)
        ,name2("")
        ,cds_start_stat (cds_stat::none)
        ,cds_end_stat (cds_stat::none)
{

}


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

void print_iso_var_map (const std::map <string, std::map <string, Isoform> > & iso_var_map){
    for (auto ext_it = iso_var_map.begin(); ext_it != iso_var_map.end(); ++ext_it){
        cout << "Chromosome: " << ext_it->first << endl;
        for (auto int_it = ext_it->second.begin(); int_it != ext_it->second.end(); ++int_it){
            cout << setw(10) << "  isoform: "
                 << setw(15) << int_it->first
                 << setw(10) << " index: "
                 << setw(3) << int_it->second.index
                 << setw(10) << " length: "
                 << setw(5) << int_it->second.length
                 << setw(15) << " total_reads: "
                 << setw(5) << int_it->second.total_reads
                 << setw(10) << " density: "
                 << setw(15) << int_it->second.density
                 << setw(10) << " rpkm: "
                 << setw(5) << int_it->second.rpkm << endl;
        }
    }
}

void print_iso_var_map_to_file (const std::map <string, std::map <string, Isoform> > & iso_var_map, const string path){
    ofstream output_stream (path);
    if (output_stream.is_open())
    {
        for (auto ext_it = iso_var_map.begin(); ext_it != iso_var_map.end(); ++ext_it){
            output_stream << "Chromosome: " << ext_it->first << endl;
            for (auto int_it = ext_it->second.begin(); int_it != ext_it->second.end(); ++int_it){
                output_stream << setw(10) << "  isoform: "
                     << setw(15) << int_it->first
                     << setw(10) << " index: "
                     << setw(3) << int_it->second.index
                     << setw(10) << " length: "
                     << setw(5) << int_it->second.length
                     << setw(15) << " total_reads: "
                     << setw(5) << int_it->second.total_reads
                     << setw(10) << " density: "
                     << setw(15) << int_it->second.density
                     << setw(10) << " rpkm: "
                     << setw(5) << int_it->second.rpkm << endl;
            }
        }
        output_stream.close();
    }
    else cout << "Unable to open output file: " << path << endl;
}



// global_annotation_map_ptr : key - chromosome name, value - multimap of annotations, sorted by not-unique key - start pose of annotation
// NOTE : forward list of annotations should be sorted by start pose with rule a<b
bool load_annotation (const string & full_path_name,
                      std::map <string, multimap <long, GffRecordPtr> > & global_annotation_map_ptr,
                      std::map <string, std::map <string, Isoform> > & iso_var_map){
    ifstream input_stream (full_path_name);
    if (!input_stream) {
        cout << "Cannot open file " << full_path_name << endl;
        return false;
    }
    string line;
    while (getline(input_stream, line)) {
        if (string_tools::include_key(line, "name")) { // to filter header lines
            continue;
        }
        Isoform current_isoform;
        try {
            current_isoform = Isoform (line);
        } catch (...){
            cout << "Skipped line [" << line << "]" << endl;
            continue;
        }

        current_isoform.index = iso_var_map[current_isoform.chrom].size()+1;
        pair <string, Isoform> internal_pair_for_iso_var_map (current_isoform.name, current_isoform);
        std::map <string, Isoform> internal_iso_var_map;
        internal_iso_var_map.insert(internal_pair_for_iso_var_map);
        pair <std::map <string, std::map <string, Isoform> >::iterator, bool> res;
        pair <string, std::map <string, Isoform> > external_pair_for_iso_var_map (current_isoform.chrom, internal_iso_var_map);
        res = iso_var_map.insert (external_pair_for_iso_var_map);
        if (res.second == false){
            res.first->second.insert(internal_pair_for_iso_var_map);
        }

        GffRecordPtr previous_annotation; // TODO maybe this is the reason of error. Try to create empty GffRecord object
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