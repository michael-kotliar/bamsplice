#include <iostream>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>
#include <map>
#include <forward_list>
#include <boost/icl/interval_map.hpp>

using namespace std;
using namespace boost::icl;

class BamRecord;
class GffRecord;
class MapElement;

typedef boost::shared_ptr<BamRecord> BamRecordPtr;
typedef boost::shared_ptr<GffRecord> GffRecordPtr;

class BamRecord {
public:
    size_t start_pose; // start position of the read
    size_t end_pose; // stop position of the read
    string read_id; // text identificator of the read
    short slices; // if the read is spliced - set here total amount of parts, that form this read, default = 1
    // ADD HERE OTHER INMPORTANT FIELDS

    // Constructor with parameters
    BamRecord (size_t start, size_t end, string read, short parts_number)
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
};

class GffRecord {
public:
    size_t start_pose; // start position of annotation
    size_t end_pose; // stop position of annotation
    string exon_id; // text-identificator of current annotation
    string isoform_id; // set the name of the isofom=rm to which current annotation belongs
    vector < BamRecordPtr > bam_records; // array of pointers to all of the reads, which belongs to this annotation
    size_t reads_count; // total number of reads, which belongs to this annotation

    // CONSTRUCTOR WITH PARAMETERS
    GffRecord (size_t start, size_t end, string exon, string isoform)
            : start_pose (start)
            , end_pose (end)
            , exon_id (exon)
            , isoform_id (isoform)
            , reads_count (0)
    {}

    // EMPTY CONSTRUCTOR
    GffRecord ()
            : start_pose (0)
            , end_pose (0)
            , exon_id ("")
            , isoform_id ("")
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




size_t last_bam_pose = 0; // only for function get_bam_record to save the last index of the bam record from input array
// This function should be upfated to the similar one when we load the real bam/sam file
bool get_bam_record (const vector <BamRecordPtr> & bam_arr, BamRecordPtr & bam_record, bool freeze = false){
    if (not freeze){
        last_bam_pose++;
    }
    if (bam_arr.empty() or last_bam_pose > bam_arr.size()){
        cout << "bam array is empty " << bam_arr.empty() << endl;
        cout << last_bam_pose << " >= " << bam_arr.size() << endl;
        return false;
    }
    if (last_bam_pose == 0){
        bam_record = bam_arr [last_bam_pose];
    } else
        bam_record = bam_arr [last_bam_pose - 1];
    return true;
}

// Updates iterator current_gtf_records_splitted_it for the interval map segment which includes the starting position of the current_bam_read
// If fails - returns false
// If we need to skip bam_record, we set freeze = false to get the new current_bam_record
// If we need to skip interval map segment we set freeze to true to avoid changing current_bam_record, when iterating over the loop,
// that get new value from bam records array
// If skip interval map segment - increment iterator current_gtf_records_splitted_it
bool find_start_segment_annotation (BamRecordPtr current_bam_record, interval_map<size_t, MapElement>::iterator & current_gtf_records_splitted_it, bool & freeze){
    if (current_bam_record->start_pose < current_gtf_records_splitted_it->first.lower()) {
        cout << "   Skip read " << current_bam_record->read_id << " [" <<
        current_bam_record->start_pose << "," <<
        current_bam_record->end_pose << "]" << endl;
        freeze = false; // Set freeze to false to change current_bam_record
        return false;
    }
    if (current_bam_record->start_pose >= current_gtf_records_splitted_it->first.upper()) {
        cout << current_bam_record->start_pose << " > " << current_gtf_records_splitted_it->first.upper() << endl;
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
                                   interval_map<size_t,MapElement>::iterator & temp_gtf_records_splitted_it,
                                   interval_map<size_t,MapElement>::iterator max_segment_annotation,
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
void print_segment_annotation (const string & title, interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it){
    cout << title << " " << "[" << current_gtf_records_splitted_it->first.lower() << ","
    << current_gtf_records_splitted_it->first.upper() << "] :";
    for (auto start_segment_annotation_it = current_gtf_records_splitted_it->second.gtf_records.begin();
         start_segment_annotation_it != current_gtf_records_splitted_it->second.gtf_records.end(); ++start_segment_annotation_it){
        cout << " " << (*start_segment_annotation_it)->exon_id;
    }
    cout << endl;
}


// Returns set of annotation pointers which are made as a result of intersection between two annotation pointers sets from two interval map segments
set<GffRecordPtr> get_intersection (interval_map<size_t, MapElement>::iterator input_1, interval_map<size_t, MapElement>::iterator input_2){
    set<GffRecordPtr> intersection;
    input_1->second.gtf_records.sort();
    input_2->second.gtf_records.sort();
    set_intersection(input_1->second.gtf_records.begin(),input_1->second.gtf_records.end(),
                     input_2->second.gtf_records.begin(),input_2->second.gtf_records.end(),
                     std::inserter(intersection, intersection.begin()));
                cout << "   Intersection : ";
                for (auto intersection_segment_annotation_it = intersection.begin();
                     intersection_segment_annotation_it != intersection.end(); ++intersection_segment_annotation_it){
                    cout << " " << (*intersection_segment_annotation_it)->exon_id;
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
bool fit_spliced_read_condition(const size_t & current_slice, const size_t & slice_number, const set <GffRecordPtr> & temp_set, BamRecordPtr current_bam_record){
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


int main() {

    // read from BAM/SAM file to bam_records_input, save as a set of pointers
    // Should be sorted according to the start position
    vector <BamRecordPtr> bam_records_input;
    BamRecordPtr r1 (new BamRecord (12,15, "1", 1));
    BamRecordPtr r2 (new BamRecord (14,21, "2", 1));
    BamRecordPtr r3 (new BamRecord (17,20, "6", 2));
    BamRecordPtr r6 (new BamRecord (36,45, "6", 2));
    BamRecordPtr r4 (new BamRecord (22,33, "3", 1));
    BamRecordPtr r5 (new BamRecord (22,45, "4", 1));
    BamRecordPtr r7 (new BamRecord (40,45, "5", 1));
    bam_records_input.push_back(r1);
    bam_records_input.push_back(r2);
    bam_records_input.push_back(r3);
    bam_records_input.push_back(r6);
    bam_records_input.push_back(r4);
    bam_records_input.push_back(r5);
    bam_records_input.push_back(r7);

    // FOR DEBUG USE ONLY
    cout << "READS" << endl;
    for (auto bam_record_it = bam_records_input.begin(); bam_record_it != bam_records_input.end(); ++bam_record_it){
        cout << (*bam_record_it)->read_id << " - [" << (*bam_record_it)->start_pose << "," << (*bam_record_it)->end_pose << "]" << endl;
    }

    // read from GFF/GFT file to gff_records_input
    forward_list <GffRecordPtr> gff_records_input;
    GffRecordPtr a1 (new GffRecord (10,20, "1" , "iso_1"));
    GffRecordPtr a2 (new GffRecord (22,33, "2" , "iso_1"));
    GffRecordPtr a3 (new GffRecord (36,47, "3" , "iso_1"));
    GffRecordPtr a4 (new GffRecord (50,58, "4" , "iso_1"));
    GffRecordPtr a5 (new GffRecord (17,33, "5" , "iso_2"));
    GffRecordPtr a6 (new GffRecord (40,45, "6" , "iso_2"));
    GffRecordPtr a7 (new GffRecord (50,58, "7" , "iso_2"));
    GffRecordPtr a8 (new GffRecord (22,30, "8" , "iso_3"));
    GffRecordPtr a9 (new GffRecord (36,58, "9" , "iso_3"));

    gff_records_input.push_front(a1);
    gff_records_input.push_front(a2);
    gff_records_input.push_front(a3);
    gff_records_input.push_front(a4);
    gff_records_input.push_front(a5);
    gff_records_input.push_front(a6);
    gff_records_input.push_front(a7);
    gff_records_input.push_front(a8);
    gff_records_input.push_front(a9);

    gff_records_input.reverse();  // we need to do reverse because we can use only push_front,
                                  // but we want to make it sorted with rule a>=b

    // FOR DEBUG USE ONLY
    cout << "ANNOTATIONS" << endl;
    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        cout << (*it)->exon_id << " " << (*it)->isoform_id << " - [";
        cout << (*it)->start_pose << "," << (*it)->end_pose << "]" << endl;
        for (auto bam_record_it = (*it)->bam_records.begin(); bam_record_it != (*it)->bam_records.end(); ++bam_record_it){
            cout << (*bam_record_it)->read_id << " ";
        }
    }

    // Making an interval map on the base of the annotation.
    // By default for each current_map_element we add the corresponding annotation pointer
    interval_map<size_t, MapElement> gtf_records_splitted;
    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        GffRecordPtr temp_ptr = *it;
        MapElement current_map_element;
        current_map_element.gtf_records.push_front(temp_ptr);
        gtf_records_splitted.add( make_pair( interval<size_t>::closed( (*it)->start_pose, (*it)->end_pose ),
                                            current_map_element )
                                );
    }

    // FOR DEBUG USE ONLY
    cout << "GENERATE INTERVAL MAP" << endl;
    for (auto temp_it = gtf_records_splitted.begin(); temp_it !=  gtf_records_splitted.end(); ++temp_it) {
        cout << "[" << temp_it->first.lower() << "," << temp_it->first.upper() << "]" << " : ";
        for (auto it = temp_it->second.gtf_records.begin();
             it != temp_it->second.gtf_records.end(); ++it) {
            cout << (*it)->exon_id << " ";
        }
        cout <<endl;
    }
    cout << endl;

    // The main iterator to iterate over interval map segment array
    interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();
    // Additinal iterator to backup current_gtf_records_splitted_it iterator when we are processing spliced read]
    // If the read is standard the backup iterator doesn't influence on the result
    interval_map<size_t, MapElement>::iterator backup_current_gtf_records_splitted_it = gtf_records_splitted.begin();

    BamRecordPtr current_bam_record; // pointer to save current bam record, which we are trying to align
    bool freeze = false; // if true - calling the get_bam_record function return's the same read as it it did it before
    size_t slice_number = 1; // if read is spliced slice_number > 1
    size_t current_slice = 1; // defines the position of cuurent read as a part of a big spliced read
    std::map <string, set <GffRecordPtr> > iso_map; // map to arrange annotation according to the isoform key
    BamRecord previous_bam_record; // temporal pointer to bam record to detect the moment when next bam record isn't a part of big scpliced read

    while (get_bam_record (bam_records_input, current_bam_record, freeze)){ // Check if I can get new record from BAM file
        // Check if gtf records array is already empty. Break the while loop
        if (current_gtf_records_splitted_it == gtf_records_splitted.end()) break;

        // FOR DEBUG USE ONLY
        cout << endl << "Process read " << current_bam_record->read_id << " [" <<
                                           current_bam_record->start_pose << "," <<
                                           current_bam_record->end_pose << "] - " <<
                                           current_bam_record->slices << " parts" << endl;

        // get the number of parts from current bam record
        slice_number = current_bam_record->slices;

        // Find start segment annotation
        // If false call get_bam_record again. freeze is true if we need to skip the read and false if we want to skip annotation
        if (not find_start_segment_annotation (current_bam_record, current_gtf_records_splitted_it, freeze))
            continue;

        // FOR DEBUG USE ONLY
        print_segment_annotation ("   start segment annotation", current_gtf_records_splitted_it);

        // Find stop segment annotation
        // If false call get_bam_record with freeze = false ===> skip current bam record and get the next one
        interval_map<size_t, MapElement>::iterator temp_gtf_records_splitted_it = current_gtf_records_splitted_it;
        if (not find_stop_segment_annotation (current_bam_record, temp_gtf_records_splitted_it, gtf_records_splitted.end(), freeze) )
            continue;

        // FOR DEBUG USE ONLY
        print_segment_annotation ("   stop segment annotation", temp_gtf_records_splitted_it);

        // if the next read isn't a part of the same big spliced read (in this case we check if read_id is equal)
        // we set current_slice equal to 1 and clear iso_map, because we don't want to have annotation which are relevant for the previous read
        // make a backup for interval map iterator.
        // We need to make a backup because the next read could be a part of spliced read
        // and in that case we will change current_gtf_records_splitted_it which we will restore from
        // its backup version when we rich the last part in spliced read
        // TODO maybe it should be the other way to check it
        if (previous_bam_record.read_id != current_bam_record->read_id){
            current_slice = 1;
            iso_map.clear();
            backup_current_gtf_records_splitted_it = current_gtf_records_splitted_it; // update the buckup for gff iterator
            cout << "iso_map is cleared" << endl;
        }
        previous_bam_record = *current_bam_record; // update previous_bam_record with current value

        // find intersection of two sets
        // we receive set of annotation which is similar between start and end interval map segments
        // in other words we receive pointers to the annotations which includes current bam read
        set<GffRecordPtr> gff_intersection = get_intersection (current_gtf_records_splitted_it, temp_gtf_records_splitted_it);

        // iteratore over the intersection set
        for (auto gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it){
            set <GffRecordPtr> temp_set; // we save only one pointer to this set, but we need to use set, because we want to add it into iso_map
            temp_set.insert ((*gff_it));
            // check if slice_number == 1 which is equal that read isn't spliced or
            // current annotation fits the condition which are required for specific part of the spliced read (start, middle or end part of spliced read)
            if (slice_number == 1 or fit_spliced_read_condition(current_slice, slice_number, temp_set, current_bam_record) ){
                // Add new element into map
                pair < std::map <string, set <GffRecordPtr> >::iterator, bool > ret;
                pair < string, set <GffRecordPtr> > input_pair ( (*gff_it)->isoform_id, temp_set);
                ret = iso_map.insert (input_pair);
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
        if (current_slice > slice_number){
            // iterating over the isoforms from iso_map
            for(auto map_iterator = iso_map.begin(); map_iterator != iso_map.end(); map_iterator++){
                // if current isoform has the same number of annotation as number of parts from the spliced read
                // we caught the correct spliced read
                // if not - go to next isoform
                if (map_iterator->second.size() == slice_number){
                    // iterating iver annotation of one specific isoform from iso_map and add to each of them
                    // pointer to a current bam record
                    for (auto gff_it = map_iterator->second.begin(); gff_it != map_iterator->second.end(); ++gff_it){
                        (*gff_it)->bam_records.push_back(current_bam_record);
                        (*gff_it)->reads_count++;
                        cout << "   into annotation " << (*gff_it)->exon_id << " added read " << current_bam_record->read_id << endl;
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

    // FOR DEBUG USE ONLY
    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        cout << "Annotation "<< (*it)->exon_id << " : ";
        for (auto bam_record_it = (*it)->bam_records.begin(); bam_record_it != (*it)->bam_records.end(); ++bam_record_it){
            cout << (*bam_record_it)->read_id << " ";
        }
        cout << endl;
    }

    return 0;
}
