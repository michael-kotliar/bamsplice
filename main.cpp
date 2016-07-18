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
    size_t start_pose;
    size_t end_pose;
    string read_id;
    short slices; // total ammount of parts if read is sliced
    // add other important fields

    BamRecord (size_t start, size_t end, string read, short parts_number)
            : start_pose (start)
            , end_pose (end)
            , read_id (read)
            , slices (parts_number)
    {}
    BamRecord ()
            : start_pose (0)
            , end_pose (0)
            , read_id ("")
            , slices (1)
    {}
};

class GffRecord {
public:
        size_t start_pose;
        size_t end_pose;
        string exon_id;
        string isoform_id;
        vector < BamRecordPtr > bam_records;
        size_t reads_count;

        GffRecord (size_t start, size_t end, string exon, string isoform)
                : start_pose (start)
                , end_pose (end)
                , exon_id (exon)
                , isoform_id (isoform)
                , reads_count (0)
        {}

        GffRecord ()
                : start_pose (0)
                , end_pose (0)
                , exon_id ("")
                , isoform_id ("")
                , reads_count (0)
        {}

};

class MapElement {
public:
    forward_list < GffRecordPtr > gtf_records;

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




size_t last_bam_pose = 0;
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

void print_segment_annotation (const string & title, interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it){
    cout << title << " " << "[" << current_gtf_records_splitted_it->first.lower() << ","
    << current_gtf_records_splitted_it->first.upper() << "] :";
    for (auto start_segment_annotation_it = current_gtf_records_splitted_it->second.gtf_records.begin();
         start_segment_annotation_it != current_gtf_records_splitted_it->second.gtf_records.end(); ++start_segment_annotation_it){
        cout << " " << (*start_segment_annotation_it)->exon_id;
    }
    cout << endl;
}


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

int main() {
    // read from BAM/SAM file to bam_records_input, save as a set of pointers
    vector <BamRecordPtr> bam_records_input;

            // FOR TESTING ONLY
            // Should be sorted according to the start position

            BamRecordPtr r1 (new BamRecord (12,15, "1", 1));
            BamRecordPtr r2 (new BamRecord (14,21, "2", 1));
            BamRecordPtr r3 (new BamRecord (17,20, "6", 2));
            BamRecordPtr r6 (new BamRecord (22,25, "6", 2));
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

    cout << "READS" << endl;
    for (auto bam_record_it = bam_records_input.begin(); bam_record_it != bam_records_input.end(); ++bam_record_it){
        cout << (*bam_record_it)->read_id << " - [" << (*bam_record_it)->start_pose << "," << (*bam_record_it)->end_pose << "]" << endl;
    }

    // read from GFF/GFT file to gff_records_input
    forward_list <GffRecordPtr> gff_records_input;

            // FOR TESTING ONLY
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

    gff_records_input.reverse();

    cout << "ANNOTATIONS" << endl;
    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        cout << (*it)->exon_id << " " << (*it)->isoform_id << " - [";
        cout << (*it)->start_pose << "," << (*it)->end_pose << "]" << endl;
        for (auto bam_record_it = (*it)->bam_records.begin(); bam_record_it != (*it)->bam_records.end(); ++bam_record_it){
            cout << (*bam_record_it)->read_id << " ";
        }
    }

    interval_map<size_t, MapElement> gtf_records_splitted;

    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        GffRecordPtr temp_ptr = *it;
        MapElement current_map_element;
        current_map_element.gtf_records.push_front(temp_ptr);
        gtf_records_splitted.add( make_pair( interval<size_t>::closed( (*it)->start_pose, (*it)->end_pose ),
                                            current_map_element )
                                );
    }

    // FOR DEBUG ONLY
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

    interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();



    BamRecordPtr current_bam_record;
    bool freeze = false;
    size_t slice_number = 1;
    size_t current_slice = 1;
    std::map <string, set <GffRecordPtr> > iso_map; // isoform is a key

    while (get_bam_record (bam_records_input, current_bam_record, freeze)){ // Check if I can get new record from BAM file
        // Check if gtf records array is already empty. Exit if empty
        cout << endl;
        if (current_gtf_records_splitted_it == gtf_records_splitted.end()) break;

                    cout << "Process read " << current_bam_record->read_id << " [" <<
                                               current_bam_record->start_pose << "," <<
                                               current_bam_record->end_pose << "] - " <<
                                               current_bam_record->slices << endl;

        slice_number = current_bam_record->slices;

        // Find start segment annotation
        if (not find_start_segment_annotation (current_bam_record, current_gtf_records_splitted_it, freeze))
            continue;

                    print_segment_annotation ("   start segment annotation", current_gtf_records_splitted_it);



        // Find stop segment annotation
        interval_map<size_t, MapElement>::iterator temp_gtf_records_splitted_it = current_gtf_records_splitted_it;
        if (not find_stop_segment_annotation (current_bam_record, temp_gtf_records_splitted_it, gtf_records_splitted.end(), freeze) )
            continue;

                    print_segment_annotation ("   stop segment annotation", temp_gtf_records_splitted_it);



        // find intersection of two sets

        set<GffRecordPtr> gff_intersection = get_intersection (current_gtf_records_splitted_it, temp_gtf_records_splitted_it);


        for (auto gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it){
            set <GffRecordPtr> temp_set;
            temp_set.insert ((*gff_it));
            pair < std::map <string, set <GffRecordPtr> >::iterator, bool > ret;
            pair < string, set <GffRecordPtr> > input_pair ( (*gff_it)->isoform_id, temp_set);
            ret = iso_map.insert (input_pair);
            if (ret.second == false) {
                cout << "Isoform " << input_pair.first << " already exists" << endl;
                cout << "Updating the set with a value " << (*gff_it)->exon_id << endl;
                ret.first->second.insert(temp_set.begin(), temp_set.end());
            }
        }

        current_slice++;

        if (current_slice > slice_number){
            cout << "   c " << current_slice << " > " << slice_number << endl;
            for(auto map_iterator = iso_map.begin(); map_iterator != iso_map.end(); map_iterator++){
                if (map_iterator->second.size() == slice_number){
                    for (auto gff_it = map_iterator->second.begin(); gff_it != map_iterator->second.end(); ++gff_it){
                        (*gff_it)->bam_records.push_back(current_bam_record);
                        (*gff_it)->reads_count++;
                        cout << "   into annotation " << (*gff_it)->exon_id << " added read " << current_bam_record->read_id << endl;
                    }
                }
            }
            current_slice = 1;
            iso_map.clear();
            cout << "map is cleared" << endl;
        }

        freeze = false;

    }
//
                    // FOR TESTING ONLY
                    for (auto it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
                        cout << "Annotation "<< (*it)->exon_id << " : ";
                        for (auto bam_record_it = (*it)->bam_records.begin(); bam_record_it != (*it)->bam_records.end(); ++bam_record_it){
                            cout << (*bam_record_it)->read_id << " ";
                        }
                        cout << endl;
                    }





        return 0;
}
