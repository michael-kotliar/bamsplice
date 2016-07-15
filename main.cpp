#include <iostream>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>
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
    // add other important fields

    BamRecord (size_t start, size_t end, string read)
            : start_pose (start)
            , end_pose (end)
            , read_id (read)
    {}
    BamRecord ()
            : start_pose (0)
            , end_pose (0)
            , read_id ("")
    {}
};

class GffRecord {
public:
        size_t start_pose;
        size_t end_pose;
        string exon_id;
        string isoform_id;
        set < BamRecordPtr > bam_records;
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
    set < GffRecordPtr > gtf_records;

    inline MapElement& operator+=(const MapElement& other_element) {
        this->gtf_records.insert(other_element.gtf_records.begin(), other_element.gtf_records.end());
        return *this;
    }

    inline bool operator==(const MapElement& other_element) const
    {
        return this->gtf_records == other_element.gtf_records;
    }

};




int main() {


    // read from BAM/SAM file to bam_records_input, save as a set of pointers
    set <BamRecordPtr> bam_records_input;

            // FOR TESTING ONLY
            BamRecordPtr r1 (new BamRecord (12,15, "1"));
            BamRecordPtr r2 (new BamRecord (14,21, "2"));
            BamRecordPtr r3 (new BamRecord (22,33, "3"));
            BamRecordPtr r4 (new BamRecord (22,45, "4"));
            BamRecordPtr r5 (new BamRecord (40,45, "5"));

            bam_records_input.insert(r1);
            bam_records_input.insert(r2);
            bam_records_input.insert(r3);
            bam_records_input.insert(r4);
            bam_records_input.insert(r5);





    // read from GFF/GFT file to gff_records_input
    set <GffRecordPtr> gff_records_input;

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

            gff_records_input.insert(a1);
            gff_records_input.insert(a2);
            gff_records_input.insert(a3);
            gff_records_input.insert(a4);
            gff_records_input.insert(a5);
            gff_records_input.insert(a6);
            gff_records_input.insert(a7);
            gff_records_input.insert(a8);
            gff_records_input.insert(a9);



    interval_map<size_t, MapElement> gtf_records_splitted;

    for (std::set<GffRecordPtr>::iterator it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        // create shared pointer from iterator and add it to set of pointers
        GffRecordPtr temp_ptr;
        temp_ptr = *it;
        MapElement current_map_element;
        current_map_element.gtf_records.insert(temp_ptr);
        gtf_records_splitted.add( make_pair( interval<size_t>::closed( (*it)->start_pose, (*it)->end_pose ),
                                            current_map_element )
                                );
    }

    // FOR DEBUG ONLY
    cout << "GENERATE INTERVAL MAP" << endl;
    for (interval_map<size_t, MapElement>::iterator temp_it = gtf_records_splitted.begin(); temp_it !=  gtf_records_splitted.end(); ++temp_it) {
        cout << "[" << temp_it->first.lower() << "," << temp_it->first.upper() << "]" << " : ";
        for (set<GffRecordPtr>::iterator it = temp_it->second.gtf_records.begin();
             it != temp_it->second.gtf_records.end(); ++it) {
            cout << (*it)->exon_id << " ";
        }
        cout <<endl;
    }
    cout << endl;

    set<BamRecordPtr>::iterator current_bam_record_it = bam_records_input.begin();
    interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();

    while (current_bam_record_it != bam_records_input.end() and current_gtf_records_splitted_it != gtf_records_splitted.end()){
        cout << "Process read " << (*current_bam_record_it)->read_id << " [" <<
                (*current_bam_record_it)->start_pose << "," <<
                (*current_bam_record_it)->end_pose << "]" << endl;
        if ((*current_bam_record_it)->start_pose < current_gtf_records_splitted_it->first.lower()) {
            cout << "   Skip read " << (*current_bam_record_it)->read_id << " [" <<
                                       (*current_bam_record_it)->start_pose << "," <<
                                       (*current_bam_record_it)->end_pose << "]" << endl;
            continue;
        }
        if ((*current_bam_record_it)->start_pose >= current_gtf_records_splitted_it->first.upper()) {
            cout << (*current_bam_record_it)->start_pose << " > " << current_gtf_records_splitted_it->first.upper() << endl;
            cout << "   Skip segment annotation : " << "[" <<  current_gtf_records_splitted_it->first.lower() << ","
                                                  << current_gtf_records_splitted_it->first.upper() << "]" << endl;
            current_gtf_records_splitted_it++;
            continue;
        }
        interval_map<size_t, MapElement>::iterator temp_gtf_records_splitted_it = current_gtf_records_splitted_it;
        cout << "   start segment annotation " << "[" <<  current_gtf_records_splitted_it->first.lower() << ","
                                               << current_gtf_records_splitted_it->first.upper() << "] :";

        for (set<GffRecordPtr>::iterator start_segment_annotation_it = current_gtf_records_splitted_it->second.gtf_records.begin();
             start_segment_annotation_it != current_gtf_records_splitted_it->second.gtf_records.end(); ++start_segment_annotation_it){
            cout << " " << (*start_segment_annotation_it)->exon_id;
        }
        cout << endl;

        while ((*current_bam_record_it)->end_pose < temp_gtf_records_splitted_it->first.lower() or
                (*current_bam_record_it)->end_pose > temp_gtf_records_splitted_it->first.upper()){
            // find our closing annotation for current read
            if (temp_gtf_records_splitted_it == gtf_records_splitted.end()){
                current_bam_record_it++;
                continue;
            }
            temp_gtf_records_splitted_it++;
        }

        cout << "   stop segment annotation " << "[" <<  temp_gtf_records_splitted_it->first.lower() << ","
        << temp_gtf_records_splitted_it->first.upper() << "] :";

        for (set<GffRecordPtr>::iterator stop_segment_annotation_it = temp_gtf_records_splitted_it->second.gtf_records.begin();
             stop_segment_annotation_it != temp_gtf_records_splitted_it->second.gtf_records.end(); ++stop_segment_annotation_it){
            cout << " " << (*stop_segment_annotation_it)->exon_id;
        }
        cout << endl;



        // find intersection of two sets

        set<GffRecordPtr> gff_intersection;
        set_intersection(current_gtf_records_splitted_it->second.gtf_records.begin(),current_gtf_records_splitted_it->second.gtf_records.end(),
                         temp_gtf_records_splitted_it->second.gtf_records.begin(),temp_gtf_records_splitted_it->second.gtf_records.end(),
                         std::inserter(gff_intersection, gff_intersection.begin()));

        cout << "   Intersection : ";
        for (set<GffRecordPtr>::iterator intersection_segment_annotation_it = gff_intersection.begin();
             intersection_segment_annotation_it != gff_intersection.end(); ++intersection_segment_annotation_it){
            cout << " " << (*intersection_segment_annotation_it)->exon_id;
        }
        cout << endl;


        // iterate over gff_intersection set and write correct read to each of the annotation
        for (set<GffRecordPtr>::iterator gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it){
            (*gff_it)->bam_records.insert(*current_bam_record_it);
            (*gff_it)->reads_count++;
            cout << "   into annotation " << (*gff_it)->exon_id << " added read " << (*current_bam_record_it)->read_id << endl;
        }

        current_bam_record_it++;


    }

            // FOR TESTING ONLY
            for (std::set<GffRecordPtr>::iterator it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
                cout << "Annotation "<< (*it)->exon_id << " : ";
                for (set<BamRecordPtr>::iterator bam_record_it = (*it)->bam_records.begin(); bam_record_it != (*it)->bam_records.end(); ++bam_record_it){
                    cout << (*bam_record_it)->read_id << " ";
                }
                cout << endl;
            }





        return 0;
}
