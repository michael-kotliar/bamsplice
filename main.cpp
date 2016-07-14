#include <iostream>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>
#include <boost/icl/interval_map.hpp>

using namespace std;
using namespace boost::icl;

class BamRecord;
class GffRecord;
typedef boost::shared_ptr<BamRecord> BamRecordPtr;
typedef boost::shared_ptr<GffRecord> GffRecordPtr;


class BamRecord {
    public:
        size_t start_pose;
        size_t end_pose;
        // add other important fields
};

class GffRecord {
public:
        size_t start_pose;
        size_t end_pose;
        string exon_id;
        string isoform_id;
        vector < BamRecordPtr > bam_records;
        size_t reads_count;
};

class MapElement {
public:
        set < GffRecordPtr > gtf_records;
        // overload operator += to push_back
    MapElement& operator+=(const MapElement& new_element) {
        this->gtf_records.insert(new_element.gtf_records.begin(), new_element.gtf_records.end());
        return *this;
    }
};





int main() {

    // read from BAM/SAM file to bam_records_input, save as a set of pointers
    set <BamRecordPtr> bam_records_input;
    // read from GFF/GFT file to gff_records_input
    set <GffRecordPtr> gff_records_input;
    interval_map<size_t,MapElement> gtf_records_splitted;

    for (std::set<GffRecordPtr>::iterator it = gff_records_input.begin(); it!=gff_records_input.end(); ++it){
        // create shared pointer from iterator and add it to set of pointers
        GffRecordPtr temp_ptr;
        temp_ptr = *it;
        MapElement current_map_element;
        current_map_element.gtf_records.insert(temp_ptr);
        gtf_records_splitted += make_pair( interval<size_t>::closed((*it)->start_pose, (*it)->end_pose), current_map_element );
    }


    set<BamRecordPtr>::iterator current_bam_record_it;
    interval_map<size_t, MapElement>::iterator current_gtf_records_splitted_it = gtf_records_splitted.begin();

    for (current_bam_record_it = bam_records_input.begin(); current_bam_record_it != bam_records_input.end(); ++current_bam_record_it){
        if ((*current_bam_record_it)->start_pose < current_gtf_records_splitted_it->first.lower()) {
            continue;
            // didn't reach even first annotation
        }

        interval_map<size_t, MapElement>::iterator temp_gtf_records_splitted_it;
        temp_gtf_records_splitted_it = current_gtf_records_splitted_it;

        while (temp_gtf_records_splitted_it->first.upper() <  (*current_bam_record_it)->end_pose){
            // find our closing annotation for current read
            temp_gtf_records_splitted_it++;
        }

        // find intersection of two sets

        set<GffRecordPtr> gff_intersection;
        set_intersection(current_gtf_records_splitted_it->second.gtf_records.begin(),current_gtf_records_splitted_it->second.gtf_records.end(),
                         temp_gtf_records_splitted_it->second.gtf_records.begin(),temp_gtf_records_splitted_it->second.gtf_records.end(),
                         std::inserter(gff_intersection, gff_intersection.begin()));

        if (gff_intersection.empty()){
            gff_intersection = current_gtf_records_splitted_it->second.gtf_records;
        }

        // iterate over gff_intersection set and write correct read to each of the annotation
        for (set<GffRecordPtr>::iterator gff_it = gff_intersection.begin(); gff_it != gff_intersection.end(); ++gff_it){
            (*gff_it)->bam_records.push_back(*current_bam_record_it);
        }


    }





        return 0;
}
