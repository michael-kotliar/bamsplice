//
// Created by kot4or on 8/8/16.
//

#include "rpkm_calculation.h"

// map to save <chromosome name, <isoform name, correspondent Isoform object> >
void print_isoform_by_name (const vector<vector<double> > & data_array,
                            std::map <string, std::map <string, Isoform> > iso_var_map,
                            string chr,
                            string isoform_name,
                            std::ostream & out){
    int index = iso_var_map[chr][isoform_name].index;
    out << "[" << chr << "] - [" << isoform_name << "]" << endl;
    for (int j = 0; j < data_array[index].size(); j++) {
        if (data_array[index][j] == 0){
//            out << setw(16) << "-----";
        } else {
            out << std::left << setw(16) << data_array[index][j] << "   [" << data_array[0][j] << "]   ";
        }
    }
    out << endl;
}



void print_weight_array(const vector<vector<double> > & weight_array, const string & title){
    cerr << endl << title << endl;
    cerr << std::setprecision(9) <<  std::fixed;
    for (int i = 0; i < weight_array.size(); i++) {
        cerr << setw(4) << i << " ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] == 0){
                cerr << setw(16) << "-----";
            } else {
                cerr << std::left << setw(16) << weight_array[i][j];
            }
        }
        cerr << endl;
    }
}

void transform_to_density (vector <vector <double> > & weight_array){
    for (int i = 1; i < weight_array.size(); i++) {
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != 0){
                weight_array[i][j] = weight_array[i][j] / weight_array[0][j];
            }
        }
    }
}

double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row, int & count){
    assert ( (row >= 0 ) and (row < weight_array.size()) );
    count = 0;
    double sum = 0;
    for (int j = 0; j < weight_array[row].size(); j++) {
        if (weight_array[row][j] != 0 ) {
            sum += weight_array[row][j];
            count++;
        }
    }
    return sum;
}

double get_average_by_row (const vector <vector <double> > & weight_array, const int & row){
    int count = 0;
    double sum = get_sum_by_row(weight_array, row, count);
    if (count == 0) return 0;
    return sum / count;
}

double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int & count){
    assert ( weight_array.size() > 1 );
    assert ( (column >= 0 ) and (column < weight_array[0].size()) );
    double sum = 0;
    count = 0;
    for (int j = 1; j < weight_array.size(); j++) {
        if (weight_array[j][column] != 0) {
            sum += weight_array[j][column];
            count++;
        }
    }
//    cout << "debug sum = " << sum << endl;
    return sum;
}

double get_average_by_column (const vector <vector <double> > & weight_array, const int & column){
    int count = 0;
    double sum = get_sum_by_column(weight_array, column, count);
    if (count == 0) return 0;
    return sum / count;
}

vector <double> get_average_density_by_all_isoforms (const vector <vector <double> > & weight_array){
    vector <double> average_density;
    cout << std::setprecision(9);
    for (int i = 1; i < weight_array.size(); i++) {
        double average_by_row = get_average_by_row(weight_array, i);
//        cout << "Average by row " << i <<" : " << average_by_row << endl;
        average_density.push_back( average_by_row );
    }
    return average_density;
}

vector <double> get_sum_density_by_all_intervals (const vector <vector <double> > & weight_array){
    assert ( weight_array.size() > 1 );
    vector <double> sum_density;
    int count;
    for (int i = 0; i < weight_array[0].size(); i++) {
        double sum_dens = get_sum_by_column(weight_array, i, count);
//        cout << "DEBUG: " << sum_dens << endl;
        sum_density.push_back( sum_dens );
    }
    return sum_density;
}

void update_isoforms_density_to_average_for_isoform (vector <vector <double> > & weight_array){
    vector <double> average_densities = get_average_density_by_all_isoforms (weight_array);
    for (int i = 1; i < weight_array.size(); i++) {
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != 0 ){
                weight_array[i][j] = average_densities[i-1];
            }
        }
    }
}

void adjust_isoforms_density_by_coef (vector <vector <double> > & weight_array, const vector <double> & original_densities_by_isoforms, const vector <double> & new_densities_by_isoforms){
    for (int i = 1; i < weight_array.size(); i++) {
//        cout << " k ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != 0){
                double coef = original_densities_by_isoforms[j]/new_densities_by_isoforms[j];
//                cout << setw(8) << coef;
                weight_array[i][j] = weight_array[i][j] * coef;
            } else {
//                cout << setw(8) << "-----";
            }
        }
    }
}



double sum_all (const vector <vector <double> > dens_matrix) {
    double sum = 0;
    for (int i = 1; i < dens_matrix.size(); i++) {
        for (int j = 0; j < dens_matrix[i].size(); j++) {
            sum+=abs(dens_matrix[i][j]);
        }
    }
    return sum;
}


void subtract_matrix (vector <vector <double> > & first, const vector <vector <double> > & second){
    assert (first.size() == second.size() && first.size() > 1);
    assert (first[0].size() == second[0].size());
    for (int i = 1; i < first.size(); i++) {
        for (int j = 0; j < first[i].size(); j++) {
            first[i][j] -= second[i][j];
        }
    }
}


int run_cycle (vector <vector <double> > & weight_array){
    double cutoff = 10e-9; // TODO put it in separate configuration file
    int cycles = 0;
    vector <vector <double> > tmp_matrix;
    tmp_matrix = weight_array;
    // Get array of original densities sum
    vector <double> original_sum_dens =  get_sum_density_by_all_intervals (weight_array);

//    print_array (original_sum_dens, "Original density sums");

    for (int i = 0; i < 2000; i++){
        // Update original density array with average values for isoforms
//        cout << endl << "Set average by row" << endl;
        update_isoforms_density_to_average_for_isoform (weight_array);
//        print_weight_array(weight_array, "Average density array");

//        cout << endl;

        // Get array of temporary densities sum
        vector <double> new_sum_dens =  get_sum_density_by_all_intervals (weight_array);
//        print_array (original_sum_dens, "Original density sums");
//        print_array (new_sum_dens, "New density sums");

//        cout << endl << "Update according to the differenceds in original and new density sums" << endl;

        // Update density array according to the coef calculated by new density sums
        adjust_isoforms_density_by_coef (weight_array, original_sum_dens, new_sum_dens);
//        print_weight_array(weight_array, "Updated density array");
//        cout << endl;

        cycles++;
//        cerr << "Cycle: " << cycles;
        subtract_matrix(tmp_matrix, weight_array);
//        print_weight_array(tmp_matrix, "Substracted matrix");
        double sum = sum_all (tmp_matrix);
//        cerr << " Sum: " << sum << endl;
        if( sum < cutoff ){
            break;
        }
//
        tmp_matrix = weight_array;
    }

    cout << "Cycles: " << cycles << endl;
    return cycles;
}


void print_array (const vector <double> & intput_array, const string & title, string tab_spacer){
    cout << title << endl;
    cout << std::setprecision(9) <<  std::fixed;
    cout << tab_spacer;
    for (int i = 0; i < intput_array.size(); i++){
        cout << std::left << setw(16) << intput_array [i];
    }
    cout << endl;
}


void calculate_totReads_density (const vector<vector<double> > & weight_array, std::map <string, Isoform> & iso_map,  const std::map <string, int> & correspondence_map){
    boost::mutex::scoped_lock scoped_lock(iso_var_map_mutex);
//    for (auto iso_it = iso_map.begin(); iso_it != iso_map.end(); ++iso_it) {
//        int index = iso_it->second.index;
//        for (int j = 0; j < weight_array[index].size(); j++) {
//            if (weight_array[index][j] != 0){
//                iso_it->second.density += weight_array[index][j] * weight_array[0][j];
//            }
//        }
//        iso_it->second.total_reads = (int)iso_it->second.density;
//        iso_it->second.density = 1000 * iso_it->second.total_reads / (double)iso_it->second.length;
//    }

        for (auto corr_map_it = correspondence_map.begin(); corr_map_it != correspondence_map.end(); ++corr_map_it) {
            double density = 0;
            double total_reads = 0;
            for (int j = 0; j < weight_array[corr_map_it->second].size(); j++){
                density += weight_array[corr_map_it->second][j] * weight_array[0][j];
            }
            total_reads = (int)density;

            iso_map[corr_map_it->first].total_reads = total_reads;
            iso_map[corr_map_it->first].density = 1000 * total_reads / (double) iso_map[corr_map_it->first].length;;

        }

}




void calculate_rpkm (std::map <string, std::map <string, Isoform> > & iso_var_map, const long & aligned) {
    for (auto chrom_it = iso_var_map.begin(); chrom_it != iso_var_map.end(); ++chrom_it) {
        for (auto iso_it = chrom_it->second.begin(); iso_it !=  chrom_it->second.end(); ++iso_it) {
            iso_it->second.rpkm = iso_it->second.density / ( (double)aligned / (double)1000000 );
        }
    }
}