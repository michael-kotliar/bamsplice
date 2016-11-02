//
// Created by kot4or on 8/8/16.
//

#include "rpkm_calculation.h"



void print_weight_array(const vector<vector<double> > & weight_array, const string & title){
    cout << endl << title << endl;
    cout << std::setprecision(9) <<  std::fixed;
    for (int i = 0; i < weight_array.size(); i++) {
        cout <<  setw(4) << i;
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] == -1){
                cout << setw(12) << "-----";
            } else {
                cout << std::left << setw(12) << weight_array[i][j];
            }
        }
        cout << endl;
    }
}

void transform_to_density (vector <vector <double> > & weight_array){
    for (int i = 1; i < weight_array.size(); i++) {
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != -1){
                weight_array[i][j] = weight_array[i][j] / weight_array[0][j];
            }
        }
    }
}

double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row, int & count){
    assert ( (row >= 0 ) and (row < weight_array.size()) );
    double sum;
    count = 0;
    for (int j = 0; j < weight_array[row].size(); j++) {
        if (weight_array[row][j] != -1) {
            sum += weight_array[row][j];
            count++;
        }
    }
    return sum;
}

double get_average_by_row (const vector <vector <double> > & weight_array, const int & row){
    int count;
    double sum = get_sum_by_row(weight_array, row, count);
    return sum / count;
}

double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int & count, int start_row){
    assert ( start_row >= 0 );
    assert ( weight_array.size() > start_row );
    assert ( (column >= 0 ) and (column < weight_array[0].size()) );
    double sum = 0;
    count = 0;
    for (int j = start_row; j < weight_array.size(); j++) {
        if (weight_array[j][column] != -1) {
            sum += weight_array[j][column];
            count++;
        }
    }
//    cout << "debug sum = " << sum << endl;
    return sum;
}

double get_average_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row){
    int count;
    double sum = get_sum_by_column(weight_array, column, count, start_row);
    return sum / count;
}

vector <double> get_average_density_by_all_isoforms (const vector <vector <double> > & weight_array, int start_row){
    vector <double> average_density;
    cout << std::setprecision(9);
    for (int i = start_row; i < weight_array.size(); i++) {
        double average_by_row = get_average_by_row(weight_array, i);
        cout << "Average by row " << i <<" : " << average_by_row << endl;
        average_density.push_back( average_by_row );
    }
    return average_density;
}

vector <double> get_sum_density_by_all_intervals (const vector <vector <double> > & weight_array, int start_row){
    assert ( weight_array.size() > start_row );
    vector <double> sum_density;
    int count;
    for (int i = 0; i < weight_array[0].size(); i++) {
        double sum_dens = get_sum_by_column(weight_array, i, count, start_row);
//        cout << "DEBUG: " << sum_dens << endl;
        sum_density.push_back( sum_dens );
    }
    return sum_density;
}

void update_isoforms_density_to_average_for_isoform (vector <vector <double> > & weight_array){
    vector <double> average_densities = get_average_density_by_all_isoforms (weight_array);
//    cout << "DEBUG" << endl;
//    print_array (average_densities, "average_densities");
    for (int i = 1; i < weight_array.size(); i++) {
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != -1){
                weight_array[i][j] = average_densities[i-1];
            }
        }
    }
}

void adjust_isoforms_density_by_coef (vector <vector <double> > & weight_array, const vector <double> & original_densities_by_isoforms, const vector <double> & new_densities_by_isoforms){
    for (int i = 1; i < weight_array.size(); i++) {
//        cout << " k ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] != -1){
                double coef = original_densities_by_isoforms[j]/new_densities_by_isoforms[j];
//                cout << setw(8) << coef;
                weight_array[i][j] = weight_array[i][j] * coef;
            } else {
//                cout << setw(8) << "-----";
            }
        }
    }
}


void run_cycle (vector <vector <double> > & weight_array){
    // Get array of original densities sum
    vector <double> original_sum_dens =  get_sum_density_by_all_intervals (weight_array);
    print_array (original_sum_dens, "Original density sums");
    for (int i = 0; i < 2000; i++){
        // Update original density array with average values for isoforms
        cout << endl << "Set average by row" << endl;
        update_isoforms_density_to_average_for_isoform (weight_array);
//        print_weight_array(weight_array, "Average density array");

        cout << endl;

        // Get array of temporary densities sum
        vector <double> new_sum_dens =  get_sum_density_by_all_intervals (weight_array);
        print_array (original_sum_dens, "Original density sums");
        print_array (new_sum_dens, "New density sums");

        cout << endl << "Update according to the differenceds in original and new density sums" << endl;

        // Update density array according to the coef calculated by new density sums
        adjust_isoforms_density_by_coef (weight_array, original_sum_dens, new_sum_dens);
//        print_weight_array(weight_array, "Updated density array");
        cout << endl;
    }


}


void print_array (const vector <double> & intput_array, const string & title, string tab_spacer){
    cout << title << endl;
    cout << std::setprecision(9) <<  std::fixed;
    cout << tab_spacer;
    for (int i = 0; i < intput_array.size(); i++){
        cout << std::left << setw(12) << intput_array [i];
    }
    cout << endl;
}


