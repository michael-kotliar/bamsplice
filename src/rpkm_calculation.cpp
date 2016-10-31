//
// Created by kot4or on 8/8/16.
//

#include "rpkm_calculation.h"



void print_weight_array(const vector<vector<double> > & weight_array, const string & title){
    cout << endl << title << endl;
    cout << std::setprecision(3) <<  std::fixed;
    for (int i = 0; i < weight_array.size(); i++) {
        cout <<  i <<") ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] == -1){
                cout << setw(8) << "-----";
            } else {
                cout << std::left << setw(8) << weight_array[i][j];
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

double get_sum_by_row (const vector <vector <double> > & weight_array, const int & row){
    assert ( (row >= 0 ) and (row < weight_array.size()) );
    double sum;
    for (int j = 0; j < weight_array[row].size(); j++) {
        if (weight_array[row][j] != -1) {
            sum += weight_array[row][j];
        }
    }
    return sum;
}

double get_average_by_row (const vector <vector <double> > & weight_array, const int & row){
    return get_sum_by_row(weight_array, row) / weight_array[row].size();
}

double get_sum_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row){
    assert ( start_row >= 0);
    assert ( weight_array.size() > start_row );
    assert ( (column >= 0 ) and (column < weight_array[0].size()) );
    double sum;
    for (int j = start_row; j < weight_array.size(); j++) {
        if (weight_array[j][column] != -1) {
            sum += weight_array[j][column];
        }
    }
    return sum;
}

double get_average_by_column (const vector <vector <double> > & weight_array, const int & column, int start_row){
    return get_sum_by_column(weight_array, column, start_row) / weight_array.size();
}