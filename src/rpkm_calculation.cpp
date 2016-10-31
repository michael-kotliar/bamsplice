//
// Created by kot4or on 8/8/16.
//

#include "rpkm_calculation.h"
#include <iomanip>


void print_array(const vector<vector<double> > &weight_array){
    cout << endl << "WEIGHT ARRAY" << endl;
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

void update_to_density (vector <vector <double> > & weight_array){

}