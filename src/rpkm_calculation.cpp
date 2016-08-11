//
// Created by kot4or on 8/8/16.
//

#include "rpkm_calculation.h"


void print_weight_array (const vector <vector <double> > & weight_array){
    cout << "WEIGHT ARRAY" << endl;
    cout << std::setprecision(3) <<  std::fixed;
    for (int i = 0; i < weight_array.size(); i++) {
        cout <<  i <<") ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            cout << weight_array[i][j] << "   ";
        }
        cout << endl;
    }
}