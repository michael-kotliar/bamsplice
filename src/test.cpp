//
// Created by kot4or on 11/1/16.
//
#include "test.h"

void print_weight_array_test(const vector<vector<double> > & weight_array, const string & title){
    cerr << endl << title << endl;
    cerr << std::setprecision(3) <<  std::fixed;
    for (int i = 0; i < weight_array.size(); i++) {
        cerr <<  i <<") ";
        for (int j = 0; j < weight_array[i].size(); j++) {
            if (weight_array[i][j] == -1){
                cerr << setw(8) << "-----";
            } else {
                cerr << std::left << setw(8) << weight_array[i][j];
            }
        }
        cerr << endl;
    }
}