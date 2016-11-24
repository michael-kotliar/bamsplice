//
// Created by kot4or on 11/1/16.
//
#include "test.h"

void print_weight_array_test( vector <vector <double> > weight_array, const string & title, const string & path){
    ofstream output_stream;
    output_stream.open (path, std::ofstream::app);
    if (output_stream.is_open()){
        output_stream << endl << title << endl;
        output_stream << std::setprecision(3) <<  std::fixed;
        for (int i = 0; i < weight_array.size(); i++) {
            output_stream <<  i <<") ";
            for (int j = 0; j < weight_array[i].size(); j++) {
                if (weight_array[i][j] == 0){
                    output_stream << setw(8) << "-----";
                } else {
                    output_stream << std::left << setw(8) << weight_array[i][j];
                }
            }
            output_stream << endl;
        }
        output_stream.close();
    }
}