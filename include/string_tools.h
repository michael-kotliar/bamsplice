//
// Created by kot4or on 8/8/16.
//

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>

#ifndef TEST_1_STRING_TOOLS_H
#define TEST_1_STRING_TOOLS_H
#endif //TEST_1_STRING_TOOLS_H

using namespace std;

namespace string_tools {

    inline void remove_spaces(vector<string> &array) {
        vector<string> array_filtered;
        for (int i = 0; i < array.size(); i++)
            if (array[i].size() > 0 and array[i] != "\t") {
                array_filtered.push_back(array[i]);
            }
        array = array_filtered;
    }

    inline vector<string> split_line(const string &line, const string & key = "\t") {
        vector<string> line_splitted;
        boost::split(line_splitted, line, boost::is_any_of(key));
        remove_spaces(line_splitted);
        return line_splitted;
    }

    inline bool include_key(const string &line, const string &key_sequence) {
        vector<string> key_splitted = split_line(key_sequence);
        for (int i = 0; i < key_splitted.size(); i++) {
            if (line.find(key_splitted[i]) != std::string::npos) {
                return true;
            }
        }
        return false;
    }

    inline bool str_to_long_ptr(boost::shared_ptr<long> &ptr, const string &value){
        long temp;
        try {
            temp = boost::lexical_cast<long>(value);
        }
        catch(...){
            cout << "Bad lexical cast of " << value << " as long" << endl;
            return false;
        }
        ptr.reset (new long (temp));
        return true;
    }

    inline bool str_to_long(long &var, const string &value){
        long temp;
        try {
            temp = boost::lexical_cast<long>(value);
        }
        catch(...){
            cout << "Bad lexical cast of " << value << " as long" << endl;
            return false;
        }
        var = temp;
        return true;
    }

    inline bool str_to_int_ptr(boost::shared_ptr<int> &ptr, const string &value){
        int temp;
        try {
            temp = boost::lexical_cast<int>(value);
        }
        catch(...){
            cout << "Bad lexical cast of " << value << " as int" << endl;
            return false;
        }
        ptr.reset (new int (temp));
        return true;
    }

    inline bool str_to_int(int &var, const string &value){
        int temp;
        try {
            temp = boost::lexical_cast<int>(value);
        }
        catch(...){
            cout << "Bad lexical cast of " << value << " as int" << endl;
            return false;
        }
        var = temp;
        return true;
    }

    inline bool str_array_to_long_ptr_array(const vector<string> &input, vector<boost::shared_ptr<long> > &output){
        for (int i = 0; i < input.size(); i++){
            boost::shared_ptr <long> temp;
            if (not str_to_long_ptr(temp, input[i])){
                return false;
            }
            output.push_back(temp);
        }
        return true;
    }


    inline bool str_array_to_long_array(const vector<string> &input, vector<long> &output){
        for (int i = 0; i < input.size(); i++){
            long temp;
            try {
                temp = boost::lexical_cast<long>(input[i]);
            }
            catch(...){
                cout << "Bad lexical cast of " << input[i] << " as long" << endl;
                return false;
            }
            output.push_back(temp);
        }
        return true;
    }

    inline void print_vector (const vector <string> & in, string title = ""){
        cout << title << endl;
        for (int i = 0; i < in.size(); i++){
            cout << i <<") " << in[i] << endl;
        }
    }

}