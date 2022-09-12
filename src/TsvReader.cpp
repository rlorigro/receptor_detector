#include "TsvReader.hpp"
#include "Sequence.hpp"

#include <fstream>

using std::ifstream;

namespace receptor_detector{

TsvReader::TsvReader(path tsv_path):
        tsv_path(tsv_path)
{}


void TsvReader::for_item_in_tsv(const function<void(Sequence& sequence)>& f) const{
    ifstream file(tsv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not open file: " + tsv_path.string());
    }

    string name;
    string sequence;

    char c;
    size_t n_delimiters = 0;
    size_t n_lines = 0;

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            Sequence s(name, sequence);

            f(s);

            name.clear();
            sequence.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                name += c;
            }
            else if (n_delimiters == 1){
                sequence += c;
            }
            else {
                throw runtime_error("ERROR: incorrect number of delimiters for sequence TSV file");
            }
        }
    }
}




void TsvReader::for_rle_item_in_tsv(const function<void(Sequence& sequence)>& f) const{
    ifstream file(tsv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not open file: " + tsv_path.string());
    }

    string name;
    string sequence;

    char c;
    size_t n_delimiters = 0;
    size_t n_lines = 0;

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            Sequence s(name, sequence);

            f(s);

            name.clear();
            sequence.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                name += c;
            }
            else if (n_delimiters == 1){
                if (not sequence.empty()){
                    if (c != sequence.back()){
                        sequence += c;
                    }
                }
                else {
                    sequence += c;
                }
            }
            else {
                throw runtime_error("ERROR: incorrect number of delimiters for sequence TSV file");
            }
        }
    }

}


}