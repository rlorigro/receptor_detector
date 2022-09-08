#pragma once

#include <stdexcept>
#include <string>

using std::runtime_error;
using std::string;


namespace receptor_detector{

class Sequence {
public:
    string name;
    string sequence;

    Sequence()=default;
    Sequence(string& name, string& sequence);
    size_t size() const;
};


char get_reverse_complement(char c);

void get_reverse_complement(const string& fc, string& rc, size_t length);

}

