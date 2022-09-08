#include "Sequence.hpp"


namespace receptor_detector{


size_t Sequence::size() const{
    return sequence.size();
}


Sequence::Sequence(string& name, string& sequence):
    name(name),
    sequence(sequence)
{}


char get_reverse_complement(char c){
    if (c == 'A'){
        return 'T';
    }
    else if (c == 'C'){
        return 'G';
    }
    else if (c == 'G'){
        return 'C';
    }
    else if (c == 'T'){
        return 'A';
    }
    else if (c == 'N'){
        return 'N';
    }
    else {
        throw runtime_error("ERROR: uncomplementable character in sequence " + string(1,c));
    }
}


void get_reverse_complement(const string& fc, string& rc, size_t length){
    auto l = int64_t(length);

    for (int64_t i=l-1; i >= 0; i--){
        rc.append(1, get_reverse_complement(fc[i]));
    }
}


}
