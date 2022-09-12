#include "Cluster.hpp"

namespace receptor_detector{

Cluster::Cluster(int32_t c):
        mass(0),
        max(0),
        c(c),
        start(0),
        stop(0),
        counter(c)
{}


void Cluster::update(int32_t i, double x){
    if (mass == 0){
        start = i;
    }

    if (x > max){
        max = x;
    }

    mass += x;
    stop = i;
    counter = c;
}


void Cluster::decrement(){
    counter--;
}

}

ostream& operator<<(ostream& o, const receptor_detector::Cluster& c){
    o << c.mass << ',' << c.max << ',' << c.start << ',' << c.stop << ',' << c.stop - c.start << ',' << c.counter;

    return o;
}
