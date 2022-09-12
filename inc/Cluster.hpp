#pragma once

#include <iostream>

using std::ostream;

namespace receptor_detector{

///
/// Very simple chaining mechanism for sequential data. Purely based on a distance threshold.
///
class Cluster{
public:
    double mass;
    double max;
    int32_t c;
    int32_t start;
    int32_t stop;
    int32_t counter;

    Cluster(int32_t k);
    void update(int32_t i, double x);
    void decrement();
};

}

ostream& operator<<(ostream& o, const receptor_detector::Cluster& c);

