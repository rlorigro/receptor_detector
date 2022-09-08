#pragma once

#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::chrono::duration;

using std::chrono::hours;
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;

using std::ostream;
using std::string;


namespace receptor_detector{

class Timer{
    steady_clock::time_point start;

public:
    Timer();
    string elapsed() const;
    void reset();
};

}

ostream& operator<<(ostream& o, const receptor_detector::Timer& t);
