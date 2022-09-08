#pragma once

#include "Filesystem.hpp"
#include "Sequence.hpp"
#include <functional>

using ghc::filesystem::path;
using std::function;
using receptor_detector::Sequence;

namespace receptor_detector {

class TsvReader {
public:
    path tsv_path;

    TsvReader()=default;
    TsvReader(path tsv_path);
    void for_item_in_tsv(const function<void(Sequence& sequence)>& f) const;
    void for_rle_item_in_tsv(const function<void(Sequence& sequence)>& f) const;
};

}