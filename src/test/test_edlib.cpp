#include "Filesystem.hpp"
#include "TsvReader.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "SvgPlot.hpp"
#include "CLI11.hpp"
#include "edlib.h"

using receptor_detector::Hasher2;
using receptor_detector::Sequence;
using receptor_detector::TsvReader;


void plot_edlib_alignment(const string& ref, const string& query, SvgPlot& plot){
    auto config = edlibNewAlignConfig(5000, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);

    EdlibAlignResult result = edlibAlign(ref.c_str(), ref.size(), query.c_str(), query.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    vector<size_t> ref_move = {1,1,0,1};
    vector<size_t> query_move = {1,0,1,1};

    size_t prev_ref_index = 0;
    size_t prev_query_index = 0;
    size_t ref_index = ref_move[result.alignment[0]];
    size_t query_index = query_move[result.alignment[0]];

    for (size_t i=0; i<size_t(result.alignmentLength); i++) {
        cerr << int(result.alignment[i]);
    }
    cerr << '\n';

    cerr << ref << '\n';
    cerr << query << '\n';

    for (size_t i=1; i<size_t(result.alignmentLength); i++){
        if (i > 0){
            if (result.alignment[i-1] != result.alignment[i]) {
                plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, ref.size()/200, "orange");

                prev_query_index = query_index;
                prev_ref_index = ref_index;
            }
        }

        query_index += query_move[result.alignment[i]];
        ref_index += ref_move[result.alignment[i]];
    }

    plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, ref.size()/200, "orange");

    edlibFreeAlignResult(result);
}


void align(){
    path this_path = __FILE__;
    path project_directory = this_path.parent_path().parent_path().parent_path();

    path tsv_ref_path = project_directory / "data" / "linker_and_receptor_segments.tsv";
    path tsv_reads_path = project_directory / "data" / "unlabelled_reads.tsv";

    Sequence linker;
    vector<Sequence> sequences;

    TsvReader ref_reader(tsv_ref_path);
    TsvReader read_reader(tsv_reads_path);

    ref_reader.for_item_in_tsv([&](Sequence& s){
        if (s.name == "LINKER"){
            linker = s;
        }
    });

    read_reader.for_item_in_tsv([&](Sequence& s){
        sequences.emplace_back(s);
    });

    path output_dir = "test";
    create_directories(output_dir);

    for (const auto& s: sequences){
        path plot_path = output_dir / (s.name + "_vs_LINKER.svg");

        cerr << plot_path << '\n';
        SvgPlot plot(plot_path, 800, 800, 0, s.size(), 0, linker.size(), true);
        plot_edlib_alignment(s.sequence, linker.sequence, plot);
    }

}


int main (int argc, char* argv[]){
    align();

    return 0;
}
