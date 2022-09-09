#include "SvgPlot.hpp"
#include "Filesystem.hpp"
#include "TsvReader.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "CLI11.hpp"
#include "edlib.h"

using receptor_detector::Hasher2;
using receptor_detector::Sequence;
using receptor_detector::TsvReader;

#include <algorithm>
#include <iostream>
#include <deque>

using std::ostream;
using std::sort;
using std::deque;


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

    size_t kernel_size = 10;

    int64_t left_trim = 0;
    deque<int16_t> d;
    for (size_t i=0; i<size_t(result.alignmentLength/2); i++) {
        left_trim += ref_move[result.alignment[i]];

        d.push_back(result.alignment[i]);

        if (d.size() > kernel_size){
            d.pop_front();
        }

//        cerr << i << '\n';
        size_t score = 0;
        for (auto x: d){
//            cerr << x << ' ';
            score += (x == 0);
        }
//        cerr << '\n';

        double identity = double(score) / double(kernel_size);

//        cerr << identity << '\n';

        if (identity > 0.7){
            break;
        }
    }

    int64_t right_trim = ref.size();
    d.clear();
    for (size_t i=size_t(result.alignmentLength-1); i>size_t(result.alignmentLength/2); i--) {
        right_trim -= ref_move[result.alignment[i]];

        d.push_back(result.alignment[i]);

        if (d.size() > kernel_size){
            d.pop_front();
        }

//        cerr << i << '\n';
        size_t score = 0;
        for (auto x: d){
//            cerr << x << ' ';
            score += (x == 0);
        }
//        cerr << '\n';

        double identity = double(score) / double(kernel_size);

//        cerr << identity << '\n';

        if (identity > 0.6){
            break;
        }
    }

    left_trim = max(int64_t(0), int64_t(left_trim - kernel_size));
    right_trim = min(int64_t(ref.size()), int64_t(right_trim + kernel_size));

    plot.add_line(left_trim, int64_t(0), left_trim, int64_t(query.size()), 1, "gray");
    plot.add_line(right_trim, int64_t(0), right_trim, int64_t(query.size()), 1, "gray");

    cerr << "editDistance: " << result.editDistance << '\n';
    cerr << "left_trim: " << left_trim << '\n';
    cerr << "right_trim: " << right_trim << '\n';

    string color = "orange";
    for (size_t i=1; i<size_t(result.alignmentLength); i++){
        if (i > 0){
            if (result.alignment[i-1] != result.alignment[i]) {
                plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, 1, color);

                prev_query_index = query_index;
                prev_ref_index = ref_index;
            }
        }

        query_index += query_move[result.alignment[i]];
        ref_index += ref_move[result.alignment[i]];
    }

    plot.add_line(prev_ref_index, prev_query_index, ref_index, query_index, 1, "orange");

    edlibFreeAlignResult(result);
}


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


Cluster::Cluster(int32_t c):
        mass(0),
        max(0),
        c(c),
        start(0),
        stop(0),
        counter(c)
{}


ostream& operator<<(ostream& o, const Cluster& c){
    o << c.mass << ',' << c.max << ',' << c.start << ',' << c.stop << ',' << c.stop - c.start << ',' << c.counter;

    return o;
}


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


void find_clusters(vector<double> x, int32_t distance, vector<Cluster>& clusters){
    clusters.clear();

    for (size_t i=0; i<x.size(); i++){
        if (x[i] > 0){
            // Make a new cluster if the old one expired or there are none
            if (not clusters.empty()){
                if (clusters.back().counter < 0){
                    clusters.emplace_back(distance);
                }
            }
            else{
                clusters.emplace_back(distance);
            }

            clusters.back().update(i, x[i]);
        }
        else{
            if (not clusters.empty()){
                clusters.back().decrement();
            }
        }

//        cerr << i << ' ' << x[i] << ' ';
//        cerr << '\n';
    }
}


int32_t get_cluster_distance(const Cluster& a, const Cluster& b){
    if (a.start > b.start){
        return a.start - b.stop;
    }
    else {
        return b.start - a.stop;
    }
}


void select_top_two(int32_t min_distance, vector<Cluster>& clusters, vector<Cluster>& top_two){
    top_two.clear();

    sort(clusters.begin(), clusters.end(), [&](const Cluster& a, const Cluster& b){
        return a.mass > b.mass;
    });

    for (auto& c: clusters) {
        cerr << c << '\n';
    }
    cerr << '\n';

    for (auto& c: clusters){
        if (top_two.size() == 2) {
            break;
        }

        if (top_two.empty()){
            top_two.emplace_back(c);
        }

        else if (get_cluster_distance(c, top_two.back()) > min_distance){
            top_two.emplace_back(c);
        }
    }
}


void classify(path output_directory, double sample_rate, size_t k, size_t n_iterations, size_t n_threads){
    create_directories(output_directory);

    path this_path = __FILE__;
    path project_directory = this_path.parent_path().parent_path().parent_path();

    path tsv_ref_path = project_directory / "data" / "linker_and_receptor_segments.tsv";
    path tsv_reads_path = project_directory / "data" / "unlabelled_reads.tsv";

    Sequence linker;
    vector<Sequence> ref_sequences;
    vector<Sequence> read_sequences;

    TsvReader ref_reader(tsv_ref_path);
    TsvReader read_reader(tsv_reads_path);

    int32_t min_length = 999999;
    ref_reader.for_item_in_tsv([&](Sequence& s){
        if (s.name == "LINKER"){
            linker = s;
        }

        ref_sequences.emplace_back(s);

        if (s.size() < min_length){
            min_length = s.size();
        }
    });

    cerr << min_length << '\n';

    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);
    hasher.hash(ref_sequences);

    vector<double> match_probabilities;
    vector<Cluster> clusters;
    vector<Cluster> top_two;

    read_reader.for_item_in_tsv([&](Sequence& s){
        cerr << s.name << '\n';

        hasher.classify_kmers(s, 0, "LINKER", match_probabilities);
        find_clusters(match_probabilities, k*2, clusters);
        select_top_two(min_length, clusters, top_two);

        for (auto& c: top_two) {
            cerr << c << '\n';
        }
        cerr << '\n';

        if (top_two.size() < 2) {
            cerr << "WARNING" << '\n';
        }

        for (auto c: top_two){
            auto diff = linker.size() - (c.stop - c.start);
            if (diff > 0){
                c.start -= (diff + k);
                c.stop += (diff + k);
            }

            auto l = s.sequence.substr(c.start, c.stop - c.start);

            cerr << l << '\n';
            cerr << linker.sequence << '\n';

            path plot_path = output_directory / (s.name + "_vs_LINKER" + to_string(c.start) + '-' + to_string(c.stop) + ".svg");
            SvgPlot plot(plot_path, 800, 800, 0, l.size(), 0, linker.size(), true);
            plot_edlib_alignment(l, linker.sequence, plot);


        }
    });
}


int main (int argc, char* argv[]){
    path output_directory;
    double sample_rate = 0.1;
    size_t k = 22;
    size_t n_iterations = 10;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-o,--output_directory",
            output_directory,
            "Path to directory where output should be written")
            ->required();

    app.add_option(
            "-r,--sample_rate",
            sample_rate,
            "Sample rate. Proportion [0-1] of k-mers to retain during in comparison")
            ->required();

    app.add_option(
            "-k,--kmer_length",
            k,
            "Length of k-mer to use for hashing")
            ->required();

    app.add_option(
            "-n,--n_iterations",
            n_iterations,
            "Number of iterations (different hash functions), each at a rate of sample_rate/n_iterations")
            ->required();

    app.add_option(
            "-t,--n_threads",
            n_threads,
            "Maximum number of threads to use")
            ->required();

    CLI11_PARSE(app, argc, argv);

    if (sample_rate < 0 or sample_rate > 1){
        throw std::runtime_error("ERROR: sample rate must be between 0 and 1.0");
    }

    classify(output_directory, sample_rate, k, n_iterations, n_threads);

    return 0;
}
