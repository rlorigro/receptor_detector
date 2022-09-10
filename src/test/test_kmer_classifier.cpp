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
using std::numeric_limits;


void plot_edlib_alignment(const string& ref, const string& query, SvgPlot& plot){
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);

    EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), ref.c_str(), ref.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    vector<size_t> ref_move = {1,0,1,1};
    vector<size_t> query_move = {1,1,0,1};

    for (size_t i=0; i<size_t(result.alignmentLength); i++) {
        cerr << int(result.alignment[i]);
    }
    cerr << '\n';

    int64_t left_trim = 0;
    int64_t right_trim = ref.size();

    // Find gaps provided by edlib local aligner. Choose longest alignment when multiple options exist.
    if (result.numLocations > 0){
        left_trim = result.startLocations[result.numLocations-1];
        right_trim = result.endLocations[result.numLocations-1];
    }

    plot.add_line(left_trim, int64_t(0), left_trim, int64_t(query.size()), 1, "gray");
    plot.add_line(right_trim, int64_t(0), right_trim, int64_t(query.size()), 1, "gray");

    cerr << "editDistance: " << result.editDistance << '\n';
    cerr << "left_trim: " << left_trim << '\n';
    cerr << "right_trim: " << right_trim << '\n';

    size_t prev_ref_index = left_trim;
    size_t prev_query_index = 0;
    size_t ref_index = prev_ref_index + ref_move[result.alignment[0]];
    size_t query_index = query_move[result.alignment[0]];

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


pair<int64_t, int64_t> get_alignment_bounds(const string& ref, const string& query){
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, nullptr, 0);

    EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), ref.c_str(), ref.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    int64_t left_trim = 0;
    int64_t right_trim = ref.size();

    // Find end gaps provided by local aligner. Choose longest alignment when multiple options exist.
    if (result.numLocations > 0){
        left_trim = result.startLocations[result.numLocations-1];
        right_trim = result.endLocations[result.numLocations-1];
    }

    cerr << "editDistance: " << result.editDistance << '\n';
    cerr << "left_trim: " << left_trim << '\n';
    cerr << "right_trim: " << right_trim << '\n';

    edlibFreeAlignResult(result);

    return {left_trim, right_trim};
}


int64_t get_edit_distance(const string& ref, const string& query){
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);

    EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), ref.c_str(), ref.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    cerr << "editDistance: " << result.editDistance << '\n';

    edlibFreeAlignResult(result);

    return result.editDistance;
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
            // Decay the cluster counter. When it reaches 0, a new cluster is made
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

    int32_t min_length = numeric_limits<int32_t>::max();
    ref_reader.for_item_in_tsv([&](Sequence& s){
        if (s.name == "LINKER"){
            linker = s;
        }

        ref_sequences.emplace_back(s);

        if (s.size() < min_length){
            min_length = s.size();
        }
    });

    Hasher2 hasher(7, 0.9999999999999999, 1, n_threads);
    hasher.hash(ref_sequences);

    vector<double> match_probabilities;
    vector<Cluster> clusters;
    vector<Cluster> top_two;

    vector <Sequence> split_sequences;
    string prefix = "split_read_";

    size_t s_index = 0;
    read_reader.for_item_in_tsv([&](Sequence& s){
        cerr << s.name << '\n';

        hasher.classify_kmers(s, 0, "LINKER", match_probabilities);
        find_clusters(match_probabilities, k*2, clusters);
        select_top_two(min_length, clusters, top_two);

        for (auto& c: top_two) {
            cerr << c << '\n';

            // Pad the region if it isn't the expected length
            auto diff = linker.size() - (c.stop - c.start);
            if (diff > 0){
                c.start -= (diff + k);
                c.stop += (diff + k);
            }

            auto l = s.sequence.substr(c.start, c.stop - c.start);

            cerr << l << '\n';
            cerr << linker.sequence << '\n';

//            path plot_path = output_directory / (s.name + "_vs_LINKER" + to_string(c.start) + '-' + to_string(c.stop) + ".svg");
//            SvgPlot plot(plot_path, 800, 800, 0, l.size(), 0, linker.size(), true);
//            plot_edlib_alignment(l, linker.sequence, plot);

            // Resolve the optimal start/stop using local alignment
            auto [left_trim, right_trim] = get_alignment_bounds(l, linker.sequence);

            // Update the bounds
            c.stop = c.start + right_trim;
            c.start += left_trim;
        }

        sort(top_two.begin(), top_two.end(), [&](const Cluster& a, const Cluster& b){
            return a.start < b.start;
        });

        for (auto& c: top_two) {
            cerr << c << '\n';
        }

        // Break out regions
        size_t a_start = 0;
        size_t a_stop = top_two[0].start;

        size_t b_start = top_two[0].stop;
        size_t b_stop = top_two[1].start;

        size_t c_start = top_two[1].stop;
        size_t c_stop = s.size();

        size_t a_length = a_stop - a_start;
        size_t b_length = b_stop - b_start;
        size_t c_length = c_stop - c_start;

        Sequence a;
        Sequence b;
        Sequence c;

        auto a_id = to_string(s_index) + "_a";
        a_id = string(6 - a_id.size(), '0') + a_id;
        a.name = prefix + a_id;
        a.sequence = s.sequence.substr(a_start, a_length);

        auto b_id = to_string(s_index) + "_b";
        b_id = string(6 - b_id.size(), '0') + b_id;
        b.name = prefix + b_id;
        b.sequence = s.sequence.substr(b_start, b_length);

        auto c_id = to_string(s_index) + "_c";
        c_id = string(6 - c_id.size(), '0') + c_id;
        c.name = prefix + c_id;
        c.sequence = s.sequence.substr(c_start, c_length);

        // Update vector of split reads
        split_sequences.emplace_back(a);
        split_sequences.emplace_back(b);
        split_sequences.emplace_back(c);

        cerr << s.sequence << '\n';
        cerr << a.sequence << string(top_two[0].stop - top_two[0].start,'_') << b.sequence << string(top_two[1].stop - top_two[1].start,'_') << c.sequence << '\n';

        s_index++;
    });

    unordered_map<string, size_t> name_to_sequence_index;
    for (size_t i=0; i < split_sequences.size(); i++){
        auto& name = split_sequences[i].name;
        name_to_sequence_index[name] = i;
    }

    // Dump all the sequences together for hashing again!
    split_sequences.insert(split_sequences.end(), ref_sequences.begin(), ref_sequences.end());

    Hasher2 rehasher(k, sample_rate, n_iterations, n_threads);
    rehasher.hash(split_sequences);

    size_t max_hits = 4;
    double min_similarity = 0;

    path hash_log_path = output_directory / "hash_results.csv";
    ofstream hash_log_csv(hash_log_path);
    if (not (hash_log_csv.good() and hash_log_csv.is_open())){
        throw runtime_error("ERROR: couldn't write to file: " + hash_log_path.string());
    }

    map<string, pair<string, int64_t> > best_matches;

    for (auto& [name, index]: name_to_sequence_index){
        best_matches[name] = {"",numeric_limits<int64_t>::max()};
    }

    rehasher.for_each_overlap(max_hits, min_similarity,[&](const string& a, const string& b, int64_t n_hashes, int64_t total_hashes){
        // Skip the reference sequences
        if (a.substr(0,11) != prefix){
            return;
        }
        // Skip the read-to-read matches
        if (b.substr(0,11) == prefix){
            return;
        }

        auto& ref = split_sequences[name_to_sequence_index[a]].sequence;
        auto& query = split_sequences[name_to_sequence_index[b]].sequence;

        auto edit_distance = get_edit_distance(ref, query);

        hash_log_csv << a << ',' << b << ',' << double(n_hashes)/double(total_hashes) << n_hashes << ',' << total_hashes << ',' << edit_distance << ',' << double(edit_distance)/ref.size() << '\n';

        cerr << a << ' ' << b << '\n';
        cerr << ref << '\n';
        cerr << query << '\n';
        cerr << '\n';

        // Pre-filled the map, guaranteed previous entry exists
        int64_t prev_edit_distance = best_matches.at(a).second;

        if (edit_distance < prev_edit_distance) {
            best_matches[a] = {b, edit_distance};
        }
    });

    path best_matches_path = output_directory / "best_matches.csv";
    ofstream best_matches_csv(best_matches_path);
    if (not (best_matches_csv.good() and best_matches_csv.is_open())){
        throw runtime_error("ERROR: couldn't write to file: " + best_matches_path.string());
    }

    string prev_name = "";
    vector<string> result;
    for (const auto& [name, item]: best_matches){
        auto& [other_name, score] = item;

        if (name.back() == 'a') {
            cerr << name << '\n';
            auto id = name.substr(0, name.find_last_of('_'));
            cerr << id << '\n';
            id = id.substr(id.find_last_of('_') + 1);
            cerr << id << '\n';
            id = to_string(stoi(id));
            cerr << id << '\n';

            cerr << "read" + id << '\t' << other_name << ':';
        }
        else if (name.back() == 'b') {
            cerr << other_name << ':';
        }
        else if (name.back() == 'c') {
            cerr << other_name << '\n';
        }
    }
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
