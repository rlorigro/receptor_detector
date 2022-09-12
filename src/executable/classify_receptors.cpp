#include "Filesystem.hpp"
#include "TsvReader.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Cluster.hpp"
#include "SvgPlot.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "edlib.h"

using receptor_detector::TsvReader;
using receptor_detector::Sequence;
using receptor_detector::Hasher2;
using receptor_detector::Cluster;
using receptor_detector::Timer;

#include <algorithm>
#include <iostream>
#include <utility>
#include <deque>

using std::numeric_limits;
using std::ostream;
using std::deque;
using std::tuple;
using std::sort;


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

    edlibFreeAlignResult(result);

    return {left_trim, right_trim};
}


int64_t get_edit_distance(const string& ref, const string& query){
    auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);

    EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), ref.c_str(), ref.size(), config);
    if (result.status != EDLIB_STATUS_OK) {
        throw runtime_error("ERROR: edlib alignment failed");
    }

    edlibFreeAlignResult(result);

    return result.editDistance;
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


void split_sequence(
        const Sequence& linker,
        const Sequence& s,
        size_t k,
        vector<Cluster>& top_two,
        vector <Sequence>& split_sequences
        ){

    for (auto& c: top_two) {

        // Pad the region if it isn't the expected length (will be trimmed later after local alignment)
        auto diff = linker.size() - (c.stop - c.start);
        if (diff > 0){
            c.start -= (diff + k);
            c.stop += (diff + k);
        }

        auto l = s.sequence.substr(c.start, c.stop - c.start);

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

    a.name = s.name + "_a";
    a.sequence = s.sequence.substr(a_start, a_length);

    b.name = s.name + "_b";
    b.sequence = s.sequence.substr(b_start, b_length);

    c.name = s.name + "_c";
    c.sequence = s.sequence.substr(c_start, c_length);

    // Update vector of split reads
    split_sequences.emplace_back(a);
    split_sequences.emplace_back(b);
    split_sequences.emplace_back(c);

//    cerr << s.sequence << '\n';
//    cerr << a.sequence << string(top_two[0].stop - top_two[0].start,'_') << b.sequence << string(top_two[1].stop - top_two[1].start,'_') << c.sequence << '\n';

}


class HashResult{
public:
    string other_name;
    int64_t edit_distance;
    int64_t n_hits;

    HashResult()=default;
    HashResult(string other_name, int64_t edit_distance);
};


HashResult::HashResult(
        string other_name,
        int64_t edit_distance):
    other_name(other_name),
    edit_distance(edit_distance),
    n_hits(0)
{}


void sort_results_by_read_id(
        const map<string, HashResult>& best_matches,
        vector <pair <string, HashResult> >& sorted_results
){

    for (const auto& [name, item]: best_matches){
        sorted_results.push_back({name, item});
    }

    sort(sorted_results.begin(), sorted_results.end(), [&](const pair <string, HashResult>& a, const pair <string, HashResult>& b){
        auto name_a = a.first;
        auto name_b = b.first;

        auto id_a = name_a.substr(4);
        auto id_b = name_b.substr(4);

        id_a = id_a.substr(0, id_a.find_last_of('_'));
        id_b = id_b.substr(0, id_b.find_last_of('_'));

        size_t i_a = stoi(id_a);
        size_t i_b = stoi(id_b);

        auto sub_id_a = name_a.back();
        auto sub_id_b = name_b.back();

        bool ordinal;
        if (i_a == i_b){
            ordinal = sub_id_a < sub_id_b;
        }
        else {
            ordinal = i_a < i_b;
        }

        return ordinal;
    });
}


void write_results_to_file(path output_directory, const map<string, HashResult>& best_matches){
    vector <pair <string, HashResult> > sorted_results;
    sort_results_by_read_id(best_matches, sorted_results);

    path best_matches_path = output_directory / "best_matches.tsv";
    ofstream best_matches_csv(best_matches_path);
    if (not (best_matches_csv.good() and best_matches_csv.is_open())){
        throw runtime_error("ERROR: couldn't write to file: " + best_matches_path.string());
    }

    for (auto& [name, item]: sorted_results){
        // Some short sequences have no hits, so just report a placeholder
        if (item.other_name.empty()){
            item.other_name = "unknown";
        }

        if (name.back() == 'a') {
            auto id = name.substr(0, name.find_last_of('_'));
            best_matches_csv << id << '\t' << item.other_name << ':';
        }
        else if (name.back() == 'b') {
            best_matches_csv << item.other_name << ':';
        }
        else if (name.back() == 'c') {
            best_matches_csv << item.other_name << '\n';
        }
    }
}


void classify(path tsv_ref_path, path tsv_reads_path, path output_directory, size_t n_threads){
    Timer t;

    create_directories(output_directory);

    Sequence linker;
    vector<Sequence> ref_sequences;
    vector<Sequence> read_sequences;

    TsvReader ref_reader(tsv_ref_path);
    TsvReader read_reader(tsv_reads_path);

    // This value is used later to tell how far apart linkers should be at a minimum
    int32_t min_length = numeric_limits<int32_t>::max();

    cerr << t << "Loading target sequences" << '\n';
    ref_reader.for_item_in_tsv([&](Sequence& s){
        if (s.name == "LINKER"){
            linker = s;
        }

        ref_sequences.emplace_back(s);

        if (s.size() < min_length){
            min_length = s.size();
        }
    });

    cerr << t << "Hashing target sequences" << '\n';

    // Do an exhaustive kmer comparison
    size_t k = 7;
    Hasher2 hasher(k, 0.9999999999999999, 1, n_threads);
    hasher.hash(ref_sequences);

    vector<double> match_probabilities;
    vector<Cluster> clusters;
    vector<Cluster> top_two;

    vector <Sequence> split_sequences;

    cerr << t << "Splitting reads" << '\n';

    // Iterate input reads, find linkers, and split into subreads
    size_t s_index = 0;
    read_reader.for_item_in_tsv([&](Sequence& s){
        hasher.classify_kmers(s, 0, "LINKER", match_probabilities);
        find_clusters(match_probabilities, k*2, clusters);
        select_top_two(min_length, clusters, top_two);
        split_sequence(linker, s, k, top_two, split_sequences);

        s_index++;
    });

    map<string, HashResult> best_matches;

    // Initialize best matches
    for (auto& s: split_sequences){
        best_matches[s.name] = {"",numeric_limits<int64_t>::max()};
    }

    // Dump all the sequences together for hashing again!
    split_sequences.insert(split_sequences.end(), ref_sequences.begin(), ref_sequences.end());

    // Build map from name to index
    unordered_map<string, size_t> name_to_sequence_index;
    for (size_t i=0; i < split_sequences.size(); i++){
        auto& name = split_sequences[i].name;
        name_to_sequence_index[name] = i;
    }

    cerr << t << "Hashing subreads" << '\n';

    // Do hashing
    k = 6;
    Hasher2 rehasher(k, 0.6, 8, n_threads);
    rehasher.hash(split_sequences);

    // Open a log for hashing results
    path hash_log_path = output_directory / "hash_results.csv";
    ofstream hash_log_csv(hash_log_path);
    if (not (hash_log_csv.good() and hash_log_csv.is_open())){
        throw runtime_error("ERROR: couldn't write to file: " + hash_log_path.string());
    }

    // Write header
    hash_log_csv << "a" << ',' << "b" << ',' << "hash_similarity" << ',' << "n_hashes" << ',' << "total_hashes" << ',' << "edit_distance" << ',' << "identity" << '\n';

    size_t max_hits = 6;

    cerr << t << "Computing alignments for subreads" << '\n';

    // Iterate hash results and compute full alignments for top n hits
    rehasher.for_each_overlap(100, 0,[&](const string& a, const string& b, int64_t n_hashes, int64_t total_hashes){
        // Skip the reference sequences
        if (a.substr(0,4) != "read"){
            return;
        }
        // Skip the read-to-read matches
        if (b.substr(0,4) == "read"){
            return;
        }

        // Pre-filled map, guaranteed entry exists
        auto& result = best_matches[a];

        // Stop considering results after a number of hits
        if (result.n_hits == max_hits){
            return;
        }

        auto& ref = split_sequences[name_to_sequence_index[a]].sequence;
        auto& query = split_sequences[name_to_sequence_index[b]].sequence;

        // Do global alignment on the substring, assuming linker splitting was good
        auto edit_distance = get_edit_distance(ref, query);

        hash_log_csv << a << ',' << b << ',' << double(n_hashes)/double(total_hashes) << ',' << n_hashes << ',' << total_hashes << ',' << edit_distance << ',' << double(edit_distance)/ref.size() << '\n';
        cerr << "prev: " << result.edit_distance << ' ' << edit_distance << ' ' << result.n_hits << '\n';

        if (edit_distance < result.edit_distance) {
            cerr << "better" << '\n';
            result.other_name = b;
            result.edit_distance = edit_distance;
        }

        result.n_hits++;
    });

    cerr << t << "Writing results to: "  << output_directory << '\n';

    write_results_to_file(output_directory, best_matches);

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path output_directory;
    path reference_path;
    path reads_path;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-o,--output_directory",
            output_directory,
            "Path to directory where output should be written")
            ->required();

    app.add_option(
            "-r,--reference",
            reference_path,
            "Path to tsv containing receptor and linker sequences")
            ->required();

    app.add_option(
            "-i,--input_reads",
            reads_path,
            "Path to tsv containing reads to be split and labeled")
            ->required();

    app.add_option(
            "-t,--n_threads",
            n_threads,
            "Maximum number of threads to use")
            ->required();

    CLI11_PARSE(app, argc, argv);

    classify(reference_path, reads_path, output_directory, n_threads);

    return 0;
}
