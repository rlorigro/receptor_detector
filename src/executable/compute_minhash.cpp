#include "Filesystem.hpp"
#include "TsvReader.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "CLI11.hpp"

using receptor_detector::Hasher2;
using receptor_detector::Sequence;
using receptor_detector::TsvReader;


int compute_minhash(path tsv_path, path output_directory, double sample_rate, size_t k, size_t n_iterations, size_t n_threads){
    create_directories(output_directory);

    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);

    vector<Sequence> sequences;

    TsvReader reader(tsv_path);

    reader.for_item_in_tsv([&](Sequence& s){
        sequences.emplace_back(s);
    });

    hasher.hash(sequences);
    hasher.write_results(output_directory);

    return 0;
}


int main (int argc, char* argv[]){
    path file_path;
    path output_directory;
    double sample_rate = 0.1;
    size_t k = 22;
    size_t n_iterations = 10;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_tsv",
            file_path,
            "Path to GFA")
            ->required();

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

    compute_minhash(file_path, output_directory, sample_rate, k, n_iterations, n_threads);

    return 0;
}
