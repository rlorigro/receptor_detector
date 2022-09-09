#include "Hasher2.hpp"


namespace receptor_detector{


Hasher2::Hasher2(size_t k, double sample_rate, size_t n_iterations, size_t n_threads):
        k(k),
        n_possible_bins(numeric_limits<uint64_t>::max()),
        n_iterations(n_iterations),
        total_sample_rate(sample_rate),
        iteration_sample_rate(sample_rate/double(n_iterations)),
        n_threads(n_threads)
{
    if (k > 32){
        throw runtime_error("ERROR: cannot perform robust 64bit hashing on kmer of length > 32");
    }

    n_bins = round(double(n_possible_bins)*iteration_sample_rate);

    cerr << "Using " << n_bins << " of " << n_possible_bins << " possible bins, for " << n_iterations
         << " iterations at a rate of " << iteration_sample_rate << '\n';
}


const vector<uint64_t> Hasher2::seeds = {
        2502369103967696952, 7135383713162725540, 13627014539775970274,
        4796741865460292034, 871560157616224954, 4803805340556337874,
        16419522238146181018, 2510544324051818521, 14730669708888401360,
        15771340054476377792, 9720085704046261736, 10056860032332123480,
        12402047041257800396, 15576524076748207768, 4544521901418279742,
        12511381447012767832, 8602441838744772705, 32484916947076159,
        13151201422946980157, 2509063378545863014, 14824159858063963905,
        885214842288208064, 15350874426397356478, 6497888167064186043,
        11929651735409705723, 6992204336913461076, 13414749134176657021,
        2625384403220687424, 1699601720342247628, 13322311195547716555,
        664080459378528451, 10579886813101962970, 2272107295286101060,
        6953264339271847999, 12499425564272938082, 8944540437674412670,
        15049164062555153936, 11492017226215766095, 8470764489325659499,
        1195888186965799934, 1451808535441337876, 15339137035898384211,
        5403968531279409728, 13269631949495200182, 370318925754183147,
        963549157246851466, 16406619304931695195, 9820707863630930290,
        9033896960103175519, 13675235632463118992, 2559561853484687358,
        2651434670380249965, 10995033789647094898, 514739756668963464,
        14659830264257792589, 3772038222307391067, 13640673517995500257,
        11203476974462366311, 2471601015170514172, 5059946906817196428,
        4640483839834766811, 4956910866326391215, 15382186762016539279,
        9006118601272042222, 17710828155857495220, 16814657398572709278,
        3676604234035115931, 9091606640466583270, 12871502770142896716,
        9421775905944046331, 12996574870094471825, 5233542693936021032,
        11739159484596970007, 8759703307818868101, 301590180242423745,
        8073837335570366087, 7136899178665934330, 14148922824375835145,
        1318395012090810332, 10251670999942663955, 14285987500822006731,
        17492932937437155077, 7102823587961144657, 10062324406689118391,
        3932036019385053742, 9321633453393523433, 15014189997857221125,
        4202779944578005508, 1715699542033256391, 6103254574080548671,
        1867805817184926607, 16558862004619463128, 4941896307182207749,
        768507597447242884, 4793176833258765644, 3559442299447860782,
        2913348424260394462, 5537057559772751678, 371726285994419264,
        10629356763957722439};


uint64_t Hasher2::hash(const BinarySequence<uint64_t>& kmer, size_t seed_index) const{
    return MurmurHash64A(kmer.sequence.data(), int(kmer.get_byte_length()), seeds[seed_index]);
}


///
/// \param sequence
/// \param i iteration of hashing to compute, corresponding to a hash function
void Hasher2::hash_sequence(const Sequence& sequence, size_t hash_index) {
    BinarySequence<uint64_t> kmer;

    // Forward iteration
    for (auto& c: sequence.sequence) {
        if (BinarySequence<uint64_t>::base_to_index.at(c) == 4){
            // Reset kmer and don't hash any region with non ACGT chars
            kmer = {};
            continue;
        }

        if (kmer.length < k) {
            kmer.push_back(c);
        } else {
            kmer.shift(c);
        }

        if (kmer.length == k){
            uint64_t h = hash(kmer, hash_index);

            if (h < n_bins){
                auto& m = bin_mutexes.at(h % bin_mutexes.size());
                auto& bin = bins.at(h % bins.size());

                m.lock();
                bin.emplace(sequence.name);
                m.unlock();
            }
        }
    }
}


///
/// \param sequence
/// \param i iteration of hashing to compute, corresponding to a hash function
void Hasher2::classify_kmers(const Sequence& sequence, size_t hash_index, const string& target_name, vector<double>& probabilities) {
    probabilities.clear();

    BinarySequence<uint64_t> kmer;

    // Forward iteration
    for (size_t i=0; i < sequence.size(); i++) {
        char c = sequence.sequence[i];

        double match_probability = 0;

        if (BinarySequence<uint64_t>::base_to_index.at(c) == 4){
            // Reset kmer and don't hash any region with non ACGT chars
            kmer = {};
            probabilities.emplace_back(match_probability);
            continue;
        }

        if (kmer.length < k) {
            kmer.push_back(c);
        } else {
            kmer.shift(c);
        }

        if (kmer.length == k){
            uint64_t h = hash(kmer, hash_index);

            if (h < n_bins){
                auto& bin = bins.at(h % bins.size());

                bool on_target = (bin.count(target_name) > 0);

                // Considering the previously hashed sequences as an exhaustive set, find the probability
                // that the current k-mer is a match, based on the number of on- vs off-target hits
                if (on_target) {
                    size_t total = bin.size();
                    match_probability = 1.0/double(total);
                }
            }
        }

        probabilities.emplace_back(match_probability);
    }
}



void Hasher2::write_hash_frequency_distribution() const{
    map <size_t, size_t> distribution;

    for (auto& b: bins){
        distribution[b.size()]++;
    }

    for (auto& [size, frequency]: distribution){
        cerr << size << '\t' << frequency << '\n';
    }
}


void Hasher2::hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index, const size_t hash_index){
    size_t i;
    while (job_index < sequences.size()){
        i = job_index.fetch_add(1);
        hash_sequence(sequences[i], hash_index);
    }
}


void Hasher2::hash(const vector<Sequence>& sequences){
    size_t max_kmers_in_sequence = 0;
    for (auto& sequence: sequences) {
        max_kmers_in_sequence += sequence.size();
    }

    cerr << max_kmers_in_sequence << " possible unique kmers in sequence" << '\n';

    // The maximum observable unique kmers is the total number of kmers in the sequence * the sample rate
    max_kmers_in_sequence = size_t(double(max_kmers_in_sequence) * total_sample_rate);

    cerr << max_kmers_in_sequence << " kmers after downsampling" << '\n';
    cerr << max_kmers_in_sequence * bins_scaling_factor << " bins allocated" << '\n';

    // Aggregate results
    for (size_t h=0; h<n_iterations; h++){
        cerr << "Beginning iteration: " << h << '\n';

        bins.clear();
        bins.resize(max_kmers_in_sequence * bins_scaling_factor);

        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        // Launch threads
        for (uint64_t t=0; t<n_threads; t++){
            try {
                threads.emplace_back(thread(
                        &Hasher2::hash_sequences,
                        this,
                        ref(sequences),
                        ref(job_index),
                        h
                ));
            } catch (const exception &e) {
                cerr << e.what() << "\n";
                exit(1);
            }
        }

        // Wait for threads to finish
        for (auto& t: threads){
            t.join();
        }

        // Iterate all hash bins for this iteration (unique hash function)
        for (auto& bin: bins){
            if (bin.size() > max_bin_size){
                continue;
            }

            vector <string> items(bin.size());

            size_t n = 0;
            for (auto& name: bin){
                items[n] = name;
                n++;
            }

            // Iterate all combinations of names found in this bin, including self hits, bc they'll be used as a
            // normalization denominator later.
            for (size_t a=0; a<items.size(); a++){
                for (size_t b=a; b<items.size(); b++){
                    overlaps[items[a]][items[b]]++;

                    // Only increment the reciprocal if it's not a self hit
                    if (a != b) {
                        overlaps[items[b]][items[a]]++;
                    }
                }
            }
        }
    }
}


void Hasher2::get_best_matches(map<string, string>& matches, double certainty_threshold) const{
    for (auto& [name, results]: overlaps){
        auto total_hashes = double(results.at(name));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        if (sorted_scores.empty()){
            continue;
        }

        string max_name = sorted_scores.rbegin()->second;
        auto max_hashes = double(sorted_scores.rbegin()->first);

//        double all_hashes = 0;
//
//        size_t i = 0;
//        // Report the top hits by % Jaccard similarity for each
//        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
//            auto score = iter->first;
//            auto other_name = iter->second;
//
//            // First item is the maximum score
//            if (i == 0){
//                max_hashes = double(score);
//                max_name = other_name;
//            }
//
//            all_hashes += double(score);
//
//            i++;
//
//            if (i == 10){
//                break;
//            }
//        }

        // Just take any top hit with greater than % threshold match, later will be used during symmetry filtering
        if (max_hashes/double(total_hashes) > certainty_threshold){
            matches[name] = max_name;
        }
    }
}


void Hasher2::get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold) const{
    map<string, string> matches;

    get_best_matches(matches, certainty_threshold);

    for (auto& [a,b]: matches){
        auto b_to_a = matches.find(b);

        if (b_to_a != matches.end()){
            if (b_to_a->second == a){
                // Only accept edges which are symmetrical
                symmetrical_matches.emplace(min(a,b), max(a,b));
            }
        }
    }
}


void Hasher2::for_each_overlap(const function<void(const pair<string,string>, int64_t weight)>& f) const{
    for (const auto& [a, result]: overlaps){
        for (const auto& [b, count]: result){
            if (a < b){
                f({a,b}, count);
            }
        }
    }
}


void Hasher2::for_each_overlap(
        size_t max_hits,
        double min_similarity,
        const function<void(const string& a, const string& b, int64_t n_hashes, int64_t total_hashes)>& f) const{

    for (auto& [name, results]: overlaps){
        // Self-hit is the total number of hashes the parent sequence had
        auto total_hashes = double(results.at(name));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        if (sorted_scores.empty()){
            continue;
        }

        string max_name = sorted_scores.rbegin()->second;
        auto max_hashes = double(sorted_scores.rbegin()->first);

        size_t i = 0;
        // Report the top hits by % similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            float similarity = float(score)/(float(total_hashes) + 1e-12f);

            if (similarity > min_similarity) {
                f(name, other_name, int64_t(score), int64_t(total_hashes));
            }

            i++;
            if (i == max_hits){
                break;
            }
        }
    }
}


void Hasher2::write_results(path output_directory) const{
    path overlaps_path = output_directory / "overlaps.csv";

    ofstream overlaps_file(overlaps_path);

    if (not overlaps_file.good() or not overlaps_file.is_open()){
        throw runtime_error("ERROR: could not write file: " + overlaps_path.string());
    }

    overlaps_file << "name" << ',' << "other_name" << ',' << "score" << ',' << "total_hashes" << ',' << "similarity" << '\n';

    for (auto& [name, results]: overlaps){
        size_t total_hashes = results.at(name);

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        size_t i = 0;

        // Report the top hits by % similarity for each (unidirectional, not Jaccard)
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            double similarity = double(score)/double(total_hashes);

            overlaps_file << name << ',' << other_name << ',' << score << ',' << total_hashes << ',' << similarity << '\n';
            i++;

            if (i == 10){
                break;
            }
        }
    }
}

}