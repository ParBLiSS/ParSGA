#include <iostream>
#include <iomanip>
#include "graph_util.hpp"
#include "csr_align.hpp"
#include "read_csr.hpp"

template<typename... Args>
void print_variables(std::ofstream& out, Args&&... args) {
    std::ostringstream oss;
    ((oss << std::forward<Args>(args) << ","), ...);
    oss << "\n";
    out << oss.str();
}

template <typename T>
void printVector(const std::string& name, const std::vector<T>& vec) {
    std::cout << name << ": ";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << static_cast<int64_t>(vec[i]);
        if (i < vec.size() - 1) {
            std::cout << " ";
        }
    }
        std::cout << '\n';
}

template <typename T>
bool areVectorsEqual(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    return (vec1.size() == vec2.size() && std::equal(vec1.begin(), vec1.end(), vec2.begin()));
}

int main(int argc, char **argv) {
    if (argc < 5) {
        throw std::runtime_error("Provided less than 4 arguments");
    }

    const bool match_cut = false;
    int match_score, subst_score, del_score, ins_score;
    match_score = 0;
    subst_score = 1;
    del_score = 1;
    ins_score = 1;
    // print with single-point precision
    std::cout << std::fixed << std::setprecision(6);
    double input_time = 0, serial_time = 0, parallel_time = 0, edge_parallel_time = 0, vertex_parallel_time = 0;

    int n, m;
    std::vector<char> labels;
    std::vector<int> offsets, edges, cut_offsets_out[ALPHABET_SIZE], cut_adjcny_out[ALPHABET_SIZE], component_indices[ALPHABET_SIZE], component_labels[ALPHABET_SIZE], thread_component_offsets[ALPHABET_SIZE];
    std::vector<std::string> reads;

    std::string input_read_length = argv[3];

    int const max_threads = omp_get_max_threads(); // remember what max threads are before setting num threads to 1
    int const num_threads = std::stoi(argv[4]);
    if (num_threads > max_threads) {
        throw std::runtime_error("Requested " + std::to_string(num_threads) + " on machine with " + std::to_string(max_threads) + " cores\n");
    }
    omp_set_num_threads(num_threads);
    std::cout << "set threads: " << omp_get_max_threads() << '\n' << std::flush;

    double input_start_time = omp_get_wtime();
    parallel_input(std::string(argv[1]), reads, n, m, labels, offsets, edges, cut_offsets_out, cut_adjcny_out, component_indices, component_labels, thread_component_offsets, input_read_length);
    input_time = omp_get_wtime() - input_start_time;
    std::cout << "input time: " << input_time << '\n' << std::flush;
    std::string output_filename(argv[2]);
    std::cout<<"\n";
    double avg_time[NUM_OPERATIONS] = {0.0};
    double stddev_time[NUM_OPERATIONS] = {0.0};
    double avg_imblance[NUM_OPERATIONS] = {0.0};
    double stddev_imblance[NUM_OPERATIONS] = {0.0};
    int const deletion_variant = 7;
    std::vector<uint8_t> parallel_score = parallel_align_CSR<uint8_t>(reads, n, m, labels, offsets, edges, cut_offsets_out, cut_adjcny_out, component_indices, component_labels, thread_component_offsets, match_score, subst_score, del_score, ins_score, avg_time, stddev_time, avg_imblance, stddev_imblance, deletion_variant);
    if(true) {
      std::cout << "AVERAGES: \n" << std::flush;
    }
    for(int op_num = 0; op_num < NUM_OPERATIONS; op_num++) {
        std::cout << static_cast<int64_t>(num_threads) << "-thread parallel " + operationToString(static_cast<Operation>(op_num)) + " time: " << avg_time[op_num] << " +/- " << stddev_time[op_num] << '\n';
        std::cout << static_cast<int64_t>(num_threads) << "-thread parallel " + operationToString(static_cast<Operation>(op_num)) + " imblance: " << avg_imblance[op_num] << " +/- " << stddev_imblance[op_num] << '\n';
    }
    printVector("score", parallel_score);
    return 0;
}