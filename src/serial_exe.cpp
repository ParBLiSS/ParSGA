#include <iostream>
#include <iomanip>
#include <numeric> // for accumulate
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

double imbalance(double const &max_element, std::vector<double>& elements) {
    return max_element / (std::accumulate(elements.begin(), elements.end(), 0.0) / elements.size());
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

int main(int argc, char **argv) {
    if (argc < 4) {
        throw std::runtime_error("Provided less than 3 arguments");
    }
    const bool match_cut = false;
    int match_score, subst_score, del_score, ins_score;
    match_score = 0;
    subst_score = 1;
    del_score = 1;
    ins_score = 1;
    // print with single-point precision
    std::cout << std::fixed << std::setprecision(6);
    double input_time = 0, serial_time = 0, preprocessing_time = 0, parallel_time = 0, edge_parallel_time = 0, vertex_parallel_time = 0;

    int n = 0;
    int m = 0;
    std::vector<char> labels;
    std::vector<int> offsets;
    std::vector<int> edges;
    std::vector<std::string> reads;

    std::string input_read_length = argv[3];
    double input_start_time = omp_get_wtime();
    serial_input(std::string(argv[1]), reads, n, m, labels, offsets, edges, input_read_length);
    input_time = omp_get_wtime() - input_start_time;
    std::cout << "input time: " << input_time << std::endl;
    std::string output_filename(argv[2]);

    int const max_threads = omp_get_max_threads();
    int const deletion_variant = 7;
    bool const edge_parallel = false;
    std::cout<<"deletion variant " << deletion_variant << ":\n";
    std::ofstream output;
    if(edge_parallel) {
        output.open(output_filename+"-"+std::to_string(deletion_variant)+"-ep.csv");
    } else {
        output.open(output_filename+"-"+std::to_string(deletion_variant)+"-vp.csv");
    }
    // want to use /dev/null as dummy output in tests
    /*
    if (!output.is_open()) {
        std::cerr << "Error opening the file: " << output_filename << std::endl;
        return 1;
    }
    */
    output << std::fixed << std::setprecision(6);
    output << "Threads,Total,Deletion,MatchMismatch,Insertion,DeletionImbalance,MatchMismatchImbalance,InsertionImbalance\n";
    //double deletion_time = 0, match_mismatch_time = 0, insertion_time = 0;
    double serial_start_time = omp_get_wtime();
    //double serial_score = serial_align_CSR<uint8_t>(reads, n, m, labels, offsets, edges, match_score, subst_score, del_score, ins_score, match_cut, deletion_time, match_mismatch_time, insertion_time);
    double serial_avg_time[NUM_OPERATIONS] = {0.0};
    double serial_stddev_time[NUM_OPERATIONS] = {0.0};
    std::vector<uint8_t> serial_score = serial_align_CSR<uint8_t>(reads, n, m, labels, offsets, edges, match_score, subst_score, del_score, ins_score, match_cut, serial_avg_time, serial_stddev_time);
    serial_time = omp_get_wtime() - serial_start_time;
    std::cout << "serial time: " << serial_time << std::endl;
    /*
    std::cout << "deletion time: " << deletion_time << std::endl;
    std::cout << "match/mismatch time: " << match_mismatch_time << std::endl;
    std::cout << "insertion time: " << insertion_time << std::endl;
    std::cout << "serial score: " << serial_score << std::endl;
    */
    printVector("serial score", serial_score);
    //print_variables(output, 1, serial_time, deletion_time, match_mismatch_time, insertion_time, 1, 1, 1);
    return 0;
}