#include <iostream>
#include <iomanip>
#include "preprocessing.hpp"
#include "csr_align.hpp"
#include <random> // for std::mt19937 and std::uniform_int_distribution

template<typename... Args>
void print_variables(std::ofstream& out, Args&&... args) {
    std::ostringstream oss;
    ((oss << std::forward<Args>(args) << ","), ...);
    oss << "\n";
    out << oss.str();
}

template<typename Func, typename... Args>
void my_omp_time(Func func, double& timing, Args&&... args) {
  double thread_start_time = omp_get_wtime();
  func(std::forward<Args>(args)...);
  timing += (omp_get_wtime() - thread_start_time);
}

double imbalance(double const &max_element, std::vector<double>& element_vector) {
    return max_element / (accumulate(element_vector.begin(), element_vector.end(), 0.0) / element_vector.size());
}

int main(int argc, char **argv) {
    if (argc < 3) {
        throw std::runtime_error("Provided less than 2 arguments");
    }
    bool match_cut;
    int match_score, subst_score, del_score, ins_score;
    double graph_loading_time;

    double start_time = omp_get_wtime();
    psgl::graphLoader g;
    std::string graph_filename = argv[1];
    if(graph_filename.substr(graph_filename.length()-3)==".vg") {
        g.loadFromVG(graph_filename);
    } else if (graph_filename.substr(graph_filename.length()-4)==".txt") {
        g.loadFromTxt(graph_filename);
    } else {
        std::cerr << "error: graph filename should end with .vg or .txt\n" <<std::endl;
    }
    graph_loading_time = omp_get_wtime() - start_time;
    std::cout << "graph-loading time (s):\t" << graph_loading_time << '\n';
    std::string output_filename(argv[2]);

    int const max_threads = omp_get_max_threads();
    bool const edge_parallel = false;
    std::ofstream output;
    output.open(output_filename+".csr");
    if (!output.is_open()) {
        std::cerr << "Error opening the file: " << output_filename << std::endl;
        return 1;
    }
    double write_graph_code_time;
    auto write_graph_code = [&]() {
        output << g.diCharGraph.numVertices << '\n';
        output << g.diCharGraph.numEdges << '\n';
        for(int i = 0; i < g.diCharGraph.numVertices; i++) {
            output << g.diCharGraph.vertex_label[i] << " "; 
        }
        output << "\n";
        for(int i = 0; i < g.diCharGraph.numVertices+1; i++) {
            output << g.diCharGraph.offsets_out[i] << " "; 
        }
        output << "\n";
        for(int i = 0; i < g.diCharGraph.numEdges; i++) {
            output << g.diCharGraph.adjcny_out[i] << " "; 
        }
        output << "\n";
    };
    my_omp_time(write_graph_code, write_graph_code_time);

    std::vector<double> match_cut_timing(omp_get_max_threads(), 0);
    double components_timing = 0;
    parlay::sequence<index_label_tuple> *tuples;
    double match_cut_time = 0;
    int **thread_component_offsets;
    double preprocessing_time = 0;
    pvector<int64_t> char_offset_out[ALPHABET_SIZE];
    int32_t* char_adj_out[ALPHABET_SIZE];
    Graph* m;
    auto preprocessing_code = [&]() {
        double thread_start_time = omp_get_wtime();
        match_cut_to_offsets_neighs(g.diCharGraph.numVertices, g.diCharGraph.numEdges, g.diCharGraph.adjcny_out, g.diCharGraph.offsets_out, g.diCharGraph.vertex_label, edge_parallel, match_cut_timing, char_offset_out, char_adj_out); // i think i can remove this part, i think it was leftover from debugging since its called inside match_cut_graphs... no, need to fix this, i need to keep access to char_offset_out and char_adj_out
        match_cut_timing = std::vector<double>(omp_get_max_threads(), 0);
        m = match_cut_graphs(g.diCharGraph.numVertices, g.diCharGraph.numEdges, g.diCharGraph.adjcny_out, g.diCharGraph.offsets_out, g.diCharGraph.vertex_label, edge_parallel, match_cut_timing);
        match_cut_time = omp_get_wtime() - thread_start_time;
        thread_start_time = omp_get_wtime();
        label_and_sort_components(m,tuples,thread_component_offsets);
        components_timing = omp_get_wtime() - thread_start_time;
    };
    my_omp_time(preprocessing_code, preprocessing_time);
    std::ofstream cut_output;
    std::ofstream component_output;
    for(char base : bases) {
        cut_output.open(output_filename+"-"+base+".csr");
        component_output.open(output_filename+"-"+base+".components");
        cut_output << g.diCharGraph.numVertices << '\n';
        cut_output << char_offset_out[base_index(base)][g.diCharGraph.numVertices]  << '\n';
        for(int i = 0; i < g.diCharGraph.numVertices+1; i++) {
            cut_output << char_offset_out[base_index(base)][i] << " "; 
        }
        cut_output << "\n";
        for(int i = 0; i < char_offset_out[base_index(base)][g.diCharGraph.numVertices]; i++) {
            cut_output << char_adj_out[base_index(base)][i] << " "; 
        }
        cut_output << "\n";
        cut_output.close();
        for(int i = 0; i < tuples[base_index(base)].size(); i++) {
            component_output << tuples[base_index(base)][i].index << " ";
        }
        component_output << '\n';
        for(int i = 0; i < tuples[base_index(base)].size(); i++) {
            component_output << tuples[base_index(base)][i].label << " ";
        }
        component_output << '\n';
        component_output.close();
        int const max_threads = omp_get_max_threads();
        for (int num_threads = 1; num_threads <= max_threads; num_threads = ((num_threads <= max_threads/ 2) || (num_threads == max_threads)) ? (num_threads *= 2) : max_threads) {
            std::ofstream thread_component_offsets_output;
            thread_component_offsets_output.open(output_filename+"-"+base+"-"+std::to_string(num_threads)+".components");
            omp_set_num_threads(num_threads);
            label_and_sort_components(m,tuples,thread_component_offsets);
            for(int i = 0; i < num_threads+1; i++) {
                thread_component_offsets_output << thread_component_offsets[base_index(base)][i] << " ";
            }
            thread_component_offsets_output << '\n';
            thread_component_offsets_output.close();
        }
    }
    delete[] m;
    delete[] tuples;
    for(char base : bases) {
        delete[] thread_component_offsets[base_index(base)];
    }
    delete[] thread_component_offsets;
    double max_cut_time = *std::max_element(match_cut_timing.begin(), match_cut_timing.end());
    double cut_imbalance = imbalance(max_cut_time, match_cut_timing);
    return 0;
}
