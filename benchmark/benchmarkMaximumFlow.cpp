#include <cassert>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>

#include <flow/MaximumFlow.hpp>
#include <flow/PushRelabel.hpp>
#include <flow/BoykovKolmogorovFlow.hpp>
#include <flow/KingRaoTarjanMaximumFlow.hpp>
#include <flow/MalhotraKumarMaheshwariFlow.hpp>
#include <io/DimacsGraphReader.hpp>

using namespace std::chrono;

template <typename FlowAlgorithm>
void run_flow_algorithm(const std::string &file_path, const std::string &name) {
    auto [G, s, t] = Koala::DimacsGraphReader().read_all(file_path);
    auto start = high_resolution_clock::now();
    auto maximum_flow = FlowAlgorithm(G, s, t);
    maximum_flow.run();
    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    std::cout << name << ": Maximum flow = " << maximum_flow.getFlowSize()
              << ", Time taken = " << elapsed.count() << " seconds\n";
}

std::map<std::string, int> ALGORITHM = {
    { "PushRelabel", 1 },
    { "BK", 2 },
    { "MKM", 3 },
    { "KRT", 4 }
};

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <algorithm> <file>" << std::endl;
        return 1;
    }

    std::string algorithm(argv[1]);
    std::string file_path(argv[2]);

    if (!std::filesystem::exists(file_path)) {
        std::cerr << "File " << file_path << " does not exist" << std::endl;
        return 1;
    }
    if (std::filesystem::is_directory(file_path)) {
        std::cerr << file_path << " is a directory" << std::endl;
        return 1;
    }

    std::cout << "\nProcessing file: " << file_path << std::endl;

    if (ALGORITHM[algorithm] == 1) {
        run_flow_algorithm<Koala::PushRelabel>(file_path, "PushRelabel");
    } else if (ALGORITHM[algorithm] == 2) {
        run_flow_algorithm<Koala::BoykovKolmogorovFlow>(file_path, "BKFlow");
    } else if (ALGORITHM[algorithm] == 3) {
        run_flow_algorithm<Koala::MalhotraKumarMaheshwariFlow>(file_path, "MKMFlow");
    } else if (ALGORITHM[algorithm] == 4) {
        run_flow_algorithm<Koala::KingRaoTarjanMaximumFlow>(file_path, "KingRaoTarjan");
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    return 0;
}
