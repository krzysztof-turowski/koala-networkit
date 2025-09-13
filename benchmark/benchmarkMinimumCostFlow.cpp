#include <cassert>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>

#include <flow/MinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/EdmondsKarpMCF.hpp>
#include <flow/minimum_cost_flow/OrlinMCF.hpp>
#include <flow/minimum_cost_flow/SuccessiveApproxMCC.hpp>
#include <io/DimacsGraphReader.hpp>

template <typename FlowAlgorithm>
void run_mcf_algorithm(const std::string &file_path, const std::string &name) {
    auto [G, costs, b] = Koala::DimacsGraphReader().read_minimum_cost_flow(file_path);
    auto start = std::chrono::high_resolution_clock::now();
    auto minimum_cost_flow = FlowAlgorithm(G, costs, b);
    minimum_cost_flow.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << name << ": Minimum cost flow flow = " << minimum_cost_flow.getMinCost()
              << ", Time taken = " << elapsed.count() << " seconds\n";
}

std::map<std::string, int> ALGORITHM = {
    { "EK", 1 },
    { "O", 2 },
    { "SA", 3 },
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
        run_mcf_algorithm<Koala::EdmondsKarpMCF>(file_path, "EdmondsKarp");
    } else if (ALGORITHM[algorithm] == 2) {
        run_mcf_algorithm<Koala::OrlinMCF>(file_path, "Orlin");
    } else if (ALGORITHM[algorithm] == 3) {
        run_mcf_algorithm<Koala::SuccessiveApproxMCC>(file_path, "SuccessiveApproximation");
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    return 0;
}
