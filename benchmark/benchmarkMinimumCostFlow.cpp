#include <cassert>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include <flow/MinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/EdmondsKarpMinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/OrlinMinimumCostFlow.hpp>
#include <flow/minimum_cost_flow/SuccessiveApproximationMinimumCostFlow.hpp>
#include <io/DimacsGraphReader.hpp>

template <typename FlowAlgorithm>
void run_mcf_algorithm(const std::string &file_path, const std::string &name) {
    auto [G, costs, b] = Koala::DimacsGraphReader().read_minimum_cost_flow(file_path);
    Koala::MCFlowNetwork network(G,
        std::unordered_map<NetworKit::Edge, std::int64_t>(costs.begin(), costs.end()),
        std::unordered_map<NetworKit::node, std::int64_t>(b.begin(), b.end()));
    auto start = std::chrono::high_resolution_clock::now();
    auto minimum_cost_flow = FlowAlgorithm(network);
    minimum_cost_flow.run();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << name << ": Minimum cost flow flow = " << minimum_cost_flow.getMinCost()
              << ", Time taken = " << elapsed.count() << " seconds\n";
}

std::map<std::string, int> ALGORITHM = {
    { "EdmondsKarp", 1 },
    { "Orlin", 2 },
    { "GoldbergTarjan", 3 },
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
        run_mcf_algorithm<Koala::EdmondsKarpMinimumCostFlow>(file_path, "EdmondsKarp");
    } else if (ALGORITHM[algorithm] == 2) {
        run_mcf_algorithm<Koala::OrlinMinimumCostFlow>(file_path, "Orlin");
    } else if (ALGORITHM[algorithm] == 3) {
        run_mcf_algorithm<Koala::SuccessiveApproximationMinimumCostFlow>(
            file_path, "SuccessiveApproximation");
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    return 0;
}
