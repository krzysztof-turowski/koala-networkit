#include <cassert>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <benchmark/utils.hpp>
#include <matching/MaximumMatching.hpp>

enum OutputFormat { weight, csv };
enum class Algorithm : uint32_t { EDMONDS, GABOW, MICALI, SCALING };

std::map<std::string, Algorithm> ALGORITHM = {
    { "edmonds", Algorithm::EDMONDS },
    { "gabow", Algorithm::GABOW },
    { "micali", Algorithm::MICALI },
    { "scaling", Algorithm::SCALING }
};

struct Arguments {
    std::string algorithmName;
    Algorithm algorithm;
    std::string filename;
    bool perfect;
    bool checkPerfect;
    Koala::BlossomMaximumMatching::InitializationStrategy initialization;
    OutputFormat outputFormat;
};

Arguments parseArguments(int argc, char**argv) {
    std::string current_exec_name = argv[0];
    std::vector<std::string> all_args;
    all_args.assign(argv + 1, argv + argc);

    // Check if required arguments are present
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " algorithm filename"
            << "[--perfect] [--check-has-perfect] [--initialization=<strategy>]"
            << "[--output=<output format>]" << std::endl;
        exit(1);
    }

    Arguments args;

    // Read algorithm and file name
    args.algorithmName = all_args[0];
    auto algorithm = ALGORITHM.find(args.algorithmName);
    if (algorithm == ALGORITHM.end()) {
        std::cout << "Unknown algorithm: " << args.algorithmName << std::endl;
        exit(1);
    }
    args.algorithm = algorithm->second;
    args.filename = all_args[1];

    // Check if the file exists
    if (!std::filesystem::exists(args.filename)) {
        std::cout << "File " << args.filename << " does not exist" << std::endl;
        exit(1);
    }

    // Parse options
    args.perfect = false;
    args.checkPerfect = false;
    args.initialization = Koala::BlossomMaximumMatching::InitializationStrategy::empty;
    args.outputFormat = weight;
    for (std::size_t i = 2; i < all_args.size(); ++i) {
        std::string option = all_args[i];

        if (option == "--perfect") {
            args.perfect = true;
        } else if (option == "--check-has-perfect") {
            args.checkPerfect = true;
        } else if (option.starts_with("--initialization=")) {
            std::string initOption = option.substr(17);
            if (initOption == "empty") {
                args.initialization = Koala::BlossomMaximumMatching::InitializationStrategy::empty;
            } else if (initOption == "greedy") {
                args.initialization = Koala::BlossomMaximumMatching::InitializationStrategy::greedy;
            } else {
                std::cout << "Unknown initialization option '" << initOption << "'" << std::endl;
                exit(1);
            }
        } else if (option.starts_with("--output=")) {
            std::string outputOption = option.substr(9);
            if (outputOption == "weight") {
                args.outputFormat = weight;
            } else if (outputOption == "csv") {
                args.outputFormat = csv;
            } else {
                std::cout << "Unknown output format option '" << outputOption << "'" << std::endl;
                exit(1);
            }
        } else {
            std::cout << "Unknown option " << option << std::endl;
            exit(1);
        }
    }

    return args;
}

bool check_matching_perfect(NetworKit::Graph& G) {
    Koala::MicaliVaziraniMatching mcm(G);
    mcm.run();
    auto matching = mcm.getMatching();

    for (auto [v, u] : matching)
        if (u == NetworKit::none)
            return false;
    return true;
}

Koala::MaximumWeightMatching* get_algorithm(NetworKit::Graph& G, const Arguments& args) {
    switch (args.algorithm) {
    case Algorithm::EDMONDS:
        return new Koala::EdmondsMaximumMatching(G, args.perfect, args.initialization);
    case Algorithm::GABOW:
        return new Koala::GabowMaximumMatching(G, args.perfect, args.initialization);
    case Algorithm::MICALI:
        return new Koala::GalilMicaliGabowMaximumMatching(G, args.perfect, args.initialization);
    case Algorithm::SCALING:
        return new Koala::GabowScalingMatching(G, args.perfect);
    default:
        throw std::logic_error("Unhandled algorithm");
    }
}

std::pair<int64_t, int> calc_matching_weight(
        const NetworKit::Graph& G, const std::map<NetworKit::node, NetworKit::node> matching) {
    int64_t weight = 0;
    int cardinality = 0;
    for (auto [u, v] : matching) {
        if (v != NetworKit::none) {
            weight += static_cast<int>(G.weight(u, v));
            cardinality++;
        }
    }
    return {weight / 2, cardinality / 2};
}

void test_algorithm(const std::string &label, NetworKit::Graph& G, const Arguments &args) {
    auto start = std::chrono::high_resolution_clock::now();
    auto algorithm = get_algorithm(G, args);
    algorithm->run();
    auto end = std::chrono::high_resolution_clock::now();

    auto elapsed_time = end - start;
    auto microseconds = elapsed_time.count() / 1000;

    auto matching = algorithm->getMatching();
    auto [weight, cardinality] = calc_matching_weight(G, matching);

    if (args.outputFormat == OutputFormat::weight) {
        std::cout << weight << std::endl;
    } else {
        std::cout
            << label << ","
            << G.upperNodeIdBound() << ","
            << G.upperEdgeIdBound() << ","
            << args.algorithmName << ","
            << (args.perfect ? "perfect" : "weight") << ","
            << (args.initialization == Koala::BlossomMaximumMatching::InitializationStrategy::empty
                ? "empty" : "greedy") << ","
            << weight << ","
            << cardinality << ","
            << microseconds
            << std::endl;
    }

    delete algorithm;
}

int main(int argc, char **argv) {
    // Parse options
    auto args = parseArguments(argc, argv);

    Koala::Benchmark::executeForEachGraph(
        args.filename, [&](const std::string &label, NetworKit::Graph G) {
        G.indexEdges(true);
        if (args.perfect && args.checkPerfect && !check_matching_perfect(G)) {
            throw std::invalid_argument("The provided graph does not contain a perfect matching");
        }
        test_algorithm(label, G, args);
    });

    return 0;
}
