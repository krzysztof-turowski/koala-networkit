#include <cassert>
#include <filesystem>
#include <iostream>
#include <chrono>

#include <matching/MaximumMatching.hpp>
#include <io/DimacsGraphReader.hpp>

enum OutputFormat { weight, csv };

struct Arguments {
    std::string algorithm;
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
            << "[--pefrect] [--check-has-perfect] [--initialization=<strategy>]"
            << "[--output=<output format>]" << std::endl;
        exit(1);
    }

    Arguments args;

    // Read algorithm and file name
    args.algorithm = all_args[0];
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
    for (int i = 2; i < all_args.size(); ++ i) {
        std::string option = all_args[i];

        if (option == "--perfect") {
            args.perfect = true;
        } else if (option == "--check-has-perfect") {
            args.checkPerfect = true;
        } else if (option.starts_with("--initialization=")) {
            std::string initOption = option.substr(17);
            if (initOption == "empty")
                args.initialization = Koala::BlossomMaximumMatching::InitializationStrategy::empty;
            else if (initOption == "greedy")
                args.initialization = Koala::BlossomMaximumMatching::InitializationStrategy::greedy;
            else {
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
    if (args.algorithm == "edmonds") {
        return new Koala::EdmondsMaximumMatching(G, args.perfect, args.initialization);
    } else if (args.algorithm == "gabow") {
        return new Koala::GabowMaximumMatching(G, args.perfect, args.initialization);
    } else if (args.algorithm == "micali") {
        return new Koala::GalilMicaliGabowMaximumMatching(G, args.perfect, args.initialization);
    } else if (args.algorithm == "scaling") {
        return new Koala::GabowScalingMatching(G, args.perfect);
    } else {
        std::cout << "Unknown algorithm: " << args.algorithm << std::endl;
        exit(1);
    }
}

std::pair<int, int> calc_matching_weight(
        const NetworKit::Graph& G, const std::map<NetworKit::node, NetworKit::node> matching) {
    int weight = 0;
    int cardinality = 0;
    for (auto [u, v] : matching) {
        if (v != NetworKit::none) {
            weight += static_cast<int>(G.weight(u, v));
            cardinality++;
        }
    }
    return {weight, cardinality};
}

void test_algorithm(NetworKit::Graph& G, Arguments args) {
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
            << args.filename << ","
            << G.upperNodeIdBound() << ","
            << G.upperEdgeIdBound() << ","
            << args.algorithm << ","
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
    
    // Read graph in
    auto G = Koala::DimacsGraphReader().read(args.filename);
    G.indexEdges(true);

    // Check if G has a perfect matching if needed
    if (args.perfect && args.checkPerfect && !check_matching_perfect(G)) {
        std::cout << "The provided graph does not contain a perfect matching" << std::endl;
        exit(1);
    }

    // Run algorithm
    test_algorithm(G, args);

    return 0;
}
