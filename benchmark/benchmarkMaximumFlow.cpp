#include <cassert>
#include <filesystem>
#include <iostream>
#include <chrono>  
#include <map>

#include <flow/MaximumFlow.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {
    if (argc == 1) {
        std::cout << "Directory name is empty" << std::endl;
        return 0;
    }
    std::string directory_path(argv[1]);
    if (!std::filesystem::exists(directory_path)) {
        std::cout << "Directory " << directory_path << " does not exist" << std::endl;
        return 0;
    }
    if (!std::filesystem::is_directory(directory_path)) {
        std::cout << directory_path << " is not a directory" << std::endl;
        return 0;
    }

    for (const auto& entry : std::filesystem::directory_iterator(directory_path)) {
        if (entry.is_regular_file()) {
            std::string path = entry.path().string();
            std::cout << "\nProcessing file: " << path << std::endl;

            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
                std::cout << "Running PushRelabel..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::PushRelabel(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "PushRelabel: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }

            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
                std::cout << "Running BKFlow..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::BKFlow(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "BKFlow: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }

            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
                std::cout << "Running MKMFlow..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::MKMFlow(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "MKMFlow: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }

            // {
            //     auto [G, s, t] = Koala::DimacsGraphReader().read_all(path);
            //     std::cout << "Running KingRaoTarjanMaximumFlow..." << std::endl;
            //     auto start = std::chrono::high_resolution_clock::now();

            //     auto maximum_flow = Koala::KingRaoTarjanMaximumFlow(G, s, t);
            //     maximum_flow.run();

            //     auto end = std::chrono::high_resolution_clock::now();
            //     std::chrono::duration<double> elapsed = end - start;
            //     std::cout << "KingRaoTarjanMaximumFlow: Maximum flow = " << maximum_flow.getFlowSize() 
            //               << ", Time taken = " << elapsed.count() << " seconds\n";
            // }
        }
    }
    return 0;
}
