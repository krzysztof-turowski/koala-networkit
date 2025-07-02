#include <cassert>
#include <filesystem>
#include <iostream>
#include <chrono>  
#include <map>

#include <flow/MaximumFlow.hpp>
#include <io/DimacsGraphReader.hpp>

int main(int argc, char **argv) {
    if (argc == 1) {
        std::cout << "File name is empty" << std::endl;
        return 0;
    }
    if (argc == 2) {
        std::cout << "Algorithm not specified" << std::endl;
        return 0;
    }

    std::string file_path(argv[1]);
    int algo_num = std::stoi(argv[2]);
    if (!std::filesystem::exists(file_path)) {
        std::cout << "File " << file_path << " does not exist" << std::endl;
        return 0;
    }
    if (std::filesystem::is_directory(file_path)) {
        std::cout << file_path << " is a directory" << std::endl;
        return 0;
    }

    std::cout << "\nProcessing file: " << file_path << std::endl;
        if(algo_num == 1)
            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(file_path);
                std::cout << "Running PushRelabel..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::PushRelabel(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "PushRelabel: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }
        if(algo_num == 2)
            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(file_path);
                std::cout << "Running BKFlow..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::BKFlow(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "BKFlow: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }
        if(algo_num == 3)
            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(file_path);
                std::cout << "Running MKMFlow..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::MKMFlow(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "MKMFlow: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }
        if(algo_num == 4)
            {
                auto [G, s, t] = Koala::DimacsGraphReader().read_all(file_path);
                std::cout << "Running KingRaoTarjanMaximumFlow..." << std::endl;
                auto start = std::chrono::high_resolution_clock::now();

                auto maximum_flow = Koala::KingRaoTarjanMaximumFlow(G, s, t);
                maximum_flow.run();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "KingRaoTarjanMaximumFlow: Maximum flow = " << maximum_flow.getFlowSize() 
                          << ", Time taken = " << elapsed.count() << " seconds\n";
            }
    return 0;
}
