#include <cstring>
#include <fstream>
#include <random>
#include <string>

std::string generate_filename(const std::string &directory) {
    const int FILENAME_LENGTH = 64;
    const std::string CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution distribution(0, static_cast<int>(CHARS.size()) - 1);
    std::string result(FILENAME_LENGTH, 0);
    std::generate_n(result.begin(), FILENAME_LENGTH, [&]() {
        return CHARS[distribution(generator)];
    });
    return directory + "/" + result;
}

bool compare_files(const std::string &leftPath, const std::string &rightPath) {
    const int BUFFER_SIZE = 8 * 1024 * 1024;
    std::ifstream leftFile(leftPath, std::ios::binary);
    std::ifstream rightFile(rightPath, std::ios::binary);
    if (!leftFile.good() || !rightFile.good()) {
        return false;
    }
    std::string leftBuffer(BUFFER_SIZE, 0), rightBuffer(BUFFER_SIZE, 0);
    while (leftFile.good() || rightFile.good()) {
        leftFile.read(leftBuffer.data(), BUFFER_SIZE);
        rightFile.read(rightBuffer.data(), BUFFER_SIZE);
        std::streamsize leftBytesCount = leftFile.gcount();
        std::streamsize rightBytesCount = rightFile.gcount();
        if (leftBytesCount != rightBytesCount) {
            return false;
        }
        if (std::memcmp(leftBuffer.data(), rightBuffer.data(), leftBytesCount) != 0) {
            return false;
        }
    }
    return true;
}

NetworKit::Graph build_graph(
        const int &N, const std::list<std::pair<int, int>> &E, bool directed = true) {
    NetworKit::Graph G(N, false, directed);
    for (const auto &[u, v] : E) {
        G.addEdge(u, v);
    }
    return G;
}

NetworKit::Graph build_graph(
        const int &N, const std::list<std::tuple<int, int, int>> &E, bool directed = true) {
    NetworKit::Graph G(N, true, directed);
    for (const auto &[u, v, w] : E) {
        G.increaseWeight(u, v, w);
        G.increaseWeight(v, u, 0);
    }
    return G;
}
