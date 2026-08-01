#pragma once

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

#include <io/DimacsGraphReader.hpp>
#include <io/G6GraphReader.hpp>

namespace Koala::Benchmark {

inline bool hasExtension(const std::string &path, const std::string &extension) {
    return path.ends_with("." + extension);
}

template <typename Callback>
void executeForEachGraph(const std::string &path, Callback callback) {
    if (hasExtension(path, "g6")) {
        std::fstream file(path, std::fstream::in);
        std::string line;
        while (file >> line) {
            callback(line, Koala::G6GraphReader().readline(line));
        }
        return;
    }

    if (hasExtension(path, "gr") || hasExtension(path, "col") || hasExtension(path, "max")) {
        callback(path, Koala::DimacsGraphReader().read(path));
        return;
    }

    throw std::invalid_argument("File type not supported: " + path);
}

}  // namespace Koala::Benchmark
