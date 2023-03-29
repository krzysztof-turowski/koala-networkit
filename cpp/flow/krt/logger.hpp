#include <utility>

#ifndef BETTER_MAX_FLOW_LOGGER_HPP
#define BETTER_MAX_FLOW_LOGGER_HPP

#define DBG if(1)

class Logger {
    std::string tag;
public:

    Logger(std::string tag) : tag(std::move(tag)) {}

    void log(std::string msg, std::pair<int, int> e, int val) {
        return log(msg, e.first, e.second, val);
    }

    void log(std::string msg, int u, int v, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " " + "(" << u << ", " << v << ") = " << val << std::endl;
    }

    void log(std::string msg, int u, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " " + "(" << u << ") = " << val << std::endl;
    }

    void log(std::string msg, int val) {
        DBG std::cerr << "[" + tag + "]" + " " + msg + " = " << val << std::endl;
    }

    void log(std::string msg) {
        DBG std::cerr << "[" + tag + "]" + " " + msg << std::endl;
    }
};

#endif //BETTER_MAX_FLOW_LOGGER_HPP
