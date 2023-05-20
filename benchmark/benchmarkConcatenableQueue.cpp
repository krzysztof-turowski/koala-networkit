#include <cassert>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>

#include <matching/PriorityQueues.hpp>

using CQ = Koala::ConcatenableQueue<std::string, int, int>;
std::map<std::string, CQ*> Q;
std::map<int, CQ::ElementRef> E;

void process_line(std::string line) {
    std::string command, name, name2, name3;
    int after, val, val2, length;
    std::stringstream ss; ss << line;
    ss >> command;
    if (command == "h" || command == "help") {
        std::cout << "h: help\n";
        std::cout << "l: list\n";
        std::cout << "n: new name\n";
        std::cout << "b: build name elements...\n";
        std::cout << "g: range name start stop\n";
        std::cout << "p: print name\n";
        std::cout << "f: find val\n";
        std::cout << "i: insert name after val\n";
        std::cout << "a: append name val\n";
        std::cout << "r: remove name val\n";
        std::cout << "m: min name\n";
        std::cout << "s: split name split nameL nameR\n";
        std::cout << "c: concat nameL nameR name\n";
        std::cout << "t: reset\n";

    } else if (command == "n" || command == "new") {
        ss >> name;

        Q[name] = new CQ(name);

    } else if (command == "b" || command == "build") {
        ss >> name;

        CQ* q = new CQ(name);

        while (ss >> val) {
            if (E.find(val) != E.end()) { std::cerr << "Element " << val << " already exists\n"; break; }

            E[val] = q->append(val, val);
        }

        Q[name] = q;

        Q[name]->print(std::cerr);
        Q[name]->print_elements(std::cerr);

    } else if (command == "g" || command == "range") {
        ss >> name >> val >> val2;

        CQ* q = new CQ(name);

        for (int i = val; i <= val2; ++ i) {
            if (E.find(i) != E.end()) { std::cerr << "Element " << i << " already exists\n"; break; }

            E[i] = q->append(i, i);
        }

        Q[name] = q;

        Q[name]->print(std::cerr);
        Q[name]->print_elements(std::cerr);

    } else if (command == "l" || command == "list") {
        for (const auto& q : Q) {
            std::cout << q.first << "|";
            q.second->print_elements(std::cout);
        }

    } else if (command == "p" || command == "print") {
        ss >> name;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }

        Q[name]->print(std::cout);
        Q[name]->print_elements(std::cout);

    } else if (command == "f" || command == "find") {
        ss >> val;

        if (E.find(val) == E.end()) { std::cerr << "No such element\n"; return; }

        auto q = E[val]->find_queue();
        q->print(std::cerr);

    } else if (command == "i" || command == "insert") {
        ss >> name >> after >> val;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }
        if (E.find(after) == E.end()) { std::cerr << "No such element: " << after << "\n"; return; }
        if (E.find(val) != E.end()) { std::cerr << "Element " << val << " already exists\n"; return; }

        E[val] = Q[name]->insert_after(E[after], val, 0);
        Q[name]->print(std::cerr);
        Q[name]->print_elements(std::cerr);

    } else if (command == "a" || command == "append") {
        ss >> name >> val;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }
        if (E.find(val) != E.end()) { std::cerr << "Element " << val << " already exists\n"; return; }

        E[val] = Q[name]->append(val, val);
        Q[name]->print(std::cerr);
        Q[name]->print_elements(std::cerr);

    } else if (command == "r" || command == "remove") {
        ss >> name >> val;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }
        if (E.find(val) == E.end()) { std::cerr << "No such element: " << val << "\n"; return; }

        Q[name]->remove(E[val]);
        E.erase(val);
        Q[name]->print(std::cerr);
        Q[name]->print_elements(std::cerr);

    } else if (command == "m" || command == "min") {
        ss >> name;

        std::cout << Q[name]->find_min().first << std::endl;
    } else if (command == "s" || command == "split") {
        ss >> name >> val >> name2 >> name3;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }
        if (E.find(val) == E.end()) { std::cerr << "No such element: " << val << "\n"; return; }

        auto [q1, q2] = Q[name]->split(E[val], name2, name3);
        delete Q[name];
        Q.erase(name);
        Q[name2] = q1;
        Q[name3] = q2;

    } else if (command == "c" || command == "concat") {
        ss >> name >> name2 >> name3;

        if (Q.find(name) == Q.end()) { std::cerr << "No such queue '" << name << "'\n"; return; }
        if (Q.find(name2) == Q.end()) { std::cerr << "No such queue '" << name2 << "'\n"; return; }

        CQ* q1 = Q[name];
        CQ* q2 = Q[name2];
        q1->concat(std::move(*q2), name3);
        Q.erase(name);
        Q.erase(name2);
        Q[name3] = q1;
        delete q2;

        Q[name3]->print_elements(std::cerr);
    } else if (command == "q" || command == "quit") {
        exit(0);
    } else if (command == "t" || command == "reset") {
        std::cout << "Reseting...\n";
        for (auto [n, q] : Q) delete q;
        Q.clear();
        E.clear();
    } else {
        std::cerr << "Unknown command '" << command << "'\n";
    }
}

int main(int argc, char **argv) {   
    std::string line;

    std::cerr << "> ";
    while (std::getline(std::cin, line)) {
        process_line(line);
        std::cerr << "> ";    
    }

    for (auto [n, q] : Q) delete q;

    return 0;
}
