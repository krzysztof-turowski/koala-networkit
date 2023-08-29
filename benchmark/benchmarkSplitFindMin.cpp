#include <bits/stdc++.h>

#include <matching/PriorityQueues.hpp>

using sfm_id = int;
using SFM = Koala::SplitFindMin<int, sfm_id, int, int>;
using SFMN = Koala::SplitFindMinNaive<int, sfm_id, int, int>;

constexpr int inf = 99;

SFM sfm{1000, inf, -1, 3};

void process_line(std::string line) {
    std::string command;
    sfm_id id, id2;
    int key, value, node, low, high;
    std::stringstream ss; ss << line;
    ss >> command;
    if (command == "h" || command == "help") {
        std::cout << "h: help\n";
        std::cout << "i: init id nodes...\n";
        std::cout << "r: range id low high\n";
        std::cout << "f: find node\n";
        std::cout << "p: print node\n";
        std::cout << "s: split node id id2\n";
        std::cout << "d: decreaseKey node key\n";
        std::cout << "m: findMin node\n";
        std::cout << "q: quit\n";

    } else if (command == "i" || command == "init") {
        ss >> id;

        std::list<int> nodes;
        while (ss >> node) nodes.push_back(node);
        
        auto L = sfm.init(nodes, id);

        std::cout << L->id << " : ";
        sfm.print(L);
        std::cout << std::endl;

    } else if (command == "r" || command == "range") {
        ss >> id >> low >> high;

        std::list<int> nodes;
        for (int i = low; i <= high; ++ i) nodes.push_back(i);
        
        auto L = sfm.init(nodes, id);

        std::cout << L->id << " : ";
        sfm.print(L);
        std::cout << std::endl;

    } else if (command == "f" || command == "find") {
        ss >> node;

        std::cout << sfm.list(node)->id << std::endl;

    } else if (command == "p" || command == "print") {
        ss >> node;

        sfm.print(sfm.list(node));
        std::cout << std::endl;

    } else if (command == "s" || command == "split") {
        ss >> node >> id >> id2;

        auto [L1, L2] = sfm.split(node, id, id2);

        std::cout << L1->id << " : ";
        sfm.print(L1);
        std::cout << "\n" << L2->id << " : ";
        sfm.print(L2);
        std::cout << std::endl;
    } else if (command == "d" || command == "decreaseKey") {
        ss >> node >> key >> value;

        sfm.decreaseKey(node, key, node);

    } else if (command == "m" || command == "findMin") {
        ss >> node;

        auto L = sfm.list(node);
        auto [k, v] = sfm.findMin(L);

        std::cout << L->id << " : " << k << " " << v << std::endl;

    } else if (command == "q" || command == "quit") {
        exit(0);
    }
}

int rand_int(int a, int b) {
    return rand() % (b - a + 1) + a;
}

int main(int argc, char **argv) {
    // std::string line;
    // std::cerr << "> ";
    // while (std::getline(std::cin, line)) {
    //     process_line(line);
    //     std::cerr << "> ";    
    // }
    // return 0;

    srand(time(NULL));

    for (int t = 0; t < 10000; ++ t) {
        int size = rand_int(500, 2000);
        int inf = 1000;
        int key_changes_per_iteration = 5;

        std::cerr << "=====================================================\n";
        std::cerr << "TEST " << t << " FOR SIZE " << size << std::endl;

        int max_i = SFM::alpha(size * key_changes_per_iteration, size);

        std::cerr << "MAX I = " << max_i << std::endl;

        SFM fast(size, inf, -1, max_i);
        SFMN naive(size, inf, -1);

        int id_counter = 0;
        std::vector<int> to_split;
        for (int i = 0; i < size; ++ i) to_split.push_back(i);

        fast.init({to_split.begin(), to_split.end()}, id_counter);
        naive.init({to_split.begin(), to_split.end()}, id_counter ++);

        to_split.pop_back();

        while (to_split.size() > 0) {
            int split = rand_int(0, to_split.size() - 1);
            int split_point = to_split[split];
            std::swap(to_split[split], to_split[to_split.size() - 1]);
            to_split.pop_back();

            int id1 = id_counter ++;
            int id2 = id_counter ++;

            // std::cerr << "------------------------------------------------------\n";
            // std::cerr << "split on " << split_point << std::endl;
            
            fast.split(split_point, id1, id2);
            naive.split(split_point, id1, id2);

            for (int i = 0; i < size; ++ i) {     
                auto fl = fast.list(i);
                auto nl = naive.list(i);

                if (fl->id != nl->id) {
                    std::cerr << "wrong id for " << i << " " << fl->id << " =/= " << nl->id << std::endl;
                    exit(1);
                }

                if (fl->min_key != nl->min_key) {
                    std::cerr << "wrong min_key for " << i << " - " << fl->min_key << " =/= " << nl->min_key << std::endl;
                    exit(1);
                }

                // if (fl->min_val != nl->min_val) {
                //     std::cerr << "wrong min_val for " << i << " - " << fl->min_val << " =/= " << nl->min_val << std::endl;
                //     exit(1);
                // }
            }

            for (int i = 0; i < key_changes_per_iteration; ++ i) {
                int x = rand_int(0, size - 1);
                int key = rand_int(0, inf - 1);

                // std::cerr << "> decreaseKey(" << x << ", " << key << ")\n";

                fast.decreaseKey(x, key, x);
                naive.decreaseKey(x, key, x);
            }

            for (int i = 0; i < size; ++ i) {     
                auto fl = fast.list(i);
                auto nl = naive.list(i);

                if (fl->id != nl->id) {
                    std::cerr << "wrong id for " << i << " " << fl->id << " =/= " << nl->id << std::endl;
                    exit(1);
                }

                if (fl->min_key != nl->min_key) {
                    std::cerr << "wrong min_key for " << i << " - " << fl->min_key << " =/= " << nl->min_key << std::endl;
                    exit(1);
                }

                // if (fl->min_val != nl->min_val) {
                //     std::cerr << "wrong min_val for " << i << " - " << fl->min_val << " =/= " << nl->min_val << std::endl;
                //     exit(1);
                // }
            }
        }
    }
    
    return 0;
}
