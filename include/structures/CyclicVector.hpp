//
// Created by mikolajtwarog on 2021-04-28.
//

#ifndef TECHNIKA_BAKER_CYCLIC_VECTOR_HPP
#define TECHNIKA_BAKER_CYCLIC_VECTOR_HPP


#include <vector>

template <typename T>
class cyclic_vector : public std::vector<T> {
public:
    T& operator[] (signed int x) {
        while (x < 0) {
            x += this->size();
        }
        x %= this->size();
        return this->at(x);
    }
};


#endif //TECHNIKA_BAKER_CYCLIC_VECTOR_HPP
