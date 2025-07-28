/**
 * SoftHeap.hpp
 * 
 *  Created on: 07.07.2025
 *      Author: Hubert Bernacki (hubert.bernacki@student.uj.edu.pl)
 */

#pragma once

#include <list>
#include <memory>
#include <optional>

template<typename T>
class SoftHeap {
public:
    SoftHeap();
    SoftHeap(const SoftHeap&);
    SoftHeap(SoftHeap&&);
    ~SoftHeap();
    SoftHeap& operator=(const SoftHeap&);
    SoftHeap& operator=(SoftHeap&&);

    void insert(T);
    void meld(SoftHeap&&);
    T extractMin();

// private:
    // class Node {
    //     std::optional<std::unique_ptr<Node>> left;
    //     std::optional<std::unique_ptr<Node>> right;
    //     T val;
    //     int rank;
    //     std::list<T>
    // };

    // class Tree {
    //     std::optional<std::unique_ptr<Node>> root;
    // };

};

template<typename T>
SoftHeap<T>::~SoftHeap(){

}

template<typename T>
SoftHeap<T>::SoftHeap() {

}

template<typename T>
SoftHeap<T>::SoftHeap(const SoftHeap&) {

}

template<typename T>
SoftHeap<T>::SoftHeap(SoftHeap&&) {

}

// template<typename T>
// SoftHeap<T>::~SoftHeap() {

// }

template<typename T>
SoftHeap<T>& SoftHeap<T>::operator=(const SoftHeap&) {

}

template<typename T>
SoftHeap<T>& SoftHeap<T>::operator=(SoftHeap<T>&&) {

}

template<typename T>
void SoftHeap<T>::insert(T) {

}

template<typename T>
void SoftHeap<T>::meld(SoftHeap&&) {

}

template<typename T>
T SoftHeap<T>::extractMin() {

}