/*
 * Heap.hpp
 *
 *  Created on: 24.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <queue>

#include <structures/heap/BinomialHeap.hpp>
#include <structures/heap/FibonacciHeap.hpp>
#include <structures/heap/PairingHeap.hpp>

namespace Koala {

template<class T>
using Heap = std::priority_queue<T>;

}  // namespace Koala
