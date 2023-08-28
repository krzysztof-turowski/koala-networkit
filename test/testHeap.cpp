#include <gtest/gtest.h>

#include <list>

#include <structures/heap/BinomialHeap.hpp>
#include <structures/heap/FibonacciHeap.hpp>
#include <structures/heap/PairingHeap.hpp>

typedef testing::Types<
    Koala::BinomialHeap<int, std::greater<int>>,
    Koala::FibonacciHeap<int, std::greater<int>>,
    Koala::PairingHeap<int, std::greater<int>>> Heaps;

template <typename Heap>
class HeapTest : public testing::Test {
 public:
    Heap heap;
};

TYPED_TEST_SUITE(HeapTest, Heaps);

TYPED_TEST(HeapTest, PushPopIncreasing) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = 0; i < MAX; i++) {
        ASSERT_EQ(this->heap.top(), i);
        this->heap.pop();
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, PushPopDecreasing) {
    const int MAX = 1000;
    for (int i = MAX - 1; i >= 0; i--) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = 0; i < MAX; i++) {
        ASSERT_EQ(this->heap.top(), i);
        this->heap.pop();
        this->heap.check();
    }
}
