#include <gtest/gtest.h>

#include <structures/Heap.hpp>

typedef testing::Types<
    Koala::BinomialHeap<int, std::greater<int>>,
    Koala::FibonacciHeap<int, std::greater<int>>,
    Koala::PairingHeap<int, std::greater<int>>> Heaps;

template <typename Heap>
class HeapTest : public testing::Test {
 public:
    Heap heap;

    typename Heap::iterator get_iterator(NetworKit::index i) {
        return typename Heap::iterator(i);
    }
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

TYPED_TEST(HeapTest, EraseIncreasing) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = 0; i < MAX; i++) {
        ASSERT_EQ(this->heap.top(), i);
        this->heap.erase(this->get_iterator(i));
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, EraseIncreasingNext) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = 1; i < MAX; i++) {
        ASSERT_EQ(this->heap.top(), 0);
        this->heap.erase(this->get_iterator(i));
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, EraseDecreasing) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = MAX - 1; i > 0; i--) {
        ASSERT_EQ(this->heap.top(), 0);
        this->heap.erase(this->get_iterator(i));
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, UpdateDecreasing) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = 0; i < MAX; i++) {
        this->heap.update(this->get_iterator(i), -i);
        ASSERT_EQ(this->heap.top(), -i);
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, UpdateDecreasing_Half) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = MAX - 1; i >= MAX / 2; i--) {
        this->heap.update(this->get_iterator(i), i - MAX / 2);
        ASSERT_EQ(this->heap.top(), 0);
        this->heap.check();
    }
}

TYPED_TEST(HeapTest, UpdateDecreasing_Other) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        this->heap.push(i);
        this->heap.check();
    }
    for (int i = MAX - 1; i >= 0; i--) {
        this->heap.update(this->get_iterator(i), i - MAX);
        ASSERT_EQ(this->heap.top(), i - MAX);
        this->heap.check();
    }
}
