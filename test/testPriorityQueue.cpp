#include <gtest/gtest.h>
#include <structures/priority_queue/VanEmdeBoasTree.hpp>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <memory>
#include <string>

class IntPriorityQueueTest :
    public testing::TestWithParam<std::shared_ptr<Koala::PriorityQueue<int>>> {
 protected:
    std::shared_ptr<Koala::PriorityQueue<int>> queue;

    void SetUp() override {
        queue = GetParam();
    }
};

INSTANTIATE_TEST_SUITE_P(
    IntPriorityQueue,
    IntPriorityQueueTest,
    testing::Values(
        std::make_shared<Koala::VanEmdeBoasTree>(1024)
        // TODO: Add X-fast Trie, Y-fast Trie, Fusion Tree
    )
);

TEST_P(IntPriorityQueueTest, PushPopIncreasing) {
    const int MAX = 1000;
    for (int i = 0; i < MAX; i++) {
        queue->push(i);
        ASSERT_FALSE(queue->empty());
    }
    for (int i = 0; i < MAX; i++) {
        ASSERT_EQ(queue->peek(), i);
        ASSERT_EQ(queue->pop(), i);
    }
    ASSERT_TRUE(queue->empty());
}

TEST_P(IntPriorityQueueTest, PushPopDecreasing) {
    const int MAX = 1000;
    for (int i = MAX - 1; i >= 0; i--) {
        queue->push(i);
        ASSERT_FALSE(queue->empty());
    }
    for (int i = 0; i < MAX; i++) {
        ASSERT_EQ(queue->peek(), i);
        ASSERT_EQ(queue->pop(), i);
    }
    ASSERT_TRUE(queue->empty());
}

TEST_P(IntPriorityQueueTest, PeekEmptyThrows) {
    ASSERT_TRUE(queue->empty());
    ASSERT_THROW(queue->peek(), std::runtime_error);
}

TEST_P(IntPriorityQueueTest, PopEmptyThrows) {
    ASSERT_TRUE(queue->empty());
    ASSERT_THROW(queue->pop(), std::runtime_error);
}

TEST_P(IntPriorityQueueTest, PushPopRandomOrder) {
    std::vector<int> values = {5, 1, 9, 3, 7};
    for (int value : values) {
        queue->push(value);
    }
    std::sort(values.begin(), values.end());
    for (int value : values) {
        ASSERT_EQ(queue->peek(), value);
        ASSERT_EQ(queue->pop(), value);
    }
    ASSERT_TRUE(queue->empty());
}

// TODO: Add generic priority queues (Brodal Queue, Weak Heap)
template <typename T>
class GenericPriorityQueueTest : public testing::Test {
 protected:
};
