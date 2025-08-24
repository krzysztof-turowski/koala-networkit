#include <gtest/gtest.h>
#include <structures/priority_queue/VanEmdeBoasTree.hpp>
#include <structures/priority_queue/XFastTrie.hpp>
#include <structures/priority_queue/YFastTrie.hpp>
#include <structures/priority_queue/WeakHeap.hpp>
#include <structures/priority_queue/SkewHeap.hpp>
#include <structures/priority_queue/RankPairingHeap.hpp>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <memory>
#include <string>

class PriorityQueueTest :
    public testing::TestWithParam<std::shared_ptr<Koala::PriorityQueue<uint32_t>>> {
 protected:
    std::shared_ptr<Koala::PriorityQueue<uint32_t>> queue;

    void SetUp() override {
        queue = GetParam();
    }
};

INSTANTIATE_TEST_SUITE_P(
    PriorityQueue,
    PriorityQueueTest,
    testing::Values(
        std::make_shared<Koala::VanEmdeBoasTree<uint32_t>>(1024),
        std::make_shared<Koala::XFastTrie<uint32_t>>(1024),
        std::make_shared<Koala::YFastTrie<uint32_t>>(1024),
        std::make_shared<Koala::WeakHeap<uint32_t>>(),
        std::make_shared<Koala::SkewHeap<uint32_t>>(),
        std::make_shared<Koala::RankPairingHeap<uint32_t>>()
    )
);

TEST_P(PriorityQueueTest, PushPopIncreasing) {
    const uint32_t MAX = 1000;
    for (uint32_t i = 0; i < MAX; i++) {
        queue->push(i);
        ASSERT_FALSE(queue->empty());
    }
    for (uint32_t i = 0; i < MAX; i++) {
        ASSERT_EQ(queue->top(), i);
        queue->pop();
    }
    ASSERT_TRUE(queue->empty());
}

TEST_P(PriorityQueueTest, PushPopDecreasing) {
    const uint32_t MAX = 1000;
    uint32_t i = MAX;
    do {
        i--;
        queue->push(i);
        ASSERT_FALSE(queue->empty());
    } while (i > 0);
    for (uint32_t i = 0; i < MAX; i++) {
        ASSERT_EQ(queue->top(), i);
        queue->pop();
    }
    ASSERT_TRUE(queue->empty());
}

TEST_P(PriorityQueueTest, TopEmptyThrows) {
    ASSERT_TRUE(queue->empty());
    ASSERT_THROW(queue->top(), std::runtime_error);
}

TEST_P(PriorityQueueTest, PopEmptyThrows) {
    ASSERT_TRUE(queue->empty());
    ASSERT_THROW(queue->pop(), std::runtime_error);
}

TEST_P(PriorityQueueTest, PushPopRandomOrder) {
    std::vector<uint32_t> values = {5, 1, 9, 3, 7};
    for (uint32_t value : values) {
        queue->push(value);
    }
    std::sort(values.begin(), values.end());
    for (uint32_t value : values) {
        ASSERT_EQ(queue->top(), value);
        queue->pop();
    }
    ASSERT_TRUE(queue->empty());
}
