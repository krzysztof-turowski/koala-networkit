#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <coloring/DRVertexQueue.hpp>

struct VertexQueueParameters {
    std::vector<std::pair<int, int>> nodes;
    std::vector<Koala::NodeValuePair> expected;
};

struct ChangeValueParameters {
    std::vector<std::pair<int, int>> nodes;
    std::vector<std::pair<int, int>> changes;
    std::vector<Koala::NodeValuePair> expected;
};

class VertexQueueTest: public testing::TestWithParam<VertexQueueParameters> {};
class VertexQueueChangeValueTest : public testing::TestWithParam<ChangeValueParameters> {};

TEST_P(VertexQueueTest, test) {
    VertexQueueParameters const& params = GetParam();
    Koala::DRVertexQueue queue(params.nodes);
    for (auto const& node : params.expected) {
        EXPECT_EQ(node, queue.pop());
    }

    EXPECT_TRUE(queue.empty());
}

INSTANTIATE_TEST_SUITE_P(
    test_example,
    VertexQueueTest,
    testing::Values(
        VertexQueueParameters{
            {{1, 2}, {2, 3}, {0, 2}, {4, 1}, {5, 0}},
            {{2, 3}, {0, 2}, {1, 2}, {4, 1}, {5, 0}}
        }
    )
);

TEST_P(VertexQueueChangeValueTest, test) {
    ChangeValueParameters const& params = GetParam();
    Koala::DRVertexQueue queue(params.nodes);

    for (auto const& change : params.changes) {
        queue.updateValue(change.first, change.second);
    }

    for (auto const& node : params.expected) {
        EXPECT_EQ(node, queue.pop());
    }
}

INSTANTIATE_TEST_SUITE_P(
    test_example,
    VertexQueueChangeValueTest,
    testing::Values(
        ChangeValueParameters{
            {{1, 2}, {2, 3}, {0, 2}, {4, 1}, {5, 0}},
            {{4, 5}},
            {{4, 5}, {2, 3}, {0, 2}, {1, 2}, {5, 0}}
        },
        ChangeValueParameters{
            {{1, 2}, {2, 3}, {0, 2}, {4, 1}, {5, 0}},
            {{4, 5}, {2, 2}},
            {{4, 5}, {0, 2}, {1, 2}, {2, 2}, {5, 0}}
        }
    )
);
