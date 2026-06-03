#include <gtest/gtest.h>

#include <type_traits>
#include <vector>

#include <structures/DynamicTree.hpp>
#include <structures/dynamic_tree/LinkCutTree.hpp>
#include <structures/dynamic_tree/NaiveDynamicTree.hpp>

namespace {

void expectSharedOperations(Koala::DynamicTree<int> &tree) {
    tree.link(1, 0, 7);
    tree.link(2, 1, 4);

    EXPECT_EQ(tree.findRoot(2), 0);

    tree.cut(2, 1);
    EXPECT_EQ(tree.findRoot(2), 2);
}

}  // namespace

static_assert(std::is_same_v<
              Koala::DefaultDynamicTree<int>, Koala::LinkCutTree<int>>);

TEST(DynamicTreeTest, link_cut_tree_supports_shared_operations) {
    Koala::LinkCutTree<int> tree;
    tree.initialize(3);

    expectSharedOperations(tree);
}

TEST(DynamicTreeTest, naive_tree_supports_shared_operations) {
    std::vector<std::vector<int>> weights(3, std::vector<int>(3, 0));
    Koala::NaiveDynamicTree<int> tree(3, weights);

    expectSharedOperations(tree);
}

TEST(LinkCutTreeTest, supports_link_cut_and_path_updates) {
    Koala::LinkCutTree<int> tree;
    tree.initialize(4);
    tree.link(1, 0, 7);
    tree.link(2, 1, 4);
    tree.link(3, 1, 9);

    EXPECT_EQ(tree.findRoot(2), 0);
    EXPECT_EQ(tree.findParent(2), 1);
    EXPECT_EQ(tree.findParent(1), 0);
    EXPECT_EQ(tree.getMinimumPathResidualCapacity(2), 4);
    EXPECT_EQ(tree.getMinimumPathResidualCapacity(3), 7);
    EXPECT_TRUE(tree.findChildren(1).contains(2));
    EXPECT_TRUE(tree.findChildren(1).contains(3));

    tree.addValue(2, -4);
    EXPECT_EQ(tree.getValue(2), 0);
    EXPECT_EQ(tree.getValue(1), 3);
    EXPECT_EQ(tree.findSaturatedEdge(2), NetworKit::Edge(2, 1));

    tree.cut(2, 1);
    EXPECT_EQ(tree.findRoot(2), 2);
    EXPECT_EQ(tree.findParent(2), NetworKit::none);
    EXPECT_FALSE(tree.findChildren(1).contains(2));
    EXPECT_EQ(tree.getValue(2), std::nullopt);
    EXPECT_EQ(tree.getMinimumPathResidualCapacity(2), std::nullopt);

    tree.addValue(3, -3);
    EXPECT_EQ(tree.findSaturatedEdge(3), NetworKit::Edge(1, 0));

    tree.initialize(2);
    EXPECT_EQ(tree.findRoot(1), 1);
    EXPECT_EQ(tree.findParent(1), NetworKit::none);
    EXPECT_EQ(tree.getValue(1), std::nullopt);
}

TEST(NaiveDynamicTreeTest, supports_electrical_flow_path_operations) {
    std::vector<std::vector<double>> weights(3, std::vector<double>(3, 0));
    weights[2][1] = 0.25;
    weights[1][2] = -0.25;
    weights[1][0] = 0.75;
    weights[0][1] = -0.75;
    Koala::NaiveDynamicTree<double> tree(3, weights);
    tree.link(1, 0, weights[1][0]);
    tree.link(2, 1, weights[2][1]);

    EXPECT_DOUBLE_EQ(tree.pathSum(2, 0), 1);
    EXPECT_EQ(tree.pathMin(2, 0), NetworKit::Edge(2, 1));

    tree.pathAdd(2, 0, 0.25);
    EXPECT_DOUBLE_EQ(tree.getWeight(2, 1), 0);
    EXPECT_DOUBLE_EQ(tree.getWeight(1, 2), 0);
    EXPECT_DOUBLE_EQ(tree.getWeight(1, 0), 0.5);
    EXPECT_DOUBLE_EQ(tree.getWeight(0, 1), -0.5);
}
