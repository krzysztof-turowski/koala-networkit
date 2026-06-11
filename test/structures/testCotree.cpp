#include <gtest/gtest.h>

#include <vector>

#include "structures/Cotree.hpp"

TEST(ConodeTest, initializes_defaults) {
    Koala::Conode node(1, 2, 3);

    EXPECT_EQ(node.first_child, 1);
    EXPECT_EQ(node.next_sibling, 2);
    EXPECT_EQ(node.previous_sibling, NetworKit::none);
    EXPECT_EQ(node.parent, 3);
    EXPECT_EQ(node.type, Koala::NodeType::UNKNOWN);
    EXPECT_EQ(node.size, 0);
}

TEST(CotreeTest, adds_nodes_and_maintains_child_sibling_links) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    tree.setRoot(root);

    tree.addChild(root, first);
    tree.addChild(root, second);

    EXPECT_EQ(tree.getRoot(), root);
    EXPECT_EQ(tree.getNode(root).first_child, second);
    EXPECT_EQ(tree.getNode(root).size, 2);
    EXPECT_EQ(tree.getNode(second).parent, root);
    EXPECT_EQ(tree.getNode(second).previous_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).parent, root);
    EXPECT_EQ(tree.getNode(first).previous_sibling, second);
    EXPECT_EQ(tree.getNode(first).next_sibling, NetworKit::none);
}

TEST(CotreeTest, removes_and_replaces_children_in_constant_time_shape) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    auto third = tree.add(Koala::NodeType::LEAF, 3);
    auto replacement = tree.add(Koala::NodeType::COMPLEMENT_NODE);
    tree.setRoot(root);

    tree.addChild(root, first);
    tree.addChild(root, second);
    tree.addChild(root, third);
    tree.removeChild(root, second);

    EXPECT_EQ(tree.getNode(root).first_child, third);
    EXPECT_EQ(tree.getNode(root).size, 2);
    EXPECT_EQ(tree.getNode(third).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).previous_sibling, third);
    EXPECT_EQ(tree.getNode(second).parent, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).next_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).previous_sibling, NetworKit::none);

    tree.replaceChild(root, third, replacement);

    EXPECT_EQ(tree.getNode(root).first_child, replacement);
    EXPECT_EQ(tree.getNode(root).size, 2);
    EXPECT_EQ(tree.getNode(replacement).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).previous_sibling, replacement);
    EXPECT_EQ(tree.getNode(third).parent, NetworKit::none);
}

TEST(CotreeTest, moves_existing_child_to_front_without_changing_size) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    auto third = tree.add(Koala::NodeType::LEAF, 3);
    tree.setRoot(root);
    tree.addChild(root, first);
    tree.addChild(root, second);
    tree.addChild(root, third);

    tree.moveChildToFront(root, first);

    EXPECT_EQ(tree.getNode(root).first_child, first);
    EXPECT_EQ(tree.getNode(root).size, 3);
    EXPECT_EQ(tree.getNode(first).next_sibling, third);
    EXPECT_EQ(tree.getNode(third).previous_sibling, first);
}

TEST(CotreeTest, clear_resets_storage_and_root) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE);
    tree.setRoot(root);

    tree.clear();

    EXPECT_EQ(tree.getRoot(), NetworKit::none);
    EXPECT_EQ(tree.upperNodeIdBound(), 0);
    EXPECT_FALSE(tree.prepared);
}
