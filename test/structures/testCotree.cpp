#include <gtest/gtest.h>

#include <vector>

#include "structures/Cotree.hpp"

TEST(ConodeTest, initializes_defaults) {
    Koala::Conode node(1, 2, 3);

    EXPECT_EQ(node.first_child, 1);
    EXPECT_EQ(node.next_sibling, 2);
    EXPECT_EQ(node.previous_sibling, NetworKit::none);
    EXPECT_EQ(node.parent, 3);
    EXPECT_EQ(node.size, 0);
    EXPECT_EQ(node.type, Koala::NodeType::UNKNOWN);
    EXPECT_EQ(node.marked, Koala::Marked::UNMARKED);
    EXPECT_EQ(node.md, 0);
    EXPECT_EQ(node.d, 0);
    EXPECT_FALSE(node.in_graph);
    EXPECT_EQ(node.number_of_vertices_in_subtree, 0);
    EXPECT_EQ(node.time_in, 0);
    EXPECT_EQ(node.time_out, 0);
    EXPECT_EQ(node.get_up[0], NetworKit::none);
}

TEST(CotreeTest, adds_nodes_and_maintains_child_sibling_links) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE, 0);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    tree.setRoot(root);

    tree.addChild(root, first);
    tree.addChild(root, second);

    EXPECT_EQ(tree.getRoot(), root);
    EXPECT_EQ(tree.getNode(root).first_child, second);
    EXPECT_EQ(tree.getNode(root).d, 2);
    EXPECT_EQ(tree.getNode(second).parent, root);
    EXPECT_EQ(tree.getNode(second).previous_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).parent, root);
    EXPECT_EQ(tree.getNode(first).previous_sibling, second);
    EXPECT_EQ(tree.getNode(first).next_sibling, NetworKit::none);
}

TEST(CotreeTest, removes_and_replaces_children_in_constant_time_shape) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE, 0);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    auto third = tree.add(Koala::NodeType::LEAF, 3);
    auto replacement = tree.add(Koala::NodeType::COMPLEMENT_NODE, 1);
    tree.setRoot(root);

    tree.addChild(root, first);
    tree.addChild(root, second);
    tree.addChild(root, third);
    tree.removeChild(root, second);

    EXPECT_EQ(tree.getNode(root).first_child, third);
    EXPECT_EQ(tree.getNode(root).d, 2);
    EXPECT_EQ(tree.getNode(third).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).previous_sibling, third);
    EXPECT_EQ(tree.getNode(second).parent, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).next_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(second).previous_sibling, NetworKit::none);

    tree.replaceChild(root, third, replacement);

    EXPECT_EQ(tree.getNode(root).first_child, replacement);
    EXPECT_EQ(tree.getNode(root).d, 2);
    EXPECT_EQ(tree.getNode(replacement).next_sibling, first);
    EXPECT_EQ(tree.getNode(first).previous_sibling, replacement);
    EXPECT_EQ(tree.getNode(third).parent, NetworKit::none);
}

TEST(CotreeTest, moves_existing_child_to_front_without_changing_size) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE, 0);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    auto third = tree.add(Koala::NodeType::LEAF, 3);
    tree.setRoot(root);
    tree.addChild(root, first);
    tree.addChild(root, second);
    tree.addChild(root, third);

    tree.moveChildToFront(root, first);

    EXPECT_EQ(tree.getNode(root).first_child, first);
    EXPECT_EQ(tree.getNode(root).d, 3);
    EXPECT_EQ(tree.getNode(first).next_sibling, third);
    EXPECT_EQ(tree.getNode(third).previous_sibling, first);
}

TEST(CotreeTest, removes_marked_prefix_and_unmarked_suffix) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE, 0);
    auto first = tree.add(Koala::NodeType::LEAF, 1);
    auto second = tree.add(Koala::NodeType::LEAF, 2);
    auto third = tree.add(Koala::NodeType::LEAF, 3);
    tree.setRoot(root);
    tree.addChild(root, first);
    tree.addChild(root, second);
    tree.addChild(root, third);
    tree.unmark(third);
    tree.unmark(second);

    auto removed = tree.removeWereMarked(root);

    ASSERT_EQ(removed.size(), 2);
    EXPECT_EQ(removed[0], third);
    EXPECT_EQ(removed[1], second);
    EXPECT_EQ(tree.getNode(root).first_child, first);
    EXPECT_EQ(tree.getNode(root).d, 1);

    auto fourth = tree.add(Koala::NodeType::LEAF, 4);
    auto fifth = tree.add(Koala::NodeType::LEAF, 5);
    tree.addChild(root, fourth);
    tree.addChild(root, fifth);
    tree.unmark(fifth);

    tree.removeWereNotMarked(root);

    EXPECT_EQ(tree.getNode(root).first_child, fifth);
    EXPECT_EQ(tree.getNode(fifth).next_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(root).d, 1);
}

TEST(CotreeTest, build_tree_encodes_habib_paul_order_as_indexed_nary_tree) {
    NetworKit::Graph reduced_graph(1);
    Koala::Cotree tree(reduced_graph);
    tree.setOrder({
        {{0, 1}, 0},
        {{1, 2}, 1},
        {{2, NetworKit::none}, 3},
    });

    tree.buildTree();

    ASSERT_TRUE(tree.prepared);
    EXPECT_EQ(tree.getRoot(), 3);
    EXPECT_GE(tree.upperNodeIdBound(), 6);
    EXPECT_EQ(tree.getNode(3).first_child, 4);
    EXPECT_EQ(tree.getNode(4).type, Koala::NodeType::COMPLEMENT_NODE);
    EXPECT_EQ(tree.getNode(4).first_child, 5);
    EXPECT_EQ(tree.getNode(5).type, Koala::NodeType::UNION_NODE);
    EXPECT_EQ(tree.getNode(5).first_child, 0);
    EXPECT_EQ(tree.getNode(0).next_sibling, 1);
    EXPECT_EQ(tree.getNode(1).parent, 5);
    EXPECT_EQ(tree.getNode(1).next_sibling, NetworKit::none);
    EXPECT_EQ(tree.getNode(5).next_sibling, 2);
}

TEST(CotreeTest, clear_resets_storage_and_root) {
    Koala::Cotree tree;
    auto root = tree.add(Koala::NodeType::UNION_NODE, 0);
    tree.setRoot(root);

    tree.clear();

    EXPECT_EQ(tree.getRoot(), NetworKit::none);
    EXPECT_EQ(tree.upperNodeIdBound(), 0);
    EXPECT_FALSE(tree.prepared);
}
