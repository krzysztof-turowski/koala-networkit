#include "techniques/separator/BalancedSeparator.hpp"

BalancedSeparator::BalancedSeparator(const NetworKit::Graph& graph)
    : graph(graph) {}

const BalancedSeparator::Partition& BalancedSeparator::getPartition() const {
    assureFinished();
    return partition;
}

const std::vector<NetworKit::node>& BalancedSeparator::getSeparator() const {
    assureFinished();
    return partition.separator;
}

const std::vector<NetworKit::node>& BalancedSeparator::getPartitionA() const {
    assureFinished();
    return partition.A;
}

const std::vector<NetworKit::node>& BalancedSeparator::getPartitionB() const {
    assureFinished();
    return partition.B;
}
