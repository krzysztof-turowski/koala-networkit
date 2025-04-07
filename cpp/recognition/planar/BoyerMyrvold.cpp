/*
 * BoyerMyrvold.cpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Implementation of Boyer-Myrvold planarity testing algorithm
 */

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <networkit/graph/Graph.hpp>
#include <recognition/planar/PlanarGraphRecognition.hpp>

namespace Koala {

class BoyerBlock {
 private:
        std::list<NetworKit::node> separatedDFSChildList;
        std::list<NetworKit::Edge> edges;
        bool flipped;
        NetworKit::node root;
        std::map<NetworKit::node, int> vertexOrientation;

     public:
        BoyerBlock(NetworKit::node r) : flipped(false), root(r) {}
        
        void addChild(NetworKit::node v) {
            separatedDFSChildList.push_back(v);
        }
        
        void addEdge(NetworKit::Edge e) {
            edges.push_back(e);
        }
        
        void flip() {
            flipped = !flipped;
            separatedDFSChildList.reverse();
            edges.reverse();
            
            // Update vertex orientations
            for (auto& pair : vertexOrientation) {
                if (pair.second == LEFT) {
                    pair.second = RIGHT;
                } else if (pair.second == RIGHT) {
                    pair.second = LEFT;
                }
            }
        }
        
        bool isFlipped() const { return flipped; }
        
        const std::list<NetworKit::node>& getChildren() const { 
            return separatedDFSChildList; 
        }
        
        const std::list<NetworKit::Edge>& getEdges() const { 
            return edges; 
        }
        
        NetworKit::node getRoot() const { return root; }
        
        void setVertexOrientation(NetworKit::node v, int orientation) {
            vertexOrientation[v] = orientation;
        }
        
        int getVertexOrientation(NetworKit::node v) const {
            auto it = vertexOrientation.find(v);
            if (it != vertexOrientation.end()) {
                return it->second;
            }
            return 0; // Not set
        }
        
        void merge(BoyerBlock& other) {
            separatedDFSChildList.splice(separatedDFSChildList.end(), 
                                    other.separatedDFSChildList);
            edges.splice(edges.end(), other.edges);
            
            // Merge vertex orientations
            for (const auto& pair : other.vertexOrientation) {
                vertexOrientation[pair.first] = pair.second;
            }
        }
};

    BoyerMyrvold::BoyerMyrvold(NetworKit::Graph &graph, bool embedding)
            : PlanarGraphRecognition(graph, embedding) {
            N = graph.numberOfNodes();
        }

    void  BoyerMyrvold::run() {
            std::vector<bool> reached(N, false);
            std::vector<NetworKit::node> dfsnum(N);
            std::vector<NetworKit::node> lowpt(N);
            std::vector<NetworKit::node> lowpt2(N);
            std::vector<NetworKit::node> parent(N, NetworKit::none);
            std::map<NetworKit::node, BoyerBlock*> blocks;
            std::map<NetworKit::node, std::list<NetworKit::node>> vertexChildListMap;
            std::map<NetworKit::node, int> vertexLowpointMap;
            std::map<NetworKit::node, std::set<NetworKit::node>> externallyActiveVertices;
            
            // Perform DFS to get numbering and lowpoints
            int dfs_count = 0;
            std::vector<NetworKit::Edge> treeEdges;
            dfs_in_reorder(graph, 0, treeEdges, dfs_count, reached, dfsnum, lowpt, lowpt2, parent);
            
            // Compute externally active vertices
            for (auto v : graph.nodeRange()) {
                for (auto w : graph.neighborRange(v)) {
                    if (dfsnum[w] < dfsnum[v] && parent[v] != w) {
                        // w is a back edge endpoint
                        NetworKit::node current = v;
                        while (current != w) {
                            externallyActiveVertices[current].insert(w);
                            current = parent[current];
                        }
                    }
                }
            }
            
            // Sort children by lowpoint values
            for (auto v : graph.nodeRange()) {
                std::list<NetworKit::node> children;
                for (auto w : graph.neighborRange(v)) {
                    if (parent[w] == v) {
                        children.push_back(w);
                    }
                }
                children.sort([&](NetworKit::node a, NetworKit::node b) {
                    return lowpt[a] < lowpt[b];
                });
                vertexChildListMap[v] = children;
            }
            
            // Process vertices in reverse DFS order
            for (int i = N - 1; i >= 0; i--) {
                NetworKit::node v = i;
                
                // Process tree edges
                for (auto w : vertexChildListMap[v]) {
                    NetworKit::Edge e(v, w);
                    BoyerBlock b(w);
                    b.addEdge(e);
                    b.addChild(w);
                    b.setVertexOrientation(v, LEFT);
                    b.setVertexOrientation(w, RIGHT);
                    blocks[w] = &b;
                }
                
                // Process back edges
                std::list<NetworKit::Edge> backEdges;
                for (auto w : graph.neighborRange(v)) {
                    if (dfsnum[w] < dfsnum[v] && parent[v] != w) {
                        backEdges.push_back(NetworKit::Edge(v, w));
                    }
                }
                
                // Walkup phase
                for (auto& e : backEdges) {
                    NetworKit::node w = e.u == v ? e.v : e.u;
                    std::stack<NetworKit::node> path;
                    NetworKit::node current = w;
                    
                    // Build path from w to v
                    while (current != v) {
                        path.push(current);
                        current = parent[current];
                    }
                    
                    // Process blocks along the path
                    std::list<BoyerBlock*> pathBlocks;
                    while (!path.empty()) {
                        NetworKit::node u = path.top();
                        path.pop();
                        
                        if (blocks.find(u) != blocks.end()) {
                            pathBlocks.push_back(blocks[u]);
                        }
                    }
                    
                    // Merge blocks if necessary
                    if (pathBlocks.size() > 1) {
                        BoyerBlock* firstBlock = pathBlocks.front();
                        for (auto it = std::next(pathBlocks.begin()); it != pathBlocks.end(); ++it) {
                            firstBlock->merge(**it);
                        }
                    }
                    
                    // Flip blocks if necessary
                    for (auto block : pathBlocks) {
                        if (block->isFlipped()) {
                            block->flip();
                        }
                    }
                }
                
                // Walkdown phase
                for (auto& e : backEdges) {
                    NetworKit::node w = e.u == v ? e.v : e.u;
                    std::stack<NetworKit::node> path;
                    NetworKit::node current = w;
                    
                    // Build path from w to v
                    while (current != v) {
                        path.push(current);
                        current = parent[current];
                    }
                    
                    // Try to embed the back edge
                    bool embedded = false;
                    std::list<BoyerBlock*> pathBlocks;
                    
                    // Collect blocks along the path
                    while (!path.empty()) {
                        NetworKit::node u = path.top();
                        path.pop();
                        
                        if (blocks.find(u) != blocks.end()) {
                            pathBlocks.push_back(blocks[u]);
                        }
                    }
                    
                    // Try to embed in each block
                    for (auto block : pathBlocks) {
                        // Check if block can be flipped to embed the edge
                        if (!block->isFlipped()) {
                            // Check if w is externally active with respect to this block
                            bool isExternallyActive = false;
                            for (auto child : block->getChildren()) {
                                if (externallyActiveVertices[child].find(w) != 
                                    externallyActiveVertices[child].end()) {
                                    isExternallyActive = true;
                                    break;
                                }
                            }
                            
                            if (!isExternallyActive) {
                                block->flip();
                                embedded = true;
                                break;
                            }
                        }
                    }
                    
                    if (!embedded) {
                        is_planar = State::NOT_PLANAR;
                        return;
                    }
                }
            }
            
            is_planar = State::PLANAR;
        }

        void  BoyerMyrvold::dfs_in_reorder(NetworKit::Graph &g, NetworKit::node v,
                        std::vector<NetworKit::Edge> &Add, int &dfs_count,
                        std::vector<bool> &reached,
                        std::vector<NetworKit::node> &dfsnum,
                        std::vector<NetworKit::node> &lowpt,
                        std::vector<NetworKit::node> &lowpt2,
                        std::vector<NetworKit::node> &parent) {
            dfsnum[v] = dfs_count++;
            lowpt2[v] = lowpt[v] = dfsnum[v];
            reached[v] = true;
            
            for (auto w : g.neighborRange(v)) {
                if (!reached[w]) {
                    Add.push_back(NetworKit::Edge(v, w));
                    parent[w] = v;
                    dfs_in_reorder(g, w, Add, dfs_count, reached, dfsnum, lowpt, lowpt2, parent);
                    lowpt[v] = std::min(lowpt[v], lowpt[w]);
                } else {
                    lowpt[v] = std::min(lowpt[v], dfsnum[w]);
                    if (dfsnum[w] < dfsnum[v] && parent[v] != w) {
                        Add.push_back({v, w});
                    }
                }
            }
            
            for (auto w : g.neighborRange(v)) {
                if (parent[w] == v) {
                    if (lowpt[w] != lowpt[v]) {
                        lowpt2[v] = std::min(lowpt2[v], lowpt[w]);
                    }
                    lowpt2[v] = std::min(lowpt2[v], lowpt2[w]);
                } else {
                    if (lowpt[v] != dfsnum[w]) {
                        lowpt2[v] = std::min(lowpt2[v], dfsnum[w]);
                    }
                }
            }
    };

} /* namespace Koala */ 