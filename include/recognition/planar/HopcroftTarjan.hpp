/*
 * IPlanarGraphRecognition.hpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include "recognition/planar/IPlanarGraphRecognition.hpp"

namespace Koala {

const int left = 1;
const int right = 2;
/**
 * @ingroup recognition
 * The Interface for class for recognition of planar graphs p.
 *
 */

class Block {
   private:
      std::list<int> Latt, Ratt;
      std::list<NetworKit::WeightedEdge> Lseg, Rseg;
   public:
      Block(NetworKit::WeightedEdge e, std::list<int>& A){
         Lseg.push_back(e);
         Latt.splice(Latt.end(), A);
      }
      ~Block(){}
      void flip() {
         std::list<int> ha;
         std::list<NetworKit::WeightedEdge> he;
         ha.splice(ha.end(), Ratt); Ratt.splice(Ratt.end(), Latt); Latt.splice(Latt.end(), ha);
         he.splice(he.end(), Rseg); Rseg.splice(Rseg.end(), Lseg); Lseg.splice(Lseg.end(), he);
      }
      int head_of_Latt() {
         return Latt.front();
      }
      bool empty_Latt() {
         return Latt.empty();
      }
      int head_of_Ratt() {
         return Ratt.front();
      }
      bool empty_Ratt() {
         return Ratt.empty();
      }
      bool left_interlace(std::stack<Block*> &s){
         if(!s.empty() && !((s.top())->empty_Latt()) && Latt.back() < (s.top())->head_of_Latt()){
            return true;
         }
         return false;
      }

      bool right_interlace(std::stack<Block*> &s){
         if(!s.empty() && !((s.top())->empty_Ratt()) && Latt.back() < (s.top())->head_of_Ratt()){
            return true;
         }
         return false;
      }

      void combine(Block *& b){
         Latt.splice(Latt.end(), b->Latt);
         Ratt.splice(Ratt.end(), b->Ratt);
         Lseg.splice(Lseg.end(), b->Lseg);
         Rseg.splice(Rseg.end(), b->Rseg);
         delete b;
      }

      bool clean(int dfsnum_w, std::vector<int> & alpha, std::vector<NetworKit::node> & dfsnum) {
         while(!Latt.empty() && Latt.front() == dfsnum_w) {
            Latt.pop_front();
         }
          while(!Ratt.empty() && Ratt.front() == dfsnum_w) {
            Ratt.pop_front();
         }
         if(!Latt.empty() || !Ratt.empty()){
            return false;
         }
         for(auto u: Lseg){
            alpha[u.weight] = left;
         }
         for(auto u: Rseg){
            alpha[u.weight] = right;
         }
         return true;
      }

      void add_to_Att(std::list<int> &Att, int dfsnum_w0, std::vector<int> & alpha) {
         if(!Ratt.empty() && head_of_Ratt() > dfsnum_w0) flip();
         Att.splice(Att.end(), Latt);
         Att.splice(Att.end(), Ratt);
         for(auto u: Lseg){
            alpha[u.weight] = left;
         }
         for(auto u: Rseg){
            alpha[u.weight] = right;
         }
      }

};

class HopcroftTarjan : public IPlanarGraphRecognition {
   public:
      HopcroftTarjan::HopcroftTarjan(NetworKit::Graph &graph, bool embedding);
      void run() override;
   private:
      NetworKit::count N;
     

      void dfs(NetworKit::Graph &g, NetworKit::node v, std::vector<bool> & reached);
      void dfs_in_make_biconnected_graph(NetworKit::Graph &g, NetworKit::node v, 
             int &dfs_count, std::vector<bool> &reached, std::vector<NetworKit::node> &dfsnum,
             std::vector<NetworKit::node> &lowpt, std::vector<NetworKit::node> &parent);
      
      void make_biconnected_graph(NetworKit::Graph &g);  

      void dfs_in_reorder(NetworKit::Graph &g, 
                           NetworKit::node v, 
                           std::vector<NetworKit::Edge> & Add,
                           int & dfs_count, 
                           std::vector<bool> & reached,
                           std::vector<NetworKit::node> & dfsnum,
                           std::vector<NetworKit::node> &lowpt, 
                           std::vector<NetworKit::node> &lowpt2, 
                           std::vector<NetworKit::node> &parent);

      void reorder(NetworKit::Graph &g, 
                std::vector<NetworKit::node> & parent, 
                std::vector<NetworKit::node> & dfsnum,
                 std::vector<NetworKit::Edge> & order
                );
      
      bool strong_planar(NetworKit::Graph &g, 
                       NetworKit::Edge e, 
                       std::list<int> & Att,
                       std::vector<int> & alpha,
                       std::vector<NetworKit::node> & dfsnum,
                       std::vector<NetworKit::node> & parent
                    );

      bool planar();
};

} /* namespace Koala */
