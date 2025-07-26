#pragma once

#include <cmath>
#include <random>
#include <iostream>
#include <map>

#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ_p.h>
#include <networkit/graph/Graph.hpp>

namespace Koala {
    constexpr double EPS = 1e-8;
    typedef NTL::ZZ_p Zp;
    typedef NTL::Vec<Zp> VecZp;
    typedef NTL::Mat<Zp> MatZp;

    inline bool eq(double a, double b) {
        return fabs(a - b) <= EPS;
    }

    inline void printGraph(const NetworKit::Graph& G) {
        std::cerr << "v\t";
        for (auto v: G.nodeRange()) {
            std::cerr << v << ' ';
        };
        std::cerr << '\n';
        for (auto [u,v]: G.edgeRange()) {
            std::cerr << "e\t" << u << ' ' << v << std::endl;
        };
    }

    inline void printMatrix(const MatZp& A) {
        std::cerr << '\n';
        for (int i = 0; i < A.NumCols(); ++i) {
        for (int j = 0; j < A.NumRows(); ++j) {
            std::cerr << A[i][j] << ' ';
        }
        std::cerr << '\n';
        }
    }

    inline void initZp(int p) {
        Zp::init(NTL::conv<NTL::ZZ>(p));
    }

    inline Zp generateRandom() {
        return Zp(rand());
    }

    inline VecZp getCol(const MatZp& A, int c) {
        VecZp col;
        col.SetLength(A.NumRows());
        for (long i = 0; i < col.length(); ++i) {
            col[i] = A[i][c];
        }
        return col;
    }

    inline void setCol(MatZp& A, int c, const VecZp& vec) {
        for (long i = 0; i < vec.length(); ++i) {
            A[i][c] = vec[i];
        }
    }


    inline VecZp getRow(const MatZp& A, int r) {
        VecZp row;
        row.SetLength(A.NumCols());
        for (long i = 0; i < row.length(); ++i) {
            row[i] = A[r][i];
        }
        return row;
    }

    inline void setRow(MatZp& A, int r, const VecZp& vec) {
        for (long i = 0; i < vec.length(); ++i) {
            A[r][i] = vec[i];
        }
    }

    inline VecZp divVec(const VecZp& vec, Zp a) {
        VecZp res;
        res.SetLength(vec.length());
        for (long i = 0; i < vec.length(); ++i) {
            res[i] = vec[i] / a;
        }
        return res;
    }

    inline VecZp segment(const VecZp & vec, long b, long e) {
        VecZp  seg;
        seg.SetLength(e-b);
        for (long i = b; i < e; ++i) {
            seg[i-b] = vec[i];
        }
        return seg;
    }

    inline VecZp zeroVec(int n) {
        VecZp vec;
        vec.SetLength(n);
        return vec;
    }

    inline MatZp zeroMat(int r, int c) {
        MatZp A;
        A.SetDims(r, c);
        return A;
    }


    inline void subBlock(MatZp& A, long r1, long c1, const MatZp& B) {
        for (long i = 0; i <  B.NumRows(); ++i) {
            for (long j = 0; j < B.NumCols(); ++j) {
                A[r1 + i][c1 + j] -= B[i][j];
            }
        }
    }

    inline std::pair<NetworKit::Graph, std::vector<int>> reindexGraph(const NetworKit::Graph& G) {
        int n = G.numberOfNodes();
        NetworKit::Graph G1(n, false, false);
        std::map<int, int> indexes;
        std::vector<int> labels(n);
        for (int i=0; auto v: G.nodeRange()) {
            indexes[v] = i;
            labels[i] = v;
            i++;
        }
        for (auto [u,v] : G.edgeRange()) { G1.addEdge(indexes[u], indexes[v]); }
        return { G1, labels };
    }

    inline std::tuple<NetworKit::Graph, std::vector<int>, MatZp> reindexGraph(const NetworKit::Graph& G, const MatZp& AG) {
        int n = G.numberOfNodes();
        auto [G1, labels] = reindexGraph(G);

        auto AG1 = zeroMat(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                AG1[i][j] = AG[labels[i]][labels[j]];
            }
        }
        
        return { G1, labels, AG1 };
    }
}