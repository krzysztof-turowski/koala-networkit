#include <graph/GraphTools.hpp>

#include "cograph_rec/Cotree.hpp"

namespace Koala {
    void Cotree::SubtreeColors(long long v) {
        if (left_son[v] != -1) {
            SubtreeColors(left_son[v]);
        }

        if (right_son[v] != -1) {
            SubtreeColors(right_son[v]);
        }

        if (type[v] == 2) {
            color[v]=0;
            number_of_colors[v]=1;
        } else {
            if(type[v]==1)
            {
                if (left_son[v] != -1) {
                    number_of_colors[v]+= number_of_colors[left_son[v]];
                }

                if (right_son[v] != -1) {
                    color[right_son[v]]+=number_of_colors[left_son[v]];
                    number_of_colors[v]+= number_of_colors[right_son[v]];
                }
            }
            else
            {
                if (left_son[v] != -1) {
                    number_of_colors[v]=number_of_colors[left_son[v]];
                }

                if (right_son[v] != -1 && number_of_colors[right_son[v]]>number_of_colors[left_son[v]]) {
                    number_of_colors[v]= number_of_colors[right_son[v]];
                }
            }

        }
    }

    void Cotree::EndOfColoring(long long v)
    {
                if (left_son[v] != -1) {
                    color[left_son[v]]+=color[v];
                    EndOfColoring(left_son[v]);
                }

                if (right_son[v] != -1) {
                    color[right_son[v]]+=color[v];
                    EndOfColoring(right_son[v]);
                }
    }

    void Cotree::Coloring() {
        long long n = graph->numberOfNodes();
        if (prepared == false) {
            BuildTree();
        }

        for(int i=0;i<2*n+2;i++)
        {
            color.push_back(0);
            number_of_colors.push_back(0);
        }

        SubtreeColors(n);
        EndOfColoring(n);
    }

    long long Cotree::GetColor(long long i)
    {
        if(prepared==false)
        {
            BuildTree();
        }
        return color[i];
    }

    bool Cotree::СheckColoring() {
        long long n=graph->numberOfNodes(),i,j;

        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                if(graph->hasEdge(i,j) && GetColor(i)== GetColor(j))
                {
                    return false;
                }
            }
        }
        return true;
    }
}