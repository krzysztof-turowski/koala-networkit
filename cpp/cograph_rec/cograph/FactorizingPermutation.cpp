
#include <list>
#include <set>
#include <unordered_map>
#include <map>
#include <graph/GraphTools.hpp>

#include "cograph_rec/FactorizingPermutation.hpp"

namespace Koala {

    FactorizingPermutation::FactorizingPermutation() {
        l = nullptr;
        r = nullptr;
        first_part = nullptr;
        last_part = nullptr;
        origin = nullptr;
    }

    FactorizingPermutation::~FactorizingPermutation() {
        for (auto v: garbage) {
            delete v;
        }

        for (auto v: E) {
            delete v;
        }

        garbage.clear();
        E.clear();
    }

    void FactorizingPermutation::Lcheck() {
        if (origin->my_part != l && l->next != nullptr) {
            if (l->next->first != l->next->last) {
                l = (l->next);
            }
        }

        while (l->previous != nullptr && (l->first == l->last || l->first == nullptr)) {
            l = l->previous;
        }

    }

    void FactorizingPermutation::Rcheck() {
        if (origin->my_part != r && r->previous != nullptr) {
            if (r->previous->first != r->previous->last) {
                r = r->previous;
            }
        }

        while (r->next != nullptr && (r->first == r->last || r->first == nullptr)) {
            r = r->next;
        }
    }

    void FactorizingPermutation::AddElementToPart(part *P, element *v) {
        P->size++;
        v->my_part = P;
        if (P->first == nullptr) {
            P->first = v;
            P->last = v;
            v->previous = nullptr;
            v->next = nullptr;
        } else {
            v->previous = P->last;
            v->next = nullptr;
            P->last->next = v;
            P->last = v;
        }
    }

    void FactorizingPermutation::EraseElement(element *v) {
        element *A, *B;
        A = v->previous;
        B = v->next;
        v->my_part->size--;
        if (A != nullptr) {
            A->next = B;
        } else {
            v->my_part->first = B;
        }

        if (B != nullptr) {
            B->previous = A;
        } else {
            v->my_part->last = A;
        }

        if (v->my_part->pivot == v) {
            v->my_part->pivot = nullptr;
        }

    }

    void FactorizingPermutation::AddPart(part *x, part *y) {
        part *A = x->next;
        x->next = y;
        y->previous = x;
        y->next = A;
        if (A != nullptr) {
            A->previous = y;
        } else {
            last_part = y;
        }
    }

    void FactorizingPermutation::ShowTheOrder() {
        part *x;
        element *y, *z;
        x = first_part;
        while (true) {
            y = x->first;
            if (y != nullptr) {
                std::cout << y->num << " ";
                while (y->next != nullptr) {
                    y = y->next;
                    std::cout << y->num << " ";
                }
            }
            std::cout << " pivot=";
            if (x->pivot == nullptr) {
                std::cout << "nullptr";
            } else {
                std::cout << x->pivot->num;
            }
            std::cout << " amount=" << x->amount << std::endl;
            std::cout << " size=" << x->size << std::endl;
            if (x != last_part) {
                x = x->next;
            } else {
                break;
            }


        }
        std::cout << std::endl;
    }

    void FactorizingPermutation::ErasePart(part *p) {

        Lcheck();
        Rcheck();

        if (p->first == nullptr) {
            garbage.push_back(p);
        }

        part *A, *B;
        A = p->previous;
        B = p->next;
        if (A == nullptr) {
            first_part = B;
        } else {
            A->next = B;
        }

        if (B == nullptr) {
            last_part = A;
        } else {
            B->previous = A;
        }

    }
}