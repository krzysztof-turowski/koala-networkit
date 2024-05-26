#include <graph/GraphTools.hpp>

#include "recognition/сograph/FactorizingPermutation.hpp"

namespace Koala {
    FactorizingPermutation::FactorizingPermutation() {
        l = nullptr;
        r = nullptr;
        first_part = nullptr;
        last_part = nullptr;
        origin = nullptr;
        last_used_part = 0;
    }

    FactorizingPermutation::~FactorizingPermutation() {
        E.clear();
        P.clear();
    }

    NetworKit::count FactorizingPermutation::newPart() {
        if (P[last_used_part].size == 0) {
            return last_used_part;
        }
        last_used_part++;
        return last_used_part;
    }


    void FactorizingPermutation::lCheck() {
        if (origin->my_part != l && l->next != nullptr) {
            if (l->next->first != l->next->last) {
                l = (l->next);
            }
        }
        while (l->previous != nullptr && (l->first == l->last || l->first == nullptr)) {
            l = l->previous;
        }
    }

    void FactorizingPermutation::rCheck() {
        if (origin->my_part != r && r->previous != nullptr) {
            if (r->previous->first != r->previous->last) {
                r = r->previous;
            }
        }

        while (r->next != nullptr && (r->first == r->last || r->first == nullptr)) {
            r = r->next;
        }
    }

    void FactorizingPermutation::addElementToPart(part &P, element &v) {
        P.size++;
        v.my_part = &P;
        if (P.first == nullptr) {
            P.first = &v;
            P.last = &v;
            v.previous = nullptr;
            v.next = nullptr;
        } else {
            v.previous = P.last;
            v.next = nullptr;
            P.last->next = &v;
            P.last = &v;
        }
    }

    void FactorizingPermutation::eraseElement(element &v) {
        element *A, *B;
        A = v.previous;
        B = v.next;
        v.my_part->size--;
        if (A != nullptr) {
            A->next = B;
        } else {
            v.my_part->first = B;
        }

        if (B != nullptr) {
            B->previous = A;
        } else {
            v.my_part->last = A;
        }

        if (v.my_part->pivot == &v) {
            v.my_part->pivot = nullptr;
        }
    }

    void FactorizingPermutation::AddPart(part &x, part &y) {
        part *A = x.next;
        x.next = &y;
        y.previous = &x;
        y.next = A;
        if (A != nullptr) {
            A->previous = &y;
        } else {
            last_part = &y;
        }
    }

    void FactorizingPermutation::erasePart(part &p) {
        lCheck();
        rCheck();

        part *A, *B;
        A = p.previous;
        B = p.next;
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
} /* namespace Koala */
