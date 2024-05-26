#pragma once

#include <list>
#include <graph/GraphTools.hpp>

#include "Part.hpp"

namespace Koala {
class FactorizingPermutation {
 public:
    part *l, *r, *first_part, *last_part;
    element *origin;
    std::vector<element> E;
    std::vector<part> P;
    NetworKit::count last_used_part;

    FactorizingPermutation();

    ~FactorizingPermutation();

    void check() const;

    void AddPart(part &l, part &r);

    void eraseElement(element &v);

    void erasePart(part &p);

    void addElementToPart(part &P, element &v);

    void rCheck();

    void lCheck();

    void ShowTheOrder();

    NetworKit::count newPart();
};
} /* namespace Koala */
