#ifndef BETTER_MAX_FLOW_EDGE_DESIGNATOR_HPP
#define BETTER_MAX_FLOW_EDGE_DESIGNATOR_HPP

class EdgeDesignator {
public:
    virtual int ce(int, int) = 0;

    virtual void response_adversary(int, int, int, int, bool) = 0;
};

#endif //BETTER_MAX_FLOW_GENERIC_EDGE_DESIGNATOR_H
