#pragma once

#include <networkit/Globals.hpp>
#include <vector>
#include <list>

namespace Koala {

template<typename Event>
class ArrayPriorityQueue {
 public:
    ArrayPriorityQueue() { reset(); }

    void scheduleEvent(int time, Event event) {
        if (time < time) return;
        if (time >= event_queue.size())
            event_queue.resize(time * 2);
        event_queue[time].push_back(event);
    }

    Event getEvent() {
        advanceToNext();
        return event_queue[time][it];
    }

    int timeNow() { return time; }

    void reset() {
        event_queue.clear();
        event_queue.resize(4);
        it = -1;
        time = 0;
    }

 private:
    NetworKit::index time;
    NetworKit::index it;
    std::vector<std::vector<Event>> event_queue;

    void advanceToNext() {
        it++;
        while (it >= event_queue[time].size()) {
            time++;
            it = 0;
        }
    }
};

} /* namespace Koala */
