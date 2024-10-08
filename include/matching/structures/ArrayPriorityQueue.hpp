#pragma once

#include <networkit/Globals.hpp>
#include <vector>
#include <list>

namespace Koala {

template<typename Event>
class ArrayPriorityQueue {
 public:
    explicit ArrayPriorityQueue(NetworKit::count bucket_size): bucket_size(bucket_size) { reset(); }

    void scheduleEvent(NetworKit::index event_time, Event event) {
        if (event_time < timeNow()) return;
        auto bucket = event_time / bucket_size;
        if (bucket == current_bucket)
            add_to_current_bucket(event_time % bucket_size, event);
        else
            add_to_bucket(bucket, event_time, event);
    }

    Event getEvent() {
        advanceToNext();
        return event_queue[time][it];
    }

    int timeNow() const {
        return current_bucket * bucket_size + time;
    }

    void reset() {
        event_queue.clear();
        buckets.clear();
        it = -1;
        time = 0;
        current_bucket = 0;
    }

 private:
    NetworKit::count bucket_size;
    NetworKit::index current_bucket;
    NetworKit::index time;
    NetworKit::index it;
    std::vector<std::vector<Event>> event_queue;
    std::vector<std::vector<std::pair<NetworKit::index, Event>>> buckets;

    void advanceToNext() {
        it++;
        while (it >= event_queue[time].size()) {
            time++;
            it = 0;
            if (time == event_queue.size()) advance_bucket();
        }
    }

    void advance_bucket() {
        event_queue.clear();
        event_queue.resize(bucket_size);
        current_bucket++;
        for (auto [t, event] : buckets[current_bucket]) {
            event_queue[t % bucket_size].push_back(event);
        }
        buckets[current_bucket].clear();
        time = 0;
    }

    void add_to_current_bucket(NetworKit::index index, Event event) {
        if (index >= event_queue.size()) {
            event_queue.resize(std::min(std::max(2 * event_queue.size(), index + 1), bucket_size));
        }
        event_queue[index].push_back(event);
    }

    void add_to_bucket(NetworKit::index bucket, NetworKit::index event_time, Event event) {
        if (bucket >= buckets.size()) {
            buckets.resize(std::max(2 * buckets.size(), bucket + 1));
        }
        buckets[bucket].emplace_back(event_time % bucket_size, event);
    }
};

} /* namespace Koala */
