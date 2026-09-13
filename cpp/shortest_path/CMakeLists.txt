koala_add_module(shortest_path
    FredericksonPlanarSSSP.cpp
    PlanarSSSP.cpp
    HenzingerPlanarSSSP.cpp
    Aingworth.cpp
    Seidel.cpp
    ShoshanZwick.cpp
    Spira.cpp
)

koala_add_module(planar
    planar/HenzingerPriorityQueue.cpp
    planar/SuitableRDivision.cpp
    planar/TopologyBasedHeap.cpp
)
