koala_add_module(perfect_graph_recognition
    PerfectGraphRecognition.cpp
    perfect/OddHoles.cpp
    perfect/Jewels.cpp
    perfect/Pyramids.cpp
    perfect/NearCleaners.cpp
)

koala_add_module(planar_graph_recognition
    PlanarGraphRecognition.cpp
)

koala_add_module(cograph_recognition
    CographRecognition.cpp
    cograph/CorneilStewartPerlCographRecognition.cpp
    cograph/BretscherCorneilHabibPaulCographRecognition.cpp
    cograph/DahlhausCographRecognition.cpp
    cograph/HabibPaulCographRecognition.cpp
    cograph/FactorizingPermutation.cpp
    cograph/Twins.cpp
)
