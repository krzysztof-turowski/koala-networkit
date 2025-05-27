koala_add_module(recognition
    PerfectGraphRecognition.cpp
    perfect/OddHoles.cpp
    perfect/Jewels.cpp
    perfect/Pyramids.cpp
    perfect/NearCleaners.cpp
    planar/PlanarGraphRecognition.cpp
    planar/HopcroftTarjan.cpp
    planar/LeftRight.cpp
)

koala_add_module(cograph
    CographRecognition.cpp
    cograph/CorneilStewartPerlCographRecognition.cpp
    cograph/BretscherCorneilHabibPaulCographRecognition.cpp
    cograph/DahlhausCographRecognition.cpp
    cograph/CoTree.cpp
)

koala_add_module(cograph_other
    CographRecognitionOther.cpp
    cograph/FactorizingPermutation.cpp
    cograph/Twins.cpp
)
