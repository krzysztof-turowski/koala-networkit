koala_add_module(flow
    MaximumFlow.cpp
    PushRelabelMaximumFlow.cpp
    KingRaoTarjanMaximumFlow.cpp
    maximum_flow/KrtEdgeDesignator.cpp
    GoldbergTarjanPushRelabelMaximumFlow.cpp
    MalhotraKumarMaheshwariFlow.cpp
    BoykovKolmogorovMaximumFlow.cpp
    minimum_cost_flow/MCFlowNetwork.cpp
    minimum_cost_flow/EdmondsKarpMinimumCostFlow.cpp
    minimum_cost_flow/OrlinMinimumCostFlow.cpp
    minimum_cost_flow/SuccessiveApproximationMinimumCostFlow.cpp

    # Electrical Flow
    ElectricalFlow.cpp
    electrical_flow/ElectricalNetwork.cpp
    electrical_flow/FlowNetwork.cpp 
    electrical_flow/NaiveLaplaceSolver.cpp
)
