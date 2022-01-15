include("staticSP.jl")

function test()
    dir = "../Instances/"
    fileName = "20_USA-road-d.BAY.gr"

    # charging data (variables are global !!! )
    include(dir * fileName)

    path, vertices = cplexSolveStaticSP()
    verifyStaticSP(path, vertices)
end