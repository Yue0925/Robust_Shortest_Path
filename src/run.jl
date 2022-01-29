include("staticSP.jl")
include("cuttingPlanes.jl")


"""
Prepare the global variables used in the project
"""
function preparation(dir::String, fileName::String)
    # charging data (variables are global !!! )
    include(dir * fileName)

    # matrix adjacecy 
    global Adjacenct = falses(n, n) 
    global arcs = size(Mat, 1)
    for i in 1:arcs
        Adjacenct[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] = true
    end

end


function test()
    dir = "../Instances/"

    fileName = "20_USA-road-d.BAY.gr"

    preparation(dir, fileName)

    # path, vertices = cuttingPlanes()
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    # println("sub-problem exact")
    # path, vertices = brunchAndCut(true, false)
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    
    println("sub-problem heurisric")
    path, vertices = brunchAndCut(false, true)
    isFeasible = verifyStaticSP(path, vertices)
    println("isFeasible? ", isFeasible)

    
    # # cplex solve
    # path, vertices = cplexSolveStaticSP()
    # # verification
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    # files = ["20_USA-road-d.BAY.gr", "20_USA-road-d.COL.gr", "20_USA-road-d.NY.gr", "40_USA-road-d.BAY.gr"]
    # for file in files #readdir(dir)
    #     println(dir * file)
    #     preparation(dir, file)
    #     # cplex solve
    #     path, vertices, isOptimal = cplexSolveStaticSP()
    #     # verification
    #     isFeasible = verifyStaticSP(path, vertices)
    #     println("isFeasible? ", isFeasible)

    #     if isFeasible != isOptimal
    #         println("Error instance : ", file)
    #         break
    #     end
    # end

end