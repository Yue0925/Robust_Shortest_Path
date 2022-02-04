include("staticSP.jl")
include("cuttingPlanes.jl")
include("heuristics.jl")


"""
Prepare the global variables used in the project
"""
function preparation(dir::String, fileName::String)
    # charging data (variables are global !!! )
    # TODO : unfortunately this include works only for small insyances
    #include(dir * fileName)

    # Open the input file
    datafile = open(dir * fileName)

    global n = parse(Int64, split(readline(datafile), " ")[3])
    global s = parse(Int64, split(readline(datafile), " ")[3])
    global t = parse(Int64, split(readline(datafile), " ")[3])
    global S = parse(Int64, split(readline(datafile), " ")[3])
    global d1 = parse(Int64, split(readline(datafile), " ")[3])
    global d2 = parse(Int64, split(readline(datafile), " ")[3])

    global p = Array{Int64, 1}(undef, 0) 
    global ph = Array{Int64, 1}(undef, 0)
    global Mat = Array{Float64, 2}(undef, 0, 4)

    data = readlines(datafile)
    close(datafile)

    isMatLine = false

    # reading tables
    for eachLine in data
        line = split(eachLine, " ")

        # reading p
        if line[1] == "p"
            for i in 3:size(line, 1)
                if i == 3 # the start
                    append!(p, parse(Int64, chop(line[i], head = 1, tail = 1)))
                else
                    append!(p, parse(Int64, chop(line[i])))
                end
            end
        end

        # reading ph
        if line[1] == "ph"
            for i in 3:size(line, 1)
                if i == 3 # the start
                    append!(ph, parse(Int64, chop(line[i], head = 1, tail = 1)))
                else
                    append!(ph, parse(Int64, chop(line[i])))
                end
            end
        end

        # reading Mat
        if line[1] == "Mat"
            isMatLine = true
            continue
        end
        if isMatLine
            Mat = vcat(Mat, [parse(Int64, line[1]) parse(Int64, line[2]) parse(Int64, line[3]) parse(Float64, chop(line[4]))])
        end
    end
    

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
    heuristicPrimal()

    # println("cuting planes exact")
    # path, vertices = cuttingPlanes(false, 2)
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)


    # println("cuting planes heurisric")
    # path, vertices = cuttingPlanes(true)
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    # println("sub-problem exact")
    # path, vertices = brunchAndCut(true, false)
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    
    # println("sub-problem heurisric")
    # path, vertices = brunchAndCut(false, true)
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

    
    # # cplex solve
    # path, vertices = cplexSolveStaticSP()
    # # verification
    # isFeasible = verifyStaticSP(path, vertices)
    # println("isFeasible? ", isFeasible)

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