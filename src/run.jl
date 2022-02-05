include("staticSP.jl")
include("cuttingPlanes.jl")
include("heuristics.jl")
include("dual.jl")



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

    files = ["20_USA-road-d.BAY.gr", "20_USA-road-d.COL.gr", "20_USA-road-d.NY.gr", "40_USA-road-d.BAY.gr"]
for fileName in files
    preparation(dir, fileName)

    # cplex solve    
    # println("static")
    # solStatic = cplexSolveStaticSP()
    # sleep(5)

    # heurisric robustness
    println("heuristic")
    heuristicPrimal()
    sleep(5)


    # dual 
    println("dual")
    solRobustDual = dualSolve()
    sleep(5)


    # println("cuting planes exact")
    # solCPExact = cuttingPlanes(false)


    # println("cuting planes exact --- Influence of init senarios")
    # for choix in [0 1 2]
    #     solCPExact = cuttingPlanes(false, choix)
    #     println("choix = ", choix, " solveTime : ", solCPExact.solveTime)
    #     sleep(5)
    # end


    # println("cuting planes heurisric")
    # solCPHeur = cuttingPlanes(true)
    # println("The gap between solCPExact and solCPHeur is ", abs(solCPExact.obj_val - solCPHeur.obj_val)/solCPExact.obj_val)

    # println("Branch and Cut exact")
    # solBCExact = brunchAndCut(true, false)


    # println("Branch and Cut heurisric")
    # solBCHeur = brunchAndCut(false, true)
    # println("The gap between solBCExact and solBCHeur is ", abs(solBCExact.obj_val - solBCHeur.obj_val)/solBCExact.obj_val)
    #sleep(5)
end
    


    


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