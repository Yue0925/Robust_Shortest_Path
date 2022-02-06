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

    # ["40_USA-road-d.BAY.gr", "40_USA-road-d.COL.gr", "40_USA-road-d.NY.gr"]
    #["20_USA-road-d.BAY.gr", "20_USA-road-d.COL.gr", "20_USA-road-d.NY.gr"]
    # ["60_USA-road-d.BAY.gr", "60_USA-road-d.COL.gr", "60_USA-road-d.NY.gr"]

    files = ["100_USA-road-d.COL.gr"]
for fileName in files
println(fileName)
preparation(dir, fileName)

    # cplex solve    
    println("static")
    solStatic = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
    try
        solStatic = cplexSolveStaticSP()
    catch
        println("CPLEX ERROR")
    end
    println(solStatic)
    sleep(5)

    # # heurisric robustness
    # println("heuristic")
    # heuristicPrimal()
    # sleep(5)


    # dual 
    # println("dual")
    # solRobustDual = dualSolve()
    # sleep(5)


    # println("cuting planes exact")
    # solCPExact = cuttingPlanes(false)
    # sleep(5)


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
    # sleep(5)


    # println("Branch and Cut heurisric")
    # solBCHeur = brunchAndCut(false, true)
    # println("The gap between solBCExact and solBCHeur is ", abs(solBCExact.obj_val - solBCHeur.obj_val)/solBCExact.obj_val)
    #sleep(5)
 end
    
    # files = sort(readdir(dir), by = x -> parse(Int64, split(x, "_")[1]))
    # println(files)

end


function generateSubHeurInfl()
    dir = "../Instances/"
    files = sort(readdir(dir), by = x -> parse(Int64, split(x, "_")[1]))
    city = 0

    dir_res = "../res/subProblemsHeur/"
    KEYS = ["CPExact" "CPHeur" "BCExact" "BCHeur"]
    instances = Dict("CPExact" => [], "CPHeur" => [], "BCExact" => [], "BCHeur" => [])
    bounds = Dict("CPExact" => [], "CPHeur" => [], "BCExact" => [], "BCHeur" => [])
    times = Dict("CPExact" => [], "CPHeur" => [], "BCExact" => [], "BCHeur" => [])

    fout = open(dir_res * "Comparaison_subproblems_heuristics.tex", "w")

    latex = raw"""\begin{sidewaystable}[!h]
    \centering
    \hspace*{-1cm}\begin{tabular}{lcccccccccccccccccc}
    \toprule
    \textbf{Instance}  & \multicolumn{3}{c}{\textbf{PC Exact}} & \multicolumn{4}{c}{\textbf{PC Heur}}  & \multicolumn{3}{c}{\textbf{BC Exact}} & \multicolumn{4}{c}{\textbf{BC Heur}}
    \\
    \cmidrule(r){2-4} \cmidrule(r){5-8} \cmidrule(r){9-11} \cmidrule(r){12-15}
     & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Gap Heur} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Gap Heur}  \\
    \midrule
"""    
    println(fout, latex)

    #files = ["1000_USA-road-d.BAY.gr", "1000_USA-road-d.COL.gr", "1000_USA-road-d.NY.gr"]

    for fileName in files
        println("---------------------")
        println(fileName)
        println("---------------------")

        countFails = 0


        # for each instances
        # reading data
        preparation(dir, fileName)

        print(fout, n, " & ")

        nn = parse(Int64, split(fileName, "_")[1])
        if n - city >3
            city = nn
        else
            city += 1
        end


        println("\n cuting planes exact \n")
        solCPExact = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        try
            solCPExact = cuttingPlanes(false)
            print(fout, round(solCPExact.solveTime, digits=2), " & ", round(solCPExact.GAP, digits=2), " & ")

            if solCPExact.isFeasible
                append!(instances["CPExact"], city)
                append!(bounds["CPExact"], round(solCPExact.best_bound, digits=2))
                append!(times["CPExact"], round(solCPExact.solveTime, digits=2))
                print(fout, round(solCPExact.best_bound, digits=2), " & ")
            else
                print(fout, "-", " & ")
            end
        catch
            println("Time out ! ")
            print(fout, "- & - & - & ")
            countFails += 1
        end


        println("\n cuting planes heurisric\n ")
        solCPHeur = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        try
            solCPHeur = cuttingPlanes(true)
            print(fout, round(solCPHeur.solveTime, digits=2), " & ", round(solCPHeur.GAP, digits=2), " & ")

            if solCPHeur.best_bound > 0.0
                print(fout, round(solCPHeur.best_bound, digits=2), " & ")

                if solCPExact.obj_val > 0.0
                    gap_heur = abs(solCPExact.obj_val - solCPHeur.obj_val)/solCPExact.obj_val
                    print(fout, round(gap_heur, digits=2), " & ")
                else
                    print(fout, " - & ")
                    gap_heur = -1
                end

                append!(instances["CPHeur"], city)
                append!(bounds["CPHeur"], round(solCPHeur.best_bound, digits=2))
                append!(times["CPHeur"], round(solCPHeur.solveTime, digits=2))

                println("The gap between solCPExact and solCPHeur is ", gap_heur)
            else
                print(fout, "- & - & ")     
            end
        catch
            println("Time out ! ")
            print(fout, "- & - & - & - & ")
            countFails +=1
        end


        println("\n Branch and Cut exact \n")
        solBCExact = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        try
            solBCExact = brunchAndCut(true, false)
            print(fout, round(solBCExact.solveTime, digits=2), " & ", round(solBCExact.GAP, digits=2), " & ")
            if solBCExact.best_bound > 0.0
                append!(instances["BCExact"], city)
                append!(bounds["BCExact"], round(solBCExact.best_bound, digits=2))
                append!(times["BCExact"], round(solBCExact.solveTime, digits=2))
                print(fout, round(solBCExact.best_bound, digits=2), " & ")
            else
                print(fout, "-", " & ")
            end
        catch
            println("Time out ! ")
            print(fout, "- & - & - & ")
            countFails +=1
        end


        println("\n Branch and Cut heurisric \n")
        solBCHeur = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        try
            solBCHeur = brunchAndCut(false, true)
            println("solBCHeur\n", solBCHeur)
            print(fout, round(solBCHeur.solveTime, digits=2), " & ", round(solBCHeur.GAP, digits=2), " & ")

            if solBCHeur.best_bound > 0.0
                print(fout, round(solBCHeur.best_bound, digits=2), " & ")

                if solBCExact.obj_val > 0.0
                    gap_heur = abs(solBCExact.obj_val - solBCHeur.obj_val)/solBCExact.obj_val
                    print(fout, round(gap_heur, digits=2), " ")
                else
                    print(fout, " - ")
                    gap_heur = -1
                end
                
                println("The gap between solBCExact and solBCHeur is ", gap_heur)
                append!(instances["BCHeur"], city)
                append!(bounds["BCHeur"], round(solBCHeur.best_bound, digits=2))
                append!(times["BCHeur"], round(solBCHeur.solveTime, digits=2))
            else
                print(fout, "- & - ")
            end
        catch
            println("Time out ! ")
            print(fout, "- & - & - & - ")
            countFails +=1
        end

        println(fout, "\\\\")

        if countFails >= 4
            break
        end

    end

    latex = raw"""\bottomrule
    \end{tabular}
    \caption{.}
    \label{tab:}
    \end{sidewaystable}"""
    println(fout, latex)

    close(fout)

    for method in KEYS
        plot(instances[method], bounds[method], label = method)
    end

    legend()
    # yscale("log")
    # xscale("log")
    title("Comparison of the best objective bounds")
    xlabel("Cities")
    ylabel("Objective bounds")
    savefig(dir_res * "Comparaison_subproblems_heuristics_bounds.png")
    close()

    for method in KEYS
        plot(instances[method], times[method], label = method)
    end

    legend()
    # yscale("log")
    # xscale("log")
    title("Comparison of the computation times")
    xlabel("Cities")
    ylabel("Times(s)")
    savefig(dir_res * "Comparaison_subproblems_heuristics_times.png")
    close()
end


"""
Test on the cutting planes algorithm.
"""
function generateSenariosInfl()
    dir = "../Instances/"
    files = sort(readdir(dir), by = x -> parse(Int64, split(x, "_")[1]))
    files = ["750_USA-road-d.BAY.gr", "750_USA-road-d.COL.gr", "750_USA-road-d.NY.gr"]

    dir_res = "../res/PCSenariosInfl/"

    fout = open(dir_res * "Comparaison_initSenarios.tex", "w")

    latex = raw"""\begin{table}[!h]
    \centering
    \hspace*{-1cm}\begin{tabular}{lccccccccc}
    \toprule
    \textbf{Instance}  & \multicolumn{3}{c}{\textbf{Defaut}} & \multicolumn{3}{c}{\textbf{Uniform}}  & \multicolumn{3}{c}{\textbf{Best Bound}}
    \\
    \cmidrule(r){2-4} \cmidrule(r){5-7} \cmidrule(r){8-10}
     & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound}  \\
    \midrule"""    
    println(fout, latex)


    for fileName in files
        println("---------------------")
        println(fileName)
        println("---------------------")

        countFails = 0

        # for each instances
        # reading data
        preparation(dir, fileName)

        print(fout, n)

        # for each choice param
        for choix in [0 1 2]

            println("\n choix ", choix, "\n")
            solCPExact = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
            try
                solCPExact = cuttingPlanes(false, choix)
                print(fout, " & ", round(solCPExact.solveTime, digits=2), " & ", round(solCPExact.GAP, digits=2))
    
                if solCPExact.best_bound > 0.0
                    print(fout, " & ", round(solCPExact.best_bound, digits=2))
                else
                    print(fout, " & -")
                end
            catch
                println("Time out ! ")
                print(fout, " & - & - & -")
                countFails += 1
            end
        end
        println(fout, " \\\\")

        if countFails >= 3
            break
        end
    end

    latex = raw"""\bottomrule
    \end{tabular}
    \caption{.}
    \label{tab:init_senarios}
    \end{table}"""
    println(fout, latex)

    close(fout)
end


function totalReults()
    dir = "../Instances/"
    files = sort(readdir(dir), by = x -> parse(Int64, split(x, "_")[1]))

    dir_res = "../res/totalReults/"
    KEYS = ["CPExact", "BCExact", "Dual", "BCHeur", "Heur"]
    times = Dict(k => [0.0] for k in KEYS)
    solved = Dict(k => [0] for k in KEYS)
    ite = 0


    fout = open(dir_res * "totalReults.tex", "w")

    latex = raw"""\begin{sidewaystable}[!h]
    \centering
    \scalebox{0.9}{
    \renewcommand{\arraystretch}{1.2}
    \begin{tabular}{lcccccccccccccccc}
    \toprule
    \textbf{Cities} & \textbf{PR}  & \multicolumn{3}{c}{\textbf{Cutting Planes}} & \multicolumn{3}{c}{\textbf{Branch\&Cut}}  & \multicolumn{3}{c}{\textbf{Dualization}} & \multicolumn{3}{c}{\textbf{B\&C Heur}} & \multicolumn{3}{c}{\textbf{Heuristic}}
    \\
    \cmidrule(r){3-5} \cmidrule(r){6-8} \cmidrule(r){9-11} \cmidrule(r){12-14} \cmidrule(r){15-17}
     & & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} & \textbf{Time(s)} & \textbf{Gap} & \textbf{Best bound} \\
    \midrule"""

    println(fout, latex)


    for fileName in files

        println("---------------------")
        println(fileName)
        println("---------------------")

        countFails = 0

        # reading data
        preparation(dir, fileName)

        print(fout, n)
        robust_val = 0.0

        # static
        println("\n static \n")
        solStatic = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        try
            solStatic = cplexSolveStaticSP()
        catch
            println("Time out ! ")
            break
        end

        # cutting planes
        println("\n cuting planes exact \n")
        solCPExact = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        last = size(solved["CPExact"], 1)
        extra = 0.0
        try
            solCPExact = cuttingPlanes(false)

            if solCPExact.best_bound > 0.0
                append!(solved["CPExact"], solved["CPExact"][last]+1)
            else
                append!(solved["CPExact"], solved["CPExact"][last])
            end
            if solCPExact.isOptimal
                robust_val = solCPExact.obj_val
            end
        catch
            println("Time out ! ")
            countFails += 1
            append!(solved["CPExact"], solved["CPExact"][last])
            extra = TimeLimit
        end
        last = size(times["CPExact"], 1)
        append!(times["CPExact"], round(solCPExact.solveTime, digits=2) + times["CPExact"][last]+extra)


        # branch cut
        println("\n Branch and Cut exact \n")
        solBCExact = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        last = size(solved["BCExact"], 1)
        extra = 0
        try
            solBCExact = brunchAndCut(true, false)

            if solBCExact.best_bound > 0.0
                append!(solved["BCExact"], solved["BCExact"][last]+1)
            else
                append!(solved["BCExact"], solved["BCExact"][last])
            end
            if solBCExact.isOptimal
                robust_val = solBCExact.obj_val
            end
        catch
            println("Time out ! ")
            countFails +=1
            append!(solved["BCExact"], solved["BCExact"][last])
            extra = TimeLimit
        end
        last = size(times["BCExact"], 1)
        append!(times["BCExact"], round(solBCExact.solveTime, digits=2) + times["BCExact"][last]+extra)


        # dualisation
        println("\n dualisation \n")
        solRobustDual = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        last = size(solved["Dual"], 1)
        extra = 0
        try
            solRobustDual = dualSolve()

            if solRobustDual.best_bound > 0.0
                append!(solved["Dual"], solved["Dual"][last]+1)
            else
                append!(solved["Dual"], solved["Dual"][last])
            end
            if solRobustDual.isOptimal
                robust_val = solRobustDual.obj_val
            end
        catch
            println("Time out ! ")
            countFails +=1
            append!(solved["Dual"], solved["Dual"][last])
            extra = TimeLimit
        end
        last = size(times["Dual"], 1)
        append!(times["Dual"], round(solRobustDual.solveTime, digits=2)+times["Dual"][last]+extra)


        # branch and cut heuristique
        println("\n Branch and Cut heurisric \n")
        solBCHeur = Solution(false, false, 0.0, 0.0, 0.0, 0.0)
        gap_heur_BC = 0.0
        last = size(solved["BCHeur"], 1)
        extra = 0
        try
            solBCHeur = brunchAndCut(false, true)

            if solBCHeur.best_bound > 0.0
                if robust_val > 0.0
                    gap_heur_BC = round(abs(robust_val - solBCHeur.obj_val)/robust_val, digits=2)
                    println("The gap between solBCExact and solBCHeur is ", gap_heur_BC)
                end

                append!(solved["BCHeur"], solved["BCHeur"][last]+1)
            else
                append!(solved["BCHeur"], solved["BCHeur"][last])
            end
        catch
            println("Time out ! ")
            countFails +=1
            append!(solved["BCHeur"], solved["BCHeur"][last])
            extra = TimeLimit
        end
        last = size(times["BCHeur"], 1)
        append!(times["BCHeur"], round(solBCHeur.solveTime, digits=2) + times["BCHeur"][last]+extra)


        # heurisric primal
        println("\n heurisric primal \n")
        solHeurPrimal = SolHeur(0.0, 0.0)
        gap_heur_primal = 0.0
        last = size(solved["Heur"], 1)
        extra = 0
        try
            solHeurPrimal = heuristicPrimal()

            if solHeurPrimal.obj_val > 0.0
                if robust_val > 0.0
                    gap_heur_primal = abs(robust_val - solHeurPrimal.obj_val)/robust_val
                    println("The gap between exact and heurisric is ", gap_heur_primal)
                end

                append!(solved["Heur"], solved["Heur"][last]+1)
            else
                append!(solved["Heur"], solved["Heur"][last])
            end
        catch
            println("Time out ! ")
            countFails +=1
            append!(solved["Heur"], solved["Heur"][last])
            extra = TimeLimit
        end
        last = size(times["Heur"], 1)
        append!(times["Heur"], round(solHeurPrimal.solveTime, digits=2)+times["Heur"][last]+extra)


        # col PR
        if robust_val>0.0
            pr = round(abs(solStatic.best_bound - robust_val)/robust_val, digits=2)
            print(fout, " & ", pr)
        else
            print(fout, " & -")
        end

        # col C P
        print(fout, " & ", round(solCPExact.solveTime, digits=2), " & ", round(solCPExact.GAP, digits=2))
        if solCPExact.best_bound > 0.0
            print(fout, " & ", round(solCPExact.best_bound, digits=2))
        else
            print(fout, " & -")
        end

        # col BC
        print(fout, " & ", round(solBCExact.solveTime, digits=2), " & ", round(solBCExact.GAP, digits=2))
        if solBCExact.best_bound > 0.0
            print(fout, " & ", round(solBCExact.best_bound, digits=2))
        else
            print(fout, " & -")
        end

        # col dual
        print(fout, " & ", round(solRobustDual.solveTime, digits=2), " & ", round(solRobustDual.GAP, digits=2))
        if solRobustDual.best_bound > 0.0
            print(fout, " & ", round(solRobustDual.best_bound, digits=2))
        else
            print(fout, " & -")
        end

        # col BC Heur
        print(fout, " & ", round(solBCHeur.solveTime, digits=2))
        print(fout, " & ", gap_heur_BC)
        if solBCHeur.best_bound > 0.0
            print(fout, " & ", round(solBCHeur.best_bound, digits=2))
        else
            print(fout, " & -")
        end


        # col Heur
        print(fout, " & ", round(solHeurPrimal.solveTime, digits=2))
        print(fout, " & ", round(gap_heur_primal, digits=2))

        if solHeurPrimal.obj_val > 0.0
            print(fout, " & ", round(solHeurPrimal.obj_val, digits=2))
        else
            print(fout, " & -")
        end

        println(fout, " \\\\")

        if countFails >= 4
            break
        end

    end

    latex = raw"""\bottomrule
    \end{tabular}
    }
    \caption{.}
    \label{tab:total}
    \end{sidewaystable}"""
    println(fout, latex)

    close(fout)


    for method in KEYS
        plot(times[method], solved[method], label = method)
    end

    legend()
    # yscale("log")
    # xscale("log")
    title("Comparison different methods' performance")
    xlabel("Time(s)")
    ylabel("The number of instances solved")
    savefig(dir_res * "totalReults.png")
    close()

end

function benchmarks()
    test()

    generateSubHeurInfl()
    generateSenariosInfl()
    totalReults()
end