# This file contains the dualization of the robust shortest path problem solved with CPLEX


function dualSolve()

    # modelization
    M = Model(CPLEX.Optimizer) 
    # set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPXPARAM_TimeLimit", TimeLimit) # seconds

    # variables
    @variable(M, y[1:n], Bin)
    @variable(M, x[1:n, 1:n], Bin)
    @variable(M, α >= 0)
    @variable(M, β >= 0)
    @variable(M, ω[1:n] >= 0)
    @variable(M, λ[1:n, 1:n] >= 0)


    # -----------------
    # dual constraints
    # -----------------
    for i in 1:n
        @constraint(M, β + ω[i] >= ph[i] * y[i])
    end

    @constraint(M, sum((p[i] * y[i] + 2 * ω[i]) for i in 1:n) + d2 * β <= S)

    for i in 1:arcs
        @constraint(M, α + λ[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] >= x[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] * round(Int, Mat[i, 3]))
    end


    # ---------------
    # prefix values
    # ---------------
    @constraint(M, x[n, n] == 0) # no bucles
    @constraint(M, λ[n, n] == 0) # no bucles

    for i in 1:n-1
        @constraint(M, x[i, i] == 0)
        @constraint(M, λ[i, i] == 0)


        for j in i+1:n
            # if arc x[ij]=1, then arc x[ji]=0
            @constraint(M, x[i, j] + x[j, i] <= 1) 

            # if i, j are not adjacent, then x[ij] = 0
            if !Adjacenct[i, j]
                @constraint(M, x[i, j] == 0)
                @constraint(M, x[j, i] == 0)
                @constraint(M, λ[j, i] == 0)
                @constraint(M, λ[i, j] == 0)

            end
        end
    end

    # -----------------------------------
    # constraints of a path from s to t
    # -----------------------------------

    # only one arc outgoing s 
    @constraint(M, sum(x[s, j] for j in 1:n) == 1)
    @constraint(M, y[s] == 1)

    # only one arc incoming t
    @constraint(M, sum(x[i, t] for i in 1:n) == 1)
    @constraint(M, y[t] == 1)

    # flot conversation
    for v in 1:n
        if v == s || v == t 
            continue
        end
        @constraint(M, sum(x[v, j] for j in 1:n) - sum(x[i, v] for i in 1:n) == 0)
    end

    # for each intermediate vertex v, at most one arc outgoing / incoming
    for v in 1:n
        if v == s || v == t 
            continue
        end
        @constraint(M, sum(x[v, j] for j in 1:n) == y[v])
        @constraint(M, sum(x[i, v] for i in 1:n) == y[v])
    end


    # objective function
    @objective(M, Min, sum( (x[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 3] +
    λ[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 4]) for l in 1:arcs) + d1 * α)

    set_silent(M) # turn off cplex output

    # solve the problem
    optimize!(M)

    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = MOI.get(M, MOI.SolveTime())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    path = Array{Tuple{Int64, Int64}, 1}()
    vertices = Array{Int64, 1}()
    obj_val = 0.0
    isFeasible = false
    best_bound = 0.0

    # display solution
    println("isOptimal ? ", isOptimal)
    println("GAP = ", GAP)

    if has_values(M)
        println("the path from ", s, " to ", t, " is :")
        for i in 1:n
            if JuMP.value(y[i]) > TOL
                append!(vertices, i)
            end
            for j in 1:n 
                if JuMP.value(x[i, j]) > TOL
                    println("(", i, ", ", j, ")")
                    append!(path, [(i, j)])
                end
            end
        end

        obj_val = objective_value(M)
        best_bound = objective_bound(M)

        println("objective value : ", obj_val)
        println("best bound : ", best_bound)
        println("solveTime : ", solveTime)
        println("nodes : ", exploredNodes)

        isFeasible = verifyRobustSP(path, vertices)
        println("isFeasible ? ", isFeasible)
    end

    return Solution(isOptimal, isFeasible, obj_val, solveTime, GAP, best_bound)
end
