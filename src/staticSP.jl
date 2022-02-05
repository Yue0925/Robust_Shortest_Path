# This file contains a method to the static shortest path problem with CPLEX

# TODO : it seems that the oriented graph is symmetric (i.e. ij and ji are identical)


TOL = 0.00001


"""
Solve the static shortest path problem by CPLEX.

ATTENTION: data should be charged before calling this function !

Return : Solution
"""
function cplexSolveStaticSP()

    # modelization
    M = Model(CPLEX.Optimizer) 
    # set_optimizer_attribute(M, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    set_optimizer_attribute(M, "CPXPARAM_TimeLimit", TimeLimit) # seconds

    # variables
    @variable(M, y[1:n], Bin)
    @variable(M, x[1:n, 1:n], Bin)

    # ---------------
    # prefix values
    # ---------------
    @constraint(M, x[n, n] == 0) # no bucles

    for i in 1:n-1
        @constraint(M, x[i, i] == 0)

        for j in i+1:n
            # if arc x[ij]=1, then arc x[ji]=0
            @constraint(M, x[i, j] + x[j, i] <= 1) 

            # if i, j are not adjacent, then x[ij] = 0
            if !Adjacenct[i, j]
                @constraint(M, x[i, j] == 0)
                @constraint(M, x[j, i] == 0)
            end
        end
    end


    # weight capacity constraint
    @constraint(M, sum(y[i] * p[i] for i in 1:n) <= S)


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
    @objective(M, Min, sum(x[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 3] for l in 1:arcs))

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

        isFeasible = verifyStaticSP(path, vertices)
        println("isFeasible ? ", isFeasible)
    end

    return Solution(isOptimal, isFeasible, obj_val, solveTime, GAP, best_bound)
end


"""
Verify whether the solution solved by CPLEX is feasible 
(i.e. - a route from s to t 
      - the total weight doesn't not exceed to the limit S ) 

Return : Bool
"""
function verifyStaticSP(path::Array{Tuple{Int64, Int64}, 1}, vertices::Array{Int64, 1})
    predecessor = [0 for _ in 1:n]

    for arc in path
        predecessor[arc[2]] = arc[1]
    end

    for v in vertices
        if v == s
            continue
        end
        if predecessor[v] == 0
            return false
        end
    end

    return sum(p[v] for v in vertices) <= S
end



"""
Verify the solution is robust feasible.
"""
function verifyRobustSP(path::Array{Tuple{Int64, Int64}, 1}, vertices::Array{Int64, 1})
    # at first it is a feasible static path
    isFeasible = verifyStaticSP(path, vertices)
    if isFeasible == false
        return false
    end

    # secondly, the total robust weight should not exceed to S
    
    # modelization
    M = Model(CPLEX.Optimizer) 

    # variables
    @variable(M, 0 <= δ2[1:n] <= 2)

    for i in 1:n
        if !(i in vertices)
            @constraint(M, δ2[i] == 0)
        end
    end

    # total augmentation limit
    @constraint(M, sum(δ2[i] for i in vertices) <= d2 )

    # objective function
    @objective(M, Max, sum((p[v] + δ2[v] * ph[v]) for v in vertices ))
    
    set_silent(M) # turn off cplex output
    optimize!(M)

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    if isOptimal
        obj_val = objective_value(M)
        return obj_val<=S
    end
    println("CPLEX sol not optimal")
    return false
end