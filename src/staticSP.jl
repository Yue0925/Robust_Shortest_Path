# This file contains a method to the static shortest path problem with CPLEX

# TODO : it seems that the oriented graph is symmetric (i.e. ij and ji are identical)


TOL = 0.00001


"""
Solve the static shortest path problem by CPLEX.

ATTENTION: data should be charged before calling this function !
"""
function cplexSolveStaticSP()

    # matrix adjacecy 
    Adjacenct = falses(n, n) 
    arcs = size(Mat, 1)
    for i in 1:arcs
        Adjacenct[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] = true
    end

    # modelization
    M = Model(CPLEX.Optimizer)

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
    
    # solve the problem
    optimize!(M)

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    # display solution
    println("isOptimal ? ", isOptimal)
    println("the path from ", s, " to ", t, " is :")
    for i in 1:n
        for j in 1:n 
            if JuMP.value(x[i, j]) > TOL
                println("(", i, ", ", j, ")")
            end
        end
    end
end