# This file contains a method to the robust shortest path problem with CPLEX using cutting planes algorithm

TOL = 0.00001


"""
Branch and Cut approach, cutting planes algorithm using callback.
"""
function brunchAndCut()
    U1_star, U2_star = initSenario()
    masterPB(U1_star, U2_star, true)
    optimize!(MP)

    # status of model
    statusMP = termination_status(MP)
    isOptimalMP = statusMP==MOI.OPTIMAL
    println("masterPB isOptimal? ", isOptimalMP)


    x = MP[:x]
    y = MP[:y]

    path = Array{Tuple{Int64, Int64}, 1}()
    vertices = Array{Int64, 1}()
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

    println("objective value : ", objective_value(MP))
    println("total weight : ", sum(p[v] for v in vertices))
    return path, vertices
end

"""
Define / modelize the master problem, by defaut we dont use callback 
"""
function masterPB(U1_star::Array{Float64,2}, U2_star::Array{Int64,1}, BranchCut=false)
    # modelization
    global MP = Model(CPLEX.Optimizer) 
    # set_optimizer_attribute(MP, "CPX_PARAM_MIPDISPLAY", 2) # MIP display
    #set_optimizer_attribute(MP, "CPXPARAM_TimeLimit", 500) # seconds
    # It is imposed to use only 1 thread in Julia with CPLEX to use the callbacks
    if BranchCut
        MOI.set(MP, MOI.NumberOfThreads(), 1) #TODO : but it seems that it's not necessary !
    end

    # variables
    @variable(MP, y[1:n], Bin)
    @variable(MP, x[1:n, 1:n], Bin)
    @variable(MP, z)

    # ---------------
    # prefix values
    # ---------------
    @constraint(MP, x[n, n] == 0) # no bucles

    for i in 1:n-1
        @constraint(MP, x[i, i] == 0)

        for j in i+1:n
            # if arc x[ij]=1, then arc x[ji]=0
            @constraint(MP, x[i, j] + x[j, i] <= 1) 

            # if i, j are not adjacent, then x[ij] = 0
            if !Adjacenct[i, j]
                @constraint(MP, x[i, j] == 0)
                @constraint(MP, x[j, i] == 0)
            end
        end
    end


    # robust distance constraint
    #@constraint(MP, sum(x[i, j] * U1_star[i, j] for i in 1:n for j in 1:n) <= z)
    @constraint(MP, sum(x[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 3] for l in 1:arcs) <= z)


    # robust weight constraint
    #@constraint(MP, sum(y[i] * U2_star[i] for i in 1:n) <= S)
    @constraint(MP, sum(y[i] * p[i] for i in 1:n) <= S)


    # -----------------------------------
    # constraints of a path from s to t
    # -----------------------------------

    # only one arc outgoing s 
    @constraint(MP, sum(x[s, j] for j in 1:n) == 1)
    @constraint(MP, y[s] == 1)

    # only one arc incoming t
    @constraint(MP, sum(x[i, t] for i in 1:n) == 1)
    @constraint(MP, y[t] == 1)

    # flot conversation
    for v in 1:n
        if v == s || v == t 
            continue
        end
        @constraint(MP, sum(x[v, j] for j in 1:n) - sum(x[i, v] for i in 1:n) == 0)
    end

    # for each intermediate vertex v, at most one arc outgoing / incoming
    for v in 1:n
        if v == s || v == t 
            continue
        end
        @constraint(MP, sum(x[v, j] for j in 1:n) == y[v])
        @constraint(MP, sum(x[i, v] for i in 1:n) == y[v])
    end

    # objective function
    @objective(MP, Min, z)


    function callback_cuttingPlanes(cb_data::CPLEX.CallbackContext) # context_id::Clong
        println("callback")
        # get current variables
        x_star = zeros(n, n)
        y_star = zeros((n))
        for i in 1:n
            y_star[i] = callback_value(cb_data, y[i])
            for j in 1:n
                x_star[i, j] = callback_value(cb_data, x[i, j])
            end
        end
        z_star = callback_value(cb_data, z)

        # solve the sub problems
        SM1 = subPB1(x_star)
        optimize!(SM1)
        z1_sub = objective_value(SM1)

        SM2 = subPB2(y_star)
        optimize!(SM2)
        z2_sub = objective_value(SM2)

        # if SP1 violates
        if abs(z_star - z1_sub) > TOL
            println("SP1 violated")
            δ1_star = value.(SM1[:δ1])
            constraint1 = @build_constraint(sum(x[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])]
            * Mat[l, 3] *(1 + δ1_star[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])]) for l in 1:arcs) <= z)

            MOI.submit(MP, MOI.LazyConstraint(cb_data), constraint1)
        end

        # if the SP2 violates
        if z2_sub - S > TOL
            println("SP2 violated")
            δ2_star = value.(SM2[:δ2])
            constraint2 = @build_constraint(sum(y[i] * (p[i] + δ2_star[i] * ph[i]) for i in 1:n ) <= S)

            MOI.submit(MP, MOI.LazyConstraint(cb_data), constraint2)
        end
    end

    if BranchCut
        MOI.set(MP, MOI.LazyConstraintCallback(), callback_cuttingPlanes)
    end
    set_silent(MP) # turn off cplex output
end


"""
Define / modelize the sub problem SP1
"""
function subPB1(x_star::Array{Float64,2})
    # modelization
    SM1 = Model(CPLEX.Optimizer) 

    # variables
    @variable(SM1, δ1[1:n, 1:n] >= 0)

    # ---------------
    # prefix values
    # ---------------
    for i in 1:n
        for j in 1:n
            # if i, j are not adjacent, then δ1[ij] = 0
            if !Adjacenct[i, j]
                @constraint(SM1, δ1[i, j] == 0)
            end
        end
    end
    # for each arc ij, set the upper bound value
    for l in 1:arcs
        @constraint(SM1, δ1[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] <= Mat[l, 4])
    end

    # total augmentation limit
    @constraint(SM1, sum(δ1[i, j] for i in 1:n for j in 1:n) <= d1)

    # objective function
    @objective(SM1, Max, sum(x_star[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 3] *(1 + δ1[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])]) for l in 1:arcs))

    set_silent(SM1) # turn off cplex output
    return SM1
end


"""
Define / modelize the sub problem SP2
"""
function subPB2(y_star::Array{Float64,1})
    # modelization
    SM2 = Model(CPLEX.Optimizer) 

    # variables
    @variable(SM2, 0 <= δ2[1:n] <= 2)

    # total augmentation limit
    @constraint(SM2, sum(δ2[i] for i in 1:n) <= d2 )

    # objective function
    @objective(SM2, Max, sum(y_star[i] * (p[i] + δ2[i] * ph[i]) for i in 1:n ))
    
    set_silent(SM2) # turn off cplex output
    return SM2
end



"""
Initialize the incertitudes senarios, by defaut d_ij^1 = d_ij 
"""
function initSenario()
    # TODO: try with other possibilities
    U1_star = zeros(n, n)
    for l in 1:arcs
        U1_star[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] = Mat[l, 3]
    end

    U2_star = [i for i in p]
    return U1_star, U2_star
end



"""
The naive cutting planes algorithm
"""
function cuttingPlanes()
    # ---------------------------
    # intialize the incertitudes
    # ---------------------------
    U1_star, U2_star = initSenario()

    ite = 1
    println("--------------")
    println("step ", ite)
    println("--------------")

    # ------------------------------------
    # step 1 : resolve the master problem
    # ------------------------------------
    masterPB(U1_star, U2_star)

    # solve the master problem
    optimize!(MP)

    # status of model
    statusMP = termination_status(MP)
    isOptimalMP = statusMP==MOI.OPTIMAL
    println("masterPB isOptimal? ", isOptimalMP)

    # get variables
    x_star = value.(MP[:x])
    z_star = value(MP[:z])
    y_star = value.(MP[:y])
    x = MP[:x]
    z = MP[:z]
    y = MP[:y]
    println("master z = ", z_star)


    # ------------------------------------
    # step 2 : resolve the sub problems
    # ------------------------------------
    SM1 = subPB1(x_star)

    # solve the sub problem related to U1
    optimize!(SM1)

    # status of model
    statusSM1 = termination_status(SM1)
    isOptimalSM1 = statusSM1==MOI.OPTIMAL
    println("subPB1 isOptimal? ", isOptimalSM1)

    # the objective value z1
    z1_sub = objective_value(SM1)
    println("z1_sub = ", z1_sub)

    SM2 = subPB2(y_star)

    # solve the sub problem related to U2
    optimize!(SM2)

    # status of model
    statusSM2 = termination_status(SM2)
    isOptimalSM2 = statusSM2==MOI.OPTIMAL
    println("subPB2 isOptimal? ", isOptimalSM2)

    # the objective value z2
    z2_sub = objective_value(SM2)
    println("z2_sub = ", z2_sub)


    # ----------------------------------------------------
    # iteratively add senario to master problem
    # until reach the optimal condition
    # ----------------------------------------------------
    while abs(z_star - z1_sub) > TOL || z2_sub - S > TOL
        if ite >300
            break
        end

        ite += 1
        println("--------------")
        println("step ", ite)
        println("--------------")

        # if the SP1 violates
        if abs(z_star - z1_sub) > TOL
            println("SP1 violated")
            δ1_star = value.(SM1[:δ1])
            @constraint(MP, sum(x[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])] * Mat[l, 3] *(1 + δ1_star[round(Int, Mat[l, 1]), round(Int, Mat[l, 2])]) for l in 1:arcs) <= z)
        end

        # if the SP2 violates
        if z2_sub - S > TOL
            println("SP2 violated")
            δ2_star = value.(SM2[:δ2])
            @constraint(MP, sum(y[i] * (p[i] + δ2_star[i] * ph[i]) for i in 1:n ) <= S)
        end

        # ------------------------------------
        # step 1 : resolve the master problem
        # ------------------------------------
        optimize!(MP)

        # status of model
        statusMP = termination_status(MP)
        isOptimalMP = statusMP==MOI.OPTIMAL
        println("masterPB isOptimal? ", isOptimalMP)

        # get variables
        x_star = value.(MP[:x])
        z_star = value(MP[:z])
        y_star = value.(MP[:y])
        x = MP[:x]
        z = MP[:z]
        y = MP[:y]
        println("master z = ", z_star)


        # ------------------------------------
        # step 2 : resolve the sub problems
        # ------------------------------------
        SM1 = subPB1(x_star)
        optimize!(SM1)

        # status of model
        statusSM1 = termination_status(SM1)
        isOptimalSM1 = statusSM1==MOI.OPTIMAL
        println("subPB1 isOptimal? ", isOptimalSM1)

        # the objective value z1
        z1_sub = objective_value(SM1)
        println("z1_sub = ", z1_sub)

        SM2 = subPB2(y_star)
        optimize!(SM2)

        # status of model
        statusSM2 = termination_status(SM2)
        isOptimalSM2 = statusSM2==MOI.OPTIMAL
        println("subPB2 isOptimal? ", isOptimalSM2)

        # the objective value z2
        z2_sub = objective_value(SM2)
        println("z2_sub = ", z2_sub)

    end

    path = Array{Tuple{Int64, Int64}, 1}()
    vertices = Array{Int64, 1}()
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

    println("objective value : ", objective_value(MP))
    println("total weight : ", sum(p[v] for v in vertices))
    return path, vertices
end