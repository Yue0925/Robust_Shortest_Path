# This file contains heurisric solutions for the robust shortest path problem 

"""
A greedy primal heurisric solution.

Idea : suppose that the distance of every arc is augmented to the maximum
i.e. augmentation of D_ij % ∀ij, and select the shortest path s to t (Dijkstra) 
without exceeding to S(so the weight of each vertex is extended to the maximum)
"""
function heuristicPrimal()

    # start a chronometer
    start = time()

    Q = [v for v in 1:n]
    prec = [0 for _ in 1:n]
    dist = Dict(v => Inf for v in Q)
    dist[s] = 0.0
    augmentedDist = zeros(n, n)
    Δ = zeros(n, n)
    for i in 1:arcs
        weight = p[round(Int, Mat[i, 2])] + 4 * ph[round(Int, Mat[i, 2])]
        if round(Int, Mat[i, 2]) == t || round(Int, Mat[i, 2]) == s
            weight = 0
        end
        augmentedDist[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] = Mat[i, 3] * (1 + Mat[i, 4]) + weight
        Δ[round(Int, Mat[i, 1]), round(Int, Mat[i, 2])] = Mat[i, 4]
    end

    while size(Q, 1) > 0
        u = reduce((x, y) -> dist[x] ≤ dist[y] ? x : y, Q)
        filter!(e -> e ≠ u, Q)

        if u == t
            break
        end

        neighbours = filter(v -> Adjacenct[u, v] , Q)
        
        for v in neighbours
            alt = dist[u] + augmentedDist[u, v]
            if alt < dist[v]
                dist[v] = alt
                prec[v] = u
            end
        end
    end

    solveTime = time() - start

    path = Array{Tuple{Int64, Int64}, 1}()
    vertices = Array{Int64, 1}()
    u = t
    while prec[u] != 0
        append!(vertices, u)
        append!(path, [(prec[u], u)])
        u = prec[u]
    end
    append!(vertices, s)

    δ1 = Dict(a => 0.0 for a in path)
    acc = 0.0
    robustDist = 0.0
    sort!(path, by = x -> augmentedDist[x[1], x[2]], rev = true)
    for a in path
        if acc + Δ[a[1], a[2]] < d1
            δ1[a] = Δ[a[1], a[2]]
            acc += Δ[a[1], a[2]]
        else
            δ1[a] = d1 - acc
            acc += δ1[a]
            break
        end
    end

    for i in 1:arcs
        a = (round(Int, Mat[i, 1]), round(Int, Mat[i, 2]))
        if a in path
            robustDist += Mat[i, 3] * (1 + δ1[a])
        end
    end

    println("robustDist : ", robustDist)
    isFeasible = verifyRobustSP(path, vertices)
    println("isFeasible ? ", isFeasible)
end