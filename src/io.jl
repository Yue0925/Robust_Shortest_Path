# This file contains functions related to reading, writing and displaying experimental results.


"""
Solution format
"""
struct Solution
    isOptimal::Bool
    isFeasible::Bool
    obj_val::Float64
    solveTime::Float64
    GAP::Float64
    best_bound::Float64
end

# time limit in seconds 60 !!!
global TimeLimit = 60
