# Robust_Shortest_Path
 A project on the Robust Shortest Path Problem for the course ECMA of MPRO 2021-2022, finished by Yue only.



# Reminders

* To compile file, using include(file) command.


# Suggestions

* Rather than charging libraries in every file, uisng libraries only once in terminal before executing files. 
* Move to the src directory : 
  
```julia
using CPLEX 
using JuMP
using StatsBase
using Random
using PyPlot

include("io.jl")
include("run.jl")
test() 
benchmarks()

```

* Package required : 
  
```julia
import Pkg; Pkg.add("StatsBase")
import Pkg; Pkg.add("PyPlot")
```