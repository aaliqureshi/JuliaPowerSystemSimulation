module TDS
# using DifferentialEquations


include("TdsGenclsSmib.jl")
# include("TdsGenclsSmib2.jl")
include("TdsGenAFSmib.jl")

using .TdsGenclsSmib
# using .TdsGenclsSmib2
using .TdsGenAFSmib

export solve_TDS



function solve_TDS(params,u_pf)

    Y, Z_dict, fault_dict, param_dict = params

    if param_dict["Gen Model"] == "CLS"
        sol = TDS_GENCLS(params,u_pf)
        return sol
    elseif param_dict["Gen Model"] == "AF"
        sol = TDS_GENAF(params,u_pf)
        return sol
    end

end #function

end #module