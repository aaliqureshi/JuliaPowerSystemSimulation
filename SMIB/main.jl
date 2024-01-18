using DifferentialEquations, Plots

dir = pwd()
include(dir*"/MyModules/Y_matrix_code.jl")
include(dir*"/SMIB/PowerFlow_SMIB.jl")
include(dir*"/SMIB/TDS_SMIB.jl")
include(dir*"/SMIB/Plotter.jl")

using .YMatrix
using .PowerFlow
using .TDS
using .Plotter

begin
    n_bus = 3
    n_lines = 2

    v_10 = 1.03
    v_2 = 1.0
    θ_2 = 0.0


    R1 = 0.015
    # R1 = 0.0
    XL1 = 0.15

    R2 = 0.0198
    # R2 = 0.0
    XL2 = 0.198

    
    Z_dict = Dict([("R1", R1),
                   ("R2", R2),
                   ("XL1", XL1),
                   ("XL2", XL2)])

    fault_dict = Dict([("Rs", 0),
                       ("XLs", 1e-2),
                       ("tf", 6.0),
                       ("tc",6.1),
                       ("tspan", 15)])

    
    param_dict = Dict([("Gen Model", "AF"),
                      ("Framework","DP Full")])
    
    # plot_dict = Dict([("Interpolate", 1),
    #                   ("Overlay",1)])

    u0 = [v_10,v_2,θ_2]

    bus_impedance = Dict((1,3) => R1+XL1*1im,
                         (2,3) => R2+XL2*1im)
    # # #                      (3,3) => 0.01+0im)

    Y = build_Y(bus_impedance,n_bus)

    params = (Y, Z_dict, fault_dict, param_dict)

end

# θ_10, v_30, θ_30, v_10,v_2,θ_2
u_pf = solve_powerflow(Y, u0)

# sol = solve_TDS(Y,Z_dict,fault_dict,u_pf)
sol = solve_TDS(params,u_pf)

# plot(sol,idxs=2)
plot_data("V1",sol,params,Interpolate=1,overlay=1) 

save_data(sol,params)