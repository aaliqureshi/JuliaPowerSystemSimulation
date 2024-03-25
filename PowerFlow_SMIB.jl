module PowerFlow
# using DifferentialEquations
using OrdinaryDiffEq

export solve_powerflow

function powerflow!(du,u,p,t)
    θ_1, v_3, θ_3 = u

    Y = p[1]
    v_1, v_2, θ_2 = p[2]

    p_1 = 0.944
    p_l3 = 0
    q_l3 = 0

    Y_abs = abs.(Y)
    Y_ang = angle.(Y)

    y_11 = Y_abs[1,1]
    y_12 = Y_abs[1,2]
    y_13 = Y_abs[1,3]
    y_21 = Y_abs[2,1]
    y_22 = Y_abs[2,2]
    y_23 = Y_abs[2,3]
    y_31 = Y_abs[3,1]
    y_32 = Y_abs[3,2]
    y_33 = Y_abs[3,3]

    ϕ_11 = Y_ang[1,1]
    ϕ_12 = Y_ang[1,2]
    ϕ_13 = Y_ang[1,3]
    ϕ_21 = Y_ang[2,1]
    ϕ_22 = Y_ang[2,2]
    ϕ_23 = Y_ang[2,3]
    ϕ_31 = Y_ang[3,1]
    ϕ_32 = Y_ang[3,2]
    ϕ_33 = Y_ang[3,3]

    # v_1 = 1.03
    # θ_2 = 0.0
    # v_2 = 1.0
    
    # v_1 = v_10
    # θ_2 = θ_2 
    # v_2 = v_2
    
    du[1] = p_1  - (v_1*v_1*y_11*cos(θ_1 - θ_1 - ϕ_11) + v_1*v_2*y_12*cos(θ_1 - θ_2 - ϕ_12) + v_1*v_3*y_13*cos(θ_1 - θ_3 - ϕ_13))
    du[2] = p_l3 - (v_3*v_1*y_31*cos(θ_3 - θ_1 - ϕ_31) + v_3*v_2*y_32*cos(θ_3 - θ_2 - ϕ_32) + v_3*v_3*y_33*cos(θ_3 - θ_3 - ϕ_33)) 
    du[3] = q_l3 - (v_3*v_1*y_31*sin(θ_3 - θ_1 - ϕ_31) + v_3*v_2*y_32*sin(θ_3 - θ_2 - ϕ_32) + v_3*v_3*y_33*sin(θ_3 - θ_3 - ϕ_33))
end


function solve_powerflow(Y::Matrix{ComplexF64}, u_init)

    M0 = zeros(3,3)
    tspan = (0.0,5)
    # θ_1,v_3,θ_3=u
    u0  =[0,1.0,0]

    p = (Y,u_init)

    pf_fun = ODEFunction(powerflow!,mass_matrix = M0)
    
    pf_prob = ODEProblem(pf_fun,u0,tspan,p)
    
    sol = solve(pf_prob,Rodas5(),reltol = 1e-8);

    println(sol.retcode)

    θ_10,v_30,θ_30 = sol.u[end]

    res = [θ_10,v_30,θ_30]

    append!(res,u_init)

    return res

end

end #module