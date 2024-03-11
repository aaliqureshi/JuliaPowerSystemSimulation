module TdsGenclsSmib

using OrdinaryDiffEq

export TDS_GENCLS

function fxinit!(du,u,p,t)
    δ,ω,e_q_prime,i_d,i_q,i_32d,i_32q = u

    Y = p[1]
    Z = p[2]
    θ_1, v_3, θ_3, v_1, v_2, θ_2 = p[3]

    R1 = Z["R1"]
    R2 = Z["R2"]
    w_L1 = Z["XL1"]
    w_L2 = Z["XL2"]
    
    x_d_prime = 0.24
    p_m = 0.944
    H = 3.5

    Ω = 2π*60
    d = 0

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

    du[1] = Ω*(ω-1)
    du[2] = (p_m - (i_d*v_1*cos(θ_1) + i_q*v_1*sin(θ_1)) - d*(ω-1))/(2H)
    du[3] = v_1*cos(θ_1) - x_d_prime*i_q + e_q_prime*sin(δ-π/2)
    du[4] = v_1*sin(θ_1) - e_q_prime*cos(δ-π/2) + i_d*x_d_prime
    du[5] = i_d - (v_1*y_11*cos(θ_1 + ϕ_11) + v_2*y_12*cos(θ_2 + ϕ_12) + v_3*y_13*cos(θ_3 + ϕ_13))
    du[6] = i_32d - (v_3*sin(θ_3) - v_2*sin(θ_2) - R2*i_32q)/(w_L2)
    du[7] = i_32q - (-v_3*cos(θ_3) + v_2*cos(θ_2) + R2*i_32d)/(w_L2)

end

function fx!(du,u,p,t)
    # δ,ω,i_d,i_q,v_1d,v_1q,v_3d,v_3q,i_13d,i_13q,i_30d,i_30q,i_32d,i_32q=u
    # δ,ω,i_13d,i_13q,i_32d,i_32q,i_30d,i_30q,i_d,i_q,v_1d,v_1q,v_3d,v_3q=u
    δ,ω,i_13d,i_13q,i_32d,i_32q,i_30d,i_30q,v_1d,v_1q,v_3d,v_3q=u

    Z, fault_dict, k, e_q_prime0 = p
    
    e_q_prime = e_q_prime0
    x_d_prime = 0.24
    p_m = 0.944
    H = 3.5
    Ω = 2π*60
    
    d = 1

    R1 = Z["R1"]
    R2 = Z["R2"]
    w_L1 = Z["XL1"]
    w_L2 = Z["XL2"]

    w_L3 = fault_dict["XLs"]
    R3 = fault_dict["Rs"]

    L1 = w_L1/Ω
    L2 = w_L2/Ω
    L3 = w_L3/Ω

    Rs = k*R3 + abs(k-1)*1e100
    Ls = k*L3 + abs(k-1)*1e100

    v_2 = 1.0
    θ_2 = 0.0

    v_2 = v_2*exp(im*θ_2)

    v_2d = real(v_2)
    v_2q = imag(v_2)


    du[1] = Ω*(ω-1)
    du[2] = (p_m - (v_1d*i_13d + v_1q*i_13q) - d*(ω-1))/(2H)

    # # ## L1 equation
    du[3] = v_1d/L1 - v_3d/L1 - R1*i_13d/L1 + Ω*i_13q
    du[4] = v_1q/L1 - v_3q/L1 - R1*i_13q/L1 - Ω*i_13d

    # # # ## L2 equation
    du[5] = v_3d/L2 - v_2d/L2 - R2*i_32d/L2 + Ω*i_32q
    du[6] = v_3q/L2 - v_2q/L2 - R2*i_32q/L2 - Ω*i_32d

    # ## R3 equation
    # du[7] = v_3d/Rs - i_30d
    # du[8] = v_3q/Rs - i_30q
    

    # # ## L3 equation
    du[7] = v_3d/Ls - Rs*i_30d/Ls + Ω*i_30q
    du[8] = v_3q/Ls - Rs*i_30q/Ls - Ω*i_30d


    du[9]  = v_1d - x_d_prime*i_13q + e_q_prime*sin(δ-π/2)
    du[10] = v_1q - e_q_prime*cos(δ-π/2) + i_13d*x_d_prime

    ## bus 1
    # du[11] = i_d - i_13d
    # du[12] = i_q - i_13q
    
    ## bus 3
    du[11] =  i_30d + i_32d - i_13d
    du[12] =  i_30q + i_32q - i_13q


end #function


function condition0(u,t,integrator,save_positions=(true,true))
    t in event_time0
end

function condition1(u,t,integrator,save_positions=(true,true))
    t in event_time1
end

function condition2(u,t,integrator,save_positions=(true,true))
    t in event_time2
end

function condition3(u,t,integrator,save_positions=(true,true))
    t in event_time3
end

function affect0!(integrator)
    @show "entering affect0"
    @show integrator.t

    SciMLBase.set_proposed_dt!(integrator,50e-6)
    integrator.opts.adaptive = false
end


function affect1!(integrator)
    @show "entering affect1"
    @show integrator.t
    
    integrator.p[3] = 1

    u_modified!(integrator,true)
end

function affect2!(integrator)
    @show "entering affect2"
    @show integrator.t
    
    integrator.p[3] = 0

    u_modified!(integrator,true)
end

function affect3!(integrator)
    @show "entering affect3"
    @show integrator.t

    integrator.opts.adaptive = true

    u_modified!(integrator,true)
end

function TDS_solve(params,u_pf,e_q_prime0)

    Y, Z_dict, fault_dict, param_dict = params

    # 1δ,2ω,3i_d,4i_q,5v_1d,6v_1q,7v_3d,8v_3q,9i_13d,10i_13q,11i_30d,12i_30q,13i_32d,14i_32q=u
    if param_dict["Framework"] == "DP Full"
        M = zeros(12,12)

        for i = 1:6
            M[i,i] = 1
        end


    elseif param_dict["Framework"] == "DP Phasor"
        M = zeros(12,12)
        M[1,1] = 1
        M[2,2] = 1
    end

    u0 = u_pf
    t_span = (0.0,fault_dict["tspan"])
    k = 0
    p = [Z_dict,fault_dict,k,e_q_prime0]

    global event_time0 = fault_dict["tf"] - (1/60)*10
    global event_time1 = fault_dict["tf"]
    global event_time2 = fault_dict["tc"]
    global event_time3 = fault_dict["tc"] + (1/60)*10

    v_1d_pf = 0
    v_1q_pf = 0

    v_3d_pf = 0
    v_3q_pf = 0

    cb0 = DiscreteCallback(condition0,affect0!)
    cb1 = DiscreteCallback(condition1,affect1!)
    cb2 = DiscreteCallback(condition2,affect2!)
    cb3 = DiscreteCallback(condition3,affect3!)

    cbs2 = CallbackSet(cb1,cb2)
    cbs = CallbackSet(cb0,cb1,cb2,cb3)
    
    
    prob0 = ODEFunction(fx!,mass_matrix = M)
    prob = ODEProblem(prob0,u0,t_span,p);
    
    tstop2 = [event_time1;event_time2]
    tstop1 = [event_time0;event_time1;event_time2;event_time3]

    if param_dict["dt"] == "Adaptive"
        @show "~~~~~ Using Adaptive Time step ~~~~~"
        sol = solve(prob, Rodas5(), callback = cbs, tstops = tstop1);
        # sol = solve(prob, Trapezoid(), callback=cbs, tstops=tstop1);
        # sol = solve(prob, ImplicitEuler(), callback=cbs, tstops=tstop1);
        # sol = solve(prob, TRBDF2(), callback=cbs, tstops=tstop1);

    elseif param_dict["dt"] == "Fixed"
        @show "~~~~~ Using Fixed Time step ~~~~~"
        dT = param_dict["step_size"]
        sol = solve(prob, Rodas5(), callback=cbs2, tstops=tstop2, dt = 50e-6, adaptive = false);
        # sol = solve(prob, Trapezoid(), callback=cbs2, tstops=tstop2, dt = dT, adaptive = false);
        # sol = solve(prob, ImplicitEuler(), callback=cbs2, tstops=tstop2, dt = dT, adaptive = false);
        # sol = solve(prob, TRBDF2(), callback=cbs2, tstops=tstop2, dt = dT, adaptive = false);
    else
        @show "Stepping rule not specified"
        return 0
    end
    
    # # sol = solve(prob, Rodas5());
    # sol = solve(prob, Rodas5(), callback=cbs, tstops=tstop1, dt = 50e-6, adaptive = false);
    # # sol = solve(prob, Rodas5(), callback=cbs, tstops=tstop1);
    
    if sol.retcode == ReturnCode.Success
        @show sol.retcode
        return sol
    else
        println("---------Integration failed ----------------- :(")
        @show sol.retcode
        return 0
    end
end


function solve_init(params,p0::Array{Float64})

    Y, Z_dict, fault_dict, param_dict = params

    M0 = zeros(7,7)
    tspan0 = (0.0,5)

    # δ,ω,e_q_prime,i_d,i_q,i_32d,i_32q
    u0_init = ones(7);

    f0 = ODEFunction(fxinit!,mass_matrix=M0)

    θ_10, v_30, θ_30, v_10, v_20, θ_20 = p0

    p0_init = (Y,Z_dict,p0)

    prob0 = ODEProblem(f0,u0_init,tspan0,p0_init)
    sol0 = solve(prob0,Rodas5(),reltol=1e-8)

    δ0,ω0,e_q_prime0,i_d0,i_q0,i_32d0,i_32q0 = sol0.u[end]

    v_1 = v_10*exp(im*θ_10)
    v_2 = v_20*exp(im*θ_20)
    v_3 = v_30*exp(im*θ_30)
    
    v_1d0 = real(v_1)
    v_1q0 = imag(v_1)
    
    v_2d0 = real(v_2)
    v_2q0 = imag(v_2)
    
    v_3d0 = real(v_3)
    v_3q0 = imag(v_3)

    i_30d0 = 0
    i_30q0 = 0

    i_13d0 = -i_d0
    i_13q0 = -i_q0

    # res = [δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_13d0,i_13q0,i_30d0,i_30q0,i_32d0,i_32q0]
    
    # res = [δ0,ω0,i_13d0,i_13q0,i_32d0,i_32q0,i_30d0,i_30q0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0]
    res = [δ0,ω0,i_13d0,i_13q0,i_32d0,i_32q0,i_30d0,i_30q0,v_1d0,v_1q0,v_3d0,v_3q0]

    # return (res,e_q_prime0)
    return (res,e_q_prime0)

end



function TDS_GENCLS(params,u_pf)
    u_init,e_q_prime0 = solve_init(params,u_pf)

    sol = TDS_solve(params, u_init, e_q_prime0)

end # function

end #module