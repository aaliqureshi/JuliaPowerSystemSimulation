module TdsGenAFSmib2
# using DifferentialEquations
using OrdinaryDiffEq

export TDS_GENAF

function fxinit!(du,u,p,t)
    # @show "fixinit start"
    δ,ω,v_f,e_q_prime,e_d_prime,e_q_2prime,e_d_2prime,i_d,i_q,i_32d,i_32q = u

    Y = p[1]
    Z = p[2]
    θ_1, v_3, θ_3, v_1, v_2, θ_2 = p[3]

    R1 = Z["R1"]
    R2 = Z["R2"]
    w_L1 = Z["XL1"]
    w_L2 = Z["XL2"]
    
    p_m = 0.944
    H = 3.5

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

    Ω = 2π*60
    d = 0

    x_d = 0.245
    x_d_prime = 0.24
    T_d_prime = 15.5
    x_d_2prime = 0.22
    T_d_2prime = 100.029

    x_q = 0.245
    x_q_prime = 0.24
    T_q_prime = 0.85
    x_q_2prime = 0.22
    T_q_2prime = 100.029


    du[1] = Ω*(ω-1)
    # du[2] = (p_m - (i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1)) - d*(ω-1))/(2H)
    du[2] = (p_m - (i_d*v_1*cos(θ_1) + i_q*v_1*sin(θ_1)) - d*(ω-1))/(2H)
    du[3] = (-e_q_prime - (x_d - x_d_prime)*i_d + v_f)/T_d_prime
    du[4] = (-e_d_prime + (x_q - x_q_prime)*i_q)/T_q_prime
    du[5] = (-e_q_2prime + e_q_prime - (x_d_prime - x_d_2prime)*i_d)/T_d_2prime
    du[6] = (-e_d_2prime + e_d_prime + (x_q_prime - x_q_2prime)*i_q)/T_q_2prime

    # du[7] = v_1*cos(θ_1) - e_d_2prime*sin(δ) - e_q_2prime*cos(δ) + i_d*x_d_2prime*cos(δ) - i_q*x_q_2prime*sin(δ)
    du[7] = v_1*cos(θ_1) - e_d_2prime*sin(δ) - e_q_2prime*cos(δ) + i_d*x_d_2prime*cos(δ) - i_q*x_q_2prime*sin(δ)
    du[8] = v_1*sin(θ_1) + e_d_2prime*cos(δ) - e_q_2prime*sin(δ) + i_d*x_d_2prime*sin(δ) + i_q*x_q_2prime*cos(δ)
    du[9] = i_d*sin(δ) + i_q*cos(δ) - (v_1*y_11*cos(θ_1 + ϕ_11) + v_2*y_12*cos(θ_2 + ϕ_12) + v_3*y_13*cos(θ_3 + ϕ_13))

    du[10] = i_32d - (v_3*sin(θ_3) - v_2*sin(θ_2) - R2*i_32q)/(w_L2)
    du[11] = i_32q - (-v_3*cos(θ_3) + v_2*cos(θ_2) + R2*i_32d)/(w_L2)

    # @show "fixinit end"
end


function fx!(du,u,p,t)
    # δ,ω,e_q_prime,e_d_prime,e_q_2prime,e_d_2prime,i_d,i_q,v_1d,v_1q,v_3d,v_3q,i_13d,i_13q,i_30d,i_30q,i_32d,i_32q = u
    δ,ω,e_q_prime,e_d_prime,e_q_2prime,e_d_2prime,i_13d,i_13q,i_32d,i_32q,i_30d,i_30q,i_d,i_q,v_1d,v_1q,v_3d,v_3q = u

    Z, fault_dict, k, v_f0 = p
    
    v_f = v_f0
    p_m = 0.944
    H = 3.5
    Ω = 2π*60
    
    d = 0

    R1 = Z["R1"]
    R2 = Z["R2"]
    w_L1 = Z["XL1"]
    w_L2 = Z["XL2"]

    w_L3 = fault_dict["XLs"]
    R3 = fault_dict["Rs"]

    L1 = w_L1/Ω
    L2 = w_L2/Ω
    L3 = w_L3/Ω

    # w_L3 = fault_dict["XLs"]
    # R3 = fault_dict["Rs"]

    Rs = k*R3 + abs(k-1)*1e100
    Ls = k*L3 + abs(k-1)*1e100
 
    x_d = 0.245
    x_d_prime = 0.24
    T_d_prime = 15.5
    x_d_2prime = 0.22
    T_d_2prime = 100.029

    x_q = 0.245
    x_q_prime = 0.24
    T_q_prime = 0.85
    x_q_2prime = 0.22
    T_q_2prime = 100.029


    v_2 = 1.0
    θ_2 = 0.0

    v_2 = v_2*exp(im*θ_2)

    v_2d = real(v_2)
    v_2q = imag(v_2)


    du[1] = Ω*(ω-1)
    du[2] = (p_m - (i_d*v_1d*sin(δ) - i_d*v_1q*cos(δ) + i_q*v_1d*cos(δ) + i_q*v_1q*sin(δ)) - d*(ω-1))/(2H)

    du[3] = (-e_q_prime - (x_d - x_d_prime)*i_d + v_f)/T_d_prime
    du[4] = (-e_d_prime + (x_q - x_q_prime)*i_q)/T_q_prime
    du[5] = (-e_q_2prime + e_q_prime - (x_d_prime - x_d_2prime)*i_d)/T_d_2prime
    du[6] = (-e_d_2prime + e_d_prime + (x_q_prime - x_q_2prime)*i_q)/T_q_2prime

    # ## L1 equation
    du[7] = v_1d/L1 - v_3d/L1 - R1*i_13d/L1 + Ω*i_13q
    du[8] = v_1q/L1 - v_3q/L1 - R1*i_13q/L1 - Ω*i_13d

    # ## L2 equation
    du[9]  = v_3d/L2 - v_2d/L2 - R2*i_32d/L2 + Ω*i_32q
    du[10]  = v_3q/L2 - v_2q/L2 - R2*i_32q/L2 - Ω*i_32d

    # # # L3 equation
    # du[11]  = (v_3d - Rs*i_30d)/Ls + Ω*i_30q
    # du[12]  = (v_3q - Rs*i_30q)/Ls - Ω*i_30d

    # # ## R3 equation
    du[11]  = v_3d/Rs - i_30d
    du[12]  = v_3q/Rs - i_30q
    

    du[13] = v_1d - e_d_2prime*sin(δ) - e_q_2prime*cos(δ) + i_d*x_d_2prime*cos(δ) - i_q*x_q_2prime*sin(δ)
    du[14] = v_1q + e_d_2prime*cos(δ) - e_q_2prime*sin(δ) + i_d*x_d_2prime*sin(δ) + i_q*x_q_2prime*cos(δ)

    ## bus 1
    du[15] = i_d*sin(δ) + i_q*cos(δ) - i_13d
    du[16] = i_q*sin(δ) - i_d*cos(δ) - i_13q


    ## bus 3
    du[17] =  i_30d + i_32d - i_13d
    du[18] =  i_30q + i_32q - i_13q

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

v1d_pre_fault = 0
v1q_pre_fault = 0
v3d_pre_fault = 0
v3q_pre_fault = 0

function affect0!(integrator)
    # @show integrator.u[7] integrator.u[8] integrator.p integrator.t
    @show "entering affect0"
    @show integrator.t

    SciMLBase.set_proposed_dt!(integrator,50e-6)
    integrator.opts.adaptive = false
end

function affect1!(integrator)
    # @show integrator.u[7] integrator.u[8] integrator.p integrator.t
    @show "entering affect1"
    @show integrator.t

    integrator.p[3] = 1

    # global v1d_pre_fault = integrator.u[9]
    # global v1q_pre_fault = integrator.u[10]
    # global v3d_pre_fault = integrator.u[11]
    # global v3q_pre_fault = integrator.u[12]

    u_modified!(integrator,true)
end

function affect2!(integrator)
    @show "entering affect2"
    @show integrator.t
    
    integrator.p[3] = 0

    

    # integrator.u[9] = v1d_pre_fault
    # integrator.u[10] = v1q_pre_fault
    # integrator.u[11] = v3d_pre_fault
    # integrator.u[12] = v3q_pre_fault



    u_modified!(integrator,true)
end

function affect3!(integrator)
    # @show integrator.u[7] integrator.u[8] integrator.p integrator.t
    @show "entering affect3"
    @show integrator.t

    # integrator.p[3] = 1

    # SciMLBase.set_proposed_dt!(integrator,50e-4)
    integrator.opts.adaptive = true

    # global v1d_pre_fault = integrator.u[9]
    # global v1q_pre_fault = integrator.u[10]
    # global v3d_pre_fault = integrator.u[11]
    # global v3q_pre_fault = integrator.u[12]

    u_modified!(integrator,true)
end


function TDS_solve(params,u_pf,v_f0)

    Y, Z_dict, fault_dict, param_dict = params

    # 1δ,2ω,3e_q_prime,4e_d_prime,5e_q_2prime,6e_d_2prime,7i_d,8i_q,9v_1d,10v_1q,11v_3d,12v_3q,13i_13d,14i_13q,15i_30d,16i_30q,17i_32d,18i_32q
    
    # 2->  1δ,2ω,3e_q_prime,4e_d_prime,5e_q_2prime,6e_d_2prime,7i_13d,8i_13q,
    #9i_32d,10i_32q,11i_30d,12i_30q,13i_d,14i_q,15v_1d,16v_1q,17v_3d,18v_3q=u
    if param_dict["Framework"] == "DP Full"
        M = zeros(18,18)

        for i in 1:10
            M[i,i] = 1
        end
        # @show M

    elseif param_dict["Framework"] == "DP Phasor"
        M = zeros(18,18)

        for i in 1:6
            M[i,i] = 1
        end

    end

    u0 = u_pf
    t_span = (0.0,fault_dict["tspan"])
    k = 0
    p = [Z_dict,fault_dict,k,v_f0]




    global event_time0 = fault_dict["tf"] - (1/60)*10
    global event_time1 = fault_dict["tf"]
    global event_time2 = fault_dict["tc"]
    global event_time3 = fault_dict["tc"] + (1/60)*10


    cb0 = DiscreteCallback(condition0,affect0!)
    cb1 = DiscreteCallback(condition1,affect1!)
    cb2 = DiscreteCallback(condition2,affect2!)
    cb3 = DiscreteCallback(condition3,affect3!)

    # cb1 = ContinuousCallback(condition1,affect1!)
    # cb2 = ContinuousCallback(condition2,affect2!)

    # cbs = CallbackSet(cb1,cb2)
    cbs = CallbackSet(cb0,cb1,cb2,cb3)
    
    
    prob0 = ODEFunction(fx!,mass_matrix = M)
    prob = ODEProblem(prob0,u0,t_span,p);
    
    # tstop1 = [event_time1;event_time2]
    tstop1 = [event_time0;event_time1;event_time2;event_time3]

    
    if param_dict["dt"] == "Adaptive"
        sol = solve(prob, Rodas5(), callback = cbs, tstops = tstop1);
    elseif param_dict["dt"] == "Fixed"
        sol = solve(prob, Rodas5(), callback=cbs, tstops=tstop1, dt = 50e-6, adaptive = false);
    else
        @show "Stepping rule not specified"
        return 0
    end

    
    # # sol = solve(prob, Rodas5());
    # # sol = solve(prob, Rodas5(), callback=cbs, tstops=tstop1, dt = 50e-6, adaptive = false);
    # # sol = solve(prob, Trapezoid(), callback=cbs, tstops=tstop1, dt = 50e-6, adaptive = false);
    # sol = solve(prob, Rodas5(), callback=cbs, tstops=tstop1);
    
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
    @show "Iniside init call"

    M0 = zeros(11,11)
    tspan0 = (0.0,5)

    # δ,ω,v_f,e_q_prime,e_d_prime,e_q_2prime,e_d_2prime,i_d,i_q,i_32d,i_32q
    u0_init = ones(11);

    f0 = ODEFunction(fxinit!,mass_matrix=M0)

    θ_10, v_30, θ_30, v_10, v_20, θ_20 = p0

    p0_init=(Y,Z_dict,p0)

    prob0 = ODEProblem(f0,u0_init,tspan0,p0_init)
    sol0 = solve(prob0,Rodas5(),reltol=1e-8)
    @show "Init solved"

    δ0,ω0,v_f0,e_q_prime0,e_d_prime0,e_q_2prime0,e_d_2prime0,i_d0,i_q0,i_32d0,i_32q0 = sol0.u[end]

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

    i_13d0 = i_d0*sin(δ0) + i_q0*cos(δ0)
    i_13q0 = i_q0*sin(δ0) - i_d0*cos(δ0)


    # res = [δ0,ω0,e_q_prime0,e_d_prime0,e_q_2prime0,e_d_2prime0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_13d0,i_13q0,i_30d0,i_30q0,i_32d0,i_32q0]

    # δ,ω,e_q_prime,e_d_prime,e_q_2prime,e_d_2prime,i_13d,i_13q,i_32d,i_32q,i_30d,i_30q,i_d,i_q,v_1d,v_1q,v_3d,v_3q
    res = [δ0,ω0,e_q_prime0,e_d_prime0,e_q_2prime0,e_d_2prime0,i_13d0,i_13q0,i_32d0,i_32q0,i_30d0,i_30q0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0]

    return (res,v_f0)

end



function TDS_GENAF(params,u_pf)
    u_init,v_f0 = solve_init(params,u_pf)
    @show "Init Complete"

    sol = TDS_solve(params, u_init, v_f0)

end # function

end #module