using DifferentialEquations
using Sundials
using Plots

function fx!(res,du,u,p,t)
    δ,ω,i_d,i_q,v_1,θ_1,v_3,θ_3=u
    k=p[1]
    
    e_q_prime=1.1087455985400045
    x_d_prime=0.24
    p_m=0.944
    H=3.5

    y_12=0.0
    y_13=6.66666667
    y_11=y_12+y_13
    y_21=0.0
    y_23=5.05050505
    y_22=y_21+y_23
    y_31=6.66666667
    y_32=5.05050505
    y_33=y_31+y_32

    ϕ_11=-1.57079633
    ϕ_12=0.0
    ϕ_13= 1.57079633
    ϕ_21=0.0
    ϕ_22=-1.57079633
    ϕ_23=1.57079633
    ϕ_31= 1.57079633
    ϕ_32=1.57079633
    ϕ_33=-1.57079633
    
    v_2=1.0
    θ_2=0.0
    
    Ω=2π*60
    d=1

    r_f=0
    x_f=1e2im
    
    y_f = 1/(r_f+x_f)

    p_fault = (v_3^2)*abs(y_f)*cos(-1*angle(y_f))
    q_fault = (v_3^2)*abs(y_f)*sin(-1*angle(y_f))

    p_l3=0 - k*p_fault
    q_l3=0 - k*q_fault

    res[1] = Ω*(ω-1) - du[1]
    res[2] = (p_m - (i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1)) - d*(ω-1))/(2H) - du[2]
    res[3] = v_1*sin(δ-θ_1) - x_d_prime*i_q
    res[4] = -e_q_prime + v_1*cos(δ-θ_1) + x_d_prime*i_d
    res[5] = i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1) - (v_1*v_1*y_11*cos(θ_1 - θ_1 -ϕ_11) + v_1*v_2*y_12*cos(θ_1 - θ_2 -ϕ_12) + v_1*v_3*y_13*cos(θ_1 - θ_3 -ϕ_13))
    res[6] = i_d*v_1*cos(δ-θ_1) - i_q*v_1*sin(δ-θ_1) - (v_1*v_1*y_11*sin(θ_1 - θ_1 -ϕ_11) + v_1*v_2*y_12*sin(θ_1 - θ_2 -ϕ_12)+ v_1*v_3*y_13*sin(θ_1 - θ_3 -ϕ_13))
    res[7] = p_l3 - (v_3*v_1*y_31*cos(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*cos(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*cos(θ_3 - θ_3 -ϕ_33))
    res[8] = q_l3 - (v_3*v_1*y_31*sin(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*sin(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*sin(θ_3 - θ_3 -ϕ_33))
end

v_10=1.03
v_20=1.0
θ_20=0.0
δ0 = 0.5243271427959093
ω0 = 1.0
i_d0 = 0.41340918211007344
i_q0 = 0.8514126244761779
θ_10 = 0.32461479155703327
v_30 = 1.004014153353707
θ_30 = 0.18725718301004352
e_q_prime0 = 1.1087455985400045

u0=[δ0,ω0,i_d0,i_q0,v_10,θ_10,v_30,θ_30];
du0=[0,0,0,0,0,0,0,0]

t_span=(0.0,5.0)
k=0
p=[k]

dvs = [true,true,false,false,false,false,false,false]
prob = DAEProblem(fx!,du0,u0,t_span,differential_vars=dvs,p);

function condition1(u,t,integrator,save_positions=(true,true))
    t == 2.0
end

function affect1!(integrator)
    @show integrator.p[1] integrator.u[7] integrator.u[8]
    integrator.p[1] = 1
    # set_proposed_dt!(integrator,1e-15)
    integrator.opts.reltol = 0.01
    integrator.opts.abstol = 1e-2
    u_modified!(integrator, true)
end

# cb1 = DiscreteCallback(condition1,affect1!)
cb1 = ContinuousCallback(condition1,affect1!)
sol = solve(prob,IDA(),callback=cb1, tstops=[2.0], dtmax = 1e-2)

plot(sol,idxs=4)
# plot(sol.t[1:end-1],sol.t[2:end]-sol.t[1:end-1],y_formatter=:scientific)