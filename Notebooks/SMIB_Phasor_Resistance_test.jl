using DifferentialEquations, Plots

function fx!(du,u,p,t)
    δ,ω,i_d,i_q,v_1,θ_1,v_3,θ_3=u
    k,e_q_prime0=p

    e_q_prime=e_q_prime0
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
    
    #initial values
    # r_f=1e-2
    # x_f = 0

    #values in ANDES
    r_f = 0
    x_f = 0.0001im

    y_f = 1/(r_f+x_f)

 

    p_fault = (v_3^2)*abs(y_f)*cos(-1*angle(y_f))
    q_fault = (v_3^2)*abs(y_f)*sin(-1*angle(y_f))

 

    p_l3=0 - k*p_fault
    q_l3=0 - k*q_fault

 

    du[1] = Ω*(ω-1)
    du[2] = (p_m - (i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1)) - d*(ω-1))/(2H)
    du[3] = v_1*sin(δ-θ_1) - x_d_prime*i_q
    du[4] = -e_q_prime + v_1*cos(δ-θ_1) + x_d_prime*i_d
    du[5] = i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1) - (v_1*v_1*y_11*cos(θ_1 - θ_1 -ϕ_11) + v_1*v_2*y_12*cos(θ_1 - θ_2 -ϕ_12) + v_1*v_3*y_13*cos(θ_1 - θ_3 -ϕ_13))
    du[6] = i_d*v_1*cos(δ-θ_1) - i_q*v_1*sin(δ-θ_1) - (v_1*v_1*y_11*sin(θ_1 - θ_1 -ϕ_11) + v_1*v_2*y_12*sin(θ_1 - θ_2 -ϕ_12)+ v_1*v_3*y_13*sin(θ_1 - θ_3 -ϕ_13))
    du[7] = p_l3 - (v_3*v_1*y_31*cos(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*cos(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*cos(θ_3 - θ_3 -ϕ_33))
    du[8] = q_l3 - (v_3*v_1*y_31*sin(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*sin(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*sin(θ_3 - θ_3 -ϕ_33))
end

 

M = zeros(8,8)
M[1,1],M[2,2]=1,1

 

t_span=(0.0,15.0)

 
begin
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
end


u0=[δ0,ω0,i_d0,i_q0,v_10,θ_10,v_30,θ_30];


f = ODEFunction(fx!,mass_matrix=M)
k=0
p=[k,e_q_prime0]
prob = ODEProblem(f,u0,t_span,p);


begin
    u3=0
    u4=0
    u5=0
    u6=0
    u7=0
    u8=0
end


function condition1(u,t,integrator,save_positions=(true,true))
    t == 10.0
end

function affect1!(integrator)
    @show "entering affect1"
    @show integrator.t
    integrator.p[1] = 1
    global u3 = integrator.u[3]
    global u4 = integrator.u[4]
    global u5 = integrator.u[5]
    global u6 = integrator.u[6]
    global u7 = integrator.u[7]
    global u8 = integrator.u[8]

    u_modified!(integrator,true)
    # set_proposed_dt!(integrator, 1e-4) # does not seem to take effect
end


function condition2(u,t,integrator,save_positions=(true,true))
    t == 10.10
end

function affect2!(integrator)
    @show "entering affect2"
    @show integrator.t
    integrator.p[1] = 0
    set_proposed_dt!(integrator, 1e-4) # does not seem to take effect

    integrator.u[5] = u5
    integrator.u[6] = u6
    integrator.u[7] = u7
    integrator.u[8] = u8
    u_modified!(integrator,true)
    # probably need to do `set_u!` to restore pre-fault voltages 
end

cb1=ContinuousCallback(condition1,affect1!)
cb2=ContinuousCallback(condition2,affect2!)

# cb1 = DiscreteCallback(condition1,affect1!)
# cb2 = DiscreteCallback(condition2,affect2!)

cbs = CallbackSet(cb1, cb2)

const tstop1 = [10.0;10.1]

sol_ph=solve(prob, Trapezoid(), callback=cbs, tstops=tstop1, dt = 500e-6, adaptive = false)
# sol_ph=solve(prob, Trapezoid(), callback=cbs, tstops=tstop1)
# sol_ph=solve(prob,Rodas5P(),callback=cbs,tstops=[10.0, 10.10], dtmax=0.02)
# sol=solve(prob,Rodas5P(),callback=cbs,tstops=[0.1, 0.15], dtmax=0.02)

# δ,ω,i_d,i_q,v_1,θ_1,v_3,θ_3

# v1_mag = sol[5,:]
# plot(sol.t,v1_mag)

plot(sol_ph,idxs=5,label="Phasor")