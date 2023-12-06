using DifferentialEquations, Plots

function powerflow!(du,u,p,t)
    θ_1,v_3,θ_3=u

    p_1=0.944
    p_l3=0
    q_l3=0

    y_11=6.66666667
    y_12=0.0
    y_13=6.66666667
    y_21=y_12
    y_22=5.05050505
    y_23=5.05050505
    y_31=y_13
    y_32=y_23
    y_33=11.71717172

    ϕ_11=-1.57079633
    ϕ_12=0.0
    ϕ_13= 1.57079633
    ϕ_21=ϕ_12
    ϕ_22=-1.57079633
    ϕ_23=1.57079633
    ϕ_31=ϕ_13
    ϕ_32=ϕ_23
    ϕ_33=-1.57079633

    v_1=1.03
    θ_2=0.0
    v_2=1.0
    
    
    du[1] = p_1 - (v_1*v_1*y_11*cos(θ_1 - θ_1 -ϕ_11) + v_1*v_2*y_12*cos(θ_1 - θ_2 -ϕ_12) + v_1*v_3*y_13*cos(θ_1 - θ_3 -ϕ_13))
    # du[2] = i_d*v_1*cos(δ-θ_1) - i_q*v_1*sin(δ-θ_1) - (v_1*v_1*y_11*sin(θ_1 - θ_1 -ϕ_11) + v_1*v_3*y_13*sin(θ_1 - θ_3 -ϕ_13))
    du[2] = p_l3 - (v_3*v_1*y_31*cos(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*cos(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*cos(θ_3 - θ_3 -ϕ_33)) 
    du[3] = q_l3 - (v_3*v_1*y_31*sin(θ_3 - θ_1 -ϕ_31) + v_3*v_2*y_32*sin(θ_3 - θ_2 -ϕ_32) + v_3*v_3*y_33*sin(θ_3 - θ_3 -ϕ_33))
    # du[5] = p_l2 - (v_2*v_2*y_22*cos(θ_2 - θ_2 -ϕ_22) + v_2*v_3*y_23*cos(θ_2 - θ_3 -ϕ_23)) 
    # du[6] = q_l2 - (v_2*v_2*y_22*sin(θ_2 - θ_2 -ϕ_22) + v_2*v_3*y_23*sin(θ_2 - θ_3 -ϕ_23)) 
end

begin
    M0_pf = zeros(3,3)
t_span_pf=(0.0,5)
# θ_2,v_3,θ_3=u
u0_pf=[0,1.0,0];
pf = ODEFunction(powerflow!,mass_matrix=M0_pf)
pf0 = ODEProblem(pf,u0_pf,t_span_pf,());
sol_pf = solve(pf0,Rodas5(),abstol=1e-8,reltol=1e-8);
θ_10,v_30,θ_30=sol_pf.u[end]
end

function fxinit!(du,u,p,t)
    δ,ω,e_q_prime,i_d,i_q=u
    θ_10,v_30,θ_30=p
    
    # e_q_prime=1.1087
    x_d_prime=0.24
    p_m=0.944
    # p_e=p_m
    H=3.5

    y_11=6.66666667
    y_12=0.0
    y_13=6.66666667
    y_21=y_12
    y_22=5.05050505
    y_23=5.05050505
    y_31=y_13
    y_32=y_23
    y_33=11.71717172

    ϕ_11=-1.57079633
    ϕ_12=0.0
    ϕ_13= 1.57079633
    ϕ_21=ϕ_12
    ϕ_22=-1.57079633
    ϕ_23=1.57079633
    ϕ_31=ϕ_13
    ϕ_32=ϕ_23
    ϕ_33=-1.57079633

    v_1=1.03
    θ_1=θ_10
    v_2=1.0
    θ_2=0.0
    v_3=v_30
    θ_3=θ_30
    Ω=2π*60
    d=1
    # q_e=0.28818

    du[1] = Ω*(ω-1)
    # du[2] = (p_m - p_e - d*(ω-1))/2H
    du[2] = -v_1*sin(δ-θ_1) + x_d_prime*i_q
    du[3] = e_q_prime -v_1*cos(δ-θ_1) - x_d_prime*i_d
    # du[4] = i_d*v_1*sin(δ-θ_1) + i_q*v_1*cos(δ-θ_1) - p_e
    du[4] = i_d*sin(δ) + i_q*cos(δ) - (v_1*y_11*cos(θ_1 + ϕ_11) + v_2*y_12*cos(θ_2 + ϕ_12) + v_3*y_13*cos(θ_3 + ϕ_13))
    # du[5] = p_e - (v_1*v_1*y_11*cos(θ_1 - θ_1 -ϕ_11) + v_1*v_3*y_13*cos(θ_1 - θ_3 -ϕ_13))
    # du[5] = i_d*v_1*cos(δ-θ_1) - i_q*v_1*sin(δ-θ_1) - q_e
    du[5] = -i_d*cos(δ) + i_q*sin(δ) - (v_1*y_11*sin(θ_1 + ϕ_11) + v_2*y_12*sin(θ_2 + ϕ_12)+ v_3*y_13*sin(θ_3 + ϕ_13))
    # du[8] = q_e - (v_1*v_1*y_11*sin(θ_1 - θ_1 -ϕ_11) + v_1*v_3*y_13*sin(θ_1 - θ_3 -ϕ_13))
    # du[9] = p_l3 - (v_3*v_2*y_32*cos(θ_3 - θ_2 -ϕ_32) + v_3*v_1*y_31*cos(θ_3 - θ_1 -ϕ_31))
    # du[10] = q_l3 - (v_3*v_2*y_32*sin(θ_3 - θ_2 -ϕ_32) + v_3*v_1*y_31*sin(θ_3 - θ_1 -ϕ_31))
end

begin
    M0 = zeros(5,5)
t_span0=(0.0,5)
u0=zeros(5);

f0 = ODEFunction(fxinit!,mass_matrix=M0)
# θ_10,v_30,θ_30=p
p0=[θ_10,v_30,θ_30]
prob0 = ODEProblem(f0,u0,t_span0,p0);

sol0 = solve(prob0);

δ0,ω0,e_q_prime0,i_d0,i_q0=sol0.u[end]
end

function fx!(du,u,p,t)
    δ,ω,i_d,i_q,v_1d,v_1q,v_3d,v_3q,i_13d,i_13q,i_30d,i_30q,i_32d,i_32q=u
    # δ,ω,i_d,i_q,v_1d,v_1q,v_3d,v_3q,i_32d,i_32q=u
    # δ,ω,i_d,i_q,v_1,θ_1,v_3,θ_3=u
    k=p[1]
    
    e_q_prime0 = 1.1087455985400045
    e_q_prime=e_q_prime0
    x_d_prime=0.24
    p_m=0.944
    H=3.5
    Ω=2π*60
    d=0

    #X_L1 = 0.15 => L = 0.00039
    # L1 = 2*pi*f*L1 

    w_L1 = 0.15
    w_L2 = 0.198

    ## change short circuit reactance to resistance

    # R_s = k*(1e2) + abs(k-1)*(1e100)

    R_s = k*(1e-2) + abs(k-1)*1e10


    # R_s = 1e-2

    L1 = w_L1/Ω
    L2 = w_L2/Ω
    # L3 = w_L3/Ω

    # v_1 = abs(v_1d + im*v_1q)ac
    # θ_1 = angle(v_1d + im*v_1q)

    v_2=1.0
    θ_2=0.0

    v_2 = v_2*exp(im*θ_2)

    v_2d = real(v_2)
    v_2q = imag(v_2)




    # y_f = 1/(r_f+x_f)

    # p_fault = (v_3^2)*abs(y_f)*cos(-1*angle(y_f))
    # q_fault = (v_3^2)*abs(y_f)*sin(-1*angle(y_f))

    # s_fault = p_fault + q_fault*im
    # i_fault = -1*conj(s_fault/(v_3d+v_3q*im)) 

    # fault_i_real = k*real(i_fault)
    # fault_i_img =  k*im(i_fault)

    # s_fault = (p_fault+q_fault*im)
    
    # i_fault = -1*(y_f*(v_3d+v_3q*im))

    # fault_real = k*real(i_fault)
    # fault_im = k*imag(i_fault)

    # 1-δ,2-ω,3-i_d,4-i_q,5-v_1d,6-v_1q,7-v_3d,8-v_3q,9-i_30d,10-i_30q,11-i_32d,12-i_32q
    

    du[1] = Ω*(ω-1)
    du[2] = (p_m - (e_q_prime*(abs(v_1d + im*v_1q))*sin(δ-(angle(v_1d + im*v_1q)))/x_d_prime) - d*(ω-1))/(2H)
    
    du[3] = (abs(v_1d + im*v_1q))*sin(δ-(angle(v_1d + im*v_1q))) - x_d_prime*i_q
    du[4] = -e_q_prime + (abs(v_1d + im*v_1q))*cos(δ-(angle(v_1d + im*v_1q))) + x_d_prime*i_d

    ## bus 1
    du[5] = i_d*sin(δ) + i_q*cos(δ) - i_13d
    du[6] = i_q*sin(δ) - i_d*cos(δ) - i_13q

    ## L1 equation
    du[7] = (v_1d - v_3d + w_L1*i_13q)/L1
    du[8] = (v_1q - v_3q - w_L1*i_13d)/L1
    
    ## L2 equation
    du[9]  = (v_3d - v_2d + w_L2*i_32q)/L2
    du[10] = (v_3q - v_2q - w_L2*i_32d)/L2

    ## L3 equation
    du[11] = v_3d/R_s - i_30d 
    du[12] = v_3q/R_s - i_30q
    
    ## bus 2
    du[13] =  -1*i_30d - (i_32d - i_13d)
    du[14] =  -1*i_30q - (i_32q - i_13q)


end #function

########## Using Sauer Pai formulation ##############
begin
v_10=1.03
θ_10 = 0.32461479155703327

v_20=1.0
θ_20=0.0

v_30 = 1.004014153353707
θ_30 = 0.18725718301004352

δ0 = 0.5243271427959093
ω0 = 1.0

i_d0 = 0.41340918211007344
i_q0 = 0.8514126244761779

e_q_prime0 = 1.1087455985400045



# v_1d0 = v_10*sin(δ0 - θ_10)
# v_1q0 = v_10*cos(δ0 - θ_10)

# v_2d0 = v_20*sin(δ0 - θ_20)
# v_2q0 = v_20*cos(δ0 - θ_20)

# v_3d0 = v_30*sin(δ0 - θ_30)
# v_3q0 = v_30*cos(δ0 - θ_30) 


v_1 = v_10*exp(im*θ_10)
v_2 = v_20*exp(im*θ_20)
v_3 = v_30*exp(im*θ_30)

v_1d0 = real(v_1)
v_1q0 = imag(v_1)

v_2d0 = real(v_2)
v_2q0 = imag(v_2)

v_3d0 = real(v_3)
v_3q0 = imag(v_3)


w_L1 = 0.15
w_L2 = 0.198

i_32d0 = (v_3q0 - v_2q0)/(w_L2) 
i_32q0 = (-v_3d0 + v_2d0)/(w_L2)

i_30d0 = 0
i_30q0 = 0

i_13d0 = i_d0*sin(δ0) + i_q0*cos(δ0)
i_13q0 = i_q0*sin(δ0) - i_d0*cos(δ0)

# w_L3 = 1e2
# i_30d0 = v_3q0/w_L3
# i_30q0 = -v_3d0/w_L3
end

begin
    # u0=[δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0];
# δ,ω,i_d,i_q,v_1d,v_1q,v_3d,v_3q,i_30d,i_30q,i_32d,i_32q
u0 = [δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_13d0,i_13q0,i_30d0,i_30q0,i_32d0,i_32q0];
# u0=[δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_30d0,i_30q0,i_32d0,i_32q0];
# u0=[δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_32d0,i_32q0];

du0=zeros(14)

t_span=(0.0,15)

k = 0
p = [k]
end

begin
    M = zeros(14,14)

M[1,1] = 1
M[2,2] = 1
M[7,9] = 1
M[8,10] = 1
M[9,13] = 1
M[10,14] = 1


const event_time1=[10.0]
const event_time2=[10.1]

v_1d_pf = 0
v_1q_pf = 0

v_3d_pf = 0
v_3q_pf = 0

# i30_d_pf = 0
# i30_q_pf = 0



function condition1(u,t,integrator,save_positions=(true,true))
    t in event_time1
end

function condition2(u,t,integrator,save_positions=(true,true))
    t in event_time2
end

function affect1!(integrator)
    @show integrator.u[7] integrator.u[8] integrator.p integrator.t
    
    global v_1d_pf = integrator.u[5]
    global v_1q_pf = integrator.u[6]

    global v_3d_pf = integrator.u[7]
    global v_3q_pf = integrator.u[8]

    # global i30_3d_pf = integrator.u[7]
    # global i30_3q_pf = integrator.u[8]
    # integrator.u[7]=0.1
    # integrator.u[8]=0
    
    integrator.p[1] = 1

    u_modified!(integrator,true)
    # set_proposed_dt!(integrator,1e-15)
end

function affect2!(integrator)
    # integrator.u[7]=0.01
    @show integrator.u[7] integrator.u[8] integrator.p integrator.t
    
    integrator.p[1] = 0

    # integrator.u[3] = i_d0_pf
    # integrator.u[4] = i_q0_pf

    # integrator.u[5] = v_1d_pf
    # integrator.u[6] = v_1q_pf

    # integrator.u[7] = v_3d_pf
    # integrator.u[8] = v_3q_pf


    # integrator.u[11] = 0
    # integrator.u[12] = 0
    
    # @show integrator.u[7] integrator.u[8] integrator.t
    u_modified!(integrator,true)
end

cb1 = DiscreteCallback(condition1,affect1!)
cb2 = DiscreteCallback(condition2,affect2!)

cbs = CallbackSet(cb1,cb2)
end

prob0 = ODEFunction(fx!,mass_matrix = M)
prob = ODEProblem(prob0,u0,t_span,p);

const tstop1 = [10.0;10.1]
const tstop2 = [1]

# sol = solve(prob, Rodas5P(), callback=cbs, tstops=tstop1, dtmax = 1e-2)
# sol = solve(prob, Trapezoid(), callback=cbs, tstops=tstop1, dt = 500e-6, adaptive = false)
# sol = solve(prob, Rodas5P(), callback=cbs, tstops=tstop1, dt = 500e-6, adaptive = false)
sol = solve(prob, Rodas5P(), callback=cbs, tstops=tstop1)

# sol = solve(prob, Rodas5P())

# δ0,ω0,i_d0,i_q0,v_1d0,v_1q0,v_3d0,v_3q0,i_13d0,i_13q0,i_30d0,i_30q0,i_32d0,i_32q0

begin
v1d = sol[5,:]
v1q = sol[6,:]

v1_mag = similar(v1d)
v1_ang = similar(v1q)

@. v1_mag = abs(v1d + im*v1q)
@. v1_ang = angle(v1d + im*v1q)

# plotlyjs()
plot(sol.t,v1_mag,label="DP")
# plot(sol.t,v1_ang)
end

plot(sol.t[2:end]-sol.t[1:end-1])
# begin
#     v1d = sol[5,:]
# v1q = sol[6,:]

# v1_mag = similar(v1d)
# v1_ang = similar(v1q)

# @. v3_mag = abs(v3d + im*v3q)
# @. v3_ang = angle(v3d + im*v3q)
# end

# length(sol.t)

plotlyjs()
plot(sol,idxs=2,label="DP")