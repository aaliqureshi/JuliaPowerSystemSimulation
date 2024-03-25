# using DifferentialEquations, Plots
using OrdinaryDiffEq, Plots

include("Y_matrix_code.jl")
include("PowerFlow_SMIB.jl")
include("TDS_SMIB.jl")
include("Plotter.jl")
include("EigenAnalysis.jl")

using .YMatrix
using .PowerFlow
using .TDS
using .Plotter: plot_data, save_data
using .EigenAnalysis

begin
    n_bus = 3
    n_lines = 2

    v_10 = 1.03
    v_2 = 1.0
    θ_2 = 0.0


    R1 = 0.015
    # R1 = 1e-18
    XL1 = 0.15

    R2 = 0.0198
    # R2 = 1e-18
    XL2 = 0.198

    
    Z_dict = Dict([("R1", R1),
                   ("R2", R2),
                   ("XL1", XL1),
                   ("XL2", XL2)])

    fault_dict = Dict([("Rs", 0),
                       ("XLs", 1e-6),
                       ("tf", 2),
                       ("tc",2.1),
                       ("tspan", 4)])

    
    param_dict = Dict([("Gen Model", "CLS"),
                      ("Framework","DP Full"),
                      ("dt", "Adaptive"),
                      ("step_size",5e-6)])
    
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

sol = solve_TDS(params,u_pf)

λ, eig, p_factor = get_eigs(params, u_pf)

@show eig
@show p_factor

# save_data(sol,params;rms=0)

# plot(sol,idxs=2)
plot_data("V1",sol,params,Interpolate=1,overlay=1)
plot_data("V3",sol,params,Interpolate=1,overlay=0)
plot_data("Omega",sol,params,Interpolate=1,overlay=0)
plot_data("theta1",sol,params,Interpolate=1,overlay=0)
plot_data("theta3",sol,params,Interpolate=1,overlay=0)

begin
    t_dp = 1:50e-6:3
    sol2 = sol(t_dp)

    idd1_dp =  sol2[3,:]
    idd2_dp =  sol2[5,:]
    iqq1_dp =  sol2[4,:]
    iqq2_dp =  sol2[6,:]
end
plot(t_dp,idd1_dp)
plot(t_dp,-1*idd2_dp)
plot(t_dp,iqq1_dp)
plot(t_dp,-1*iqq2_dp)


begin
    # load and plot andes data

    dir = pwd()
    # file1 = dir*"/Andes data/Algeb.npz"
    # file1 = dir*"/Andes data/Andes.npz"
    # file1 = dir*"/Andes data/NewAlgeb.npz"


    # file1 = dir*"/Andes data/State_full.npz"
    # file1 = dir*"/Andes data/State_stime.npz"
    # file1 = dir*"/Andes data/State_fixed.npz"
    # file1 = dir*"/Andes data/State_stable.npz"
    # file1 = dir*"/Andes data/State_New.npz"
    # file1 = dir*"/Andes data/State_New2.npz"
    file1 = dir*"/Andes data/State_5us.npz"
    
    using NPZ

    data_algeb = npzread(file1)

    if occursin("State",file1)
        t_algeb = data_algeb["data"][:,1]
        omega_algeb = data_algeb["data"][:,8]
        idd1_algeb = data_algeb["data"][:,2]
        idd2_algeb = data_algeb["data"][:,3]
        iqq1_algeb = data_algeb["data"][:,4]
        iqq2_algeb = data_algeb["data"][:,5]
        igend_algeb = data_algeb["data"][:,21]

        v1_algeb = data_algeb["data"][:,13]
        a1_algeb = data_algeb["data"][:,10]

        v3_algeb = data_algeb["data"][:,15]
        a3_algeb = data_algeb["data"][:,12]
    else
        t_algeb = data_algeb["data"][:,1]
        omega_algeb = data_algeb["data"][:,4]
        idd1_algeb = data_algeb["data"][:,15]
        idd2_algeb = data_algeb["data"][:,16]
        iqq1_algeb = data_algeb["data"][:,17]
        iqq2_algeb = data_algeb["data"][:,18]
        igend_algeb = data_algeb["data"][:,21]
    end



end

begin
    t_dp = 1:50e-6:15
    sol2 = sol(t_dp)

    omega_dp = sol2[2,:]
    idd1_dp =  sol2[3,:]
    idd2_dp =  sol2[5,:]
    iqq1_dp =  sol2[4,:]
    iqq2_dp =  sol2[6,:]

    v1d = sol2[9,:]
    v1q = sol2[10,:]

    v3d = sol2[11,:]
    v3q = sol2[12,:]

    v1_dp = abs.(v1d+1im*v1q)
    a1_dp = angle.(v1d+1im*v1q)
    v3_dp = abs.(v3d+1im*v3q)
    a3_dp = angle.(v3d+1im*v3q)
    
end


begin
    plot(t_dp,omega_dp,label="DP")
    plot!(t_algeb,omega_algeb,label="Andes")
    # vline!([6.1],label="fault clear time",color=:black)
    title!("omega")
end

begin
    plot(t_dp,idd1_dp,label="DP")
    plot!(t_algeb,idd1_algeb,label="Andes")
    # vline!([6.1],label="fault clear time",color=:black)
    title!("id (line 1)")
end

begin
    plot(t_dp,-1*idd2_dp,label="DP")
    plot!(t_algeb,idd2_algeb,label="Andes")
    title!("id (line 2)")
end

begin
    plot(t_dp,iqq1_dp,label="DP")
    plot!(t_algeb,iqq1_algeb,label="Andes")
    title!("iq (line 1)")
end

begin
    plot(t_dp,-1*iqq2_dp,label="DP")
    plot!(t_algeb,iqq2_algeb,label="Andes")
    title!("iq (line 2)")
end


# begin
#     plot(t_dp,idd1_dp,label="DP")
#     plot!(t_algeb,igend_algeb,label="Andes")   
# end


begin
    plot(t_dp,v1_dp,label="DP")
    plot!(t_algeb,v1_algeb,label="Andes")  
    title!("V1") 
end

begin
    plot(t_dp,a1_dp,label="DP")
    plot!(t_algeb,a1_algeb,label="Andes")   
    title!("θ1")
end


begin
    plot(t_dp,v3_dp,label="DP")
    plot!(t_algeb,v3_algeb,label="Andes")  
    title!("V3") 
end

begin
    plot(t_dp,a3_dp,label="DP")
    plot!(t_algeb,a3_algeb,label="Andes")   
    title!("θ3")
end



begin
    plotlyjs()
    t = 1:50e-6:15
    sol2 = sol(t)

    # vd = sol2[15,:]
    # vq = sol2[16,:]
    vd = sol2[3,:]
    vq = sol2[4,:]

    v_mag = abs.(vd + im*vq)
    v_mag = clamp.(v_mag,-5,5)

    v_ang = angle.(vd + im*vq)
    v_ang = clamp.(v_ang,-5,5)

    ω = 2π*60

    data = Array{Float64}(undef,length(t))

    # @. data = V*cos(ω*t + θ)
    # @. data = 110*sqrt(2)*V*cos(ω*t + θ)

    # @. data = 110*sqrt(2)*real(v_mag*exp(1im*v_ang)*exp(1im*ω*t))
    @. data = 110*sqrt(2)*v_mag*cos(ω*t + v_ang)

    plot(t,data)
end

begin
    plotlyjs()
    t = 1:50e-6:15
    sol2 = sol(t)

    vd = sol2[3,:]
    vq = sol2[4,:]

    v_mag = abs.(vd + im*vq)
    v_mag = clamp.(v_mag,-5,5)

    v_ang = angle.(vd + im*vq)
    v_ang = clamp.(v_ang,-5,5)

    # ω = 2π*60

    data = Array{Float64}(undef,length(t))

    @. data = v_mag*cos(v_ang)

    plot(t,data)
end


begin
    t = sol.t

    id1 = 10
    id2 = length(t)
    # # id2 = 500

    plot(t[id1:id2-1],t[id1+1:id2] - t[id1:id2-1],label="")
    # plot(t[id1+1:id2] - t[id1:id2-1])  
    # dt = (t[2:end] - t[1:end-1])
    # xx = range(1,15,length(log_dt))
    
    # plot(xx,log_dt,label="")
    title!("DiffEq time step - $(param_dict["Gen Model"])($(param_dict["Framework"]))")
    ylabel!("dt")
    xlabel!("time")
end

λ, eig, p_factor = get_eigs(params, u_pf)

@show eig
@show p_factor

max_eig = maximum(real.(eig))
min_eig = minimum(real.(eig))

stiffness = max_eig/min_eig


# save_data(sol,params;rms=0)