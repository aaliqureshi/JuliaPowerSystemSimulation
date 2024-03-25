module Plotter
using Plots
using CSV, DataFrames
plotlyjs()

export plot_data
export save_data

function get_data(var,sol,params)

    params_dict = params[end]

    if var == "Omega"
        data = sol[2,:]
        return data
    end

    if params_dict["Gen Model"] == "AF" && var == "V1"
        vd = sol[15,:]
        vq = sol[16,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "AF" && var == "V3"
        vd = sol[17,:]
        vq = sol[18,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "AF" && var == "theta1"
        vd = sol[15,:]
        vq = sol[16,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    elseif params_dict["Gen Model"] == "AF" && var == "theta3"
        vd = sol[17,:]
        vq = sol[18,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    end


    if params_dict["Gen Model"] == "CLS" && var == "V1"
        vd = sol[9,:]
        vq = sol[10,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "CLS" && var == "V3"
        vd = sol[11,:]
        vq = sol[12,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "CLS" && var == "theta1"
        vd = sol[9,:]
        vq = sol[10,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    elseif params_dict["Gen Model"] == "CLS" && var == "theta3"
        vd = sol[11,:]
        vq = sol[12,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    end



end #get_data


function plot_data(var,sol,params;Interpolate=0,overlay=0)
    fault_dict, params_dict = params[end-1:end]
    t_end = fault_dict["tspan"]

    if Interpolate == 1
        t = 1:50e-6:t_end
        sol = sol(t)
    else
        t = sol.t
    end

    if var ∈ ["V1","V3","theta1","theta3","Omega"]
        data = get_data(var,sol,params)
    else
        return error("Variable does not exist")
    end

    label = params_dict["Framework"]*"("*params_dict["Gen Model"]*')'
    if overlay == 1
        plot!(t,data,label="", linewidth=1.5, legend=:bottomleft,tickfontsize=10)
        title!("$(var) - ($(label))")
        xlabel!("time")
    else
        plot(t,data,label="", linewidth=1.5, legend=:bottomleft,tickfontsize=10)
        title!("$(var) - ($(label))")
        xlabel!("time")
    end

end   #plot_data


function get_instant(data_rms,idx,t)
    idx1 = idx[1]
    idx2 = idx[2]
    # t = 1:50e-6:15
    ω = 2*π*60

    V = data_rms[:,idx1]
    θ = data_rms[:,idx2]

    data = Array{Float64}(undef,length(t))

    # @. data = V*cos(ω*t + θ)
    # @. data = 110*sqrt(2)*V*cos(ω*t + θ)

    @. data = 110*sqrt(2)*real(V*exp(1im*θ)*exp(1im*ω*t))
    # data = V


    # data = sqrt(2)*data

    # # convert back from pu
    # data = 110*data 

    return data
end


function save_data(sol,params;rms=1)

    fault_dict, params_dict = params[end-1:end]
    t_end = fault_dict["tspan"]

    folder = "./Working Codes/Cursor/Data Analysis/DP Data/"
    

    if rms == 1
        t = 1:50e-6:t_end

        dt = zeros(length(t))
        dt[1:length(sol.t)] = sol.t

        sol = sol(t)

        filename =  params_dict["Framework"]*"_"*params_dict["Gen Model"]*"_RMS"*".csv"
    
        
    
        data = Array{Float64}(undef,length(t),5)
    
        idx = 1
        for var ∈ ["V1","V3","theta1","theta3","Omega"]
            data[:,idx] = get_data(var,sol,params)
            idx+=1
        end
    
        df = DataFrame(t = sol.t, v1 = data[:,1], v3 = data[:,2], a1 = data[:,3], a3 = data[:,4], Omega = data[:,5], dt = dt)
        CSV.write(filename, df)
    end

    if rms == 0

        t = 1:50e-6:t_end
        dt = zeros(length(t))
        dt[1:length(sol.t)] = sol.t
        sol = sol(t)

        filename =  params_dict["Framework"]*"_"*params_dict["Gen Model"]*"_Time"*".csv"
    
        data_rms = Array{Float64}(undef,length(t),5)
        data_time = Array{Float64}(undef,length(t),2)
    
        idx = 1
        for var ∈ ["V1","V3","theta1","theta3","Omega"]
            data_rms[:,idx] = get_data(var,sol,params)
            idx+=1
        end

        idxs = [(1,3),(2,4)]

        ii = 1
        for idx in idxs
            data_time[:,ii] = get_instant(data_rms,idx,t)
            ii+=1
        end


        df = DataFrame(t = sol.t, v1 = data_time[:,1], v3 = data_time[:,2],Omega=data_rms[:,5], dt = dt)
        CSV.write(filename, df)

    end

    filename
end


end #module