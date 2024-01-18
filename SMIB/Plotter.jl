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
        vd = sol[9,:]
        vq = sol[10,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "AF" && var == "V3"
        vd = sol[11,:]
        vq = sol[12,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "AF" && var == "theta1"
        vd = sol[9,:]
        vq = sol[10,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    elseif params_dict["Gen Model"] == "AF" && var == "theta3"
        vd = sol[11,:]
        vq = sol[12,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    end


    if params_dict["Gen Model"] == "CLS" && var == "V1"
        vd = sol[5,:]
        vq = sol[6,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "CLS" && var == "V3"
        vd = sol[7,:]
        vq = sol[8,:]
        v_mag = abs.(vd + im*vq)
        v_mag = clamp.(v_mag,-5,5)
        return v_mag
    elseif params_dict["Gen Model"] == "CLS" && var == "theta1"
        vd = sol[5,:]
        vq = sol[6,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    elseif params_dict["Gen Model"] == "CLS" && var == "theta3"
        vd = sol[7,:]
        vq = sol[8,:]
        v_ang = angle.(vd + im*vq)
        v_ang = clamp.(v_ang,-5,5)
        return v_ang
    end



end #get_data


function plot_data(var,sol,params;Interpolate=0,overlay=0)
    params_dict = params[end]

    if Interpolate == 1
        t = 1:50e-6:15
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
        plot!(t,data,label=label, linewidth=1, legend=:bottomleft)
        title!(var)
        xlabel!("time")
    else
        plot(t,data,label=label, linewidth=1, legend=:bottomleft)
        title!(var)
        xlabel!("time")
    end

end   #plot_data



function save_data(sol,params)

    params_dict = params[end]

    t = 1:50e-6:15
    sol = sol(t)

    filename =  params_dict["Framework"]*"_"*params_dict["Gen Model"]*".csv"

    data = Array{Float64}(undef,length(t),5)

    idx = 1
    for var ∈ ["V1","V3","theta1","theta3","Omega"]
        data[:,idx] = get_data(var,sol,params)
        idx+=1
    end

    df = DataFrame(t = sol.t, v1 = data[:,1], v3 = data[:,2], a1 = data[:,3], a3 = data[:,4], Omega = data[:,5])
    CSV.write(filename, df)
    
    filename
end


end #module