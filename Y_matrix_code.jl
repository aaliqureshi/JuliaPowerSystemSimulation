module YMatrix

export build_Y


function build_Y(bus_impedance,n_bus)
    Y = zeros(ComplexF64,n_bus,n_bus)

    for ((bus1,bus2), impedance) in bus_impedance
        if bus1 != bus2
            Y[bus1,bus2] = -1/(impedance)
            Y[bus2,bus1] = -1/(impedance)
        else
            Y[bus1,bus2] = 1/(impedance)
        end
    end

    for idx = 1:n_bus
        Y[idx,idx] = sum(-1*Y[idx,:]) + 2*Y[idx,idx]
    end
    
    return Y

end

end #Module