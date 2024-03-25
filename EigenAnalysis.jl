module EigenAnalysis
using ForwardDiff, LinearAlgebra

include("TdsGenAFSmib.jl")
include("TdsGenclsSmib.jl")
# include("TdsGenclsSmib2.jl")


using .TdsGenAFSmib:solve_init as AF_init
using .TdsGenAFSmib:fx! as AF_fx!

using .TdsGenclsSmib:solve_init as CLS_init
using .TdsGenclsSmib:fx! as CLS_fx!
# using .TdsGenclsSmib2:solve_init as CLS_init
# using .TdsGenclsSmib2:fx! as CLS_fx!


export get_eigs


function fx_wrapper(fx!, u, p)
    du = similar(u)
    fx!(du, u, p, 0.0)
    return du
end


function State_Matrix(J;nx)
    fx = J[1:nx,1:nx]
    fy = J[1:nx,nx+1:end]
    gx = J[nx+1:end,1:nx]
    gy = J[nx+1:end,nx+1:end]

    As = fx - fy*inv(gy)*gx
    return As
    # return fx
end

function Eigz(J;nx)
    As = State_Matrix(J;nx=nx)
    λ = eigen(As).values
    return λ
end




function participation_factor(A)
    """
    the output matrix is relationship of k^{th} state variable with i^{th} eigenvalue

    formula taken from Sauer Pai
    """
    λ,Rv = eigen(A)
    _,Lv = eigen(A')
    w = Matrix(Lv)
    v = Matrix(Rv)
    nx = length(λ)
    PF = Matrix{Float64}(undef,nx,nx)
    for k=1:nx
        for i=1:nx
            nm = abs(v[k,i])*abs(w[k,i])
            dn = sum(abs.(v[1:nx,i]).*abs.(w[1:nx,i]))
            PF[k,i] = nm/dn
        end
    end
    return λ, PF

end #function




function get_eigs(params, u_pf)
    Y, Z_dict, fault_dict, param_dict = params

    if param_dict["Gen Model"] == "CLS"
        u_init, e_q_prime0 = CLS_init(params,u_pf)

        k = 0
        p = [Z_dict,fault_dict,k,e_q_prime0]

        J_dp = ForwardDiff.jacobian(u -> fx_wrapper(CLS_fx!, u, p), u_init)

        if param_dict["Framework"] == "DP Full"
            nx = 6
        elseif param_dict["Framework"] == "DP Phasor"
            nx = 2
        end

        λ = Eigz(J_dp,nx=nx)

        eig, p_factor = participation_factor(State_Matrix(J_dp,nx=nx))

        return λ, eig, p_factor

    elseif param_dict["Gen Model"] == "AF"
        u_init, v_f0 = AF_init(params,u_pf)

        k = 0
        p = [Z_dict,fault_dict,k,v_f0]

        J_dp = ForwardDiff.jacobian(u -> fx_wrapper(AF_fx!, u, p), u_init)

        if param_dict["Framework"] == "DP Full"
            nx = 10
        elseif param_dict["Framework"] == "DP Phasor"
            nx = 6
        end

        λ = Eigz(J_dp,nx=nx)

        eig, p_factor = participation_factor(State_Matrix(J_dp,nx=nx))

        return λ, eig, p_factor

    end

end #get_eigs

end #module