## File created 20.10.2026
## T-flow runs for the theta = 15 degrees phase diagram (J1 = 1, J2 = J1 * tan(theta)).
## J = [Jx, Jy, Jz] below is the SPIN ANISOTROPY (Jy = 1), applied to all pairs; the neighbour-shell
## strengths are the scalars [J1, J2] handed to getPyrochlore. Same convention as the theta = 5.749 runs.
using Pkg
Pkg.activate("../")

include("../src/PMFRG_xyz.jl") # This should not be needed because pmfrg_xyz will be in the env
using .PMFRG_xyz # this will become 'using PMFRG_xyz'

using SpinFRGLattices
using SpinFRGLattices.Pyrochlore
using JLD2
using DiffEqCallbacks
using ArgParse

struct save_struct{T}
    Chi_x::Array{T,2}
    Chi_y::Array{T,2}
    Chi_z::Array{T,2}
    Sigma_x::Array{T,2}
    Sigma_y::Array{T,2}
    Sigma_z::Array{T,2}
    f_int::Vector{T}
end

function main()
   args = parse_commandline()
   
   NLen  = args["nlen"]
   J = [args["Jx"], args["Jy"], args["Jz"]]
   
   theta_deg = args["theta-deg"]

   J1 = 1.0
   J2 = J1 * tan(deg2rad(theta_deg))      # theta in degrees -> radians; J2 = 0.2679
   
   
   System = getPyrochlore(NLen, [J1, J2])
   
   anisotropy = zeros(System.Npairs, 3)
   for n in 1:System.Npairs
       anisotropy[n, :] = J
   end
   
   Par = PMFRG_xyz.Params(
       System,
       N=args["N"],
       accuracy=1e-10,
       temp_max=args["tmax"],
       temp_min=args["tmin"]
   )
   
   mkpath("checkpoints_theta15")
   mkpath("data_theta15")
  
   if args["test-only"]
       exit()
   end

   sol, saved_values = PMFRG_xyz.SolveFRG(
       Par, anisotropy,
       SavedValues(PMFRG_xyz._getFloatType(Par), 
                   save_struct{PMFRG_xyz._getFloatType(Par)}),
       (State, t, integrator) -> save_all(State, t, Par, integrator)
   )
   name = filename(theta_deg,J,NLen,N)
   save_object("data/" * name, [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
   
end

include("./filename.jl")

function save_all(State, t, Par, integrator)
	chi_x = PMFRG_xyz.getChi_x(State, exp(t), Par, Par.NumericalParams.N)
	chi_y = PMFRG_xyz.getChi_y(State, exp(t), Par, Par.NumericalParams.N)
	chi_z = PMFRG_xyz.getChi_z(State, exp(t), Par, Par.NumericalParams.N)

	println("t: $(exp(t))")

	return save_struct(
		copy(chi_x), copy(chi_y), copy(chi_z),
		copy(State.x[2]), copy(State.x[3]), copy(State.x[4]), copy(State.x[1])
	)
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--nlen"
        required = true
        arg_type = Int64
        "--Jx"
        arg_type = Float64
        required = true
        "--Jy"
        arg_type = Float64
        required = true
        "--Jz"
        arg_type = Float64
        required = true
        "--theta-deg"
        arg_type = Float64
        default = 0
        "--N"
        arg_type= Int64
        default = 32
        "--tmin"
        arg_type = Float64
        default = 0.02
        "--tmax"
        arg_type = Float64
        default = 100.0
        "--test-only"
        action = :store_true
    end

    return parse_args(s)
end



main()
