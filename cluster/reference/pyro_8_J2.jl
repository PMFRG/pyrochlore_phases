## File created 20.10.2026
## T-flow runs for the theta = 15 degrees phase diagram (J1 = 1, J2 = J1 * tan(theta)).
## J = [Jx, Jy, Jz] below is the SPIN ANISOTROPY (Jy = 1), applied to all pairs; the neighbour-shell
## strengths are the scalars [J1, J2] handed to getPyrochlore. Same convention as the theta = 5.749 runs.

using Pkg
Pkg.activate("../")

include("../src/PMFRG_xyz.jl")
using .PMFRG_xyz

using SpinFRGLattices
using SpinFRGLattices.Pyrochlore
using JLD2

using DiffEqCallbacks

struct save_struct{T}
    Chi_x::Array{T,2}
    Chi_y::Array{T,2}
    Chi_z::Array{T,2}
    Sigma_x::Array{T,2}
    Sigma_y::Array{T,2}
    Sigma_z::Array{T,2}
    f_int::Vector{T}
end

counter = Ref(0)
function save_all(State, t, Par, integrator, checkfile="")
	chi_x = PMFRG_xyz.getChi_x(State, exp(t), Par, Par.NumericalParams.N)
	chi_y = PMFRG_xyz.getChi_y(State, exp(t), Par, Par.NumericalParams.N)
	chi_z = PMFRG_xyz.getChi_z(State, exp(t), Par, Par.NumericalParams.N)

	println("t: $(exp(t))")
	counter[] += 1

	if(checkfile != "")
		jldopen(checkfile, "a+") do file
			file["checkpoint/$(counter[])/chix"] = copy(chi_x)
			file["checkpoint/$(counter[])/chiy"] = copy(chi_y)
			file["checkpoint/$(counter[])/chiz"] = copy(chi_z)
			file["checkpoint/$(counter[])/sigx"] = copy(State.x[2])
			file["checkpoint/$(counter[])/sigy"] = copy(State.x[3])
			file["checkpoint/$(counter[])/sigz"] = copy(State.x[4])
		end
	end

	return save_struct(
		copy(chi_x), copy(chi_y), copy(chi_z),
		copy(State.x[2]), copy(State.x[3]), copy(State.x[4]), copy(State.x[1])
	)
end

NLen = 8

theta_deg = 15
J1 = 1.0
J2 = J1 * tan(deg2rad(theta_deg))      # theta in degrees -> radians; J2 = 0.2679

# Spin anisotropies [Jx, 1, Jz] to run (only Jz <= Jx; the Jx < Jz half is the mirror image).
# First scan, the phase boundary is NOT known yet. At theta = 5.749 (J2 = 0.1) AIAO Z orders for Jz <~ -0.13,
# at J2 = 0 for Jz <~ -0.18; extrapolating, the boundary should have moved towards Jz ~ 0 or above here.
# So every Jx column brackets Jz = -0.1 ... 0.4, with wider steps where the boundary is least certain.
J_array = [
# First run
	# Jx = 0 (Jz <= 0 only)
#	[0.0, 1.0, -0.2],
#	[0.0, 1.0, -0.1],
#	[0.0, 1.0, 0.0],
#	# Jx = 0.25
#	[0.25, 1.0, -0.1],
#	[0.25, 1.0, 0.0],
#	[0.25, 1.0, 0.1],
#	[0.25, 1.0, 0.2],
#	# Jx = 0.5
#	[0.5, 1.0, -0.1],
#	[0.5, 1.0, 0.0],
#	[0.5, 1.0, 0.15],
#	[0.5, 1.0, 0.3],
#	# Jx = 0.75
#	[0.75, 1.0, -0.1],
#	[0.75, 1.0, 0.05],
#	[0.75, 1.0, 0.2],
#	[0.75, 1.0, 0.4],
	# the tip on the diagonal (J2 = 0: Jx = Jz ~ -0.13)
#	[-0.1, 1.0, -0.1],
# Test run
#	[0.0, 1.0, -1.0]
# Second run
#	[0.0, 1.0, -0.3],
#	[0.125, 1.0, -0.3],
#	[0.25, 1.0, -0.3],
#	[0.375, 1.0, -0.3],
#	[0.5, 1.0, -0.4],
#	[0.675, 1.0, -0.4],
#	[0.75, 1.0, -0.4],
#	[-0.1, 1.0, -0.3]
# Third run: line along Jx=0.25 in T-flow
	[0.25, 1.0, 1.0],
	[0.25, 1.0, 0.8],
	[0.25, 1.0, 0.6],
	[0.25, 1.0, 0.3]
]
J_index = parse(Int, ARGS[1])
J = J_array[J_index]

System = getPyrochlore(NLen, [J1, J2])

anisotropy = zeros(System.Npairs, 3)
for n in 1:System.Npairs
    anisotropy[n, :] = J
end

Par = PMFRG_xyz.Params(
    System,
    N=32,
    accuracy=1e-10,
    temp_max=100.0,
    temp_min=0.02
)

mkpath("checkpoints_theta15")
mkpath("data_theta15")
name = "pyro_theta=$(theta_deg)_J=[$(J[1]), $(J[2]), $(J[3])]_NLen=$(NLen)_N=32.jld2"

sol, saved_values = PMFRG_xyz.SolveFRG(
    Par, anisotropy,
    SavedValues(PMFRG_xyz._getFloatType(Par), save_struct{PMFRG_xyz._getFloatType(Par)}),
    (State, t, integrator) -> save_all(State, t, Par, integrator, "checkpoints_theta15/" * name)
)
save_object("data_theta15/" * name, [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
