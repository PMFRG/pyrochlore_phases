# Structure factor chi(q) of one data point at one temperature (or Lambda), in a plane of reciprocal space.
#
# Run by hand:  julia --project=.. structure_factor.jl      (edit the settings below first)
# Output: <name>_chi<channel>.png

using JLD2, CairoMakie
using PMFRGEvaluation                 # also gives SpinFRGLattices, getFourier
using SpinFRGLattices.Pyrochlore

# ---------------- settings ----------------
file    = raw"..\data_J2_temperature\data_theta15_phase_diagram\pyro_theta=15.0_J=[0.0, 1.0, -0.3]_NLen=8_N=32.jld2"
NLen    = 8
channel = :z
T       = 0.4                         # flow value to show; the nearest saved step is used
plane   = (x, y) -> (x, x, y)         # (kx, ky, kz) as a function of the plot axes; here the (h,h,l) plane
                                      # (x, y) -> (x, y, 0) gives the (h,k,0) plane
extent  = 6π
# ------------------------------------------

struct save_struct{T}                 # layouts of the saved data (JLD2 looks these types up by name)
    Chi_x::Array{T,2}; Chi_y::Array{T,2}; Chi_z::Array{T,2}
    Sigma_x::Array{T,2}; Sigma_y::Array{T,2}; Sigma_z::Array{T,2}
    f_int::Vector{T}
end
struct Observables{T}
    Chi_x::Vector{T}; Chi_y::Vector{T}; Chi_z::Vector{T}
end

# 1. read the file, pick the step closest to T
steps  = load_object(file)
params = [s[2] for s in steps]
i      = argmin(abs.(log.(params) .- log(T)))
chi    = getproperty(steps[i][1], Symbol(:Chi_, channel))
chi_R  = ndims(chi) == 2 ? chi[:, 1] : chi          # static real-space chi (one value per inequivalent pair)

# 2. Fourier transform: chi(q) as a function
lattice = LatticeInfo(getPyrochlore(NLen), Pyrochlore)
chi_q   = getFourier(chi_R, lattice)

# 3. evaluate on a grid in the chosen plane and plot
k      = range(-extent, extent, length = 200)
values = [chi_q(plane(x, y)...) for x in k, y in k]
n      = floor(Int, extent / π)
ticks  = (collect(-n:n) .* π, string.(-n:n, "π"))

fig = Figure(size = (650, 550))
ax = Axis(fig[1, 1]; aspect = 1, xticks = ticks, yticks = ticks,
          title = "χ_$channel(q) at T = $(round(params[i]; sigdigits = 3)), NLen = $NLen")
hm = heatmap!(ax, k, k, values; colormap = :magma)
Colorbar(fig[1, 2], hm)
name = replace(basename(file), ".jld2" => "")
save("$(name)_chi$(channel).png", fig)
println("saved $(name)_chi$(channel).png   (max chi in plane = $(maximum(values)))")
