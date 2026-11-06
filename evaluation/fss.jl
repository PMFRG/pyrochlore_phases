# Finite-size scaling of one data point: xi/L along the flow for several NLen, crossing = critical T (or Lambda).
#
# Run by hand:  julia --project=.. -t 4 fss.jl      (edit the settings below first)
# Output: <name>_fss.png, and one row in fss_results.csv

using Pkg
Pkg.activate("../.")

using JLD2, CairoMakie
using PMFRGEvaluation                 # also gives SpinFRGLattices, StaticArrays, maxCorrelationLength, getFourier
using SpinFRGLattices.Pyrochlore

# ---------------- settings ----------------
file    = joinpath(@__DIR__, "..", "data", "pyro_theta=15.0_J=[0.25, 1.0, -0.3]_NLen=6_N=32.jld2")
nlens   = [6, 8, 10]                  # the file name is reused with NLen=6 -> NLen=8 etc.
channel = :z                          # susceptibility chi_x, chi_y or chi_z
Q       = SA[0.0, 0.0, 0.0] * π       # wave vector at which xi is measured ((0,0,0): all-in-all-out order)
# ------------------------------------------

# Layouts of the saved data (JLD2 looks these types up by name)
struct save_struct{T}
    Chi_x::Array{T,2}; Chi_y::Array{T,2}; Chi_z::Array{T,2}
    Sigma_x::Array{T,2}; Sigma_y::Array{T,2}; Sigma_z::Array{T,2}
    f_int::Vector{T}
end

# Read one file -> flow parameters (high to low) and real-space static chi of every step
function load_flow(path, channel)
    steps  = load_object(path) # loads the data
    params = [s[2] for s in steps]
    chis   = map(steps) do s
        chi = getproperty(s[1], Symbol(:Chi_, channel))
        ndims(chi) == 2 ? chi[:, 1] : chi
    end
    params[1] < params[end] && (reverse!(params); reverse!(chis))
    return params, chis
end

# xi/L at every step: second-moment correlation length (maxCorrelationLength) divided by NLen
function xi_over_L(chis, NLen, Q)
    lattice = LatticeInfo(getPyrochlore(NLen), Pyrochlore)
    xi = zeros(length(chis))
    Threads.@threads for i in eachindex(chis)
        xi[i] = maxCorrelationLength(Q, getFourier(chis[i], lattice), NLen) / NLen
    end
    return xi
end

# Linear interpolation of y(x) at xnew, using log10(x); x is ordered high to low
function interp(x, y, xnew)
    lx = reverse(log10.(x)); ly = reverse(y)
    return [begin
                j = clamp(searchsortedlast(lx, log10(xn)), 1, length(lx) - 1)
                t = (log10(xn) - lx[j]) / (lx[j+1] - lx[j])
                (1 - t) * ly[j] + t * ly[j+1]
            end for xn in xnew]
end

# First crossing (highest T) of two xi/L curves given on the same grid x (high to low).
# Differences below `tol` count as "equal", so the noisy start of the flow cannot fake a crossing.
function first_crossing(x, y1, y2; tol = 1e-3)
    d = y1 .- y2
    last = 0                                            # last index where the curves were clearly apart
    for i in eachindex(d)
        abs(d[i]) <= tol * max(maximum(abs, y1), maximum(abs, y2)) && continue
        if last > 0 && sign(d[i]) != sign(d[last])
            for j in last+1:i                           # locate the sign change linearly
                if d[j-1] * d[j] <= 0 && d[j-1] != d[j]
                    t = d[j-1] / (d[j-1] - d[j])
                    return 10^((1 - t) * log10(x[j-1]) + t * log10(x[j])), (1 - t) * y1[j-1] + t * y1[j]
                end
            end
        end
        last = i
    end
    return nothing
end

# ---------------- 1. load and compute xi/L for every NLen ----------------
data = Dict{Int,Any}()
for n in nlens
    path = replace(file, r"NLen=\d+" => "NLen=$n")
    isfile(path) || (println("missing: ", path); continue)
    params, chis = load_flow(path, channel)
    data[n] = (params = params, xiL = xi_over_L(chis, n, Q))
end
ns = sort(collect(keys(data)))
length(ns) >= 2 || error("need at least two NLen")

# ---------------- 2. crossings of all NLen pairs, on the mesh of the smallest NLen ----------------
lo = maximum(minimum(data[n].params) for n in ns)       # window where all NLen have data
hi = minimum(maximum(data[n].params) for n in ns)
grid = filter(p -> lo <= p <= hi, data[ns[1]].params)
curves = Dict(n => interp(data[n].params, data[n].xiL, grid) for n in ns)

crossings = [first_crossing(grid, curves[a], curves[b]) for a in ns for b in ns if a < b]
crossings = filter(!isnothing, crossings)
Tc  = isempty(crossings) ? NaN : 10^(sum(log10(c[1]) for c in crossings) / length(crossings))
xiL = isempty(crossings) ? NaN : sum(c[2] for c in crossings) / length(crossings)
println(isempty(crossings) ? "no crossing -> disordered" : "critical parameter = $Tc   (xi/L = $xiL)")

# ---------------- 3. plot ----------------
name = replace(basename(file), r"_NLen=\d+" => "", ".jld2" => "")
flowname = startswith(basename(file), "Lpyro") ? "Λ" : "T"
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; xlabel = flowname, ylabel = "ξ / L", xscale = log10, xreversed = true,
          title = "$name\nχ_$channel at Q = $(round.(Q ./ π; digits = 3))π")
for n in ns
    lines!(ax, data[n].params, data[n].xiL; label = "NLen = $n")
end
isnan(Tc) || vlines!(ax, Tc; color = :black, linestyle = :dash)
axislegend(ax; position = :lt)
save(joinpath(@__DIR__, "$(name)_fss.png"), fig)

# ---------------- 4. csv ----------------
# phase: chi_z peak at Q = 0 -> AIAO Z, chi_x peak at Q = 0 -> AIAO X, other crossing -> ordered
phase = isnan(Tc) ? "disordered" : (all(iszero, Q) && channel == :z) ? "AIAO Z" :
        (all(iszero, Q) && channel == :x) ? "AIAO X" : "ordered"
J = parse.(Float64, split(match(r"J=\[([^\]]*)\]", file)[1], ","))
csv = joinpath(@__DIR__, "fss_results.csv")
header = "name,channel,Jx,Jy,Jz,phase,Tc,xiL_at_crossing"
key = replace(name, "," => ";") * "," * string(channel) * ","
row = key * join([J..., phase, Tc, xiL], ",")
old = isfile(csv) ? readlines(csv)[2:end] : String[]
old = filter(l -> !startswith(l, key), old)             # replace an older result of the same set and channel
write(csv, join([header; old; row], "\n") * "\n")
println("phase: $phase")