#!/usr/bin/env julia
include("./filename.jl")
using ArgParse

function main(dry_run=true)
    NLens = (6,8,10)
    N=32

    # For N == 32, scales as N^4
    # This could be also a function of N and NLen instead of a Dictionary
    walltime_strings = Dict(32 => Dict(6 => "1-00:00:00",
                        8 => "3-00:00:00",
                        10 => "7-00:00:00"))

    # For N == 32, scales as N^3 
    # This could be also a function of N and NLen instead of a Dictionary
    mem_per_cpu_strings = Dict(32 => Dict( 6 => "700MB",
                              8 => "1000MB",
                              10 => "1500MB",
                             ))
    
    # List to populate (or to read from file)
    Js = [  (0.25,1.0,0.3), ]
   
    theta_degs = [ 15.0, 5.749 ] 

    for theta_deg in theta_degs
        for NLen in NLens
            for J in Js
                let (Jx, Jy, Jz) = J

                    run_if_results_not_exist(NLen,Jx,Jy,Jz,theta_deg,N,
                                         walltime=walltime_strings[N][NLen],
                                         mem_per_cpu=mem_per_cpu_strings[N][NLen],
                                         dry_run
                                    )
                end
            end
        end
    end
end

function run_if_results_not_exist(NLen,Jx,Jy,Jz,theta_deg,N, dry_run; walltime, mem_per_cpu)

    J = [Jx,Jy,Jz]
    name = filename(theta_deg,J,NLen,N)

    command_array = if dry_run ; ["echo"] else [] end

    if ! isfile("data/$name")
        command_string = "sbatch --time $walltime --mem $mem_per_cpu pyro_J2.sh --nlen $NLen --Jx $(J[1]) --Jy $(J[2]) --Jz $(J[3]) --theta-deg $theta_deg --N $N"

        command_array = vcat(command_array, split(command_string))
        run(`$command_array`)
    else 
        println("File data/$name exists, not running again.")
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dry-run"
        action = :store_true
    end

    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()
    main(args["dry_run"])
end 

