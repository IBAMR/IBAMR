#!/usr/bin/env julia

using LinearAlgebra
using Printf
using ImplicitIBCAVReference

function usage()
    println("usage: run_ex0_gmres_convergence_julia.jl [options]")
    println("")
    println("options:")
    println("  --summary-path PATH   Output TSV summary path")
    println("  --residual-dir DIR    Output residual-history directory")
    println("  --depths CSV          Depth list, e.g. 2,3,4")
    println("  --k-values CSV        Stiffness list, e.g. 1e2,1e4")
    println("  --coarse-n N          Coarsest grid cells (default: 4)")
    println("  --nu1 N               Number of pre-sweeps (default: 1)")
    println("  --nu2 N               Number of post-sweeps (default: 1)")
end

function parse_csv_int(raw::String)
    vals = Int[]
    for token in split(raw, ',')
        t = strip(token)
        isempty(t) && continue
        push!(vals, parse(Int, t))
    end
    return vals
end

function parse_csv_tokens(raw::String)
    vals = String[]
    for token in split(raw, ',')
        t = strip(token)
        isempty(t) && continue
        push!(vals, t)
    end
    return vals
end

function parse_args(args)
    opts = Dict{String, String}(
        "summary_path" => "/tmp/julia_summary.tsv",
        "residual_dir" => "/tmp/julia_residual_histories",
        "depths" => "2,3,4",
        "k_values" => "1e2,1e4,1e6,1e8",
        "coarse_n" => "4",
        "nu1" => "1",
        "nu2" => "1",
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--summary-path"
            i += 1
            i > length(args) && error("missing value for --summary-path")
            opts["summary_path"] = args[i]
        elseif arg == "--residual-dir"
            i += 1
            i > length(args) && error("missing value for --residual-dir")
            opts["residual_dir"] = args[i]
        elseif arg == "--depths"
            i += 1
            i > length(args) && error("missing value for --depths")
            opts["depths"] = args[i]
        elseif arg == "--k-values"
            i += 1
            i > length(args) && error("missing value for --k-values")
            opts["k_values"] = args[i]
        elseif arg == "--coarse-n"
            i += 1
            i > length(args) && error("missing value for --coarse-n")
            opts["coarse_n"] = args[i]
        elseif arg == "--nu1"
            i += 1
            i > length(args) && error("missing value for --nu1")
            opts["nu1"] = args[i]
        elseif arg == "--nu2"
            i += 1
            i > length(args) && error("missing value for --nu2")
            opts["nu2"] = args[i]
        elseif arg == "--help" || arg == "-h"
            usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return opts
end

function write_history(path::String, resvec::Vector{Float64})
    open(path, "w") do io
        for (idx, value) in enumerate(resvec)
            @printf(io, "%d\t%.17e\n", idx - 1, value)
        end
    end
end

function run()
    opts = parse_args(ARGS)

    summary_path = opts["summary_path"]
    residual_dir = opts["residual_dir"]
    depths = parse_csv_int(opts["depths"])
    k_tokens = parse_csv_tokens(opts["k_values"])
    coarse_n = parse(Int, opts["coarse_n"])
    nu1 = parse(Int, opts["nu1"])
    nu2 = parse(Int, opts["nu2"])

    isempty(depths) && error("no depths provided")
    isempty(k_tokens) && error("no k-values provided")

    mkpath(dirname(summary_path))
    mkpath(residual_dir)

    open(summary_path, "w") do io
        println(io,
                "case_id\tdepth\tfinest_n\tcoarse_n\tK\tjulia_converged\tjulia_iterations\tjulia_relres\tjulia_rhs_norm\tjulia_dt\tjulia_rho\tjulia_mu\tjulia_lag_ds\tjulia_num_curve_points\tjulia_nu1\tjulia_nu2\tjulia_history_file")

        for depth in depths
            depth < 1 && error("depth must be positive")
            finest_n = coarse_n * (2^(depth - 1))
            for k_token in k_tokens
                k_value = parse(Float64, k_token)
                case_id = @sprintf("depth%d_N%d_K%s", depth, finest_n, k_token)
                result = solve_smoke_example(N = finest_n,
                                             depth = depth,
                                             K = k_value,
                                             nu1 = nu1,
                                             nu2 = nu2,
                                             tol = 1.0e-8,
                                             maxiter = 150)

                rhs_norm = norm(result.problem.rhs)
                dt = result.problem.geometry.dt
                rho = result.problem.rho
                mu = result.problem.mu
                lag_ds = result.problem.geometry.ds
                num_curve_points = result.problem.geometry.approx
                history_file = joinpath(residual_dir, case_id * ".tsv")
                write_history(history_file, result.gmres.resvec)

                converged = result.gmres.converged ? 1 : 0
                @printf(io,
                        "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%d\t%d\t%d\t%s\n",
                        case_id,
                        depth,
                        finest_n,
                        coarse_n,
                        k_token,
                        converged,
                        result.gmres.iterations,
                        result.gmres.relres,
                        rhs_norm,
                        dt,
                        rho,
                        mu,
                        lag_ds,
                        num_curve_points,
                        nu1,
                        nu2,
                        history_file)
            end
        end
    end

    println("julia_summary=" * summary_path)
end

run()
