#=
Analyze 11 years of ERA5 RHI data (2010-2020) from external drive.
Computes fractal exponent β for each daily time step and analyzes seasonal trends.

To run with parallelization on M1 MacBook Air:

    julia -t auto --project=. -e 'include("verlauf.jl")'

Or in the REPL after starting Julia with:

    julia -t auto --project=.
    include("verlauf.jl")

The `-t auto` flag uses all available cores (typically 8 on M1).
=#

using Dates
using DataFrames
using JLD2
using Base.Threads

include("func.jl")
include("setup.jl")

println("Running with $(nthreads()) thread(s)")

## Load existing results or compute from scratch

results_file = "beta_results_2010_2020.jld2"
if isfile(results_file)
    println("Loading previously computed results from $results_file...")
    data = load(results_file)
    all_betas = data["all_betas"]
    all_dates = data["all_dates"]
    all_r2 = data["all_r2"]
    println("✓ Loaded $(length(all_betas)) data points")
else
    ## Process all monthly files from external drive

    data_dir = "/Volumes/ISSR/Data_raw"
    nc_files = filter(f -> endswith(f, ".nc") && startswith(f, "pl"), readdir(data_dir))
    sort!(nc_files)

    all_betas = []
    all_dates = []
    all_r2 = []

    println("Processing $(length(nc_files)) monthly files with $(nthreads()) thread(s)...")

    for nc_file in nc_files
        filepath = joinpath(data_dir, nc_file)
        
        try
            ds = NCDataset(filepath)
            q = ds["q"][461:670, 71:280, :, :]
            t = ds["t"][461:670, 71:280, :, :]
            close(ds)
            
            # Extract year-month from filename
            date_str = nc_file[3:8]
            year = parse(Int, date_str[1:4])
            month = parse(Int, date_str[5:6])
            
            rhi_300 = rhi_calc(q[:, :, 3, :], t[:, :, 3, :], 30000)
            n_times = size(rhi_300, 3)
            
            # Pre-allocate arrays for this month
            month_betas = Vector{Float64}(undef, n_times)
            month_r2 = Vector{Float64}(undef, n_times)
            month_dates = Vector{Date}(undef, n_times)
            
            # Compute β for each daily time step in parallel
            @threads for t_idx in 1:n_times
                field = rhi_300[:, :, t_idx]
                β, _, r2, _, _ = measure(field; window = true)
                month_betas[t_idx] = β
                month_r2[t_idx] = r2
                month_dates[t_idx] = Date(year, month, 1) + Day(t_idx - 1)
            end
            
            # Append month results to global arrays
            append!(all_betas, month_betas)
            append!(all_r2, month_r2)
            append!(all_dates, month_dates)
            
            println("  ✓ $nc_file ($n_times days)")
            
        catch e
            println("  ✗ Error in $nc_file: $(typeof(e).name)")
        end
    end

    ## Save computed results to file (so you don't have to recompute)

    jldsave(results_file; all_betas, all_dates, all_r2)
    println("\n✓ Results saved to $results_file")
end

## Plot β over time (seasonal trends)

if !isempty(all_betas)
    p_time = plot(all_dates, all_betas;
                  xlabel = "Date",
                  ylabel = "β (Fractal Exponent)",
                  title = "Fractal Dimension Over 11 Years - Seasonal Trends",
                  legend = false,
                  markersize = 3,
                  markerstrokewidth = 0,
                  alpha = 0.6)
    display(p_time)
    savefig(p_time, "beta_timeseries.png")
    
    # Monthly statistics plot
    month_groups = [month(d) for d in all_dates]
    monthly_means = Float64[]
    monthly_stds = Float64[]
    
    for m in 1:12
        month_betas = [all_betas[i] for i in 1:length(all_betas) if month_groups[i] == m]
        if !isempty(month_betas)
            push!(monthly_means, mean(month_betas))
            push!(monthly_stds, std(month_betas))
        else
            push!(monthly_means, NaN)
            push!(monthly_stds, NaN)
        end
    end
    
    p_month = plot(1:12, monthly_means;
                   yerror = monthly_stds,
                   xlabel = "Month",
                   ylabel = "β (mean ± std)",
                   title = "Seasonal Distribution of β",
                   legend = false,
                   markersize = 6,
                   markerstrokewidth = 0,
                   linewidth = 2)
    display(p_month)
    savefig(p_month, "beta_seasonal_boxplot.png")
    
    println("\n✓ Plots saved:")
    println("  - beta_timeseries.png")
    println("  - beta_seasonal_boxplot.png")
    
    # Print statistics
    println("\n" * "="^50)
    println("SEASONAL ANALYSIS (2010-2020)")
    println("="^50)
    for m in 1:12
        month_betas = [all_betas[i] for i in 1:length(all_betas) if month_groups[i] == m]
        if !isempty(month_betas)
            month_name = Dates.monthname(m)
            @printf("%2d %-10s: μ = %.3f, σ = %.3f, n = %4d\n", 
                    m, month_name, mean(month_betas), std(month_betas), length(month_betas))
        end
    end
    println("="^50)
end


## fuck

# Running mean with window size = 30 days (approximately 1 month)
window_size = 25

running_mean = [mean(all_betas[max(1, i-window_size):i]) for i in 1:length(all_betas)]

# Plot it
plot(all_dates, running_mean;
     xlabel = "Date",
     ylabel = "β (30-day Running Mean)",
     title = "Seasonal Trends - Running Average",
     legend = false,
     linewidth = 2)