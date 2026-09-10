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

results_file = "beta_results_monthly_2010_2020.jld2"
if isfile(results_file)
    println("Loading previously computed monthly results from $results_file...")
    data = load(results_file)
    monthly_betas = data["monthly_betas"]
    monthly_dates = data["monthly_dates"]
    monthly_r2 = data["monthly_r2"]
    println("✓ Loaded $(length(monthly_betas)) monthly data points")
else
    ## Process all monthly files from external drive

    data_dir = "/Volumes/ISSR/Data_raw"
    nc_files = filter(f -> endswith(f, ".nc") && startswith(f, "pl"), readdir(data_dir))
    sort!(nc_files)

    monthly_betas = Float64[]
    monthly_dates = Date[]
    monthly_r2 = Float64[]

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
            
            # Pre-allocate arrays for this month's daily spectra and get k from first field
            month_spectra = Matrix{Float64}(undef, 105, n_times)
            field_0 = rhi_300[:, :, 1]
            _, _, _, k_month, S_0 = measure(field_0; window = true)
            month_spectra[:, 1] = S_0
            
            # Compute spectra for remaining daily time steps in parallel
            @threads for t_idx in 2:n_times
                field = rhi_300[:, :, t_idx]
                _, _, _, _, S = measure(field; window = true)
                month_spectra[:, t_idx] = S
            end
            
            # Average spectra for this month
            S_month_mean = vec(mean(month_spectra; dims = 2))
            
            # Fit β to the monthly-averaged spectrum
            β_month, _, r2_month = fit_beta(k_month, S_month_mean)
            
            push!(monthly_betas, β_month)
            push!(monthly_r2, r2_month)
            push!(monthly_dates, Date(year, month, 1))
            
            println("  ✓ $nc_file: β = $(round(β_month; digits=3)), R² = $(round(r2_month; digits=4))")
            
        catch e
            println("  ✗ Error in $nc_file: $(typeof(e).name)")
        end
    end

    ## Save computed monthly results to file

    jldsave(results_file; monthly_betas, monthly_dates, monthly_r2)
    println("\n✓ Monthly results saved to $results_file")
    println("✓ $(length(monthly_betas)) months averaged and fitted")
end

## Plot β over time (monthly averages)

jahre = Date(year(first(monthly_dates)), 1, 1):Year(1):Date(year(last(monthly_dates)), 1, 1)


if !isempty(monthly_betas)
    p_time = plot(monthly_dates, monthly_betas;
                  xlabel = "Date",
                  ylabel = "β (Fractal Exponent)",
                  title = "Monthly-Averaged Fractal Dimension (2010-2020)",
                  legend = false,
                  marker = :circle,
                  markersize = 5,
                  markerstrokewidth = 0,
                  markerstrokecolor = :false,
                  linewidth = 2,
                  xticks = (jahre, Dates.year.(jahre)),
                  size = (1200, 400))
    display(p_time)
    savefig(p_time, "beta_monthly_timeseries.png")
    
    # Seasonal plot by month
    month_groups = [month(d) for d in monthly_dates]
    seasonal_means = Float64[]
    seasonal_stds = Float64[]
    seasonal_n = Int[]
    
    for m in 1:12
        m_betas = [monthly_betas[i] for i in 1:length(monthly_betas) if month_groups[i] == m]
        if !isempty(m_betas)
            push!(seasonal_means, mean(m_betas))
            push!(seasonal_stds, std(m_betas))
            push!(seasonal_n, length(m_betas))
        else
            push!(seasonal_means, NaN)
            push!(seasonal_stds, NaN)
            push!(seasonal_n, 0)
        end
    end
    
    p_season = plot(1:12, seasonal_means;
                    yerror = seasonal_stds,
                    xlabel = "Month",
                    ylabel = "β (mean ± std)",
                    title = "Seasonal Pattern - Monthly Spectral Averages",
                    legend = false,
                    markersize = 6,
                    markerstrokewidth = 0,
                    linewidth = 2,
                    xticks = 1:12)
    display(p_season)
    savefig(p_season, "beta_seasonal_pattern.png")
    
    println("\n✓ Plots saved:")
    println("  - beta_monthly_timeseries.png")
    println("  - beta_seasonal_pattern.png")
    
    # Print statistics
    println("\n" * "="^50)
    println("MONTHLY-AVERAGED ANALYSIS (2010-2020)")
    println("="^50)
    for m in 1:12
        idx_for_month = findall(x -> x == m, month_groups)
        if !isempty(idx_for_month)
            m_betas = monthly_betas[idx_for_month]
            m_r2 = monthly_r2[idx_for_month]
            month_name = Dates.monthname(m)
            @printf("%2d %-10s: μ = %.3f, σ = %.3f, n = %2d months, mean R² = %.4f\n", 
                    m, month_name, mean(m_betas), std(m_betas), length(m_betas), mean(m_r2))
        end
    end
    println("="^50)
    println("\nNote: Each β is fitted to the monthly-averaged spectrum (proper method)")
end

## lol

mask = monthly_dates .>= Date(2013, 1, 1)
monthly_dates_filtered = monthly_dates[mask]
monthly_betas_filtered = monthly_betas[mask]

jahre = Date(2013, 1, 1):Year(1):Date(year(last(monthly_dates_filtered)), 1, 1)

p_time = plot(monthly_dates_filtered, monthly_betas_filtered;
                xlabel = "Date",
                ylabel = "β (Fractal Exponent)",
                title = "Monthly-Averaged Fractal Dimension (2010-2020)",
                legend = false,
                marker = :circle,
                markersize = 5,
                markerstrokewidth = 0,
                markerstrokecolor = :false,
                linewidth = 2,
                xticks = (jahre, Dates.year.(jahre)),
                size = (1200, 400))
display(p_time)
savefig(p_time, "beta_monthly_timeseries2013.png")
