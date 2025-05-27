#####################################################################
#               Production Function Estimation Project              #
#                   Markus Trunschke and Ken Judd                   #
#                      - Estimation Functions -                     #
#####################################################################
# First version: 05/07/2024                                         #
# Author: Markus Trunschke                                          #
#####################################################################
# Notes: - This file contains the all functions for estimating
#          production functions with the control function approach  #
#          developed in Ackerberg, Caves, Frazer (2015). Their      #
#          approach is based on Olley and Pakes (1996), and         #
#          Levinsohn and Petrin (2003)                              #
#####################################################################

## Wrapper for acfprodest! that copies data
function acfprodest(; data::DataFrame,
                      time::Symbol,
                      id::Symbol,
                      output::Symbol,
                      flexible_inputs::Union{Symbol,Vector{Symbol}},
                      fixed_inputs::Union{Symbol,Vector{Symbol}},
                      proxy_variable::Symbol,
                      fes_poly_degree::Int = 3,
                      ω_lom_degree::Int = 3,
                      opts::Dict = Dict()
                      )

    # Copy data
    est_data = copy(data)

    # Run real acfprodest command
    acf_returns = acfprodest!(data = est_data,
                              time = time,
                              id = id,
                              output = output,
                              flexible_inputs = flexible_inputs,
                              fixed_inputs = fixed_inputs,
                              proxy_variable = proxy_variable,
                              fes_poly_degree = fes_poly_degree,
                              ω_lom_degree = ω_lom_degree,
                              opts = opts
                              )

    # Return results
    return acf_returns
end

## Wrapper of the estimation procedure
function acfprodest!(;data::DataFrame,
    time::Symbol,
    id::Symbol,
    output::Symbol,
    flexible_inputs::Union{Symbol,Vector{Symbol}},
    fixed_inputs::Union{Symbol,Vector{Symbol}},
    proxy_variable::Symbol,
    poly_names::Vector{Symbol} = [:NotDefinedByUser],
    lag_poly_names::Vector{Symbol} = [:NotDefinedByUser],
    fes_poly_degree::Int = 3,
    ω_lom_degree::Int = 3,
    opts::Dict = Dict()
    )
    
    # User input cleaning
    fixed_inputs, flexible_inputs, opts = acf_user_input_cleaner!(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, opts = opts)

    # Check for syntax errors
    # error_check_fnc(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, opts = opts)

    # Prepare data
    data, all_input_symbols, poly_names, lag_poly_names = data_prep!(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree)

    # Run point estimation
    fes_returns, ses_returns = acfprodest_estimation(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, opts = opts)
    
    # Run SE calculation
    res_tab = acf_SE_stats(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, point_est = ses_returns, opts = opts)

    # Print results summary
    acf_print_all_res(data = data, results_df = res_tab, id = id, opts = opts)

    # Return stuff
    return fes_returns, ses_returns, res_tab
end

## Main production function estimation function 
function acfprodest_estimation(;data::DataFrame,
                                time::Symbol,
                                id::Symbol,
                                output::Symbol,
                                flexible_inputs::Vector{Symbol},
                                fixed_inputs::Vector{Symbol},
                                proxy_variable::Symbol,
                                poly_names::Vector{Symbol},
                                lag_poly_names::Vector{Symbol},
                                fes_poly_degree::Int,
                                ω_lom_degree::Int,
                                opts::Dict = Dict()
                               )
                               
    ################################# First Stage ##############################
    # Regress output on inputs and non-parametric approximation of omega using
    # some demand function as the proxy
    fes_returns = acfprodest_fes!(data = data, id = id, time = time, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree, opts = opts)

    ################################# Second Stage #############################
    # Use estimated phihats and estimate the GMM part of the model
    ses_returns =  acfprodest_ses!(data = data, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, ω_lom_degree = ω_lom_degree, opts = opts)
    
    # Return results
    if opts["call_from_boot_runs"] == true
        return ses_returns
    else
        return fes_returns , ses_returns
    end
end

## Wrapper of first stage of the estimation procedure
function acfprodest_fes(;data::DataFrame,
    time::Symbol,
    id::Symbol,
    output::Symbol,
    flexible_inputs::Union{Symbol,Vector{Symbol}},
    fixed_inputs::Union{Symbol,Vector{Symbol}},
    proxy_variable::Symbol,
    poly_names::Vector{Symbol} = [:NotDefinedByUser],
    lag_poly_names::Vector{Symbol} = [:NotDefinedByUser],
    fes_poly_degree::Int = 3,
    ω_lom_degree::Int = 3,
    opts::Dict = Dict()
   )

    # Copy data
    est_data = copy(data)

    # Actual first stage estimation
    fes_returns = acfprodest_fes!(data = est_data,
                                  time = time,
                                  id = id,
                                  output = output,
                                  flexible_inputs = flexible_inputs,
                                  fixed_inputs = fixed_inputs,
                                  proxy_variable = proxy_variable,
                                  poly_names = poly_names,
                                  lag_poly_names = lag_poly_names,
                                  fes_poly_degree = fes_poly_degree,
                                  ω_lom_degree = ω_lom_degree,
                                  opts = opts
                                 )

    return fes_returns, est_data
end

## First stage of the estimation procedure
function acfprodest_fes!(;data::DataFrame,
                         time::Symbol,
                         id::Symbol,
                         output::Symbol,
                         flexible_inputs::Union{Symbol,Vector{Symbol}},
                         fixed_inputs::Union{Symbol,Vector{Symbol}},
                         proxy_variable::Symbol,
                         poly_names::Vector{Symbol} = [:NotDefinedByUser],
                         lag_poly_names::Vector{Symbol} = [:NotDefinedByUser],
                         fes_poly_degree::Int = 3,
                         ω_lom_degree::Int = 3,
                         opts::Dict = Dict()
                        )

    # Calculate polynomes of inputs if not already done before (e.g. if user calls this function by iteself)
    if poly_names == [:NotDefinedByUser]
        poly_names = polynom_series!(;data = data, var_names = [flexible_inputs..., fixed_inputs..., proxy_variable], degree = fes_poly_degree)
        panel_lag!(data = data, id = id, time = time, variable = poly_names, lag_prefix = "lag_", lags = 1, drop_missings = false, force = false)
        lag_poly_names = Symbol.("lag_" .* string.(poly_names))
    end

    # Calcualte first stage parameters (aka run OLS)
    α = fastOLS(Y = data[!, output], X = Matrix{Float64}(data[!, vcat(:constant, poly_names)]))
    
    # Calculate some predictions
    data.Φ_hat .= vec(Matrix(data[!, vcat(:constant, poly_names)])*α)
    data.lag_Φ_hat .= vec(Matrix(data[!, vcat(:constant, lag_poly_names)])*α)

    # Print results
    if opts["fes_print_results"] == true && opts["call_from_boot_runs"] != true
        print_fes_res(α, vcat(:constant, poly_names))
    end

    # Return data including estimates of phi for second stage
    return α
end

## Wrapper of second stage of the estimation procedure
function acfprodest_ses(;data::DataFrame, 
                        flexible_inputs::Vector{Symbol},
                        fixed_inputs::Vector{Symbol},
                        ω_lom_degree::Int = 3,
                        opts::Dict = opts
                        )

    # copy data
    est_data = copy(data)

    # Run real estimation command
    ses_returns = acfprodest_ses!(data = est_data, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, ω_lom_degree = ω_lom_degree, opts = opts)

    return ses_returns 
end

# ## Second stage of the estimation procedure
function acfprodest_ses!(;data::DataFrame, 
    flexible_inputs::Union{Symbol,Vector{Symbol}},
    fixed_inputs::Union{Symbol,Vector{Symbol}},
    ω_lom_degree::Int = 3,
    opts::Dict
    )

    # Define weights as an identity matrix. Inefficient GMM estimator has arbitrary weights
    weights = I

    # Optimize
    opt = optimize(par -> acfprodest_ses_criterion(data = data, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, 
    β = par, weight = weights, ω_lom_degree = ω_lom_degree), opts["lower_bound"], opts["upper_bound"], opts["startvals"], opts["Optimizer"], opts["Optimizer_Options"])

    # Throw an error message if GMM did not converge
    con_flag = Optim.converged(opt)
    if con_flag == false
        error("GMM did not converge")
    else
        β_hat = Optim.minimizer(opt)
    end

    # Calculate ω lom parameters
    ρ_hat = prod_reg!(data = data, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, ω_lom_degree = ω_lom_degree, β = β_hat)

    # Print results
    if opts["ses_print_results"] == true && opts["call_from_boot_runs"] != true
        print_res_ses(β_hat = β_hat, ρ_hat = ρ_hat, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs)
    end

    return [β_hat..., ρ_hat...]
end

## Second stage GMM criterion function
function acfprodest_ses_criterion(;data::DataFrame, 
                                   flexible_inputs::Union{Symbol,Vector{Symbol}}, 
                                   fixed_inputs::Union{Symbol,Vector{Symbol}}, 
                                   β::Vector{<:Number}, 
                                   weight::Union{Array{<:Number},LinearAlgebra.UniformScaling{Bool}},
                                   ω_lom_degree::Int
                                   )
                                   
    # Get ξ and ρ (outsourced in another program b/c it is used in multiple locations)
    prod_reg!(data = data, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, ω_lom_degree = ω_lom_degree, β = β)

    #### Calculate moments
    # First moment condition: E[xi(b0,b1,b2,b3,b4)] = 0
    # (This is not mentioned in the papers but it does not make a difference leaving it out
    # Laura Grigolon suggested in her class a different moment condition: E[omega(b0,b1,b2,b3,b4)] = 0
    # and it allows to get the constant beta_0)

    # My moments work better if you want to get β₀ but ACF define them differently
    # G = vec((1/size(data)[1]).*(data.ξ'*hcat(Array(select(data, "lag_".*string.(flexible_inputs))), Array(select(data, fixed_inputs)))))

    # # Calculate criterion
    # criterion = G'*weight*G

    # ACF define their criterion function as follows (see replication code)
    # /* Moment condition - analogue of middle two moments of equation (28) from text - note that instruments are current capital and lagged labor*/
    if size(flexible_inputs) != (0,) 
        inputs_data = Array(data[!, [Symbol.(:lag_, flexible_inputs...), fixed_inputs...]])
    else
        inputs_data = Array(data[!, fixed_inputs])
    end
    
    m = data.ξ_hat .* inputs_data

    return (size(m)[1] * mean(m, dims = 1) * inv(cov(m)) * mean(m, dims = 1)')[1]
end

## OLS regression to get ω development process
function prod_reg!(;data::DataFrame, flexible_inputs::Union{Symbol,Vector{Symbol}}, fixed_inputs::Union{Symbol,Vector{Symbol}}, β::Vector{<:Number}, ω_lom_degree::Int)

    # Get ω + ln(β_0)
    data.ω_ln_b0_hat     = data.Φ_hat .- Array(data[!, [flexible_inputs..., fixed_inputs...]]) * β
    data.lag_ω_ln_b0_hat = data.lag_Φ_hat .- Array(data[!, Symbol.(:lag_, [flexible_inputs..., fixed_inputs...])]) * β

    # Get ξ by regressing omega onto omega_lag, digi and no_digi
    # Calculate productivity development parameters
    ρ_hat = fastOLS(Y = data.ω_ln_b0_hat, X = Array(data[!, [:constant, :lag_ω_ln_b0_hat]]))

    # Define residual being the difference b/w omega, lagged omega and the innovation decisions
    data.ξ_hat .= data.ω_ln_b0_hat .- data.constant .* ρ_hat[1] .- data.lag_ω_ln_b0_hat .* ρ_hat[2]

    return ρ_hat
end

## Bootstrap function
function bootstrap_acf_prodest(;data::DataFrame,
                                time::Symbol,
                                id::Symbol,
                                output::Symbol,
                                flexible_inputs::Vector{Symbol},
                                fixed_inputs::Vector{Symbol},
                                proxy_variable::Symbol,
                                poly_names::Vector{Symbol},
                                lag_poly_names::Vector{Symbol},
                                fes_poly_degree::Int,
                                ω_lom_degree::Int,
                                point_est::Array{<:Number},
                                opts::Dict)
    
    # Set some values
    rep = 1
    ntry = 1

    # Preallocate cache
    res_vec = Array{Float64}(undef, opts["boot_reps"], length(point_est))

    # Run bootstrap iterations
    Threads.@threads for rep = 1:opts["boot_reps"]

        successfull_rep = false # Define exit flag

        while successfull_rep == false

            # Print status
            print(".")

            # Draw bootstrap sample
            boot_data = draw_sample(data = data, id = id)

            # Estimate model parameters with bootstrap sample
            try
                res_vec[rep,:] .= acfprodest_estimation(data = boot_data, time = time, id = :unique_id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, opts = opts)
            catch
                if ntry < opts["maxboottries"]
                    ntry += 1
                    continue # Continue without increasing counter to restart iteration
                else
                    throw("Error: Maximal number of bootstrap iteration tries reached at bootstrap repetition "*string(rep))
                end
            end
            successfull_rep = true # Set exit flag to true
        end
    end

    return res_vec
end

## Function calculating some SE related statistics
function acf_SE_stats(;data::DataFrame,
                   time::Symbol,
                   id::Symbol,
                   output::Symbol,
                   flexible_inputs::Vector{Symbol},
                   fixed_inputs::Vector{Symbol},
                   proxy_variable::Symbol,
                   poly_names::Vector{Symbol},
                   lag_poly_names::Vector{Symbol},
                   fes_poly_degree::Int,
                   ω_lom_degree::Int,
                   point_est::Array{<:Number},
                   opts::Dict)
    if opts["no_SE_est"] == false
        # Run Bootstrap repetitions
        opts["call_from_boot_runs"] = true

        res_vec = bootstrap_acf_prodest(data = data, time = time, id = id, output = output, flexible_inputs = flexible_inputs, fixed_inputs = fixed_inputs, proxy_variable = proxy_variable, poly_names = poly_names, lag_poly_names = lag_poly_names, fes_poly_degree = fes_poly_degree, ω_lom_degree = ω_lom_degree, point_est = point_est, opts = opts)
        
        opts["call_from_boot_runs"] = false

        # Calculate variances
        var_vec = vec(var(res_vec, dims = 1))

        # Calculate some quantities
        se_vec  = sqrt.(var_vec)
    
        ## Inference (see Cameron and Trivado (2005) ch.11)
        # Calculate t-statistic
        tstat_vec = point_est ./ se_vec
    
        # Confidence intervals (using N(0,1))
        conf_low_vec  = point_est - 1.96.*vec(se_vec)
        conf_high_vec = point_est + 1.96.*vec(se_vec)
    
        # P-values (using N(0,1))
        pvalues_vec = 2 .*(1 .- cdf.(Normal(),abs.(tstat_vec)))

        # Generate output table
        # Generate symbols for ω lom
        lom_sym_vec = ["ln(β₀)", "ω"]

        if ω_lom_degree > 1
            for j = 2:ω_lom_degree
                lom_sym_vec = hcat(lom_sym_vec, "ω"*string(superscript_this!(string(j))))
            end
        end

        # Combine all
        all_var_vec = [flexible_inputs..., fixed_inputs...,lom_sym_vec...]

        return DataFrame(Variable = vec(all_var_vec), Coefs = point_est,
                     SE = vec(se_vec), t_stat = vec(tstat_vec),
                     P_Value = vec(pvalues_vec), Conf_low = vec(conf_low_vec),
                     Conf_high = vec(conf_high_vec))
    else
        # Generate output table
        # Generate symbols for ω lom
        lom_sym_vec = ["ln(β₀)", "ω"]

        if ω_lom_degree > 1
            for j = 2:ω_lom_degree
                lom_sym_vec = hcat(lom_sym_vec, "ω"*string(superscript_this!(string(j))))
            end
        end

        # Combine all
        all_var_vec = [flexible_inputs..., fixed_inputs...,lom_sym_vec...]

        return DataFrame(Variable = vec(all_var_vec), Coefs = point_est)
    end
end

## Error check function
function error_check_fnc(; kwargs...)

    # Get a vector with elements of all input names
    all_inputs = vcat(flexible_inputs, fixed_inputs, proxy_variable)

    # Check validity of starting values
    if prod_fnc_constant == true
        if startvals != [] && size(startvals,1) != (size(vcat(flexible_inputs, fixed_inputs),1) + 1)
            error("Not enough or too many starting values!")
        end
    elseif prod_fnc_constant == false
        if startvals != [] && size(startvals,1) != size(vcat(flexible_inputs, fixed_inputs),1)
            error("Not enough or too many starting values!")
        end
    end

    # Check if order of taylor approximation is valid
    if taylor_order < 1 || taylor_order > 4
        error("Taylor approximation order in the first stage is not valid. The order can be between 1 and 4")
    end

    # Check that one variable is not entered in multiple spots (e.g. in state and free vars)
    if size(unique(all_inputs),1) != size(all_inputs,1)
        error("Variable entries are not unique: Each variable can only entered in one spot!")
    end

    ## check if variables are in dataset
    for var in all_inputs
        if var in names(data) == false # Check if column name for l_α already exists and drop it if true
            error("Variable: " * string(var) * " is not present in DataFrame")
    
        end
    end
    
    # Check if one variable inputs are only one variable and if they are present in dataset
    if issubset([revenue],names(data)) == false
        if typeof(revenue) != String
            error("revenue variable needs to be a string containing just one variable name!")
        else
            error("The revenue variable is not present in the specified dataframe!")
        end
    end
    if issubset([time],names(data)) == false && time_dummies == 1
        if typeof(time) != String
            error("The time variable needs to be a string containing just one variable name!")
        else
            error("The time variable is not present in the specified dataframe!")
        end
    end
    if issubset([industries],names(data)) == false && industry_dummies == 1
        if typeof(time) != String
            error("The industries variable needs to be a string containing just one variable name!")
        else
            error("The industries variable is not present in the specified dataframe!")
        end
    end
end

## Basic data preparation function
function data_prep!(; data::DataFrame, output::Symbol, flexible_inputs::Union{Symbol,Array{Symbol}}, fixed_inputs::Union{Symbol,Array{Symbol}}, proxy_variable::Symbol, poly_names::Vector{Symbol}, lag_poly_names::Vector{Symbol}, id::Symbol, time::Symbol, fes_poly_degree::Int)
    
    ## Select necessary variables from data frame
    # Flatten all arguments. Some might be arrays of symbols and others might just be symbols. The following iterates over all sublists and flattens them
    all_var_symbols = vec(hcat(output, flexible_inputs..., fixed_inputs... , proxy_variable, id, time))
    all_input_symbols = vec(hcat(flexible_inputs..., fixed_inputs..., proxy_variable))
    
    # Select the data
    # select!(data, all_var_symbols)

    # Add a constant to the data
    data.constant = ones(size(data)[1])

    # Drop missings
    dropmissing!(data, unique(all_var_symbols))

    # Preallocate Φ_hat and lag_Φ_hat
    data.Φ_hat .= Vector{Float64}(undef, size(data)[1])
    data.lag_Φ_hat .= Vector{Float64}(undef, size(data)[1])

    # Preallocate polynomials in data
    if poly_names == [:NotDefinedByUser]
        poly_names = polynom_series!(;data = data, var_names = [flexible_inputs..., fixed_inputs..., proxy_variable], degree = fes_poly_degree)
        panel_lag!(data = data, id = id, time = time, variable = poly_names, lag_prefix = "lag_", lags = 1, drop_missings = false, force = false)
        lag_poly_names = Symbol.("lag_" .* string.(poly_names))
    end

    dropmissing!(data, Symbol.("lag_" .* string.(poly_names))) # Drop missings in lagged polynomials because we need them to calculate ω later.

    # Return the data
    return data, all_input_symbols, poly_names, lag_poly_names
end

## ACF estimation user inputs cleaner
function acf_user_input_cleaner!(; data::DataFrame,
                                    time::Symbol,
                                    id::Symbol,
                                    output::Symbol,
                                    flexible_inputs::Union{Symbol,Vector{Symbol}},
                                    fixed_inputs::Union{Symbol,Vector{Symbol}},
                                    proxy_variable::Symbol,
                                    fes_poly_degree::Int,
                                    ω_lom_degree::Int,
                                    opts::Dict
                                )
    
    if typeof(flexible_inputs) == Symbol
        flexible_inputs = [flexible_inputs]
    end

    if typeof(fixed_inputs) == Symbol
        fixed_inputs = [fixed_inputs]
    end

    # Check if opts Dict type is any
    if typeof(opts) != Dict{String,Any}
        # Change type to Dict{Symbol,Any}
        opts = Dict{String,Any}(opts)
    end

    if "startvals" ∉ keys(opts)
        # If no startvals are given, use zeros as default
        opts["startvals"] = vec(fastOLS(Y = data[!, output], X = Matrix{Float64}(data[!, vcat(flexible_inputs..., fixed_inputs...)])))
        opts["startvals"][opts["startvals"] .< 0] .= 0.5 # Replace negative starting values with 0.5
    end

    if "boot_reps" ∉ keys(opts)
        opts["boot_reps"] = 200
    end

    if "Optimizer" ∉ keys(opts)
        opts["Optimizer"] = NelderMead()
    end

    if "fes_print_results" ∉ keys(opts)
        opts["fes_print_results"] = true
    end
    if "ses_print_results" ∉ keys(opts)
        opts["ses_print_results"] = true
    end
    if "print_results" ∉ keys(opts)
        opts["print_results"] = true
    end

    if "call_from_boot_runs" ∉ keys(opts)
        opts["call_from_boot_runs"] = false
    end

    if "maxboottries" ∉ keys(opts)
        opts["maxboottries"] = 10
    end

    if "no_SE_est" ∉ keys(opts)
        opts["no_SE_est"] = false
    end

    if "lower_bound" ∉ keys(opts)
        opts["lower_bound"] = fill(-Inf, length(opts["startvals"]))

    end
    
    if "upper_bound" ∉ keys(opts)
        opts["upper_bound"] = fill(Inf, length(opts["startvals"]))
    end

    if "Optimizer_Options" ∉ keys(opts)
        opts["Optimizer_Options"] = Optim.Options(show_trace = false)
    end

    return fixed_inputs, flexible_inputs, opts
end

function print_fes_res(α_hat::Array{<:Number}, var::Vector{Symbol})

    fes_res_tab = hcat(var, α_hat)
    header_par = (["Variable", "Estimate"])

    println("First stage parameters")
    pretty_table(fes_res_tab, header = header_par, formatters =  ft_printf("%5.5f"), limit_printing = false)

    return nothing
end

function print_res_ses(; β_hat::Array{<:Number}, ρ_hat::Array{<:Number}, fixed_inputs::Array{Symbol}, flexible_inputs::Array{Symbol})
    
    ses_lom_res_tab = Array{Union{String, Char, <:Number}}(undef, size(ρ_hat)[1], 2)

    for j = 1:size(ses_lom_res_tab)[1]
        if j > 2
            ses_lom_res_tab[j,1] = "ω"*string(superscript_this!(string(j-1)))
        elseif j == 2
            ses_lom_res_tab[j,1] = "ω"
        elseif j == 1
            ses_lom_res_tab[j,1] = "ln(β₀)"
        end

        ses_lom_res_tab[j,2] = ρ_hat[j]
    end

    ses_prod_res_tab = hcat(vcat(flexible_inputs..., fixed_inputs...), β_hat)
    ses_res_header = (["Variable", "Estimate"])

    println("Second stage parameters")
    println("Production function")
    pretty_table(ses_prod_res_tab, header = ses_res_header, formatters =  ft_printf("%5.5f"), limit_printing = false)

    println("Productivity development")
    pretty_table(ses_lom_res_tab, header = ses_res_header, formatters =  ft_printf("%5.5f"), limit_printing = false)
end

function acf_print_all_res(; data::DataFrame, results_df::DataFrame, id::Symbol, opts::Dict)

    # Show complete set of results
    if opts["print_results"] == true
        println("")
        println("Number of observations: "*string(length(data[!, id])))
        pretty_table(results_df, show_subheader = false, formatters =  ft_printf("%5.5f"), limit_printing = false)
    end

    return results_df
end