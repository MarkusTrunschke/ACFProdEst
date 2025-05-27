# Rad in dependencies
using ACFProdEst, DataFrames, StatFiles, CSV

# Read in the data
path_to_data = joinpath(pwd(),"ACFProdEst", "data", "sim_data_ACF.dta") # Define path to data in string. I use the joinpath function because the paths can have different folder deviders for Mac and Windows (/ vs \). pwd() gives the current working directory.
df = DataFrame(load(path_to_data)) # Load the data into a DataFrame

# Run estimation with defaul options (minimal example)
acfprodest(data = df,
        time = :year,
        id = :ID,
        output = :y,
        flexible_inputs = [:l],
        fixed_inputs = [:k],
        proxy_variable = :m,
        fes_poly_degree = 1,
        ω_lom_degree = 1,
        opts = Dict("boot_reps" => 10))



using Optim # Need to load this package to set the Optimizer option in the opts dictionary

# Run estimation with custom options
opts = Dict("startvals" => [0.5, 0.5], # Initial values for output elasticities
            "boot_reps" => 10, # Number of bootstrap replications to perform
            "Optimizer" => Fminbox(NelderMead()), # Algorithm to use for optimization. Fminbox is used to set bounds on the output elasticities.
            "lower_bound" => [0.0, 0.0], # Set lower bound for output elasticities (needs Fminbox optimizer)
            "upper_bound" => [1.0, 1.0], # Set upper bound for output elasticities (needs Fminbox optimizer)
            "fes_print_results" => true, # Show estimates of first stage right after they are computed
            "ses_print_results" => true, # Show estimates of second stage right after they are computed
            "maxboottries" => 10, # Number of times the code tries to estimate the model in a bootstrap iteration before giving up. (After each unsuccessful try it selectes a new sample and tries again)
            "Optimizer_Options" => Optim.Options(f_reltol = 1e-2,
                                                 x_reltol = 1e-2,
                                                 g_tol = 1e-2, # √(Σ(yᵢ-ȳ)²)/n ≤ 1.0e-13 (only sets g_abstol, not outer_g_abstol)
                                                 allow_f_increases = true,
                                                 show_trace = false,
                                                 extended_trace = false,
                                                 show_every = 1,
                                                 time_limit = NaN,
                                                 store_trace = false,
                                                 iterations = 100000 # Optimization options for the Optim package
                                                 ))

acfprodest(data = df,
        time = :year,
        id = :ID,
        output = :y,
        flexible_inputs = [:l],
        fixed_inputs = [:k],
        proxy_variable = :m,
        fes_poly_degree = 1,
        ω_lom_degree = 1,
        opts = opts)