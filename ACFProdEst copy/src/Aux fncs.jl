#####################################################################
#               Production Function Estimation Project              #
#                   Markus Trunschke and Ken Judd                   #
#                    - Auxiliary functions file -                   #
#####################################################################
# First version: 25/06/2024                                         #
# Author: Markus Trunschke                                          #
#####################################################################
# Notes: - This file contains auxiliary functions that are used     #
#          throughout the package.                                  #
#####################################################################

## Generate a vector of polynomials of the given scalar. Skip the ^0 polynomial. Don't need that. (Thanks GPT)
function polynomial_vector(x::Number, degree::Int)
    return [x^i for i in 1:degree]
end

## Save simulation data and all parameters of the data generation
function save_fnc(;kwargs...)

    # Unpack some stuff
    save_path = kwargs[:save_path]
    save_file = kwargs[:save_file]

    # Check if path exists up to last part and create it if not
    if isdir(joinpath(splitpath(save_path)[begin:end-1])) && !isdir(save_path)
        mkdir(save_path)
        println("Path to save the data did not exist but its parentfolder did. I created the missing folder.")
    end

    # Check if file already exists and increase counter if not
    save_name = save_file
    if kwargs[:save_replace] == false
        exit_flag = false
        j = 1
        while exit_flag == false
            if any(isfile(joinpath(save_path, save_file*"_"*string(j)*".jld2")) || isfile(joinpath(save_path, save_file*"_"*string(j)*".csv"))) == false
                save_name = save_file*"_"*string(j)
                exit_flag = true # Jump out of loop if file does not exist
            end
            j += 1
        end
    end

    # Save everything
    if kwargs[:save_as_jld] == true
        jldsave(joinpath(save_path, save_name*".jld2"); kwargs...)
    end
    if kwargs[:save_as_csv] == true # Save a CSV version of it as well
        CSV.write(joinpath(save_path, save_name*".csv"), kwargs[:sim_data])
    end
end

## Function calculating the within variance of a vector
function withinvar(input_vec::Vector)
    dm_vec = input_vec .- mean(input_vec)

    return var(dm_vec)
end

## Polynomial series generating function
function polynom_series!(;data::DataFrame, var_names::Union{Vector{Symbol},Symbol}, degree::Int)

    # Check if user put in an invalid degree
    if degree < 1
        throw(error("Polynomial series degree must be at least 1"))
    end

    # Convert to a vector if only one symbol put in to make my life easier
    if typeof(var_names) == Symbol
        var_names = [var_names]
    end

    # Preallocate ooutput
    poly_var_names = Set(Symbol[])

    # Recursively define polynomail calc fnc
    function generate_polynomials(data, var_names, degree, prefix, idx)

        # Stop if degree is 0
        if degree == 0
            return
        end

        # Calculate poynomials
        for i in idx:length(var_names)
            new_prefix = prefix == "" ? string(var_names[i]) : string(prefix, '⋅', var_names[i])
            varn = Symbol(new_prefix)
            push!(poly_var_names, varn)
            if varn ∈ names(data)
                data = select(data, Not(varn))
            end

            if prefix == ""
                data[:, varn] = data[:, var_names[i]]
            else
                prev_var = Symbol(prefix)
                data[:, varn] = data[:, prev_var] .* data[:, var_names[i]]
            end

            generate_polynomials(data, var_names, degree - 1, new_prefix, i)
        end
    end

    # Get vector of o
    for var in var_names
        push!(poly_var_names, var)
    end

    # Run the just-defined fnc
    generate_polynomials(data, var_names, degree, "", 1)
    
    return collect(poly_var_names)
end

## If there are already prepared columns in the dataframe and their values just need to be updated, jump in here. This is the fast version with no dynamic allocations at runtime.
function polynomial_fnc_fast!(poly_mat::Union{Array{<:Number}, SubArray{<:Number}}, degree::Int; par_cal::Bool = false)
    # Compute polynomial columns (each column of the matrix represents the i's polynomial of the first column)
    if par_cal == false
        for i in 2:degree
            poly_mat[:,i] .= @view(poly_mat[:,1]) .^ i
        end
    else
        Threads.@threads for i in 2:degree
            poly_mat[:,i] .= @view(poly_mat[:,1]) .^ i
        end
    end

    return poly_mat 
end

## Function that generates lagged values in a panel
# Panel lag function with return df
function panel_lag(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Symbol,Vector{Symbol}}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false, force::Bool = false)
    
    # Clean input
    if typeof(variable) == Symbol
        variable = [variable]
    end

    # Sort data and create new df. This is the only difference to panel_lag!()
    new_df = copy(sort(data, [id, time]))

    # Drop lag columns if they already exist if force is true. Throw error if not.
    if any(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))
        if force == true
            data = data[!, Not(filter(in(Symbol.(names(data))), Symbol.(lag_prefix, variable)))]
        else
            throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
        end
    end

    # Do the actual lagging
    lagging_that_panel!(data = new_df, id = id, time = time, variable = variable, lag_prefix = lag_prefix, lags = lags, drop_missings = drop_missings)
    
    # Return resulting Dataframe
    return new_df
end

# Panel lag function manipulating the original df
function panel_lag!(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Array{Symbol},Symbol}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false, force::Bool = false)
    
    # Clean input
    if typeof(variable) == Symbol
        variable = [variable]
    end

    if any(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))
        if force == true
            select!(data, Not(filter(in(Symbol.(names(data))), Symbol.(lag_prefix, variable))))
            # df2 = df2[!, Not(Symbol.(lag_prefix, variable))]
        else
            throw("Specified name for lag of variable already present in specified dataframe. Either set force = true, choose difference lag variable name, or rename the column.")
        end
    end

    # Sort data
    sort!(data, [id, time])

    # # Do the actual lagging
    data = lagging_that_panel!(data = data, id = id, time = time, variable = variable, lag_prefix = lag_prefix, lags = lags, drop_missings = drop_missings)

    return data
end

function lagging_that_panel!(;data::DataFrame, id::Symbol, time::Symbol, variable::Union{Symbol,Vector{Symbol}}, lag_prefix::String = "lag_", lags::Int = 1, drop_missings::Bool = false)

    # Generate lagged values per id and select all but the original variable (causes problems in join). The ShiftedArrays.lag function names the lagged column itself with variable_lag
    df_lag = select(data, [id, time, variable...])
    
    lag_variables::Vector{Symbol} = []

    # Lag all variables
    for lag_var in [time, variable...]
        transform!(groupby(df_lag, id), lag_var => ShiftedArrays.lag)

        push!(lag_variables, Symbol(string(lag_var) .* "_lag"))
    end

    # Drop missings in lagged variables we just generated
    dropmissing!(df_lag, lag_variables)
    
    # Check if lag is actually only the expected lag apart
    for var in lag_variables[2:end]
        df_lag[!, var] = ifelse.(df_lag[!,time] .- lags .== df_lag[!,Symbol(string(time)*"_lag")], df_lag[!,var], missing)
    end

    select!(df_lag, [time, id, lag_variables[2:end]...]) # Drop lagged time variable from df

    # Combine lagged variable with original data and sort it.
    sort!(leftjoin!(data, df_lag, on = [id, time]), [id, time])

    # Drop missings in lagged variable we just generated if user wants to
    if drop_missings == true
        dropmissing!(data, lag_variables[2:end])
    end
    
    # Rename variable to user-specified name
    for var in variable
        rename!(data, Symbol(string(var)*"_lag") => Symbol(lag_prefix*string(var)))
    end

    # Return result
    return data
end

## Fast OLS function
function fastOLS(; Y::Union{Matrix, Vector}, X::Union{Matrix{<:Number}, Vector{<:Number}}, multicolcheck::Bool = true, force::Bool = false) # Y::Vector{<:Number}, X::Matrix{<:Number}

    x_len::Int = 1
    if typeof(X) <: Matrix
        x_len = size(X)[2]
    end

    # Check for perfect multicoliniarity in regressors
    if multicolcheck == true && rank(X) != x_len
        if force == false
            throw("ERROR: Regressors are perfectly multicolliniar")
        else
            for j = eachindex(eachcol(X))
                if rank(X[:,Not(j)]) == x_len - 1
                    X = copy(X)
                    X = X[:,Not(j)]

                    x_len = size(X)[2]

                    println("Dropped regressor "*string(j)*" because of multicolliniarity!")
                    break
                end
            end
        end
    end

    coefs = Array{Float64}(undef, x_len, 1)
    cache = Array{Float64}(undef, x_len, x_len)
    
    # OLS
    mul!(cache, X', X) # X'X
    mul!(coefs, X', Y) # X'*Y
    ldiv!(cholesky!(Hermitian(cache)), coefs) # inv(X'*X)*X*y, see https://discourse.julialang.org/t/memory-allocation-left-division-operator/103111/3 for an explination why this is the fastest way (even though not optimal for ill-conditioned matrices)
    
    return coefs
end

## Function that returns the superscript of the corresponding input character (b/c Julia does not have a simple function for that)
# Define a dictionary for superscript characters
const superscript_map = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
    'a' => 'ᵃ', 'b' => 'ᵇ', 'c' => 'ᶜ', 'd' => 'ᵈ', 'e' => 'ᵉ',
    'f' => 'ᶠ', 'g' => 'ᵍ', 'h' => 'ʰ', 'i' => 'ⁱ', 'j' => 'ʲ',
    'k' => 'ᵏ', 'l' => 'ˡ', 'm' => 'ᵐ', 'n' => 'ⁿ', 'o' => 'ᵒ',
    'p' => 'ᵖ', 'r' => 'ʳ', 's' => 'ˢ', 't' => 'ᵗ', 'u' => 'ᵘ',
    'v' => 'ᵛ', 'w' => 'ʷ', 'x' => 'ˣ', 'y' => 'ʸ', 'z' => 'ᶻ',
    '+' => '⁺', '-' => '⁻', '=' => '⁼', '(' => '⁽', ')' => '⁾'
)

function superscript_this!(c::String) # Need to use a string as input because I don't understand Chars in Julia. Char(5) returns a different unicode than string(5). And the superscript of Char(5) does not  work
    # Return the superscript character if it exists in the map, else return the original character
    return get(superscript_map, c[1], c[1])
end

## Function drawing a sample of firms
function draw_sample(;data::DataFrame, id::Symbol)

    id_list = unique(data[!, id])
    
    ## Randomly chose observation
    # Generate vector with sample of IDs
    boot_choose = DataFrame(id => sample(id_list, length(id_list); replace=true))

    # Add unique identifier
    boot_choose.unique_id = range(1,length = length(boot_choose[!, id]))

    # Select obs in dataset that match the drawn numbers
    sample_data = innerjoin(data, boot_choose, on = id)

    return sample_data
end