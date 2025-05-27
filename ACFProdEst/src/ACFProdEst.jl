#####################################################################
#                  ACF Production Function Estimation               #
#                          Markus Trunschke                         #
#                    - Main Project Module File -                   #
#####################################################################
# First version: 27/06/2025                                         #
# Author: Markus Trunschke                                          #
#####################################################################
# Notes: - This file contains the main module for estimating        #
#          production functions using the ACF approach.              #
#####################################################################
module ACFProdEst

    # Read in dependencies
    using DataFrames, Revise, LinearAlgebra, Optim, ShiftedArrays, PrettyTables, Statistics, Distributions #, JLD2, CSV, NaNMath BenchmarkTools, SMTPClient, Latexify, LaTeXStrings, ProgressMeter, Plots

    # Include other files
    include("ACF Estimation Fncs.jl")
    include("Aux fncs.jl")

    # Export functions
    export acfprodest, acfprodest!, acfprodest_fes, acfprodest_fes!, acfprodest_ses, acfprodest_ses!
end
