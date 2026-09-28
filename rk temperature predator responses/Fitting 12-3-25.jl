# ============================================================================
# PARAMECIUM LOGISTIC GROWTH - BAYESIAN FITTING ONLY
# ============================================================================
# This script performs Bayesian fitting of logistic growth curves
# Outputs: Fitted parameters and curve data for visualization
# ============================================================================

# ----------------------------------------------------------------------------
# 1. SETUP: ENVIRONMENT & PACKAGES
# ----------------------------------------------------------------------------
dir = "/Users/thisisme710/Desktop/Research/Paramecium logistic growth/Fall 2025/Data Sheets"
cd(dir)

using MCMCChains, Turing, DataFrames, CSV, Distributions, Random
using Statistics
using Plots
using StatsPlots


Turing.setprogress!(false)

trace_dir = "/Users/thisisme710/Desktop/Research/Paramecium logistic growth/Fall 2025/plots/Trace plots"

isdir(trace_dir) || mkpath(trace_dir)   

# ----------------------------------------------------------------------------
# 2. DATA LOADING
# ----------------------------------------------------------------------------
df = CSV.read("Paramecium_logistic_growth_data.csv", DataFrame)
df4 = CSV.read("July_Phenotype_Data.csv", DataFrame)

# ----------------------------------------------------------------------------
# 3. INITIAL DATA PREPARATION
# ----------------------------------------------------------------------------
dishes_to_use = unique(df.DishID)
num_studies = length(dishes_to_use)

# Initialize result vectors
fitted_k = zeros(Float64, num_studies)
fitted_r = zeros(Float64, num_studies)
env_temp = zeros(Int64, num_studies)
copepods = zeros(Int64, num_studies)
athome = zeros(Int64, num_studies)
tech_rep = zeros(Int64, num_studies)

# Extract population names
names_df = combine(groupby(df, :DishID), :Bp => first => :names)
names = names_df.names
pop_names = unique(names)

# ----------------------------------------------------------------------------
# 4. BAYESIAN LOGISTIC GROWTH FITTING
# ----------------------------------------------------------------------------
 






@model function logistic_growth(times,density,starting_d) 
    #r0 ~ truncated(Normal(1,20),0,15) # maximum growth rate
    #r0 ~ truncated(Normal(1,10), lower = 0) # maximum growth rate
    #k ~ truncated(Normal(60000,60000),0,6*maximum(density)) # carrying capacity
    #k ~ truncated(Normal(20000,20000), lower =0) # carrying capacity
    r0 ~ truncated(Normal(1,20),0, Inf)
    k ~ truncated(Normal(60000,60000),0, Inf)
σ ~ InverseGamma(2,3)



	for i in 1:length(times)
        density[i] ~ truncated(Normal(k./(1+((k-starting_d)./starting_d).*exp(-r0.*times[i]))),0,maximum(density))
	end
end

# Fit models for each dish
for i = 1:num_studies
    println("Fitting dish $i of $num_studies")
    indices = findall(df.DishID .== i)
    days = df.Day_of_experiment[indices]
    density = df.Abundance[indices]
    temps = df.Temp[indices]
    cops = df.COP_01[indices]
    reps = df.Rep[indices]
    homes = df.Home[indices]
    
    model_logistic = logistic_growth(days, density, density[1])
    
    chain_logistic = sample(
        model_logistic,
        NUTS(20000, 0.9),
        MCMCThreads(),

        20000,
        4
    )
    
 #save trace plots 
  plot_trace= plot(chain_logistic)
  savefig(plot_trace, joinpath(trace_dir, "Traceplot_Dish_$(i).png"))

    fitted_params = DataFrame(summarystats(chain_logistic))
    fitted_k[i] = fitted_params[2, 2]
    fitted_r[i] = fitted_params[1, 2]
    env_temp[i] = temps[1]
    copepods[i] = cops[1]
    tech_rep[i] = reps[1]
    athome[i] = homes[1]
end

# ----------------------------------------------------------------------------
# 5. GENERATE FITTED CURVES DATA
# ----------------------------------------------------------------------------
fitted_curves_data = DataFrame()

for i = 1:num_studies
    indices = findall(df.DishID .== i)
    days = df.Day_of_experiment[indices]
    density = df.Abundance[indices]
    
    # Generate predicted values from fitted parameters
    n0 = density[1]
    days_predict = LinRange(0, maximum(days), 100)
    predicted_abundance = fitted_k[i] ./ 
        (1.0 .+ (fitted_k[i] / n0 - 1.0) .* exp.(-fitted_r[i] .* days_predict))
    
    # Create DataFrame for this dish
    dish_data = DataFrame(
        DishID = fill(i, length(days_predict)),
        Day = days_predict,
        Predicted_Abundance = predicted_abundance,
        Fitted_r = fill(fitted_r[i], length(days_predict)),
        Fitted_K = fill(fitted_k[i], length(days_predict)),
        Temperature = fill(env_temp[i], length(days_predict)),
        Copepods = fill(copepods[i], length(days_predict)),
        Population = fill(names[i], length(days_predict))
    )
    
    append!(fitted_curves_data, dish_data)
end

# Raw observed data for overlay
raw_data_for_plots = DataFrame(
    DishID = df.DishID,
    Day = df.Day_of_experiment,
    Observed_Abundance = df.Abundance,
    Temperature = df.Temp,
    Copepods = df.COP_01,
    Population = df.Bp
)

# ----------------------------------------------------------------------------
# 6. PHENOTYPE DATA PREPARATION
# ----------------------------------------------------------------------------
cell_length = df4.mean_major
cell_width = df4.mean_minor
cell_vol = 1.333333 .* pi .* cell_width.^2 .* cell_length

df4[!, :cell_vol] .= cell_vol
df4[!, :aspect_ratio] = df4.mean_major ./ df4.mean_minor
df4[!, :speed] = df4.net_speed
df4[!, :turning] = df4.sd_turning
df4[!, :major] = df4.mean_major
df4[!, :minor] = df4.mean_minor

# Aggregate phenotype data by population
gd = groupby(df4, [:Dish_pop])
pheno_means = combine(gd,
    :cell_vol => mean => :cell_vol,
    :aspect_ratio => mean => :aspect_ratio,
    :speed => mean => :speed,
    :turning => mean => :sd_turning,
    :major => mean => :mean_major,
    :minor => mean => :mean_minor
)
pheno_means = DataFrame(pheno_means)

# ----------------------------------------------------------------------------
# 7. CREATE FULL DATASET (df2)
# ----------------------------------------------------------------------------
df2 = DataFrame(
    ID = 1:num_studies,
    fitted_k = fitted_k,
    fitted_r = fitted_r,
    env_temp = env_temp,
    copepods = copepods,
    athome = athome,
    tech_rep = tech_rep,
    names = names,
    tech_rep2 = collect(1:num_studies),
    cell_vol = Vector{Union{Missing, Float64}}(missing, num_studies),
    aspect_ratio = Vector{Union{Missing, Float64}}(missing, num_studies),
    speed = Vector{Union{Missing, Float64}}(missing, num_studies),
    mean_major = Vector{Union{Missing, Float64}}(missing, num_studies),
    mean_minor = Vector{Union{Missing, Float64}}(missing, num_studies),
    sd_turning = Vector{Union{Missing, Float64}}(missing, num_studies)
)

# Match phenotype data to fitted data
for i = 1:nrow(df2)
    name_to_use = df2.names[i]
    row_to_use = findfirst(pheno_means.Dish_pop .== name_to_use)
    if row_to_use === nothing
        @warn "No phenotype match for $(name_to_use)"
        continue
    end
    
    df2.cell_vol[i] = pheno_means.cell_vol[row_to_use]
    df2.aspect_ratio[i] = pheno_means.aspect_ratio[row_to_use]
    df2.speed[i] = pheno_means.speed[row_to_use]
    df2.sd_turning[i] = pheno_means.sd_turning[row_to_use]
    df2.mean_major[i] = pheno_means.mean_major[row_to_use]
    df2.mean_minor[i] = pheno_means.mean_minor[row_to_use]
end

# Remove problematic populations (26_COP_04)
deleteat!(df2, [85, 86, 87, 88, 89, 90])

# Create home subset
df2_home = df2[findall(df2.athome .== 1), :]

# ----------------------------------------------------------------------------
# 8. CREATE HOME/AWAY COMPARISON DATASET (df3)
# ----------------------------------------------------------------------------
num_rows = 16

# Initialize comparison vectors
k_at_home = zeros(Float64, 1, num_rows)
r_at_home = zeros(Float64, 1, num_rows)
k_away = zeros(Float64, 1, num_rows)
r_away = zeros(Float64, 1, num_rows)
home_home = zeros(Int64, 1, num_rows)
home_away = zeros(Int64, 1, num_rows)
env_temp_home = zeros(Int64, 1, num_rows)
env_temp_away = zeros(Int64, 1, num_rows)
copepods_home = zeros(Int64, 1, num_rows)
copepods_away = zeros(Int64, 1, num_rows)
speed_home = zeros(Float64, 1, num_rows)
speed_away = zeros(Float64, 1, num_rows)
aspect_ratio_home = zeros(Float64, 1, num_rows)
aspect_ratio_away = zeros(Float64, 1, num_rows)
cellvol_home = zeros(Float64, 1, num_rows)
cellvol_away = zeros(Float64, 1, num_rows)
minor_home = zeros(Float64, 1, num_rows)
minor_away = zeros(Float64, 1, num_rows)
major_home = zeros(Float64, 1, num_rows)
major_away = zeros(Float64, 1, num_rows)
turning_home = zeros(Float64, 1, num_rows)
turning_away = zeros(Float64, 1, num_rows)
k_diff = zeros(Float64, 1, num_rows)
r_diff = zeros(Float64, 1, num_rows)
zero_vec = zeros(Float64, 1, num_rows)

# Calculate home/away comparisons
for i = 1:16
    indices_home = findall(df2.names .== pop_names[i] .&& df2.athome .== 1)
    indices_away = findall(df2.names .== pop_names[i] .&& df2.athome .== 0)
    
    if isempty(indices_home) || isempty(indices_away)
        @warn "No home/away data for $(pop_names[i])"
        continue
    end
    
    # HOME
    k_at_home[i] = mean(df2.fitted_k[indices_home])
    r_at_home[i] = mean(df2.fitted_r[indices_home])
    env_temp_home[i] = df2.env_temp[indices_home][1]
    copepods_home[i] = df2.copepods[indices_home][1]
    cellvol_home[i] = df2.cell_vol[indices_home][1]
    aspect_ratio_home[i] = df2.aspect_ratio[indices_home][1]
    speed_home[i] = df2.speed[indices_home][1]
    turning_home[i] = df2.sd_turning[indices_home][1]
    major_home[i] = df2.mean_major[indices_home][1]
    minor_home[i] = df2.mean_minor[indices_home][1]
    home_home[i] = 1
    
    # AWAY
    k_away[i] = mean(df2.fitted_k[indices_away])
    r_away[i] = mean(df2.fitted_r[indices_away])
    env_temp_away[i] = df2.env_temp[indices_away][1]
    copepods_away[i] = df2.copepods[indices_away][1]
    cellvol_away[i] = df2.cell_vol[indices_away][1]
    aspect_ratio_away[i] = df2.aspect_ratio[indices_away][1]
    speed_away[i] = df2.speed[indices_away][1]
    turning_away[i] = df2.sd_turning[indices_away][1]
    major_away[i] = df2.mean_major[indices_away][1]
    minor_away[i] = df2.mean_minor[indices_away][1]
    home_away[i] = 0
    
    # DIFFERENCES
    k_diff[i] = k_away[i] - k_at_home[i]
    r_diff[i] = r_away[i] - r_at_home[i]
end

# Construct df3
df3 = DataFrame(ID = 1:32)
df3.mean_k = vcat(vec(k_at_home), vec(k_away))
df3.mean_r = vcat(vec(r_at_home), vec(r_away))
df3.cellvol = vcat(vec(cellvol_home), vec(cellvol_away))
df3.aspect_ratio = vcat(vec(aspect_ratio_home), vec(aspect_ratio_away))
df3.speed = vcat(vec(speed_home), vec(speed_away))
df3.turning = vcat(vec(turning_home), vec(turning_away))
df3.major = vcat(vec(major_home), vec(major_away))
df3.minor = vcat(vec(minor_home), vec(minor_away))
df3.home = vcat(vec(home_home), vec(home_away))
df3.copepods = vcat(vec(copepods_home), vec(copepods_away))
df3.env_temp = vcat(vec(env_temp_home), vec(env_temp_away))
names16 = collect(first(pop_names, num_rows))
df3.names = vcat(names16, names16)
df3.r_diffs = vcat(vec(r_diff), vec(zero_vec))
df3.k_diffs = vcat(vec(k_diff), vec(zero_vec))

# Remove 26_COP_04 from df3
deleteat!(df3, [15, 31])

# Create home subset
df3_home = df3[findall(df3.home .== 1), :]

# ----------------------------------------------------------------------------
# 9. SAVE ALL OUTPUTS
# ----------------------------------------------------------------------------
CSV.write("rK_full_fitted_dataset.csv", df2)
CSV.write("rK_means_fitted_dataset.csv", df3)
CSV.write("rK_d2_home.csv", df2_home)
CSV.write("rK_d3_home.csv", df3_home)
CSV.write("fitted_curves_predicted.csv", fitted_curves_data)
CSV.write("fitted_curves_observed.csv", raw_data_for_plots)

