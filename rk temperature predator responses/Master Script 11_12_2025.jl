dir = "/Users/thisisme710/Desktop/Research/Paramecium logistic growth/Fall 2025/Data Sheets"
cd(dir)

using MCMCChains
using Pkg
using Measures   # ← add Measures here
using Colors  # optional, for RGB()
using Pkg
using CategoricalArrays
using Tables
using StatsPlots
using XLSX
using Turing, MCMCChains 
using DataFrames, CSV, Distributions, Random
using MixedModels, GLM, Statistics
using LsqFit
using Turing
Turing.setprogress!(false)

df = CSV.read("Paramecium_logistic_growth_data.csv",DataFrame)
# Drop rows where ID is in a list of IDs
#ids_to_drop = [85,86,87,88,89,90]
#filter!(:DishID => id -> !(id in ids_to_drop), df)

dishes_to_use = unique(df.DishID)
num_studies = length(dishes_to_use)
fitted_k   = zeros(Float64, num_studies)
fitted_r   = zeros(Float64, num_studies)
env_temp   = zeros(Int64,   num_studies)
copepods   = zeros(Int64,   num_studies)
athome     = zeros(Int64,   num_studies)
tech_rep   = zeros(Int64,   num_studies)
estimated_k = zeros(Float64, num_studies)
estimated_r = zeros(Float64, num_studies)
names_df = combine(groupby(df, :DishID), :Bp => first => :names)
names = names_df.names
pop_names = unique(names)

#=########################### logistic fits using least squares
#for i = 1:num_studies
    println(i)
    indices = findall(df.DishID .== i)
    days = df.Day_of_experiment[indices]
    density = df.Abundance[indices]
    row_name = df.Bp[indices]
    env_temp[i] = df.Temp[indices][1]
    copepods[i] = df.COP_01[indices][1]

    # start with some empirical estimates of r and k
    estimated_k[i] = mean(density[end-1:end])
    estimated_r[i] = log(density[2]/density[1])/(days[2]-days[1])
    if estimated_r[i] <= 0
        estimated_r[i] = log(density[3]/density[2])/(days[3]-days[2])
    end
    if estimated_r[i] <= 0
        estimated_r[i] = log(density[4]/density[3])/(days[4]-days[3])
    end
    if isnan(estimated_r[i]) || estimated_r[i] <= 0
        estimated_r[i] = 0.1
    end
    if estimated_k[i] == 0
        estimated_k[i] = 10
    end

    n0 = density[1]
    @. model(x,p) = p[2] ./ 
                (1.0 .+ (p[2] / n0 - 1.0) .* exp.(-p[1] .* x))

    p0 = [estimated_r[i], estimated_k[i]]
    lower_p = [0.0; 0.0]
    upper_p = [30.0; 2.0*maximum(density)]
    fit = curve_fit(model, days, density, p0; lower = lower_p, upper = upper_p)
    param = fit.param
    # or coef(fit)
    cis = confidence_interval(fit, 0.05)

    fitted_k[i] = param[2]
    fitted_r[i] = param[1]

    xrange = LinRange(0, maximum(days), 20)
    y_fit = model(xrange, fit.param)

    #p1 = scatter(days,density,markercolor=:black)
    #plot!(p1,xrange, y_fit, color=:blue, label="", linewidth=2)
    #xlabel!("Days")
    #ylabel!("Density (cells per mL)")

    #display(p1)

   # push!(names,row_name[1])

end =#

############################ logistic fits using bayesian approach
# logistic growth function

@model function logistic_growth(times,density,starting_d) 
    r0 ~ truncated(Normal(1,10),0,3) # maximum growth rate
    k ~ truncated(Normal(20000,20000),0,2*maximum(density)) # carrying capacity
    σ ~ InverseGamma(2,3)
	for i in 1:length(times)
        density[i] ~ truncated(Normal(k./(1+((k-starting_d)./starting_d).*exp(-r0.*times[i]))),0,maximum(density))
	end
end

for i = 1:num_studies
    println(i)
    indices = findall(df.DishID .== i)
    days = df.Day_of_experiment[indices]
    density = df.Abundance[indices]
    temps = df.Temp[indices]
    cops = df.COP_01[indices]
    row_name = df.Bp[indices]
    reps = df.Rep[indices]
    model_logistic = logistic_growth(days,density,density[1])

    # call the fitting for type 2
chain_logistic = sample(
    model_logistic,
    NUTS(1500, 0.9),
    MCMCThreads(),     # ← correct spelling
    3000,
    4
)
    

    fitted_params = DataFrame(summarystats(chain_logistic))
    fitted_k[i] = fitted_params[2,2]
    fitted_r[i] = fitted_params[1,2]
    env_temp[i] = temps[1]
    copepods[i] = cops[1]
    tech_rep[i] = reps[1]

end

# pull home and dish data
for i = 1:num_studies
    println(i)
    indices = findall(df.DishID .== i)
    homes = df.Home[indices]
    row_name = df.Bp[indices]
    athome[i] = homes[1]
end

# group results into a data frame
df2 = DataFrame(
    ID        = 1:num_studies,
    fitted_k  = fitted_k,
    fitted_r  = fitted_r,
    env_temp  = env_temp,
    copepods  = copepods,
    athome    = athome,
    tech_rep  = tech_rep,
    names     = names,
    tech_rep2 = collect(1:num_studies),
    cell_vol     = Vector{Union{Missing, Float64}}(missing, num_studies),
    aspect_ratio = Vector{Union{Missing, Float64}}(missing, num_studies),
    speed        = Vector{Union{Missing, Float64}}(missing, num_studies),
    mean_major = Vector{Union{Missing, Float64}}(missing, num_studies),
    mean_minor = Vector{Union{Missing, Float64}}(missing, num_studies),
    sd_turning = Vector{Union{Missing, Float64}}(missing, num_studies)
)

# grab and pair up the cell volume data
df4 = CSV.read("July_Phenotype_Data.csv", DataFrame)
# grab lengths and calculate the volume
cell_length = df4.mean_major
cell_width = df4.mean_minor
cell_vol = 1.333333 .* pi .* cell_width.^2 .* cell_length
# add it to the df
df4[!,:cell_vol] .= cell_vol
df4[!, :aspect_ratio] = df4.mean_major ./ df4.mean_minor
df4[!, :speed]        = df4.net_speed  
df4[!, :turning]        = df4.sd_turning   
df4[!, :major]        = df4.mean_major
df4[!, :minor]        = df4.mean_minor  
   

# set up to take the average by group
gd = groupby(df4, [:Dish_pop])
pheno_means = combine(gd,
    :cell_vol     => mean => :cell_vol,
    :aspect_ratio => mean => :aspect_ratio,
    :speed        => mean => :speed,
    :turning      => mean => :sd_turning,
    :major        => mean => :mean_major,
    :minor        => mean => :mean_minor
)
    #change to df
pheno_means = DataFrame(pheno_means)
# open an empty vector and cycle through to match means with the right rows
for i = 1:nrow(df2)
    name_to_use = df2.names[i]
    row_to_use  = findfirst(pheno_means.Dish_pop .== name_to_use)
    if row_to_use === nothing
        @warn "No phenotype match for $(name_to_use)"
        continue
    end

    df2.cell_vol[i]     = pheno_means.cell_vol[row_to_use]
    df2.aspect_ratio[i] = pheno_means.aspect_ratio[row_to_use]
    df2.speed[i]        = pheno_means.speed[row_to_use]
    df2.sd_turning[i]   = pheno_means.sd_turning[row_to_use]
    df2.mean_major[i]   = pheno_means.mean_major[row_to_use]
    df2.mean_minor[i]   = pheno_means.mean_minor[row_to_use]
end

# save the dataset
CSV.write("rK_full_fitted_dataset.csv", df2)

# just check if we are getting the same pattern
p1 = scatter(df2[!, :fitted_r], df2[!, :fitted_k], legend=false)
print = p1
############## now build the home versus away DataSet


num_rows = 16

k_at_home = zeros(Float64,1,num_rows)
r_at_home = zeros(Float64,1,num_rows)
k_away = zeros(Float64,1,num_rows)
r_away = zeros(Float64,1,num_rows)
home_home = zeros(Int64,1,num_rows)
home_away = zeros(Int64,1,num_rows)
env_temp_home = zeros(Int64,1,num_rows)
env_temp_away = zeros(Int64,1,num_rows)
copepods_home = zeros(Int64,1,num_rows)
copepods_away = zeros(Int64,1,num_rows)
speed_home = zeros(Float64,1,num_rows)
speed_away = zeros(Float64,1,num_rows)
aspect_ratio_home = zeros(Float64,1,num_rows)
aspect_ratio_away = zeros(Float64,1,num_rows)
cellvol_home = zeros(Float64,1,num_rows)
cellvol_away = zeros(Float64,1,num_rows)
minor_home = zeros(Float64,1,num_rows)
minor_away = zeros(Float64,1,num_rows)
major_home = zeros(Float64,1,num_rows)
major_away = zeros(Float64,1,num_rows)
turning_home = zeros(Float64,1,num_rows)
turning_away = zeros(Float64,1,num_rows)
k_diff = zeros(Float64,1,num_rows)
r_diff = zeros(Float64,1,num_rows)
zero_vec = zeros(Float64,1,num_rows)
names_home = String[]
names_away = String[]


for name in pop_names
    n_home = count(df2.names .== name .&& df2.athome .== 1)
    n_away = count(df2.names .== name .&& df2.athome .== 0)
    println("$(name): home=$(n_home), away=$(n_away)")
end

####.  26_COP_04: home=0, away=0 ######

for i = 1:16
    indices_home = findall(df2.names .== pop_names[i] .&& df2.athome .== 1)
    indices_away = findall(df2.names .== pop_names[i] .&& df2.athome .== 0)

     if isempty(indices_home) || isempty(indices_away)
        @warn "No home/away data for $(pop_names[i])"
        continue
    end
    
    # --- HOME ---
    k_at_home[i]     = mean(df2.fitted_k[indices_home])
    r_at_home[i]     = mean(df2.fitted_r[indices_home])
    env_temp_home[i] = df2.env_temp[indices_home][1]
    copepods_home[i] = df2.copepods[indices_home][1]
    cellvol_home[i]  = df2.cell_vol[indices_home][1]
    aspect_ratio_home[i]  = df2.aspect_ratio[indices_home][1]  
    speed_home[i]    = df2.speed[indices_home][1]
    turning_home[i] = df2.sd_turning[indices_home][1]
    major_home[i]   = df2.mean_major[indices_home][1]
    minor_home[i]   = df2.mean_minor[indices_home][1]

    home_home[i] = 1
    push!(names_home, pop_names[i])

    # --- AWAY ---
    k_away[i]     = mean(df2.fitted_k[indices_away])
    r_away[i]     = mean(df2.fitted_r[indices_away])
    env_temp_away[i] = df2.env_temp[indices_away][1]
    copepods_away[i] = df2.copepods[indices_away][1]
    cellvol_away[i]  = df2.cell_vol[indices_away][1]      
    aspect_ratio_away[i]  = df2.aspect_ratio[indices_away][1]   
    speed_away[i]    = df2.speed[indices_away][1]
    turning_away[i] = df2.sd_turning[indices_away][1]
    major_away[i]   = df2.mean_major[indices_away][1]
    minor_away[i]   = df2.mean_minor[indices_away][1]

    home_away[i] = 0
    push!(names_away, pop_names[i])

    # --- DIFFS ---
    k_diff[i] = k_away[i] - k_at_home[i]
    r_diff[i] = r_away[i] - r_at_home[i]
end

df3 = DataFrame(ID=1:32)

df3.mean_k       = vcat(vec(k_at_home),        vec(k_away))
df3.mean_r       = vcat(vec(r_at_home),        vec(r_away))
df3.cellvol      = vcat(vec(cellvol_home),     vec(cellvol_away))
df3.aspect_ratio = vcat(vec(aspect_ratio_home),vec(aspect_ratio_away))
df3.speed        = vcat(vec(speed_home),       vec(speed_away))
df3.turning      = vcat(vec(turning_home),     vec(turning_away))
df3.major        = vcat(vec(major_home),       vec(major_away))
df3.minor        = vcat(vec(minor_home),       vec(minor_away))
df3.home         = vcat(vec(home_home),        vec(home_away))
df3.copepods     = vcat(vec(copepods_home),    vec(copepods_away))
df3.env_temp     = vcat(vec(env_temp_home),    vec(env_temp_away))
names16 = collect(first(pop_names, num_rows))  # num_rows = 16 above
df3.names = vcat(names16, names16)             # 32 names to match df3
df3.r_diffs      = vcat(vec(r_diff),           vec(zero_vec))
df3.k_diffs      = vcat(vec(k_diff),           vec(zero_vec))

CSV.write("rK_means_fitted_dataset.csv", df3)

# delete the 26COP04's because they all pooped out
deleteat!(df2, [85,86,87,88,89,90])
CSV.write("rK_full_reduced.csv", df2)
deleteat!(df3, [15,31])
CSV.write("rK_means_reduced.csv", df3)

df2_home = df2[findall(df2.athome .== 1),:]
CSV.write("rK_d2_home.csv", df2_home)
df3_home = df3[findall(df3.home .== 1),:]
CSV.write("rK_d3_home.csv", df3_home)




df3_home = CSV.read("rK_d3_home.csv", DataFrame)



# ------------------------

# ADDITIVE models ( + )
# ------------------------
lm_r_add   = lm(@formula(mean_r       ~ env_temp + copepods), df3_home)
lm_k_add   = lm(@formula(mean_k       ~ env_temp + copepods), df3_home)
lm_cv_add  = lm(@formula(cellvol      ~ env_temp + copepods), df3_home)
lm_sp_add  = lm(@formula(speed        ~ env_temp + copepods), df3_home)
lm_ar_add  = lm(@formula(aspect_ratio ~ env_temp + copepods), df3_home)
lm_tr_add  = lm(@formula(turning      ~ env_temp + copepods), df3_home)
lm_mj_add  = lm(@formula(major        ~ env_temp + copepods), df3_home)
lm_mn_add  = lm(@formula(minor        ~ env_temp + copepods), df3_home)

# ------------------------
# INTERACTION models ( * )
# ------------------------
lm_r_int   = lm(@formula(mean_r       ~ env_temp * copepods), df3_home)
lm_k_int   = lm(@formula(mean_k       ~ env_temp * copepods), df3_home)
lm_cv_int  = lm(@formula(cellvol      ~ env_temp * copepods), df3_home)
lm_sp_int  = lm(@formula(speed        ~ env_temp * copepods), df3_home)
lm_ar_int  = lm(@formula(aspect_ratio ~ env_temp * copepods), df3_home)
lm_tr_int  = lm(@formula(turning      ~ env_temp * copepods), df3_home)
lm_mj_int  = lm(@formula(major        ~ env_temp * copepods), df3_home)
lm_mn_int  = lm(@formula(minor        ~ env_temp * copepods), df3_home)

# ------------------------
# NO-COP versions (temperature only)
# ------------------------
lm_r_temp  = lm(@formula(mean_r       ~ env_temp), df3_home)
lm_k_temp  = lm(@formula(mean_k       ~ env_temp), df3_home)
lm_cv_temp = lm(@formula(cellvol      ~ env_temp), df3_home)
lm_sp_temp = lm(@formula(speed        ~ env_temp), df3_home)
lm_ar_temp = lm(@formula(aspect_ratio ~ env_temp), df3_home)
lm_tr_temp = lm(@formula(turning      ~ env_temp), df3_home)
lm_mj_temp = lm(@formula(major        ~ env_temp), df3_home)
lm_mn_temp = lm(@formula(minor        ~ env_temp), df3_home) 


# some quick plots
@df df3_home boxplot(:copepods, :cellvol, fillalpha=0.75, linewidth=2)
@df df3_home boxplot(:env_temp, :cellvol, fillalpha=0.75, linewidth=2)












# ---- cellvol : µm³ → mm³ ------------------------------------------------
#df3_home[!, :cellvol] .= df3_home[!, :cellvol] ./ 1_000_000   

# ---- (optional) standardise ALL numeric predictors ----------------------
numeric_cols = [:cellvol, :aspect_ratio, :turning, :speed, :major, :minor,
                :env_temp, :copepods]

for col in numeric_cols
    μ = mean(skipmissing(df3_home[!, col]))
    σ = std(skipmissing(df3_home[!, col]))
    if σ == 0
        @warn "Column $col has zero variance – skipping standardisation"
        continue
    end
    df3_home[!, col] .= (df3_home[!, col] .- μ) ./ σ
end

# Function to extract model results into a DataFrame
function extract_model_results(model, trait_name, response_var, model_formula)
    try
        coef_table = coeftable(model)
        n_terms = length(coef_table.rownms)
        
        results = DataFrame(
            Trait = fill(trait_name, n_terms),
            Response = fill(response_var, n_terms),
            Model = fill(string(model_formula), n_terms),
            Term = coef_table.rownms,
            Coef = coef_table.cols[1],
            StdError = coef_table.cols[2],
            T = coef_table.cols[3],
            Pr_t = coef_table.cols[4],
            Lower95 = confint(model)[:, 1],
            Upper95 = confint(model)[:, 2]
        )
        
        return results
    catch e
        println("Warning: Could not extract results for $model_formula - $e")
        return DataFrame()
    end
end

# ------------------------
# DEFINE ALL MODELS (exactly what you asked for)
# ------------------------

# ---------- 1. speed ----------
lm_r_envcope_speed   = lm(@formula(mean_r ~ env_temp + copepods + speed), df3_home)
lm_r_env_speed       = lm(@formula(mean_r ~ env_temp + speed), df3_home)
lm_K_envcope_speed   = lm(@formula(mean_k ~ env_temp + copepods + speed), df3_home)
lm_K_env_speed       = lm(@formula(mean_k ~ env_temp + speed), df3_home)
lm_speed_envcope     = lm(@formula(speed ~ env_temp + copepods), df3_home)
lm_speed_env         = lm(@formula(speed ~ env_temp), df3_home)
lm_speed_cope        = lm(@formula(speed ~ copepods), df3_home)
lm_r_speed           = lm(@formula(mean_r ~ speed), df3_home)
lm_K_speed           = lm(@formula(mean_k ~ speed), df3_home)

# ---------- 2. turning ----------
lm_r_envcope_turning = lm(@formula(mean_r ~ env_temp + copepods + turning), df3_home)
lm_r_env_turning     = lm(@formula(mean_r ~ env_temp + turning), df3_home)
lm_K_envcope_turning = lm(@formula(mean_k ~ env_temp + copepods + turning), df3_home)
lm_K_env_turning     = lm(@formula(mean_k ~ env_temp + turning), df3_home)
lm_turning_envcope   = lm(@formula(turning ~ env_temp + copepods), df3_home)
lm_turning_env       = lm(@formula(turning ~ env_temp), df3_home)
lm_turning_cope      = lm(@formula(turning ~ copepods), df3_home)
lm_r_turning         = lm(@formula(mean_r ~ turning), df3_home)
lm_K_turning         = lm(@formula(mean_k ~ turning), df3_home)

# ---------- 3. cellvol ----------
lm_r_envcope_cellvol = lm(@formula(mean_r ~ env_temp + copepods + cellvol), df3_home)
lm_r_env_cellvol     = lm(@formula(mean_r ~ env_temp + cellvol), df3_home)
lm_K_envcope_cellvol = lm(@formula(mean_k ~ env_temp + copepods + cellvol), df3_home)
lm_K_env_cellvol     = lm(@formula(mean_k ~ env_temp + cellvol), df3_home)
lm_cellvol_envcope   = lm(@formula(cellvol ~ env_temp + copepods), df3_home)
lm_cellvol_env       = lm(@formula(cellvol ~ env_temp), df3_home)
lm_cellvol_cope      = lm(@formula(cellvol ~ copepods), df3_home)
lm_r_cellvol         = lm(@formula(mean_r ~ cellvol), df3_home)
lm_K_cellvol         = lm(@formula(mean_k ~ cellvol), df3_home)

# ---------- 4. aspect_ratio ----------
lm_r_envcope_aspect  = lm(@formula(mean_r ~ env_temp + copepods + aspect_ratio), df3_home)
lm_r_env_aspect      = lm(@formula(mean_r ~ env_temp + aspect_ratio), df3_home)
lm_K_envcope_aspect  = lm(@formula(mean_k ~ env_temp + copepods + aspect_ratio), df3_home)
lm_K_env_aspect      = lm(@formula(mean_k ~ env_temp + aspect_ratio), df3_home)
lm_aspect_envcope    = lm(@formula(aspect_ratio ~ env_temp + copepods), df3_home)
lm_aspect_env        = lm(@formula(aspect_ratio ~ env_temp), df3_home)
lm_aspect_cope       = lm(@formula(aspect_ratio ~ copepods), df3_home)
lm_r_aspect          = lm(@formula(mean_r ~ aspect_ratio), df3_home)
lm_K_aspect          = lm(@formula(mean_k ~ aspect_ratio), df3_home)

# ---------- 5. Multivariate (main effects) ----------
lm_r_multi = lm(@formula(mean_r ~ cellvol + aspect_ratio + speed + turning), df3_home)
lm_K_multi = lm(@formula(mean_k ~ cellvol + aspect_ratio + speed + turning), df3_home)

# ---------- 6. Full pairwise interactions ----------
lm_r_interact = lm(@formula(mean_r ~ cellvol*aspect_ratio + cellvol*speed + cellvol*turning + 
                            aspect_ratio*speed + aspect_ratio*turning + speed*turning), df3_home)
lm_K_interact = lm(@formula(mean_k ~ cellvol*aspect_ratio + cellvol*speed + cellvol*turning + 
                            aspect_ratio*speed + aspect_ratio*turning + speed*turning), df3_home)

# ------------------------
# EXTRACT RESULTS
# ------------------------

all_results = DataFrame()

# ----- speed -----
append!(all_results, extract_model_results(lm_r_envcope_speed,   "speed", "r", "r ~ env_temp + copepods + speed"))
append!(all_results, extract_model_results(lm_r_env_speed,       "speed", "r", "r ~ env_temp + speed"))
append!(all_results, extract_model_results(lm_K_envcope_speed,   "speed", "K", "K ~ env_temp + copepods + speed"))
append!(all_results, extract_model_results(lm_K_env_speed,       "speed", "K", "K ~ env_temp + speed"))
append!(all_results, extract_model_results(lm_speed_envcope,     "speed", "speed", "speed ~ env_temp + copepods"))
append!(all_results, extract_model_results(lm_speed_env,         "speed", "speed", "speed ~ env_temp"))
append!(all_results, extract_model_results(lm_speed_cope,        "speed", "speed", "speed ~ copepods"))
append!(all_results, extract_model_results(lm_r_speed,           "speed", "r", "r ~ speed"))
append!(all_results, extract_model_results(lm_K_speed,           "speed", "K", "K ~ speed"))

# ----- turning -----
append!(all_results, extract_model_results(lm_r_envcope_turning, "turning", "r", "r ~ env_temp + copepods + turning"))
append!(all_results, extract_model_results(lm_r_env_turning,     "turning", "r", "r ~ env_temp + turning"))
append!(all_results, extract_model_results(lm_K_envcope_turning, "turning", "K", "K ~ env_temp + copepods + turning"))
append!(all_results, extract_model_results(lm_K_env_turning,     "turning", "K", "K ~ env_temp + turning"))
append!(all_results, extract_model_results(lm_turning_envcope,   "turning", "turning", "turning ~ env_temp + copepods"))
append!(all_results, extract_model_results(lm_turning_env,       "turning", "turning", "turning ~ env_temp"))
append!(all_results, extract_model_results(lm_turning_cope,      "turning", "turning", "turning ~ copepods"))
append!(all_results, extract_model_results(lm_r_turning,         "turning", "r", "r ~ turning"))
append!(all_results, extract_model_results(lm_K_turning,         "turning", "K", "K ~ turning"))

# ----- cellvol -----
append!(all_results, extract_model_results(lm_r_envcope_cellvol, "cellvol", "r", "r ~ env_temp + copepods + cellvol"))
append!(all_results, extract_model_results(lm_r_env_cellvol,     "cellvol", "r", "r ~ env_temp + cellvol"))
append!(all_results, extract_model_results(lm_K_envcope_cellvol, "cellvol", "K", "K ~ env_temp + copepods + cellvol"))
append!(all_results, extract_model_results(lm_K_env_cellvol,     "cellvol", "K", "K ~ env_temp + cellvol"))
append!(all_results, extract_model_results(lm_cellvol_envcope,   "cellvol", "cellvol", "cellvol ~ env_temp + copepods"))
append!(all_results, extract_model_results(lm_cellvol_env,       "cellvol", "cellvol", "cellvol ~ env_temp"))
append!(all_results, extract_model_results(lm_cellvol_cope,      "cellvol", "cellvol", "cellvol ~ copepods"))
append!(all_results, extract_model_results(lm_r_cellvol,         "cellvol", "r", "r ~ cellvol"))
append!(all_results, extract_model_results(lm_K_cellvol,         "cellvol", "K", "K ~ cellvol"))

# ----- aspect_ratio -----
append!(all_results, extract_model_results(lm_r_envcope_aspect,  "aspect_ratio", "r", "r ~ env_temp + copepods + aspect_ratio"))
append!(all_results, extract_model_results(lm_r_env_aspect,      "aspect_ratio", "r", "r ~ env_temp + aspect_ratio"))
append!(all_results, extract_model_results(lm_K_envcope_aspect,  "aspect_ratio", "K", "K ~ env_temp + copepods + aspect_ratio"))
append!(all_results, extract_model_results(lm_K_env_aspect,      "aspect_ratio", "K", "K ~ env_temp + aspect_ratio"))
append!(all_results, extract_model_results(lm_aspect_envcope,    "aspect_ratio", "aspect_ratio", "aspect_ratio ~ env_temp + copepods"))
append!(all_results, extract_model_results(lm_aspect_env,        "aspect_ratio", "aspect_ratio", "aspect_ratio ~ env_temp"))
append!(all_results, extract_model_results(lm_aspect_cope,       "aspect_ratio", "aspect_ratio", "aspect_ratio ~ copepods"))
append!(all_results, extract_model_results(lm_r_aspect,          "aspect_ratio", "r", "r ~ aspect_ratio"))
append!(all_results, extract_model_results(lm_K_aspect,          "aspect_ratio", "K", "K ~ aspect_ratio"))

# ----- Multivariate (main effects) -----
append!(all_results, extract_model_results(lm_r_multi, "multivariate", "r", "r ~ cellvol + aspect_ratio + speed + turning"))
append!(all_results, extract_model_results(lm_K_multi, "multivariate", "K", "K ~ cellvol + aspect_ratio + speed + turning"))

# ----- Full pairwise interactions -----
append!(all_results, extract_model_results(lm_r_interact, "interactions", "r", 
        "r ~ cellvol*aspect_ratio + cellvol*speed + cellvol*turning + aspect_ratio*speed + aspect_ratio*turning + speed*turning"))
append!(all_results, extract_model_results(lm_K_interact, "interactions", "K", 
        "K ~ cellvol*aspect_ratio + cellvol*speed + cellvol*turning + aspect_ratio*speed + aspect_ratio*turning + speed*turning"))

# ------------------------
# CLEAN UP AND EXPORT
# ------------------------

rename!(all_results, 
    :StdError => Symbol("Std. Error"),
    :Pr_t => Symbol("Pr(>|t|)"),
    :Lower95 => Symbol("Lower 95%"),
    :Upper95 => Symbol("Upper 95%")
)

XLSX.writetable("regression_results.xlsx", all_results, overwrite=true)

println("Results exported to regression_results.xlsx")
println("Total rows: ", nrow(all_results))


savefig("p2 /Users/thisisme710/Desktop/Research/Paramecium logistic growth/Fall 2025/plots")
println("Plot saved!")

plot_8pheno, p7s,p2, plot_r_final

