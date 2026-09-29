# ============================================================================
# PARAMECIUM LOGISTIC GROWTH - REORGANIZED ANALYSIS SCRIPT
# ============================================================================

using CSV
using DataFrames
using GLM
using StatsPlots
using Plots
using Statistics
using MixedModels
using MultivariateStats
using Biplots
using Tables

# use the dataset with all replicates, not means
df2_home = CSV.read("rK_d2_home.csv", DataFrame)
#df2_home_labeled = copy(df2_home)
df2_home.copepods_label = ifelse.(df2_home.copepods .== 0, "No Copepods", "Copepods")
df2_home.temp_label = ifelse.(df2_home.env_temp .== 20, "20°C", "26°C")

# ============================================================================
# Using mixed models with the population as a random effect
# ============================================================================
# the base model is just the design of the experiment
# life history parameters r and K
mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * copepods + (1|names)), df2_home)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * copepods + (1|names)), df2_home)
# cell morphology and behavior
mm_cellvol_int = fit(MixedModel, @formula(cell_vol ~ env_temp * copepods + (1|names)), df2_home)
mm_speed_int = fit(MixedModel, @formula(speed ~ env_temp * copepods + (1|names)), df2_home)
mm_aspectratio_int = fit(MixedModel, @formula(aspect_ratio ~ env_temp * copepods + (1|names)), df2_home)
mm_turning_int = fit(MixedModel, @formula(sd_turning ~ env_temp * copepods + (1|names)), df2_home)
mm_major_int = fit(MixedModel, @formula(mean_major ~ env_temp * copepods + (1|names)), df2_home)
mm_minor_int = fit(MixedModel, @formula(mean_minor ~ env_temp * copepods + (1|names)), df2_home)

# pull out model tables and put them together in a dataframe
main_mm_results = DataFrame()

# append model by model, adding the dependent variable name as we go
append!(main_mm_results,hcat(DataFrame(depvar = fill("r", nrow(ct))),
    DataFrame(coeftable(mm_r_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("k", nrow(ct))),
    DataFrame(coeftable(mm_k_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("cell volume", nrow(ct))),
    DataFrame(coeftable(mm_cellvol_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("speed", nrow(ct))),
    DataFrame(coeftable(mm_speed_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("aspect ratio", nrow(ct))),
    DataFrame(coeftable(mm_aspectratio_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("sd turning", nrow(ct))),
    DataFrame(coeftable(mm_turning_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("major axis", nrow(ct))),
    DataFrame(coeftable(mm_major_int))))
append!(main_mm_results,hcat(DataFrame(depvar = fill("minor axis", nrow(ct))),
    DataFrame(coeftable(mm_minor_int))))
    
CSV.write("main_mm_results.csv", main_mm_results)

# ============================================================================
# make some plotting functions
# ============================================================================

# Global Plot Settings
default(titlefont=10, guidefont=12, tickfont=10, legendfont=15, margin=3mm)
cop_palette = [:yellow, :brown]
temp_palette = [:blue, :purple]

# Helper: Extract p-value
function get_pvalue(model, term_index)
    return coeftable(model).cols[4][term_index]
end

# Helper to make treatment boxplots
function make_effect_plot(y_sym, x_sym, pal, y_lab)
    
    @df df2_home_labeled boxplot(cols(x_sym), cols(y_sym),
        group=cols(x_sym), palette=pal, fillalpha=0.75, linewidth=1.5,
        #title="$title_str\np=$(round(pval, digits=3))", 
        ylabel=y_lab, xlabel="Treatment", legend=false)
    @df df2_home_labeled dotplot!(cols(x_sym), cols(y_sym),
            marker=(:black, 0.4),
            mode=:density)

end

# ----------------------------------------------------------------------------
# FIGURE SET 1: r and K BOXPLOTS
# ----------------------------------------------------------------------------
# make the subpanels
# 1. r vs Copepods
p1_1 = make_effect_plot(:fitted_r, :copepods_label, cop_palette, "r")
# 2. K vs Copepods
p1_2 = make_effect_plot(:fitted_k, :copepods_label, cop_palette, "K")
# 3. r vs Temperature
p1_3 = make_effect_plot(:fitted_r, :temp_label, temp_palette, "r")
# 4. K vs Temperature
p1_4 = make_effect_plot(:fitted_k, :temp_label, temp_palette, "K")

# put them together
fig_set_1 = plot(p1_1, p1_3, p1_2, p1_4, layout=(2,2), size=(800, 800))
savefig(fig_set_1, "Figure1_MainEffects.png")

# ----------------------------------------------------------------------------
# FIGURE SET 2: PHENOTYPE BOXPLOTS (8 panels - 4×2 grid)
# ----------------------------------------------------------------------------
println("Generating Figure Set 3: Phenotype Boxplots...")

# Define traits and units
traits = [
    (:cell_vol, "Cell Volume (µm³)"),
    (:aspect_ratio, "Aspect Ratio"),
    (:speed, "Speed (µm/s)"),
    (:sd_turning, "Turning (rad)")
]

# We will generate pairs: (Trait vs Cop) and (Trait vs Temp) for each trait
plots_set3 = []

for (trait, label) in traits
    # vs Copepods
    p_cop = make_effect_plot(trait, :copepods_label, cop_palette, label)
    push!(plots_set3, p_cop)
    
    # vs Temperature
    p_temp = make_effect_plot(trait, :temp_label, temp_palette, label)
    push!(plots_set3, p_temp)
end

# Layout: 4 rows (traits), 2 columns (treatments)
# The loop pushes [Row1-Col1, Row1-Col2, Row2-Col1...] which matches plot(... layout=(4,2))
fig_set_3 = plot(plots_set3..., layout=(4,2), size=(1000, 1200), margin=10mm)
savefig(fig_set_3, "Figure3_Phenotypes.png")

# ----------------------------------------------------------------------------
# PCA on phenotypes
# ----------------------------------------------------------------------------

# select columns from dataframe
df_pca = df2_home[:, [:cell_vol, :speed, :aspect_ratio, :sd_turning]]
subset_names = ["Cell volume", "Speed", "Aspect ratio", "SD of turning"]
X = Matrix(df_pca)
cor(X)
# scale and center
Xz = (X .- mean(X, dims=1)) ./ std(X, dims=1)
cor(Xz)

# run the pca
M = fit(PCA, Xz'; maxoutdim=5)
principalvars(M)
loadings = projection(M)
scores = predict(M, Xz')
df2_home[!,:PCA1] .= scores[1,:]
df2_home[!,:PCA2] .= scores[2,:]

# ----------------------------------------------------------------------------
# make the PCA biplot
# ----------------------------------------------------------------------------

# Scale loadings so arrows are visible
scale = maximum(abs.(vcat(pc1, pc2))) * 0.7

s1 = scatter(df2_home.PCA1[findall(df2_home.copepods .== 0)],df2_home.PCA2[findall(df2_home.copepods .== 0)],
    group=df2_home.temp_label[findall(df2_home.copepods .== 0)],
    shape=:circle,
    color=[:blue :purple],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8,
    ylim = (-3, 2))
scatter!(df2_home.PCA1[findall(df2_home.copepods .== 1)],df2_home.PCA2[findall(df2_home.copepods .== 1)],
    group=df2_home.temp_label[findall(df2_home.copepods .== 1)],
    shape=:pent,
    color=[:yellow :brown],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8)

    s1

# first panel, colored by copepod treatment
s1 = @df df2_home scatter(:PCA1(findall(copepods .== 0)),:PCA2(findall(copepods .== 0)),
    group=:copepods_label(findall(copepods .== 0)),
    color=[:yellow],
    shape=[:square],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8)

@df df2_home scatter(:PCA1,:PCA2,
    group=:temp_label,
    shape=[:square :star],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8)


# second panel, colored by temperature treatment
s2 = @df df2_home scatter(:PCA1,:PCA2,
    group=:temp_label,
    color=[:blue :purple],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8)

current()

for i in 1:size(loadings, 1)
    x = loadings[i, 1] * scale
    y = loadings[i, 2] * scale
    plot!(
    [0, x],
    [0, y],
    arrow=true,
    color=:black,
    label=""
    )
annotate!(x, y, string(subset_names[i]))
end

fig_set_pca = plot(s1, s2, layout=(1,2), size=(800,400))
savefig(fig_set_pca, "Figure3_PCA.png")

# ----------------------------------------------------------------------------
# Test if PCA scores predict fitted r and K
# ----------------------------------------------------------------------------

# first, do the treatments affect the PCAs?
mm_r_int = fit(MixedModel, @formula(PCA1 ~ env_temp * copepods + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(PCA1 ~ env_temp * copepods + (1|names)), df2_home)


mm_r_int = fit(MixedModel, @formula(PCA2 ~ env_temp * copepods + (1|names)), df2_home)

# second, do the PCAs affect r and K?
mm_r_int = fit(MixedModel, @formula(fitted_r ~ PCA1 * PCA2 + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(fitted_r ~ PCA1 + PCA2 + (1|names)), df2_home)

mm_r_int = fit(MixedModel, @formula(fitted_k ~ PCA1 * PCA2 + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(fitted_k ~ PCA1 + PCA2 + (1|names)), df2_home)

# ----------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------



mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ cell_vol + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ speed + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ aspect_ratio + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ sd_turning + (1|names)), df2_home)

mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ cell_vol + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ speed + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ aspect_ratio + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(log(fitted_r) ~ sd_turning + (1|names)), df2_home)



mm_r_int = fit(MixedModel, @formula(fitted_k ~ cell_vol + aspect_ratio + speed + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(fitted_r ~ cell_vol + speed + (1|names)), df2_home)
mm_r_int = fit(MixedModel, @formula(fitted_r ~ aspect_ratio + (1|names)), df2_home)


mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * copepods + (1|names)), df2_home)
# cell morphology and behavior
mm_cellvol_int = fit(MixedModel, @formula(cell_vol ~ env_temp * copepods + (1|names)), df2_home)
mm_speed_int = fit(MixedModel, @formula(speed ~ env_temp * copepods + (1|names)), df2_home)
mm_aspectratio_int = fit(MixedModel, @formula(aspect_ratio ~ env_temp * copepods + (1|names)), df2_home)
mm_turning_int = fit(MixedModel, @formula(sd_turning ~ env_temp * copepods + (1|names)), df2_home)
mm_major_int = fit(MixedModel, @formula(mean_major ~ env_temp * copepods + (1|names)), df2_home)
mm_minor_int = fit(MixedModel, @formula(mean_minor ~ env_temp * copepods + (1|names)), df2_home)



