# ============================================================================
# PARAMECIUM LOGISTIC GROWTH - REORGANIZED ANALYSIS SCRIPT
# ============================================================================

# Kristi is learning Github
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
using LaTeXStrings


dir = "C:/Users/93652672/Documents/GitHub/Paramecium-temperature-predator-ORCC-/rk temperature predator responses/"
cd(dir)
pwd() # check your pathway
readdir()

# use the dataset with all replicates, not means
dfA = CSV.read("rK_full_fitted_dataset.csv",DataFrame)
dfA.copepods_label = ifelse.(dfA.copepods .== 0, "No Copepods", "Copepods")
dfA.temp_label = ifelse.(dfA.env_temp .== 20, "20°C", "26°C")
dfA.home_label = ifelse.(dfA.athome .== 1, "Home", "Away")

# note the 26COP04's have already been deleted from this data set

# ============================================================================
# First plot the r's and K's across all treatments, home and away
# ============================================================================

# these are traditional reaction norm plots
# to match up colors with previous plots, going to do this for each treatment

# open an empty plot with dimensions and axis labels
p4 = plot(xlims=[18,28],
    xlabel = L"Temperature ($^{o}C$)",
    ylabel = L"Intrinsic rate of growth ($r$)",
    size=(500,500),
    xticks = [20,26],
    dpi=1000)
p5 = plot(xlims=[18,28],
    xlabel = L"Temperature ($^{o}C$)",
    ylabel = L"Carrying capacity ($K$)",
    size=(500,500),
    xticks = [20,26],
    dpi=1000)

function means_and_std(temp_2_plot)
    grp = groupby(temp_2_plot, [:athome]) # taking the mean across tech reps and pop reps
    rx_norm_means = combine(grp, [:fitted_k, :fitted_r] .=> mean,
        [:fitted_k, :fitted_r] .=> std) # gives us the mean and the std by groups in grp
    return rx_norm_means
end
# start with 20C and no copepods
    temp_2_plot = dfA[findall(dfA.env_temp .== 20 .&& dfA.copepods .== 0),:]
    rx_norm_means = means_and_std(temp_2_plot)

    yerror = rx_norm_means[:,5] # designate what data is the error
    plot!(p4,[19.5,25.5],rx_norm_means[:,3];yerror,
        label="",color=gradient[2])
    scatter!(p4,[19.5,25.5],rx_norm_means[:,3],color=gradient[2],
        label="20C_nopreds",markershape=collect([:circle, :rect]),markersize = 8)

    yerror = rx_norm_means[:,4] # designate what data is the error
    plot!(p5,[19.5,25.5],rx_norm_means[:,2];yerror,
        label="",color=gradient[2])
    scatter!(p5,[19.5,25.5],rx_norm_means[:,2],color=gradient[2],
        label="",markershape=collect([:circle, :rect]),markersize = 8)

# now 20C and copepods
    temp_2_plot = dfA[findall(dfA.env_temp .== 20 .&& dfA.copepods .== 1),:]
    rx_norm_means = means_and_std(temp_2_plot)

    yerror = rx_norm_means[:,5] # designate what data is the error
    plot!(p4,[19.75,25.75],rx_norm_means[:,3];yerror,
        label="",color=gradient[4])
    scatter!([19.75,25.75],rx_norm_means[:,3],color=gradient[4],
        label="20C_preds",markershape=collect([:circle, :rect]),markersize = 8)

    yerror = rx_norm_means[:,4] # designate what data is the error
    plot!(p5,[19.75,25.75],rx_norm_means[:,2];yerror,
        label="",color=gradient[4])
    scatter!(p5,[19.75,25.75],rx_norm_means[:,2],color=gradient[4],
        label="",markershape=collect([:circle, :rect]),markersize = 8)

#now 26C and no copepods
    temp_2_plot = dfA[findall(dfA.env_temp .== 26 .&& dfA.copepods .== 0),:]
    rx_norm_means = means_and_std(temp_2_plot)

    yerror = rx_norm_means[:,5] # designate what data is the error
    plot!(p4,[20.25,26.25],rx_norm_means[:,3];yerror,
        label="",color=gradient[9])
    scatter!([20.25,26.25],rx_norm_means[:,3],color=gradient[9],   
        label="26C_nopreds",markershape=collect([:rect, :circle]),markersize = 8)

    yerror = rx_norm_means[:,4] # designate what data is the error
    plot!(p5,[20.25,26.25],rx_norm_means[:,2];yerror,
        label="",color=gradient[9])
    scatter!(p5,[20.25,26.25],rx_norm_means[:,2],color=gradient[9],
        label="",markershape=collect([:rect, :circle]),markersize = 8)

#now 26C and copepods
    temp_2_plot = dfA[findall(dfA.env_temp .== 26 .&& dfA.copepods .== 1),:]
    rx_norm_means = means_and_std(temp_2_plot)

    yerror = rx_norm_means[:,5] # designate what data is the error
    plot!(p4,[20.5,26.5],rx_norm_means[:,3];yerror,
        label="",color=gradient[7])
    scatter!([20.5,26.5],rx_norm_means[:,3],color=gradient[7],   
        label="26C_preds",markershape=collect([:rect, :circle]),markersize = 8)

    yerror = rx_norm_means[:,4] # designate what data is the error
    plot!(p5,[20.25,26.25],rx_norm_means[:,2];yerror,
        label="",color=gradient[7])
    scatter!([20.25,26.25],rx_norm_means[:,2],color=gradient[7],   
        label="",markershape=collect([:rect, :circle]),markersize = 8)

figure5 = plot(p4,p5, layout=(1,2), size=(600,300), dpi=800)
savefig(figure5,"reaction_norms.png")
        
# ============================================================================
# First, is there a difference between r and K at homes
# ============================================================================
# Use mixed models with the population as a random effect
# the base model is just the design of the experiment
# life history parameters r and K

df_home = dfA[findall(dfA.athome .== 1),:]

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * copepods + (1|names)), df_home)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * copepods + (1|names)), df_home)

# ============================================================================
# Second, is there a difference between r and K between home and away
# ============================================================================
# because copepod was a prior treatment, maybe it's best to separate these out by copepod treatment
gbye_copepods = dfA[findall(dfA.copepods .== 1),:]

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * athome + (1|names)), gbye_copepods)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * athome + (1|names)), gbye_copepods)

gbye_nocopepods = dfA[findall(dfA.copepods .== 0),:]

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * athome + (1|names)), gbye_nocopepods)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * athome + (1|names)), gbye_nocopepods)







##### figure 1 - r versus K ##########
# data by treatment with all replicate subsamples
# 26cop04's are deleted
gradient = cgrad(:berlin,9,categorical = true) # create a color gradient
p1 = plot(ylims=[0,2e4],
    xlabel = "Intrinsic rate of growth, r",
    ylabel = "Carrying capacity, K",
    size=(700,700),
    dpi=1000,
    title="Showing replicate subsamples")
scatter!(p1,dfA.fitted_r[findall(dfA.env_temp .== 20 .&& dfA.copepods .== 0)],
    dfA.fitted_k[findall(dfA.env_temp .== 20 .&& dfA.copepods .== 0)],
    markerstrokecolor=gradient[2],markercolor=gradient[2],
    markersize=6,
    markershape=:square,
    label="20C_nopreds")


    

# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================



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



