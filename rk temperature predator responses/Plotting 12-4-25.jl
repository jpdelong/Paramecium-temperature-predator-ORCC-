# ============================================================================
# PARAMECIUM LOGISTIC GROWTH - REORGANIZED PLOTTING SCRIPT
# ============================================================================
# Figure Organization:
# Set 1: Main Effects (r/K vs Treatments) [2x2]
# Set 2: Growth Curves (Fits vs Observed) [2x2]
# Set 3: Phenotypes (Traits vs Treatments) [4x2]
# Set 4: Correlations (r/K vs Traits) [4x2]
# ============================================================================

# ----------------------------------------------------------------------------
# 1. SETUP & DATA LOADING
# ----------------------------------------------------------------------------
#dir = "/Users/thisisme710/Desktop/Research/Paramecium logistic growth/Fall 2025/Data Sheets"

dir = "/Users/93652672/Documents/GitHub/Paramecium-temperature-predator-ORCC-/rk temperature predator responses"
cd(dir)

using CSV, DataFrames, GLM, StatsPlots, Plots, Statistics, Measures

# Load Data
df3_home = CSV.read("rK_d3_home.csv", DataFrame)
fitted_curves = CSV.read("fitted_curves_predicted.csv", DataFrame)
observed_data = CSV.read("fitted_curves_observed.csv", DataFrame)

# Pre-processing for Labels
df3_home_labeled = copy(df3_home)
df3_home_labeled.copepods_label = ifelse.(df3_home_labeled.copepods .== 0, "No Copepods", "Copepods")
df3_home_labeled.temp_label = ifelse.(df3_home_labeled.env_temp .== 20, "20°C", "26°C")

# Helper: Extract p-value
function get_pvalue(model, term_index=2)
    try
        return coeftable(model).cols[4][term_index]
    catch
        return NaN
    end
end

# Global Plot Settings
default(titlefont=10, guidefont=9, tickfont=8, legendfont=7, margin=3mm)
cop_palette = [:yellow, :brown]
temp_palette = [:blue, :purple]

# ----------------------------------------------------------------------------
# FIGURE SET 1: MAIN EFFECTS (4 panels - 2×2 grid)
# ----------------------------------------------------------------------------
println("Generating Figure Set 1: Main Effects...")

# Helper to make treatment boxplots
function make_effect_plot(y_sym, x_sym, group_sym, pal, title_str, y_lab)
    formula = Term(y_sym) ~ Term(x_sym)
    model = lm(formula, df3_home)
    pval = get_pvalue(model)
    
    @df df3_home_labeled boxplot(cols(group_sym), cols(y_sym),
        group=cols(group_sym), palette=pal, fillalpha=0.75, linewidth=1.5,
        title="$title_str\np=$(round(pval, digits=2))", 
        ylabel=y_lab, legend=false)
end

# 1. r vs Copepods
p1_1 = make_effect_plot(:mean_r, :copepods, :copepods_label, cop_palette, "r vs Copepods", "r")
# 2. K vs Copepods
p1_2 = make_effect_plot(:mean_k, :copepods, :copepods_label, cop_palette, "K vs Copepods", "K")
# 3. r vs Temperature
p1_3 = make_effect_plot(:mean_r, :env_temp, :temp_label, temp_palette, "r vs Temperature", "r")
# 4. K vs Temperature
p1_4 = make_effect_plot(:mean_k, :env_temp, :temp_label, temp_palette, "K vs Temperature", "K")

fig_set_1 = plot(p1_1, p1_2, p1_3, p1_4, layout=(2,2), size=(800, 800))




# ----------------------------------------------------------------------------
# FIGURE SET 2: LOGISTIC GROWTH CURVES (4 panels - 2×2 grid)
# ----------------------------------------------------------------------------
println("Generating Figure Set 2: Logistic Growth Curves...")

function plot_growth_panel(temp, cop, title_str)
    # Filter data
    pred_sub = filter(row -> row.Temperature == temp && row.Copepods == cop, fitted_curves)
    obs_sub = filter(row -> row.Temperature == temp && row.Copepods == cop, observed_data)
    
    # Calculate stats
    if nrow(pred_sub) > 0
        unique_dishes = unique(pred_sub.DishID)
        avg_r = mean([pred_sub[pred_sub.DishID .== d, :Fitted_r][1] for d in unique_dishes])
        avg_k = mean([pred_sub[pred_sub.DishID .== d, :Fitted_K][1] for d in unique_dishes])
        subtitle = "Avg r=$(round(avg_r, digits=2)), Avg K=$(round(avg_k, digits=0))"
    else
        subtitle = "No Data"
        unique_dishes = []
    end

    p = plot(title="$title_str\n$subtitle", xlabel="Day", ylabel="Abundance", legend=false)
    
    colors = [:green, :blue, :purple, :red, :orange, :cyan, :magenta, :brown]
    
    for (i, dish) in enumerate(unique_dishes)
        d_pred = filter(row -> row.DishID == dish, pred_sub)
        d_obs = filter(row -> row.DishID == dish, obs_sub)
        c = colors[mod1(i, length(colors))]
        
        plot!(p, d_pred.Day, d_pred.Predicted_Abundance, color=c, lw=2, alpha=0.8)
        scatter!(p, d_obs.Day, d_obs.Observed_Abundance, color=c, ms=3, alpha=0.8)
    end
    return p
end

p2_1 = plot_growth_panel(20, 1, "20°C with Predators")
p2_2 = plot_growth_panel(20, 0, "20°C without Predators")
p2_3 = plot_growth_panel(26, 1, "26°C with Predators")
p2_4 = plot_growth_panel(26, 0, "26°C without Predators")

fig_set_2 = plot(p2_1, p2_2, p2_3, p2_4, layout=(2,2), size=(900, 900))

# ----------------------------------------------------------------------------
# FIGURE SET 3: PHENOTYPE BOXPLOTS (8 panels - 4×2 grid)
# ----------------------------------------------------------------------------
println("Generating Figure Set 3: Phenotype Boxplots...")

# Define traits and units
traits = [
    (:cellvol, "Cell Volume (µm³)"),
    (:aspect_ratio, "Aspect Ratio"),
    (:speed, "Speed (µm/s)"),
    (:turning, "Turning (rad)")
]

# We will generate pairs: (Trait vs Cop) and (Trait vs Temp) for each trait
plots_set3 = []

for (trait, label) in traits
    # vs Copepods
    p_cop = make_effect_plot(trait, :copepods, :copepods_label, cop_palette, "$label vs Cop", label)
    push!(plots_set3, p_cop)
    
    # vs Temperature
    p_temp = make_effect_plot(trait, :env_temp, :temp_label, temp_palette, "$label vs Temp", label)
    push!(plots_set3, p_temp)
end

# Layout: 4 rows (traits), 2 columns (treatments)
# The loop pushes [Row1-Col1, Row1-Col2, Row2-Col1...] which matches plot(... layout=(4,2))
fig_set_3 = plot(plots_set3..., layout=(4,2), size=(1000, 1200), margin=4mm)

# ----------------------------------------------------------------------------
# FIGURE SET 4: SCATTER PLOTS (8 panels - 4×2 grid)
# ----------------------------------------------------------------------------
println("Generating Figure Set 4: Scatter Plots...")

function make_scatter(x_sym, y_sym, x_lab, y_lab)
    # Linear Model
    formula = Term(y_sym) ~ Term(x_sym)
    model = lm(formula, df3_home)
    pval = get_pvalue(model)
    r2_val = r2(model)

    # Plot
    p = scatter(df3_home[!, x_sym], df3_home[!, y_sym], 
        xlabel=x_lab, ylabel=y_lab, legend=false, 
        markerstrokewidth=0, alpha=0.6, color=:black)
    
    # Add regression line
    plot!(p, df3_home[!, x_sym], predict(model), color=:red, linewidth=2)
    title!(p, "$y_lab vs $x_lab\np=$(round(pval, digits=3)), R²=$(round(r2_val, digits=2))")
    return p
end

# We want 4 rows (traits) x 2 columns (r and K)
plots_set4 = []

# Order: CellVol, Aspect, Turning, Speed (to match Set 3 order for consistency)
# Note: User request list order was specific, but row-by-row trait comparison is usually clearer.
# Below generates: Row 1 (r vs CV, K vs CV), Row 2 (r vs AR, K vs AR)...

traits_scatter = [
    (:cellvol, "Cell Vol"),
    (:aspect_ratio, "Aspect Ratio"),
    (:turning, "Turning"),
    (:speed, "Speed")
]

for (trait, label) in traits_scatter
    # Column 1: r vs Trait
    push!(plots_set4, make_scatter(trait, :mean_r, label, "r"))
    
    # Column 2: K vs Trait
    push!(plots_set4, make_scatter(trait, :mean_k, label, "K"))
end

fig_set_4 = plot(plots_set4..., layout=(4,2), size=(1000, 1200), margin=4mm)

# ----------------------------------------------------------------------------
# 5. SAVE FILES
# ----------------------------------------------------------------------------
println("Saving figures to PDF...")

savefig(fig_set_1, "Figure1_MainEffects.pdf")
savefig(fig_set_2, "Figure2_GrowthCurves.pdf")
savefig(fig_set_3, "Figure3_Phenotypes.pdf")
savefig(fig_set_4, "Figure4_ScatterPlots.pdf")

println("\n✓ Done! Created 4 Figure Sets:")
println("  1. Figure1_MainEffects.pdf (r/K Boxplots)")
println("  2. Figure2_GrowthCurves.pdf (Logistic Fits)")
println("  3. Figure3_Phenotypes.pdf (Traits vs Treatments)")
println("  4. Figure4_ScatterPlots.pdf (r/K vs Traits)")