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

# ============================================================================
# 1. Read in the data
# ============================================================================

# use the dataset with all replicates, not means
dfA = CSV.read("rK_full_fitted_dataset.csv",DataFrame)
dfA.copepods_label = ifelse.(dfA.copepods .== 0, "No Copepods", "Copepods")
dfA.temp_label = ifelse.(dfA.env_temp .== 20, "20°C", "26°C")
dfA.home_label = ifelse.(dfA.athome .== 1, "Home", "Away")

pop_names = unique(dfA.names)

# note the 26COP04's have already been deleted from this data set

# ============================================================================
# 2. Make a PCA on phenotypes
# ============================================================================

# select columns from dataframe
df_pca = dfA[:, [:cell_vol, :speed, :aspect_ratio, :sd_turning]]
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
dfA[!,:PCA1] .= scores[1,:]
dfA[!,:PCA2] .= scores[2,:]

# ============================================================================
# 3. Make the home dataset
# ============================================================================

df_home = dfA[findall(dfA.athome .== 1),:]

# ============================================================================
# 4. Make the PCA biplot on the home data
# ============================================================================

gradient = cgrad(:berlin,9,categorical = true) # create a color gradient

# Scale loadings so arrows are visible
scale = maximum(abs.(vcat(df_home.PCA1, df_home.PCA2))) * 0.7

s1 = scatter(df_home.PCA1[findall(df_home.copepods .== 0)],df_home.PCA2[findall(df_home.copepods .== 0)],
    group=df_home.temp_label[findall(df_home.copepods .== 0)],
    shape=:circle,
    color=[gradient[2] gradient[9]],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8,
    label = ["20°C no predators" "26°C no predators"],
    ylim = (-3, 2))
scatter!(df_home.PCA1[findall(df_home.copepods .== 1)],df_home.PCA2[findall(df_home.copepods .== 1)],
    group=df_home.temp_label[findall(df_home.copepods .== 1)],
    shape=:pent,
    color=[gradient[4] gradient[7]],
    xlabel="PC1",
    ylabel="PC2",
    markersize = 8,
    legendfontsize = 8,
    label = ["20°C predators" "26°C predators"])

current()

for i in 1:size(loadings, 1)
    x = loadings[i, 1] * scale
    y = loadings[i, 2] * scale
    plot!([0, x],[0, y],arrow=true,color=:black,label="")

    annotate!(x, y, string(subset_names[i]))
end

#fig_set_pca = plot(s1, s2, layout=(1,2), size=(800,400))
fig_set_pca = plot(s1, size=(400,400),
    legend=:bottomright)

savefig(fig_set_pca, "Figure2_PCA.png")

# ============================================================================
# 5. Plot the r's and K's across all treatments, home and away
# ============================================================================

# these are traditional reaction norm plots
# to match up colors with previous plots, going to do this for each treatment

# open an empty plot with dimensions and axis labels
p4 = plot(xlims=[18,28],
    xlabel = L"Temperature ($^{o}C$)",
    ylabel = L"Intrinsic rate of growth ($r$)",
    size=(600,500),
    xticks = [20,26],
    dpi=600,
    margin = 5Plots.mm,
    legend=:outertopright)
p5 = plot(xlims=[18,28],
    xlabel = L"Temperature ($^{o}C$)",
    ylabel = L"Carrying capacity ($K$)",
    size=(600,500),
    xticks = [20,26],
    dpi=600,
    margin = 5Plots.mm)

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

#figure5 = plot(p5,p4, layout=(1,2), size=(700,300), dpi=600)
figure5 = plot(p5,p4, layout=@layout([a{0.37w} b{0.63w}]), size=(700,300), dpi=600)

savefig(figure5,"reaction_norms.png")
        
# ============================================================================
# 6. Test for effect of treatments on r and K at homes
# ============================================================================
# Use mixed models with the population as a random effect
# this model is just the design of the experiment
# life history parameters r and K

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * copepods + (1|names)), df_home)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * copepods + (1|names)), df_home)

# ============================================================================
# 7. Test for a difference between r and K between home and away
# ============================================================================

# because copepod was a prior treatment, maybe it's best to separate these out by copepod treatment
gbye_copepods = dfA[findall(dfA.copepods .== 1),:]

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * athome + (1|names)), gbye_copepods)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * athome + (1|names)), gbye_copepods)

gbye_nocopepods = dfA[findall(dfA.copepods .== 0),:]

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * athome + (1|names)), gbye_nocopepods)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * athome + (1|names)), gbye_nocopepods)

# go ahead and put everything in together
mm_k_int = fit(MixedModel, @formula(fitted_r ~ env_temp * athome * copepods + (1|names)), dfA)
mm_k_int = fit(MixedModel, @formula(fitted_k ~ env_temp * athome * copepods + (1|names)), dfA)

# ============================================================================
# 8. Test if PCA scores predict fitted r and K
# ============================================================================

# first, do the treatments affect the PCAs?
mm_r_int = fit(MixedModel, @formula(PCA1 ~ env_temp * copepods + (1|names)), df_home)
mm_r_int = fit(MixedModel, @formula(PCA2 ~ env_temp * copepods + (1|names)), df_home)

# second, do the PCAs affect r and K?
mm_r_int = fit(MixedModel, @formula(fitted_r ~ PCA1 * PCA2 + (1|names)), df_home)
mm_r_int = fit(MixedModel, @formula(fitted_r ~ PCA1 + PCA2 + (1|names)), df_home)

mm_r_int = fit(MixedModel, @formula(fitted_k ~ PCA1 * PCA2 + (1|names)), df_home)
mm_r_int = fit(MixedModel, @formula(fitted_k ~ PCA1 + PCA2 + (1|names)), df_home)

mm_r_int = fit(MixedModel, @formula(fitted_r ~ env_temp * copepods + (1|names)), df_home)

# ============================================================================
# 9. CREATE HOME/AWAY COMPARISON DATASET (df3)
# ============================================================================
num_rows = 15

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
pc1_home = zeros(Float64, 1, num_rows)
pc1_away = zeros(Float64, 1, num_rows)
pc2_home = zeros(Float64, 1, num_rows)
pc2_away = zeros(Float64, 1, num_rows)
k_diff = zeros(Float64, 1, num_rows)
r_diff = zeros(Float64, 1, num_rows)
zero_vec = zeros(Float64, 1, num_rows)

# Calculate home/away comparisons
for i = 1:15
    indices_home = findall(dfA.names .== pop_names[i] .&& dfA.athome .== 1)
    indices_away = findall(dfA.names .== pop_names[i] .&& dfA.athome .== 0)
    
    if isempty(indices_home) || isempty(indices_away)
        @warn "No home/away data for $(pop_names[i])"
        continue
    end
    
    # HOME
    k_at_home[i] = mean(dfA.fitted_k[indices_home])
    r_at_home[i] = mean(dfA.fitted_r[indices_home])
    env_temp_home[i] = dfA.env_temp[indices_home][1]
    copepods_home[i] = dfA.copepods[indices_home][1]
    pc1_home[i] = dfA.cell_vol[indices_home][1]
    pc2_home[i] = dfA.aspect_ratio[indices_home][1]
    home_home[i] = 1
    
    # AWAY
    k_away[i] = mean(dfA.fitted_k[indices_away])
    r_away[i] = mean(dfA.fitted_r[indices_away])
    env_temp_away[i] = dfA.env_temp[indices_away][1]
    copepods_away[i] = dfA.copepods[indices_away][1]
    pc1_away[i] = dfA.cell_vol[indices_away][1]
    pc2_away[i] = dfA.aspect_ratio[indices_away][1]
    home_away[i] = 0
    
    # DIFFERENCES
    k_diff[i] = k_away[i] - k_at_home[i]
    r_diff[i] = r_away[i] - r_at_home[i]
end

# Construct df3
df3 = DataFrame(ID = 1:30)
df3.mean_k = vcat(vec(k_at_home), vec(k_away))
df3.mean_r = vcat(vec(r_at_home), vec(r_away))
df3.pc1 = vcat(vec(pc1_home), vec(pc1_away))
df3.pc2 = vcat(vec(pc2_home), vec(pc2_away))
df3.home = vcat(vec(home_home), vec(home_away))
df3.copepods = vcat(vec(copepods_home), vec(copepods_away))
df3.env_temp = vcat(vec(env_temp_home), vec(env_temp_away))
names15 = collect(first(pop_names, num_rows))
df3.names = vcat(names15, names15)
df3.r_diffs = vcat(vec(r_diff), vec(zero_vec))
df3.k_diffs = vcat(vec(k_diff), vec(zero_vec))

# Create home subset
df3_home = df3[findall(df3.home .== 1), :]

# ============================================================================
# 10. Make r and K home-away change vector diagram
# ============================================================================

# get overall means
grand_mean_r = mean(df3_home.mean_r)
grand_mean_k = mean(df3_home.mean_k)

# get means by groups
gd = groupby(df3_home, [:env_temp,:copepods])
df_means = combine(gd, [:mean_r, :mean_k] .=> mean)
df_means.r_diff = df_means.mean_r_mean .- grand_mean_r
df_means.k_diff = df_means.mean_k_mean .- grand_mean_k

#### plot averages of home versus away connected by lines
p2 = plot(ylims=[0,2e4],
    xlabel = "Intrinsic rate of growth, r",
    ylabel = "Carrying capacity, K",
    size=(500,500),
    dpi=600,
    margin = 5Plots.mm)
for i = 1:15
    plot!(p2,[df3.mean_r[i],df3.mean_r[i+15]],[df3.mean_k[i],df3.mean_k[i+14]],
        arrow=true,color=:black)
    scatter!(p2,[df3.mean_r[i]],[df3.mean_k[i]],legend=:none,
        markerstrokecolor=gradient[3],markercolor=gradient[3],)
    #scatter!(p2,[dfB.mean_r[i+15]],[dfB.mean_k[i+15]],legend=:none,
        #markerstrokecolor=gradient[2],markercolor=gradient[2],)
end
scatter!([grand_mean_r],[grand_mean_k])

########### to make arrow plot from zeros, need another df
zero_vec_4 = zeros(Float64,1,4) # make a vector of zeros
# subset out each treatment  
r_20_no_cop = df3_home.r_diffs[df3_home.env_temp .== 20 .&& df3_home.copepods .== 0]
k_20_no_cop = df3_home.k_diffs[df3_home.env_temp .== 20 .&& df3_home.copepods .== 0]
r_26_no_cop = df3_home.r_diffs[df3_home.env_temp .== 26 .&& df3_home.copepods .== 0]
k_26_no_cop = df3_home.k_diffs[df3_home.env_temp .== 26 .&& df3_home.copepods .== 0]
r_20_cop = df3_home.r_diffs[df3_home.env_temp .== 20 .&& df3_home.copepods .== 1]
k_20_cop = df3_home.k_diffs[df3_home.env_temp .== 20 .&& df3_home.copepods .== 1]
r_26_cop = df3_home.r_diffs[df3_home.env_temp .== 26 .&& df3_home.copepods .== 1]
k_26_cop = df3_home.k_diffs[df3_home.env_temp .== 26 .&& df3_home.copepods .== 1]

# open an empty plot with dimensions and axis labels
p3 = plot(ylims=[-2e4,2e4],
    xlabel = "Change in r",
    ylabel = "Change in K",
    size=(500,500),
    dpi=600,
    margin = 5Plots.mm)
# cycle through the replicates
for i = 1:4 # 4 reps for everything but 26cop
    plot!(p3,[0,r_20_no_cop[i]],[0,k_20_no_cop[i]],
    label=false,arrow=true,
        color=gradient[2])
    plot!(p3,[0,r_26_no_cop[i]],[0,k_26_no_cop[i]],
    label=false,arrow=true,
        color=gradient[4])
    plot!(p3,[0,r_20_cop[i]],[0,k_20_cop[i]],
    label=false,arrow=true,
        color=gradient[7])
end
# cycle through for 26COP which only has 3 reps
for i = 1:3
    plot!(p3,[0,r_26_cop[i]],[0,k_26_cop[i]],
    label=false,arrow=true,
        color=gradient[9])
end

# plotting one more set of lines here to get the legend
plot!(p3,[0,r_20_no_cop[1]],[0,k_20_no_cop[1]],
    label="20C_nopreds",arrow=true,
    color=gradient[2])
plot!(p3,[0,r_26_no_cop[1]],[0,k_26_no_cop[1]],
    label="20C_preds",arrow=true,
    color=gradient[4])
plot!(p3,[0,r_20_cop[1]],[0,k_20_cop[1]],
    label="26C_nopreds",arrow=true,
    color=gradient[7])
plot!(p3,[0,r_26_cop[1]],[0,k_26_cop[1]],
    label="26C_preds",arrow=true,
    color=gradient[9])

# plot again based on the 
plot!(p3,[0,df_means.r_diff[1]],[0,df_means.k_diff[1]],
    label="",arrow=true,linewidth=4,
    color=gradient[2]) #label="20C_nopreds"
plot!(p3,[0,df_means.r_diff[2]],[0,df_means.k_diff[2]],
    label="",arrow=true,linewidth=4,
    color=gradient[4]) #20C_preds
plot!(p3,[0,df_means.r_diff[3]],[0,df_means.k_diff[3]],
    label="",arrow=true,linewidth=4,
    color=gradient[7]) #26C_nopreds
plot!(p3,[0,df_means.r_diff[4]],[0,df_means.k_diff[4]],
    label="",arrow=true,linewidth=4,
    color=gradient[9]) #26C_preds

plot!(p3,[-0.2,0.6], -15000 .* [-0.2,0.6],color=:black,
    line=:dash,width=2,label="Dividing line")

# put the three figures in a layout
figure2 = plot(p2,p3, layout=(1,2), size=(600,300), dpi=1000)
savefig(figure2,"vectors.png")

# ============================================================================
# now fit models to test for treatment and trait effects on the home-away differences
# ============================================================================

diffs_mdl1 = lm(@formula(r_diffs ~ env_temp * copepods), df3_home)
diffs_mdl1 = lm(@formula(r_diffs ~ env_temp + copepods), df3_home)


diffs_mdl1 = lm(@formula(k_diffs ~ env_temp * copepods), df3_home)
diffs_mdl1 = lm(@formula(k_diffs ~ env_temp + copepods), df3_home)


diffs_mdl1 = lm(@formula(r_diffs ~ pc1 + pc2), df3_home)

diffs_mdl1 = lm(@formula(r_diffs ~ env_temp + copepods + pc1), df3_home)
diffs_mdl1 = lm(@formula(r_diffs ~ env_temp + pc1), df3_home)
diffs_mdl1 = lm(@formula(r_diffs ~ pc1), df3_home)




diffs_mdl1 = lm(@formula(k_diffs ~ pc1 + pc2), df3_home)

diffs_mdl1 = lm(@formula(k_diffs ~ env_temp + copepods + pc1 + pc2), df3_home)
diffs_mdl1 = lm(@formula(k_diffs ~ env_temp + pc1 + pc2), df3_home)

diffs_mdl1
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# extra below here, parasitize or delete
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================

# ============================================================================
# third, is there an effect of treatments on traits
# ============================================================================

# cell morphology and behavior
mm_cellvol_int = fit(MixedModel, @formula(cell_vol ~ env_temp * copepods + (1|names)), df_home)
mm_speed_int = fit(MixedModel, @formula(speed ~ env_temp * copepods + (1|names)), df_home)
mm_aspectratio_int = fit(MixedModel, @formula(aspect_ratio ~ env_temp * copepods + (1|names)), df_home)
mm_turning_int = fit(MixedModel, @formula(sd_turning ~ env_temp * copepods + (1|names)), df_home)
mm_major_int = fit(MixedModel, @formula(mean_major ~ env_temp * copepods + (1|names)), df_home)
mm_minor_int = fit(MixedModel, @formula(mean_minor ~ env_temp * copepods + (1|names)), df_home)

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




#= first panel, colored by copepod treatment
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
=#