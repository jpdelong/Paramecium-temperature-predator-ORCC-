


using MixedModels

fitted_model_plasticity = fit(MixedModel, @formula(mean_r ~ Conjugation + (1|Population)), newdf)


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
