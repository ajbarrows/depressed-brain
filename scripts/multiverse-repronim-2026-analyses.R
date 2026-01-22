# ---
# title: "multiverse analysis - repronim 2026"
# output:
#   html_document:
#     toc: true
#     toc_depth: 3
# ---

proj_root = "/shared/hackathon/working-area/"

pacman::p_load("dplyr", "lavaan", "semptools")

data <- read.csv(paste0(proj_root, ""))

process(data = data_full, cov = c("sample_meancen", "ses_meancen", "DEM_2"), y = "MASQ_GD", x = "ela_meancen", m = c("age_onset_meancen", "pam_meancen"), boot = 10000, model = 6, total = 1)

# left_AMYG_uni <- '
# 
# # latent true scores
# left_AMYGlv_BL =~ 1*left_AMYG_BL
# left_AMYGlv_Y2 =~ 1*left_AMYG_Y2
# left_AMYGlv_Y4 =~ 1*left_AMYG_Y4
# left_AMYGlv_Y6 =~ 1*left_AMYG_Y6
# 
# # perfect autoregression (LCS backbone)
# left_AMYGlv_Y2 ~ 1*left_AMYGlv_BL
# left_AMYGlv_Y4 ~ 1*left_AMYGlv_Y2
# left_AMYGlv_Y6 ~ 1*left_AMYGlv_Y4
# 
# # latent change scores
# dleft_AMYG1 =~ 1*left_AMYGlv_Y2
# dleft_AMYG2 =~ 1*left_AMYGlv_Y4
# dleft_AMYG3 =~ 1*left_AMYGlv_Y6
# 
# # residual variances (constrained)
# left_AMYG_BL ~~ res_left_AMYGvar*left_AMYG_BL
# left_AMYG_Y2 ~~ res_left_AMYGvar*left_AMYG_Y2
# left_AMYG_Y4 ~~ res_left_AMYGvar*left_AMYG_Y4
# left_AMYG_Y6 ~~ res_left_AMYGvar*left_AMYG_Y6
# 
# # self-feedback (inertia)
# dleft_AMYG1 ~ left_AMYGselfFB*left_AMYGlv_BL
# dleft_AMYG2 ~ left_AMYGselfFB*left_AMYGlv_Y2
# dleft_AMYG3 ~ left_AMYGselfFB*left_AMYGlv_Y4
# 
# # growth factors
# ileft_AMYG =~ 1*left_AMYGlv_BL
# sleft_AMYG =~ 1*dleft_AMYG1 + 1*dleft_AMYG2 + 1*dleft_AMYG3
# 
# # growth means and variances
# ileft_AMYG ~ 1
# ileft_AMYG ~~ ileft_AMYG
# sleft_AMYG ~ 1
# sleft_AMYG ~~ sleft_AMYG
# 
# # intercept–slope covariance
# ileft_AMYG ~~ sleft_AMYG
# '
# 
# fit_left_AMYG <- lavaan(left_AMYG_uni, data=data, estimator='mlr',missing='fiml')
# summary(fit_left_AMYG, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# right_AMYG_uni <- '
# 
# # latent true scores
# right_AMYGlv_BL =~ 1*right_AMYG_BL
# right_AMYGlv_Y2 =~ 1*right_AMYG_Y2
# right_AMYGlv_Y4 =~ 1*right_AMYG_Y4
# right_AMYGlv_Y6 =~ 1*right_AMYG_Y6
# 
# # perfect autoregression (LCS backbone)
# right_AMYGlv_Y2 ~ 1*right_AMYGlv_BL
# right_AMYGlv_Y4 ~ 1*right_AMYGlv_Y2
# right_AMYGlv_Y6 ~ 1*right_AMYGlv_Y4
# 
# # latent change scores
# dright_AMYG1 =~ 1*right_AMYGlv_Y2
# dright_AMYG2 =~ 1*right_AMYGlv_Y4
# dright_AMYG3 =~ 1*right_AMYGlv_Y6
# 
# # residual variances (constrained)
# right_AMYG_BL ~~ res_right_AMYGvar*right_AMYG_BL
# right_AMYG_Y2 ~~ res_right_AMYGvar*right_AMYG_Y2
# right_AMYG_Y4 ~~ res_right_AMYGvar*right_AMYG_Y4
# right_AMYG_Y6 ~~ res_right_AMYGvar*right_AMYG_Y6
# 
# # self-feedback (inertia)
# dright_AMYG1 ~ right_AMYGselfFB*right_AMYGlv_BL
# dright_AMYG2 ~ right_AMYGselfFB*right_AMYGlv_Y2
# dright_AMYG3 ~ right_AMYGselfFB*right_AMYGlv_Y4
# 
# # growth factors
# iright_AMYG =~ 1*right_AMYGlv_BL
# sright_AMYG =~ 1*dright_AMYG1 + 1*dright_AMYG2 + 1*dright_AMYG3
# 
# # growth means and variances
# iright_AMYG ~ 1
# iright_AMYG ~~ iright_AMYG
# sright_AMYG ~ 1
# sright_AMYG ~~ sright_AMYG
# 
# # intercept–slope covariance
# iright_AMYG ~~ sright_AMYG
# '
# 
# fit_right_AMYG <- lavaan(right_AMYG_uni, data=data, estimator='mlr',missing='fiml')
# summary(fit_right_AMYG, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# left_HIPPO_uni <- '
# 
# # latent true scores
# left_HIPPOlv_BL =~ 1*left_HIPPO_BL
# left_HIPPOlv_Y2 =~ 1*left_HIPPO_Y2
# left_HIPPOlv_Y4 =~ 1*left_HIPPO_Y4
# left_HIPPOlv_Y6 =~ 1*left_HIPPO_Y6
# 
# # perfect autoregression (LCS backbone)
# left_HIPPOlv_Y2 ~ 1*left_HIPPOlv_BL
# left_HIPPOlv_Y4 ~ 1*left_HIPPOlv_Y2
# left_HIPPOlv_Y6 ~ 1*left_HIPPOlv_Y4
# 
# # latent change scores
# dleft_HIPPO1 =~ 1*left_HIPPOlv_Y2
# dleft_HIPPO2 =~ 1*left_HIPPOlv_Y4
# dleft_HIPPO3 =~ 1*left_HIPPOlv_Y6
# 
# # residual variances (constrained)
# left_HIPPO_BL ~~ res_left_HIPPOvar*left_HIPPO_BL
# left_HIPPO_Y2 ~~ res_left_HIPPOvar*left_HIPPO_Y2
# left_HIPPO_Y4 ~~ res_left_HIPPOvar*left_HIPPO_Y4
# left_HIPPO_Y6 ~~ res_left_HIPPOvar*left_HIPPO_Y6
# 
# # self-feedback (inertia)
# dleft_HIPPO1 ~ left_HIPPOselfFB*left_HIPPOlv_BL
# dleft_HIPPO2 ~ left_HIPPOselfFB*left_HIPPOlv_Y2
# dleft_HIPPO3 ~ left_HIPPOselfFB*left_HIPPOlv_Y4
# 
# # growth factors
# ileft_HIPPO =~ 1*left_HIPPOlv_BL
# sleft_HIPPO =~ 1*dleft_HIPPO1 + 1*dleft_HIPPO2 + 1*dleft_HIPPO3
# 
# # growth means and variances
# ileft_HIPPO ~ 1
# ileft_HIPPO ~~ ileft_HIPPO
# sleft_HIPPO ~ 1
# sleft_HIPPO ~~ sleft_HIPPO
# 
# # intercept–slope covariance
# ileft_HIPPO ~~ sleft_HIPPO
# '
# fit_left_HIPPO <- lavaan(left_HIPPO_uni, data=data, estimator='mlr',missing='fiml')
# summary(fit_left_HIPPO, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# right_HIPPO_uni <- '
# 
# # latent true scores
# right_HIPPOlv_BL =~ 1*right_HIPPO_BL
# right_HIPPOlv_Y2 =~ 1*right_HIPPO_Y2
# right_HIPPOlv_Y4 =~ 1*right_HIPPO_Y4
# right_HIPPOlv_Y6 =~ 1*right_HIPPO_Y6
# 
# # perfect autoregression (LCS backbone)
# right_HIPPOlv_Y2 ~ 1*right_HIPPOlv_BL
# right_HIPPOlv_Y4 ~ 1*right_HIPPOlv_Y2
# right_HIPPOlv_Y6 ~ 1*right_HIPPOlv_Y4
# 
# # latent change scores
# dright_HIPPO1 =~ 1*right_HIPPOlv_Y2
# dright_HIPPO2 =~ 1*right_HIPPOlv_Y4
# dright_HIPPO3 =~ 1*right_HIPPOlv_Y6
# 
# # residual variances (constrained)
# right_HIPPO_BL ~~ res_right_HIPPOvar*right_HIPPO_BL
# right_HIPPO_Y2 ~~ res_right_HIPPOvar*right_HIPPO_Y2
# right_HIPPO_Y4 ~~ res_right_HIPPOvar*right_HIPPO_Y4
# right_HIPPO_Y6 ~~ res_right_HIPPOvar*right_HIPPO_Y6
# 
# # self-feedback (inertia)
# dright_HIPPO1 ~ right_HIPPOselfFB*right_HIPPOlv_BL
# dright_HIPPO2 ~ right_HIPPOselfFB*right_HIPPOlv_Y2
# dright_HIPPO3 ~ right_HIPPOselfFB*right_HIPPOlv_Y4
# 
# # growth factors
# iright_HIPPO =~ 1*right_HIPPOlv_BL
# sright_HIPPO =~ 1*dright_HIPPO1 + 1*dright_HIPPO2 + 1*dright_HIPPO3
# 
# # growth means and variances
# iright_HIPPO ~ 1
# iright_HIPPO ~~ iright_HIPPO
# sright_HIPPO ~ 1
# sright_HIPPO ~~ sright_HIPPO
# 
# # intercept–slope covariance
# iright_HIPPO ~~ sright_HIPPO
# '
# 
# fit_right_HIPPO <- lavaan(right_HIPPO_uni, data=data, estimator='mlr',missing='fiml')
# summary(fit_right_HIPPO, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# DEP_uni <- '
# 
# # latent true scores
# DEPlv_BL =~ 1*DEP_BL
# DEPlv_Y2 =~ 1*DEP_Y2
# DEPlv_Y4 =~ 1*DEP_Y4
# DEPlv_Y6 =~ 1*DEP_Y6
# 
# # perfect autoregression
# DEPlv_Y2 ~ 1*DEPlv_BL
# DEPlv_Y4 ~ 1*DEPlv_Y2
# DEPlv_Y6 ~ 1*DEPlv_Y4
# 
# # latent change scores
# dDEP1 =~ 1*DEPlv_Y2
# dDEP2 =~ 1*DEPlv_Y4
# dDEP3 =~ 1*DEPlv_Y6
# 
# # residual variances (constrained)
# DEP_BL ~~ resDEPvar*DEP_BL
# DEP_Y2 ~~ resDEPvar*DEP_Y2
# DEP_Y4 ~~ resDEPvar*DEP_Y4
# DEP_Y6 ~~ resDEPvar*DEP_Y6
# 
# # self-feedback
# dDEP1 ~ DEPselfFB*DEPlv_BL
# dDEP2 ~ DEPselfFB*DEPlv_Y2
# dDEP3 ~ DEPselfFB*DEPlv_Y4
# 
# # growth factors
# iDEP =~ 1*DEPlv_BL
# sDEP =~ 1*dDEP1 + 1*dDEP2 + 1*dDEP3
# 
# # growth means and variances
# iDEP ~ 1
# iDEP ~~ iDEP
# sDEP ~ 1
# sDEP ~~ sDEP
# 
# # intercept–slope covariance
# iDEP ~~ sDEP
# '
# 
# fitDEP <- lavaan(DEP_uni, data=data, estimator='mlr',missing='fiml')
# summary(fitDEP, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# 
# BDCS<-'
# 
# left_AMYGlv_BL=~1*left_AMYG_BL        
# left_AMYGlv_Y2=~1*left_AMYG_Y2        
# left_AMYGlv_Y4=~1*left_AMYG_Y4        
# left_AMYGlv_Y6=~1*left_AMYG_Y6        
# 
# DEPlv_BL=~1*DEP_BL        
# DEPlv_Y2=~1*DEP_Y2        
# DEPlv_Y4=~1*DEP_Y4        
# DEPlv_Y6=~1*DEP_Y6        
# 
# ##### The following parameters capture the core assumptions of the LCS and should not generally be modified
# 
# left_AMYGlv_Y2 ~ 1*left_AMYGlv_BL     # This parameter regresses COG_T2 perfectly on COG_T1
# left_AMYGlv_Y4 ~ 1*left_AMYGlv_Y2     # This parameter regresses COG_T3 perfectly on COG_T2
# left_AMYGlv_Y6 ~ 1*left_AMYGlv_Y4     # This parameter regresses COG_T4 perfectly on COG_T3
# 
# DEPlv_Y2 ~ 1*DEPlv_BL     # This parameter regresses NEU_T2 perfectly on NEU_T1
# DEPlv_Y4 ~ 1*DEPlv_Y2     # This parameter regresses NEU_T3 perfectly on NEU_T2
# DEPlv_Y6 ~ 1*DEPlv_Y4     # This parameter regresses NEU_T4 perfectly on NEU_T3
# 
# dleft_AMYG1 =~ 1*left_AMYGlv_Y2       # This defines the change score as measured perfectly by scores on COG_T2
# dleft_AMYG2 =~ 1*left_AMYGlv_Y4       # This defines the change score as measured perfectly by scores on COG_T3
# dleft_AMYG3 =~ 1*left_AMYGlv_Y6       # This defines the change score as measured perfectly by scores on COG_T4
# 
# dDEP1 =~ 1*DEPlv_Y2       # This defines the change score as measured perfectly by scores on NEU_T2
# dDEP2 =~ 1*DEPlv_Y4       # This defines the change score as measured perfectly by scores on NEU_T3
# dDEP3 =~ 1*DEPlv_Y6       # This defines the change score as measured perfectly by scores on NEU_T4
# 
# left_AMYG_BL~~res_left_AMYGvar*left_AMYG_BL          # This estimates the COG residual variances 
# left_AMYG_Y2~~res_left_AMYGvar*left_AMYG_Y2          # This estimates the COG residual variances 
# left_AMYG_Y4~~res_left_AMYGvar*left_AMYG_Y4          # This estimates the COG residual variances 
# left_AMYG_Y6~~res_left_AMYGvar*left_AMYG_Y6          # This estimates the COG residual variances 
# 
# DEP_BL~~resDEPvar*DEP_BL          # This estimates the NEU residual variances 
# DEP_Y2~~resDEPvar*DEP_Y2          # This estimates the NEU residual variances 
# DEP_Y4~~resDEPvar*DEP_Y4          # This estimates the NEU residual variances 
# DEP_Y6~~resDEPvar*DEP_Y6          # This estimates the NEU residual variances 
# 
# #Dynamics
# 
# dDEP1~DEPselfFB*DEPlv_BL       # This estimates the NEU self-feedback parameter (equality constrained across timepoints)
# dDEP2~DEPselfFB*DEPlv_Y2       # This estimates the NEU self-feedback parameter (equality constrained across timepoints) 
# dDEP3~DEPselfFB*DEPlv_Y4       # This estimates the NEU self-feedback parameter (equality constrained across timepoints)
# 
# dleft_AMYG1~left_AMYGselfFB*left_AMYGlv_BL       # This estimates the COG self-feedback parameter (equality constrained across timepoints)
# dleft_AMYG2~left_AMYGselfFB*left_AMYGlv_Y2       # This estimates the COG self-feedback parameter (equality constrained across timepoints) 
# dleft_AMYG3~left_AMYGselfFB*left_AMYGlv_Y4       # This estimates the COG self-feedback parameter (equality constrained across timepoints)
# 
# dDEP1~leftAMYG_to_DEP*left_AMYGlv_BL         # This estimates the COG to NEU coupling parameter 
# dDEP2~leftAMYG_to_DEP*left_AMYGlv_Y2         # This estimates the COG to NEU coupling parameter 
# dDEP3~leftAMYG_to_DEP*left_AMYGlv_Y4         # This estimates the COG to NEU coupling parameter 
# 
# dleft_AMYG1~DEP_to_left_AMYG*DEPlv_BL        # This estimates the NEU to COG coupling parameter 
# dleft_AMYG2~DEP_to_left_AMYG*DEPlv_Y2        # This estimates the NEU to COG coupling parameter 
# dleft_AMYG3~DEP_to_left_AMYG*DEPlv_Y4        # This estimates the NEU to COG coupling parameter 
# 
# 
# ileft_AMYG=~1*left_AMYGlv_T1                   # This defines the COG intercept measurement model
# sleft_AMYG=~1*dleft_AMYG1+1*dleft_AMYG2+1*dleft_AMYG3      # This defines the COG slope measurement model
# ileft_AMYG~1                             # This estimates the COG intercept intercept (mean)
# ileft_AMYG~~ileft_AMYG                         # This estimates the COG intercept variance
# sleft_AMYG~1                             # This estimates the COG slope intercept
# sleft_AMYG~~sleft_AMYG                         # This estimates the COG slope variance
# 
# iDEP=~1*DEPlv_T1                   # This defines the NEU slope measurement model
# sDEP=~1*dDEP1+1*dDEP2+1*dDEP3      # This defines the NEU slope measurement model
# iDEP~1                             # This estimates the NEU intercept intercept (mean)
# iDEP~~iDEP                         # This estimates the NEU intercept variance
# sDEP~1                             # This estimates the NEU slope intercept
# sDEP~~sDEP                         # This estimates the NEU slope variance
# 
# iDEP~~sDEP                      # This estimates the iNEU sNEU covariance
# iDEP~~sleft_AMYG                      # This estimates the iNEU sCOG covariance
# iDEP~~ileft_AMYG                      # This estimates the iNEU iCOG covariance
# ileft_AMYG~~sleft_AMYG                      # This estimates the iCOG sCOG covariance        
# ileft_AMYG~~sDEP                      # This estimates the iCOG sNEU covariance
# sleft_AMYG~~sDEP                      # This estimates the sCOG sNEU covariance  
# 
# '
# 
# fitBDCS <- lavaan(BDCS, data=simdatBDCS, estimator='mlr',missing='fiml')
# summary(fitBDCS, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE) 
# 
# 
