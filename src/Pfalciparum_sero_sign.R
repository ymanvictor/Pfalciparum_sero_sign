
###############################################################

# R Script for Analysis of Microarray data for:
# 
# Yman V, Tuju J, White M T, Kamuyu G, Mwai K, Kibinge N, Asghar M,
# Sundling C, Sondén K, Murungi L, Kiboi D, Kimathi R, Chege T, Chepsat E, 
# Kiyuka P, Nyamako L, Osier F H A, Färnert A
# 
# Distinct kinetics in antibody responses to 111 Plasmodium falciparum 
# antigens identifies novel serological markers of recent malaria exposure
# 
# Author: Victor Yman
# Stockholm 2021-10-04

###############################################################

# load packages

packages <- c(#"renv", 
              "ggrepel", "foreach",
              "cowplot", "viridis", "RColorBrewer",
              "nlme", "contrast", 
              "pROC", "ROCR",
              "Boruta",
              "MESS", "MASS", 
              "tidyverse",
              "randomForest", "doParallel", 
              "betareg"
)


invisible(lapply(packages, library, character.only = TRUE))


# source heatmap.3 function

source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# load antibody and meta data

full_data <- read_csv(file = "data/ab_data.csv") # for variable descriptions see Suplementary File S1: Supplementary Table S11.

# set factor levels of factor variables

full_data <- full_data %>% mutate(
  antigen = factor(antigen, levels= unique(antigen)),
  t_id    = factor(t_id, levels= c("Acute", "10 days", "1 month", 
                                   "3 months", "6 months", "12 months", 
                                   "Neg. control")), 
  exposure = factor(exposure, levels = c("Primary infected", 
                                         "Previously exposed", 
                                         "Neg. control"))
)

# antigen names

antigen_name_vect <- levels(factor(full_data$antigen)) 


###############################################################
###############################################################
###############################################################

#                           Part 1. 
 
#                       Descriptive plots


###############################################################
###############################################################
###############################################################

# Custom color schemes for plotting

# Colours scheme exposure status

exposure_colors <- c(
  "#FF80FF",
  rgb(t(col2rgb("grey40")/255)),
  rgb(t(col2rgb("darkred")/255))
)

# Colour scheme time-point id

t_cols <- c(brewer.pal(6, "Set2"), "#8B0000")

###############################################################

#   Plot Ab reactivity heatmap (Fig. 1)

###############################################################

#calculate mean reactivity by antigen

full_data <- full_data %>% group_by(antigen) %>% 
  
  mutate(
    
    ant_mean = mean(MFI)
  
    ) %>%
  
  ungroup() %>% 
  
  arrange(t_id, desc(ant_mean))

#create ab reactivity matrix ordered by mean reactivity across antigens

ab_mtrx <- full_data %>% select(antigen, sample_id, MFI)

ab_mtrx <- as.data.frame(ab_mtrx %>% pivot_wider(names_from = sample_id, values_from = MFI))

rownames(ab_mtrx) <- ab_mtrx$antigen

ab_mtrx <- ab_mtrx[2:283]


# Extract data for ordering of heat map and heat map color side bars

samples <- as.data.frame(colnames(ab_mtrx)) 

colnames(samples) <- "sample_id"

# calculate mean reactivity by column of matrix (i.e. by sample_id)

col_mean <- foreach(n = 1:ncol(ab_mtrx), .combine = 'c') %do% {
  mean(ab_mtrx[, n])
}

# bind sample_id and mean reactivity

samples <- cbind(samples, col_mean)

#create meta data df. 

meta_subset <- full_data %>% 
  
  select(sample_id, exposure, t_id, t_d, ant_mean) %>% 
  
  distinct(sample_id, .keep_all = T)

meta <- samples %>% 
  
  left_join(meta_subset, by = "sample_id") %>% 
  
  arrange(exposure, t_id, desc(col_mean))

column_order <- meta$sample_id

# create exposure vector

exposure <- meta$exposure

# create time-point vector

t_id <- meta$t_id

# Order data matrix columns according to meta data (sample_ids ordered by i) exposure, ii) time-point, 
# and iii) mean reactivity across all antigens)

ab_mtrx <- ab_mtrx %>% select(all_of(column_order))

# Create colour vectors for for heat map sidebars

exposure_color_vect <- exposure_colors[exposure]
t_id_color_vect <- t_cols[t_id]

clab = cbind(t_id_color_vect, exposure_color_vect)
colnames(clab)=c("Time-point", "Exposure Status")


pdf(file="results/Fig_1.pdf")

main_title=""
par(cex.main=1, mar = c(1.1, 1.1, 4.7, .5))
heatmap.3((as.matrix(ab_mtrx)),
          na.rm = T, 
          scale="none", 
          dendrogram="none",
          margins=c(6, 8),
          Rowv=F, 
          Colv=F, 
          ColSideColors=clab, 
          symbreaks=F, 
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none", 
          main=main_title, 
          labCol=FALSE, 
          cexRow=0.4,
          col=(viridis(400, alpha=1, begin = 0.1, end = 1, option = "B")),
          ColSideColorsSize=3, 
          KeyValueName="Ab intensity",
          keysize = 1, 
          colsep=c(91, 241), 
          sepwidth = 0.01, 
          breaks = seq(200, 60000, length.out = 401))

legend("bottom",
       legend=c("Primary infected",
                "Previously exposed",
                "","",
                "Acute","10 Days","1 Month",
                "3 Months","6 Months","12 Months", 
                "Neg. Control"),
       fill=c(exposure_colors[1:2],
              "white", "white",
              t_cols[1:7]), 
       border=FALSE, bty="n", 
       y.intersp = 0.9, cex=0.5, ncol = 6)

dev.off()

#################################################################

# Plot mean Ab reactivity and breadth over time (Fig. 2ab)

#################################################################

# Create mean reactivity and antibody breadth variables


full_data <- full_data %>% 
  
  group_by(sample_id) %>% 
  
  mutate(
    mean_MFI = exp(mean(log_MFI))
  )


full_data <- full_data %>% 
  
  group_by(study_idn, t_id) %>%
  
  mutate(
    breadth = sum(sero_pos)
  )


# plot mean ab reactivity over time (Across 111 antigens)

plot_mean <- full_data %>% 
  
  ggplot(aes( x = t_id, y = mean_MFI, fill = exposure))+
  
  geom_boxplot()+
  
  scale_y_log10(limits = c(200, 24000))+
  
  theme_bw()+
  
  theme(legend.position = "bottom", legend.title = element_blank(), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(1, "cm"), 
        plot.title = element_text(size = 28))+
  
  xlab("Time-point")+
  
  ylab("Mean normalised MFI")+
  
  scale_fill_manual(values = exposure_colors)+
  
  ggtitle("a")

# plot breadth over time (111 antigens)

plot_breadth <- full_data %>%
  
  distinct(antigen, .keep_all = T) %>%
  
  ggplot(aes( x = t_id, y = breadth, fill = exposure))+
  
  geom_boxplot()+
  
  ylab("Antibody breadth")+xlab("Time-point")+
  
  theme_bw()+
  
  theme(legend.position = "bottom", legend.title = element_blank(), 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(1, "cm"), 
        plot.title = element_text(size = 28))+
  
  scale_fill_manual(values = exposure_colors)+
  
  ggtitle("b")



# generate plot grid and save plot (Fig. 2ab)

ggsave(plot_grid(plot_mean, plot_breadth, ncol = 2), file = "results/Fig_2ab.pdf", scale = 1, width = 18, height = 9)


###############################################################
###############################################################





###############################################################
###############################################################
###############################################################



#                          Part 2. 

# Differential response analysis (previously exposed vs. 
# primary infected) using linear mixed effects models



###############################################################
###############################################################
###############################################################





###############################################################
###############################################################



#                         FUNCTIONS



###############################################################
###############################################################

# Function 1: Model (lme) and contrasts function 
# a) Fits a mixed effects model
# b) Calculates the linear combinations of coefficients of interest.

fun_lme_contrasts <- function(df) {
  
  # a
  
  model = lme(
    
    log(MFI)~t_id*exposure, random = ~1|study_idn, data = df
    
  )
  
  # b
  conts = contrast(
    
    model, 
    
    list(t_id = levels(df$t_id),
         exposure="Primary infected"
    ), 
    
    list(t_id = levels(df$t_id), 
         exposure = "Previously exposed"
    )
  )
  
  return(conts)
  
}

# Function 2: Extract linear combinations and corresponding p-values

fun_lincom <- function(contrasts) {
  
  contrs  = unlist(contrasts$Contrast)
  
  pval   = unlist(contrasts$Pvalue)
  
  lincoms = as.data.frame(cbind(contrs, pval))
  
  return(lincoms)
  
}

###############################################################
###############################################################



# Prepare data to fit models

# drop neg. control data and create nested (by antibody response) data frame

ab_nest <- full_data %>% 
  
  filter(exposure != "Neg. control") %>%
  
  mutate(
    # update factor levels for subset data without neg. controls
    exposure = factor(exposure, 
                      levels = c("Primary infected", "Previously exposed")
    ),
    t_id     = factor(t_id, 
                      levels = c("Acute", "10 days", "1 month", 
                                 "3 months", "6 months", "12 months"
                      )
    )
  ) %>%
  
  group_by(antigen) %>% 
  
  nest()

# fit models to nested data return linear combination of coefficients

ab_nest <- ab_nest %>% 
  
  mutate(
    
    contrasts = map(data, fun_lme_contrasts),
    lincoms   = map(contrasts, fun_lincom)
    
  )

# create separate data frame for estimated fold-differences

fold_diffs <- ab_nest %>% 
  
  select(antigen, lincoms) %>% 
  
  unnest(cols = c (lincoms))

fold_diffs$t_id <- factor(
  
  rep(
    c(
      "Acute", "10 days", "1 month",
      "3 months", "6 months", "12 months"
    ), 
    nrow(fold_diffs)/6
  ), 
  
  levels = c(
    "Acute", "10 days", "1 month",
    "3 months", "6 months", "12 months"
  )
)

# calculate log2(fold-differences) and FDR adjusted p-values

fold_diffs <- fold_diffs %>%
  
  mutate(
    
    fold_diff  = log2(exp(-contrs)), 
    
    qval       = p.adjust(pval, method = "BH"), 
    
    diff_signif = ifelse(qval<0.05 , 1, 0)
    
  )

###############################################################

# Voclano plots (Fig. 2c)

###############################################################

# create "in-plot labels" for top 20 antibody ressponses with greatest differences

diff_ab_names <- fold_diffs %>%
  
  filter(diff_signif==1) %>% 
  
  group_by(t_id) %>% 
  
  slice_max(fold_diff, n = 20)

# create volcano-plot

x_limits <- c(NA, -0.5)
y_limits <- c(0.5, NA)

plot_volcano_exposure_lme <- ggplot(
  
  fold_diffs, 
  aes( x= fold_diff,                                    
       y = -log10(qval)
  )
)+
  
  geom_point(aes(colour = factor(diff_signif)), size = 1.2)+
  
  geom_hline(aes(yintercept = -log10(0.05)), lty = "dotted")+
  
  geom_vline(aes(xintercept = 1), lty = "dotted")+
  
  geom_vline(aes(xintercept = -1), lty = "dotted")+
  
  facet_wrap(~t_id, scales = "fixed", labeller = label_value, nrow=1)+
  
  scale_x_continuous(limits = c(-3.5, 3.5))+
  
  scale_color_manual(values = c("black", "red", "blue"))+
  
  geom_label_repel(
    data = diff_ab_names, aes(label = antigen),
    size = 2, 
    label.size = 0.4,
    box.padding = 0, 
    point.padding = 0, 
    segment.alpha = 0.4,
    max.iter = 10000,
    max.overlaps = 15,
    #force = 1,
    xlim = x_limits, 
    ylim = y_limits,
    nudge_x = -0.2,
    # nudge_x = 0.02,
    direction = "y", 
    hjust = 1, 
    segment.size = 0.5)  +
  
  theme_bw()+
  
  theme(
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20))+
  
  xlab(expression("log"[2]*" (Fold difference)"))+
  ylab(expression("-log"[10]*" (FDR-adjusted P-value)"))

ggsave(plot_volcano_exposure_lme, file = "results/Fig_2C.pdf", width = 18, height = 4.5)


###############################################################
###############################################################






###############################################################
###############################################################
###############################################################




#                          Part 3
# 
# Classifying recent infection based on threshold Ab level
# towards individual antigens



###############################################################
###############################################################
###############################################################


# prepare data to fit simple threshold Ab level classifiers

full_nest <- full_data %>%
  
  group_by(antigen) %>%
  
  nest()

# fit classifiers and extract results

set.seed(6)

ab_lev_classifiers <- full_nest %>%
  
  mutate(
    
      roc_obj = data %>% 
        map(
            function(x) roc(t_bin~log_MFI, data = x, ci = T, ci.method = "delong")
      
    ),

    auc = map_dbl(roc_obj, "auc"),

    auc_ci = map(roc_obj, "ci"),

    auc_lb = map_dbl(auc_ci, 1),

    auc_ub = map_dbl(auc_ci, 3),

    roc_sens = map(roc_obj, "sensitivities"),

    roc_spec = map(roc_obj, "specificities")
    
    
  ) %>%
  
  arrange(desc(auc)) %>%
  
  select(-data, -roc_obj, -auc_ci)

###############################################################

# plot classifier results (Fig. 3)

###############################################################

pdf("results/Fig_3.pdf")

par(mfrow=c(1,1))
par(mar=c(3,3,4,1.5))
par(mgp=c(1.75, 0.6,0))

line_seq <- c(0.2, 0.4, 0.6, 0.8)

plot(
  x=c(0,1), 
  y=c(0,1),
  
  type='l', 
  lty="dashed",
  
  xlim=c(0,1.002), ylim=c(0,1.002),
  
  xaxs='i', yaxs='i', 
  xaxt='n', yaxt='n',
  
  xlab="False positive rate (100 - specificity)", 
  ylab="True positive rate (sensitivity)",
  
  main=paste("
             ROC curve for classifying recent", 
             "infection: within 90 days"
  ),
  
  cex.lab=1.1, 
  cex.axis=0.8, 
  cex.main=1.2)


for(i in 1:4)
{
  points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
  points(x=rep(1-line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1),
     labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 )

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1),
     labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 )


line_width <-  c(
  
  rep(2, 10), 
  
  rep(0.05, 
      nrow(ab_lev_classifiers)-10)
)

top_10_cols <- c(
  
  rainbow(10), 
  
  rep("lightgrey", 
      nrow(ab_lev_classifiers)-10
  )
)

for (i in 1:nrow(ab_lev_classifiers)) {
  {
    points( 
      
      x = 1-ab_lev_classifiers$roc_spec[[i]],
      
      y = ab_lev_classifiers$roc_sens[[i]],
      
      type = 's', 
      
      lwd =line_width[i], 
      
      col = top_10_cols[i])
  }
}


legend(
  x="bottomright",
  bty='n', 
  cex=0.7,
  fill=top_10_cols,
  border=top_10_cols,
  legend = c(
    paste0(ab_lev_classifiers$antigen[1:10],
           
           "; AUC = ",
           
           as.character(round(ab_lev_classifiers$auc[1:10], 3)), 
           
           " (95% CI:", 
           
           as.character(round(ab_lev_classifiers$auc_lb[1:10], 3)), 
           
           " - ", 
           
           as.character(round(ab_lev_classifiers$auc_ub[1:10], 3)), 
           
           ")"
    )
  )
  
)

dev.off()
###############################################################

# create results summary 

table_s2 <- ab_lev_classifiers %>%
  
  unnest(cols = c(roc_sens, roc_spec)) %>% 
  
  mutate(
    
    diff = sqrt((roc_sens - roc_spec)^2)
    
  ) %>% 
  
  group_by(antigen) %>%
  
  filter(diff == min(diff)) %>%
  
  select(-diff)


###############################################################
###############################################################






###############################################################
###############################################################
###############################################################



#                          Part 4.
#
#  Classifying recent infection based on combination of 2 to 5 
#  antibody responses
   
   

###############################################################
###############################################################
###############################################################



###############################################################
###############################################################




#                         FUNCTIONS



###############################################################
###############################################################

# Function 3. Performs cross-validated random forest classification
#             using paralell backends

###############################################################

fun_RF_classifier <- function(ab_combo) 
  
{
  
  N_tree = 500
  N_rep = 200
  N_plot = 1000
  train_prop = 0.6667
  
  ab_wide_local = ab_wide[, 3:ncol(ab_wide)]
  
  print(ab_combo)
  
  
  vals =  foreach(n = 1:N_rep, .combine = 'rbind',
                  
                  .packages = c("ROCR", "MASS", "MESS", 
                                "tibble", "randomForest")
                  
                  ) %dopar%
    
  {
    
    ############################################
    ## Prepare testing and training data
    
    index_train = sample( nrow(ab_wide_local), train_prop*nrow(ab_wide_local) )
    index_test = setdiff( 1:nrow(ab_wide_local), index_train )
    
    data_train = ab_wide_local[index_train,]
    data_test = ab_wide_local[index_test,]
    
    
    ############################################
    ## Fit RF classifier and create prediction object
    
    tryCatch(
      {
        RF = randomForest(
          
          as.formula(c("as.factor(t_bin)~", paste(ab_combo, collapse = "+"))),
          
          data=as.data.frame(data_train), 
          
          importance=TRUE, ntree=N_tree )
        
        RF_pred_obj = predict(RF, newdata=data_test, predict.all=TRUE)
        
        RF_votes = rowSums(RF_pred_obj$individual=="recent")/N_tree
        
        RF_pred = prediction(RF_votes, data_test$t_bin)
        
        RF_perf_xy = performance(RF_pred, "tpr", "fpr")
        
        RF_auc = performance(RF_pred, "auc")@y.values[[1]]
        
      }, error=function(e){ NULL }
    )
    
    
    ##############################################
    ## Update x and y values	
    
    each_vals = cbind(RF_perf_xy@x.values[[1]], 
                      RF_perf_xy@y.values[[1]], 
                      rep(RF_auc, length(RF_perf_xy@x.values[[1]]))
                      )
    
  }
  
  lowess_xy = lowess(vals[,1], vals[, 2], f = 0.2)
  
  RF_Xval_x = lowess_xy$x[seq(from = 1, to = length(lowess_xy$x), length=N_plot)]
  RF_Xval_y = lowess_xy$y[seq(from = 1, to = length(lowess_xy$y), length=N_plot)]
  
  RF_Xval_auc = vals[sort(sample(1:nrow(vals), size = N_plot)), 3]
  
  RF_Xval = as.data.frame(cbind(RF_Xval_x, RF_Xval_y, RF_Xval_auc))
  
  RF_Xval = RF_Xval[sort(sample(1:nrow(RF_Xval), 200)),]
  
}

###############################################################
###############################################################


# prepare data for RF classification

ab_wide <- full_data %>% 
  
  ungroup() %>%
  
  select(study_idn, idn, sample_id, antigen, t_id, t_bin, log_MFI) %>%
  
  pivot_wider(names_from = antigen, values_from = log_MFI)


###############################################################

#  Feature selection using Boruta algorithm

###############################################################


# boruta

set.seed(6)

boruta_results <- Boruta(
  
  factor(t_bin)~.,
  
  data = ab_wide[,5:ncol(ab_wide)], 
  
  doTrace = 2, maxRuns = 200   # Note: For demo run maxRuns is set to 200. For full run set maxRuns to 10000.
  
)

plotImpHistory(boruta_results)

boruta_signif <- names(
  
  boruta_results$finalDecision[
    
    boruta_results$finalDecision %in% c("Confirmed")
    ]
)

print(boruta_signif)

getConfirmedFormula(boruta_results)

###############################################################

# Plot of boruta results (Fig. 4a)

###############################################################

pdf("results/Fig_4a.pdf")
par(mfrow=c(1,1))
par(mar=c(7,5,3,1))

plot(
  boruta_results,
  
  cex.axis=.25, 
  
  las=2, xlab="", 
  
  main="Variable Importance"
  
)

dev.off()

###############################################################
###############################################################

# create objets containing all 2 to 5 way combinations of Ab responses

# 2 antigens

ab_combo_n2 <- as_tibble(
  
  t(
    combn(boruta_signif, 2)
  ), 
  
  stringsAsFactors = F, 
  
  .name_repair = "unique"
  
)

ab_combo_n2$combo <- seq(1:nrow(ab_combo_n2))

# 3 antigens

ab_combo_n3 <- as_tibble(
  
  t(
    
    combn(
      
      boruta_signif, 3
    )
  ),
  
  stringsAsFactors = F,
  
  .name_repair = "unique"
  
)

ab_combo_n3$combo <- seq(1:nrow(ab_combo_n3))

# 4 antigens

ab_combo_n4 <- as_tibble(
  
  t(
    
    combn(
      
      boruta_signif, 4
      
    )
    
  ),
  
  stringsAsFactors = F,
  
  .name_repair = "unique"
  
)

ab_combo_n4$combo <- seq(1:nrow(ab_combo_n4))

# 5 antigens

ab_combo_n5 <- as_tibble(
  
  t(
    combn(
      
      boruta_signif, 5
      
    )
  ), 
  
  stringsAsFactors = F,
  
  .name_repair = "unique"
  
)

ab_combo_n5$combo <- seq(1:nrow(ab_combo_n5))

###############################################################

# initiate parallel backends for doPar function

cl <- makePSOCKcluster(4)
registerDoParallel(cl)

###############################################################

# run cross-validated random forest classification analysis 
# (using Function 3) and extract results
# 
# Note:
# For demo run, analysis is performed only on the first 10 combinations of antibody responses.
# 
# The full exhaustive analysis of all potential combinations of responses takes a long time to run even
# using a large number cores.

###############################################################

# RF classifier fitted to all 2-way combinations

ab_combo_n2 <- ab_combo_n2[1:10,] # SIC. for demo run only - selects first 10 combinations. silence for exhaustive analysis.

ab_combo_n2 <- ab_combo_n2 %>% 
  
  group_by(combo) %>% 
  
  nest() %>%
  
  mutate(
    
    RF_results = map(data, fun_RF_classifier),
    
    auc = RF_results %>% 
      
      map_dbl(function(x)
        
        MESS::auc(
          unique(x$RF_Xval_x), 
          unique(x$RF_Xval_y)
        )
        
      ), 
    
    auc_lb = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.025)),
    
    auc_ub = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.975)) 
  )

# RF classifier fitted to all 3-way combinations 

ab_combo_n3 <- ab_combo_n3[1:10,] # SIC. for demo run only - selects first 10 combinations. silence for exhaustive analysis.

ab_combo_n3 <- ab_combo_n3 %>% 
  
  group_by(combo) %>% 
  
  nest() %>%
  
  mutate(
    
    RF_results = map(data, fun_RF_classifier),
    
    auc = RF_results %>%
      
      map_dbl(function(x) 
        
        MESS::auc(
          unique(x$RF_Xval_x),
          unique(x$RF_Xval_y)
        )
        
      ),
    
    auc_l = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.025)),
    
    auc_ub = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.975)) 
  )

# RF classifier fitted to all 4-way combinations 

ab_combo_n4 <- ab_combo_n4[1:10,] # SIC. for demo run only - selects first 10 combinations. silence for exhaustive analysis

ab_combo_n4 <- ab_combo_n4 %>% 
  
  group_by(combo) %>% 
  nest() %>%
  
  mutate(
    
    RF_results = map(data, fun_RF_classifier),
    
    auc = RF_results %>% 
      
      map_dbl(function(x) 
        
        MESS::auc(
          unique(x$RF_Xval_x), 
          unique(x$RF_Xval_y)
        )
        
      ),
    
    auc_lb = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.025)),
    
    auc_ub = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.975))
    
  )

# RF classifier fitted to all 5-way combinations

ab_combo_n5 <- ab_combo_n5[1:10,] # SIC. for demo run only - selects first 10 combinations. silence for exhaustive analysis

ab_combo_n5 <- ab_combo_n5 %>% 
  
  group_by(combo) %>% 
  nest() %>% 
  
  mutate(
    
    RF_results = map(data, fun_RF_classifier),
    
    auc = RF_results %>% 
      
      map_dbl(function(x) 
        
        MESS::auc(
          unique(x$RF_Xval_x), 
          unique(x$RF_Xval_y)
        )
        
      ), 
    
    auc_lb = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.025)),
    
    
    auc_ub = RF_results %>% 
      
      map_dbl(function(x) quantile(x$RF_Xval_auc, probs = 0.975)) 
  )

# shut down parallel workers

stopCluster(cl)

###############################################################
###############################################################

ROC_ab_combo_n2 <- ab_combo_n2 %>%
  
  ungroup() %>%
  
  slice_max(auc, n = 100) %>% 
  
  unnest(cols = c(data, RF_results)) %>%
  
  dplyr::select(- RF_Xval_auc)

colnames(ROC_ab_combo_n2) <- c("combo", "ab1", "ab2", "fpr", "sens", "AUC", "AUC_lb", "AUC_ub")

ROC_ab_combo_n2 <- ROC_ab_combo_n2 %>% 
  
  mutate(
    
    sens = 1-sens, #edited from  1-sens
    fpr  = 1-fpr
    
  )

top_combos_n2 <- ROC_ab_combo_n2 %>%
  
  distinct(combo, .keep_all = T) %>% 
  
  select(-fpr, -sens)

#################################################################

# 3 way

ROC_ab_combo_n3 <- ab_combo_n3 %>% 
  
  ungroup() %>%
  
  slice_max(auc, n = 100) %>% 
  
  unnest(cols = c(data, RF_results)) %>%
  
  dplyr::select(- RF_Xval_auc)

colnames(ROC_ab_combo_n3) <- c("combo", "ab1", "ab2", "ab3",
                               "fpr", "sens", "AUC", "AUC_lb", "AUC_ub")

ROC_ab_combo_n3 <- ROC_ab_combo_n3 %>% 
  
  mutate(
    
    sens = 1-sens, #edited from  1-sens
    fpr  = 1-fpr
    
  )

top_combos_n3 <- ROC_ab_combo_n3 %>%
  
  distinct(combo, .keep_all = T) %>% 
  
  select(-fpr, -sens)

#################################################################

# 4 way

ROC_ab_combo_n4 <- ab_combo_n4 %>% 
  
  ungroup() %>%
  
  slice_max(auc, n = 100) %>% 
  
  unnest(cols = c(data, RF_results)) %>%
  
  dplyr::select(- RF_Xval_auc)

colnames(ROC_ab_combo_n4) <- c("combo", "ab1", "ab2", "ab3", "ab4",
                               "fpr", "sens", "AUC", "AUC_lb", "AUC_ub")

ROC_ab_combo_n4 <- ROC_ab_combo_n4 %>% 
  
  mutate(
    
    sens = 1-sens, #edited from  1-sens
    
    fpr  = 1-fpr
    
  )

top_combos_n4 <- ROC_ab_combo_n4 %>%
  
  distinct(combo, .keep_all = T) %>% 
  
  select(-fpr, -sens)

#################################################################

# 5 way

ROC_ab_combo_n5 <- ab_combo_n5 %>% 
  
  ungroup() %>%
  
  slice_max(auc, n = 100) %>% 
  
  unnest(cols = c(data, RF_results)) %>%
  
  dplyr::select(- RF_Xval_auc)

colnames(ROC_ab_combo_n5) <- c("combo", "ab1", "ab2", "ab3", "ab4", "ab5",
                               "fpr", "sens", "AUC", "AUC_lb", "AUC_ub")

ROC_ab_combo_n5_test <- ROC_ab_combo_n5 %>% 
  
  mutate(
    
    fpr  = 1-fpr
    
  )

top_combos_n5 <- ROC_ab_combo_n5 %>%
  
  distinct(combo, .keep_all = T) %>%
  
  select(-fpr, -sens)


#################################################################
#################################################################



#################################################################

# Plotting the RF classifier results (Fig. 4b)

#################################################################


line_seq <- c(0.2, 0.4, 0.6, 0.8)

pdf("results/Fig_4B.pdf")

par(mfrow=c(1,1))
par(mar=c(3,3,4,1.5))
par(mgp=c(1.75, 0.6,0))


plot(
  
  x=c(0,1), y=c(0,1), 
  type='l', lty="dashed",
  xlim=c(0,1.002), ylim=c(0,1.002),
  xaxs='i', yaxs='i', 
  xaxt='n', yaxt='n',
  xlab="False positive rate (100 - specificity)", 
  ylab="True positive rate (sensitivity)",
  main=expression(
    paste("ROC curve for classifying recent ", 
          italic("P. falciparum "), 
          "infections (5 responses)")
  ),
  cex.lab=1.1, cex.axis=0.8, cex.main=1.2)


for(i in 1:4)
{
  points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
  points(x=rep(1-line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}


axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1),
     labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 )

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1),
     labels=c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.0 )


top_10_cols <- c(rainbow(10))

for (i in 1:10) {
  points(
    
    x = ROC_ab_combo_n5[which(ROC_ab_combo_n5$combo == top_combos_n5[[i,1]]),]$fpr,
    y = ROC_ab_combo_n5[which(ROC_ab_combo_n5$combo == top_combos_n5[[i,1]]),]$sens,
    
    type = 's', 
    lwd =1, 
    col = top_10_cols[i])
}

legend(x="bottomright",
       bty='n', cex=0.7,
       fill=top_10_cols,
       border=top_10_cols,
       legend = c(
         paste0(top_combos_n5$ab1[1:10], 
                ", ", 
                top_combos_n5$ab2[1:10],
                ", ", 
                top_combos_n5$ab3[1:10],
                ", ", 
                top_combos_n5$ab4[1:10],
                ", ", 
                top_combos_n5$ab5[1:10],
                "; AUC = ",
                as.character(round(top_combos_n5$AUC[1:10], 3)),
                " (95% CI: ", 
                as.character(round(top_combos_n5$AUC_lb[1:10], 3)),
                " - ",
                as.character(round(top_combos_n5$AUC_ub[1:10], 3)), 
                ")")
       )
       
)

dev.off()

##############################################################
##############################################################





##############################################################
##############################################################
##############################################################



#                          Part 5.
#
#           Analysing Ab kinetic model results



##############################################################
##############################################################
##############################################################


# load antibody kinetic model output: 1-year reduction in antibody reactivity

one_Y_r <- read_csv(file="data/ab_reduction.csv") # for variable descriptions see Suplementary File S1: Supplementary Table S7.

# transpose data

antigen_names_tmp <- one_Y_r$antigen

one_Y_r <- as.data.frame(t(one_Y_r[,2:65]))

colnames(one_Y_r) <- antigen_names_tmp

one_Y_r$study_idn <- rownames(one_Y_r)


# gather by antigen

one_Y_r_gather <- one_Y_r %>% 
  
  gather(key = "antigen", value = "reduc_1y", AARP:TLP)

# calculate median reduction by antigen

one_Y_r_gather <- one_Y_r_gather %>% 
  
  group_by(antigen) %>% 
  
  mutate(
    
    reduc_1y = 1-reduc_1y,
    med_reduc= median(reduc_1y, na.rm = T),
    col = 1, 
    col = if_else(med_reduc>0.25, 2, col),
    col = if_else(med_reduc>0.5, 3, col), 
    col = if_else(med_reduc>0.75, 4, col)
    
  )


# box-plot 1y reduction - all antigens

antigen_names_ordered <- one_Y_r_gather %>% ungroup() %>%
  
  arrange(med_reduc) %>% 
  
  dplyr::select(antigen) %>% 
  
  distinct() %>% 
  
  map(as.character)

# create character vector with names of top 10 individually most infomative antibody responses

top_10_abs <- as.character(table_s2[1:10,1]$antigen) 

# create color vector for plotting

top_col_vect <- ifelse(antigen_names_ordered$antigen %in% top_10_abs, "red","black") # SIC object name


###############################################################

# Plot reduction in antibody levels over time (Fig. 6)

###############################################################


plot_reduct_1y <- one_Y_r_gather %>% 
  
  ggplot(
    aes(
      y = reduc_1y, 
      x = reorder(antigen, med_reduc), 
      fill = factor(col)
    )
  ) +
  
  geom_boxplot() +
  
  xlab("Antigen")+
  ylab("Reduction in antibody reactivity after 1 year (%)")+
  
  theme_bw()+
  
  theme(
    
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(2, "cm"), 
    legend.title = element_text(size=24),
    legend.text = element_text(size = 16),
    
    axis.title = element_text(size = 28),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(
      angle = 90, colour = top_col_vect, size = 16,
      vjust = 0.5, hjust = 0
    )
    
  )+
  
  scale_fill_manual(
    
    name = "Median reduction (%)", 
    labels = c("<25%", "25-50%", "50-75%", ">75%"),
    values = viridis(4)
    
  )

ggsave(plot_reduct_1y, file = "results/Fig_6.pdf", height = 18, width = 26, scale = 1)

###############################################################
###############################################################

# calculate antigen-specific mean reduction in ab levels
# and bind data with threshold ab level classifier resutls (table_s2)

reduc_summary_auc <- one_Y_r_gather %>% 
  
  group_by(antigen) %>%
  
  summarise(
    
    n       = n(),
    mean    = mean(reduc_1y), 
    std     = sd(reduc_1y),
    cv      = std/mean*100,
    ste     = std/sqrt(n), 
    mean_ub = mean+1.96*ste, 
    mean_lb = mean-1.96*ste
    
  ) %>%
  
  left_join (table_s2, by = "antigen")

top_ant_summary <- reduc_summary_auc %>% 
  
  filter (antigen %in% top_10_abs)

###############################################################

# Plot mean 1-year reduction in ab levels vs.
# threshold ab level classifier AUC (Fig. 7)

###############################################################

plot_mean_reduct_auc <-  reduc_summary_auc %>%
  
  ggplot(aes(x = auc, y = mean*100))+
  
  # facet_wrap (~exposure_group)+
  geom_errorbarh(
    aes(y = mean*100, xmin = auc_lb, xmax=auc_ub), 
    color="darkblue", alpha=0.4
  )+
  
  geom_errorbar(
    aes(x = auc, ymin = mean_lb*100, ymax = mean_ub*100), 
    color="darkblue", 
    alpha=0.4)+
  
  geom_point(size = 3, shape = 21, fill = "darkblue", alpha = 0.7)+
  
  geom_point(data = top_ant_summary, 
             aes(y = mean*100, x = auc), 
             col="black", size = 4)+
  
  geom_point(data = top_ant_summary, 
             aes(y = mean*100, x = auc, col = antigen),
             size = 3)+
  
  geom_label_repel(
    data = top_ant_summary, 
    aes(label = antigen),
    box.padding = 1.3, 
    point.padding = 0.6, 
    segment.alpha = 0.4,
    force = 1, max.iter = 10000,
    max.overlaps = 30
  )+
  
  xlab("Classifier AUC")+
  ylab("Mean reduction in antibody reactivity after 1 year (%)")+
  
  scale_colour_manual(values = (top_10_cols), limits = top_10_abs)+
  
  ylim(0,100)+
  xlim(0.32, 1)+
  
  theme_bw()+
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
    
  )

ggsave(plot_mean_reduct_auc, file = "results/Fig_7.pdf", width = 8, height = 8)

###############################################################
###############################################################



###############################################################

#  Beta regression models (fractional response)

###############################################################

# load antibody kinetic model output: peak antibody levels

peak_ab <- read_csv(file = "data/ab_peak.csv") # for variable descriptions see Suplementary File S1: Supplementary Table S6.

# transpose data

peak_ab <- as.data.frame(t(peak_ab[,2:65]))
colnames(peak_ab) <- antigen_names_tmp
peak_ab$study_idn <- rownames(peak_ab)

peak_ab_gather <- peak_ab %>% 
  
  gather(key = "antigen", value = "peak_MFI", AARP:TLP)

# subset meta data to include in analysis

beta_meta_data <- full_data %>% 
  
  filter(t_id != "Neg. control") %>% ungroup()%>%
  
  dplyr::distinct(study_idn, .keep_all = T) %>%
  
  select(study_idn, exposure:years_since_endemic)

# join data to analyse

beta_input_data <- one_Y_r_gather %>% 
  
  select(study_idn, antigen, reduc_1y) %>% 
  
  left_join(peak_ab_gather, by = c("study_idn", "antigen")) %>%
  
  left_join(beta_meta_data, by = "study_idn") %>%
  
  mutate(
    
    exposure = factor(exposure, levels = c("Primary infected", "Previously exposed")) # reset factor levels for subset data
  
    )

# nest data frame

beta_input_nest <- beta_input_data %>% 
  
  group_by(antigen) %>%
  
  nest()

# fit beta regression models

beta_input_nest <- beta_input_nest %>%
  
  mutate(
    
    model_1 = data %>% map( 
      
      ~betareg(reduc_1y ~ log(peak_MFI)+exposure, data = .x, link = 'logit')
      
    ), 
    
    coefs_1 = model_1 %>% map(broom::tidy)
    
  )


beta_model_1 <- beta_input_nest %>% 
  
  select(antigen, coefs_1) %>% 
  
  unnest(cols = coefs_1) %>%
  
  ungroup() %>% 
  
  mutate(
    
    q.value = p.adjust(p.value, method = "BH"), 
    
    term = factor(term, labels = c("Intercept", 
                                   "Phi (precision parameter)",
                                   "Previously exposed", 
                                   "log (peak MFI)"))
    
  ) %>%
  
  select(-component, -statistic)


###############################################################
###############################################################


###############################################################
###############################################################
###############################################################

# END

###############################################################
###############################################################
###############################################################
