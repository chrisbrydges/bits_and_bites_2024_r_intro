# Same code as 3_complete_code, except much shorter and far fewer comments
#
library(readxl)
library(janitor)
library(rstatix)
library(mixOmics)
library(EnhancedVolcano)
library(ggpubr)
library(tidyverse)

# Read the Excel file
data <- read_excel("original_data.xlsx")

# Identified compounds only
df <- data[1:116, 1:111]  # Columns 1 to 111 contain metabolite name and intensity
chem <- data[1:116, c(1, 112:117)]  # Chemical metadata

# Data preparation
df <- df %>% 
  column_to_rownames(var = "Metabolite Name") %>%
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>%
  filter(Factors != "QC") %>% 
  mutate(Factors = factor(Factors, levels = c("MECFS", "Control")),
         across(xylose:X1.5.anhydroglucitol, ~as.numeric(.x)))
  
chem <- chem %>% 
  filter(`Metabolite Name` != "Factors")

colnames(df)[3:117] <- chem$`Metabolite Name`


# UNIVARIATE ANALYSES
# t-tests
univariate_results <- df_transform %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity") %>% 
  group_by(Metabolite) %>% 
  t_test(Intensity ~ Factors, p.adjust.method = "none") %>% 
  mutate(FDR = p.adjust(p, method = "fdr")) %>% 
  select(Metabolite, p, FDR)

# fold changes
fold_changes <- df %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity") %>% 
  group_by(Factors, Metabolite) %>% 
  summarise(M = median(Intensity)) %>% 
  pivot_wider(id_cols = Metabolite, names_from = Factors, values_from = M) %>% 
  mutate(fc = MECFS/Control,
         log2fc = log2(fc)) %>% 
  select(Metabolite, fc, log2fc)

# Join together and with chemical metadata
univariate_results <- inner_join(univariate_results, fold_changes, by = "Metabolite") %>% 
  inner_join(chem, ., by = c("Metabolite Name" = "Metabolite"))

# Prepare files for metaboanalyst pathway analysis
all_keggs <- univariate_results %>% 
  filter(!is.na(KEGG)) %>% 
  select(KEGG)

sig_keggs <- univariate_results %>% 
  filter(!is.na(KEGG) & p < 0.05) %>% 
  select(KEGG)

write.table(all_keggs, file = "all_keggs.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sig_keggs, file = "sig_keggs.csv", row.names = FALSE, col.names = FALSE, sep = ",")



# MULTIVARIATE ANALYSIS


plsda_res <- plsda(X = df[, 3:117], Y = df$Factors, ncomp = 10)


# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(plsda_res, 
          comp = 1:2, 
          group = df$Factors, 
          pch = df$Factors,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, 
          title = "PLSDA with confidence ellipses")

# Undergo performance evaluation in order to tune the number of components to use
perf_plsda <- perf(plsda_res, 
                   validation = "Mfold", # can be changed to "loo"
                   folds = 5, 
                   nrepeat = 50, # use repeated cross-validation
                   progressBar = TRUE, 
                   auc = TRUE) # include AUC values

# Plot the outcome of performance evaluation across all ten components
plot(perf_plsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")


# Let's look to see how many components we should keep in our PLSDA model
perf_plsda$choice.ncomp

# 2 seems to be the unanimous decision.

optimal.ncomp <- 2

# Let's create the final model 
final_plsda <- plsda(X = df[, 3:117], 
                     Y = df$Factors, 
                     ncomp = optimal.ncomp)

# The final model object has all the information in it needed to run performance
# analyses and create plots

perf_final_plsda <- perf(final_plsda, 
                         validation = "Mfold", # can be changed to "loo"
                         folds = 5, 
                         nrepeat = 50, # use repeated cross-validation
                         progressBar = TRUE, 
                         auc = TRUE)

# Plot samples from final model
plotIndiv(final_plsda, 
          comp = 1:2, 
          group = df$Factors, 
          pch = df$Factors, # colour by class label
          ellipse = TRUE, # include 95% confidence ellipse 
          legend = TRUE, 
          title = "PLS-DA on MECFS, comp 1 & 2")

# Plot variables from final model, but only show those with correlations with 
# components 1 and 2 > 0.5 (this really reduces the number of metabolites on the
# plot)
plotVar(final_plsda, 
        comp = 1:2, 
        cutoff = 0.5, 
        cex = 3)


# Calculate VIP scores for each metabolite and component, then join 
# them to the univariate_results and save the whole analysis as "results.csv"

vip_scores <- vip(final_plsda) %>% 
  rowname_to_column(var = "Metabolite Name")

results <- inner_join(univariate_results, vip_scores, by = "Metabolite Name")

write.csv(results, "results.csv", row.names = FALSE)




# VISUALIZATIONS

# Volcano plot
p1 <- EnhancedVolcano(univariate_results, # the data frame with our results in it
                      lab = univariate_results$`Metabolite Name`, # where to find the labels for interesting data points
                      x = "log2fc", # log2 fold change values
                      y = "p", # p values (can be FDR corrected if you want)
                      pCutoff = 0.05, # The cut-off for an 'interesting' p value
                      FCcutoff = 0, # the cut-off for an 'interesting' fold change
                      xlim = c(-1, 1), # how wide do you want the x-axis (min, max)
                      ylim = c(0, -log10(10e-5)), # how tall do you want the y-axis (min, max)
                      title = "MECFS vs. Controls", # what is the title
                      subtitle = NULL, # get rid of subtitle
                      caption = NULL, # get rid of caption
                      legendPosition = "bottom") # put legend at bottom

ggsave("VolcanoPlot.png", p1, width = 6, height = 5)

# Create boxplots for significant metabolites only
sig_metabolites <- univariate_results %>% 
  filter(p < 0.05) %>% 
  pull(`Metabolite Name`)

significant_df <- df %>% 
  select(ID, Factors, all_of(sig_metabolites))

# Make a new directory called "SignificantBoxplots"
dir.create(paste0(getwd(), "/SignificantBoxplots/"))



for (i in 3:ncol(significant_df)) { 
  
  # Create a two column data frame with Factors as one column, and the current metabolite as the other
  temp_df <- significant_df %>% 
    select(2, all_of(i)) 
  
  # Save the current metabolite name, and rename the second column of temp_df as "Intensity"
  current_metabolite <- colnames(temp_df)[2]
  colnames(temp_df)[2] <- "Intensity"
  
  # Use temp_df as the data frame in the plot, with "Intensity" as the y axis value
  temp_p <- ggplot(temp_df, aes(x = Factors, y = Intensity, fill = Factors)) +
    geom_boxplot() + 
    labs(title = current_metabolite, x = "Group", y = "Intensity") + 
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 16))
  # Save the plot in the boxplots folder with the current metabolite as the file name
  ggsave(paste0(getwd(), "/SignificantBoxplots/", current_metabolite, ".png"), temp_p, width = 6, height = 5)
}



# PCA plot
pca_res <- pca(df[, 3:117], center = TRUE, scale = TRUE)
p3 <- plotIndiv(pca_res, group = df$Factors, legend = TRUE, ellipse = TRUE, pch = df$Factors)$graph
ggsave("PCA_plot.png", p3, width = 6, height = 5)


# PLSDA plot
p4 <- plotIndiv(final_plsda, 
                comp = 1:2, 
                group = df$Factors, 
                pch = df$Factors, # colour by class label
                ellipse = TRUE, # include 95% confidence ellipse 
                legend = TRUE, 
                title = "PLS-DA on MECFS, comp 1 & 2")$graph
ggsave("PLSDA_plotIndiv.png", p4, width = 6, height = 5)

