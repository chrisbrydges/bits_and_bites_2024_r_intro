# Load necessary libraries
#
# Libraries are add-ons that contain functions that we can use, so we don't have 
# to manually create the functions ourselves.
#
# The order that libraries are called can affect whether your code will run. In 
# my experience, I find that my code will run if I call `tidyverse` last.
library(readxl)
library(janitor)
library(rstatix)
library(mixOmics)
library(EnhancedVolcano)
library(ggpubr)
library(tidyverse)

# Read the Excel file - it should be in the folder that you've made a project in
data <- read_excel("original_data.xlsx")

# We should have 264 observations of 117 columns. We can have a look at our 
# loaded data file by typing View(data) into the console pane.  


# We need to split our data file into two, as the chemical metadata is in the 
# same data frame as the GC-MS data. We'll refer to the GC-MS data as `df` 
# (short for data frame) and the chemical metadata as `chem`.

# One way of subsetting the data is to manually select the column numbers and/or
# row numbers that you want to keep/get rid of. The format in R is as follows:
# data[first_row_to_keep:last_row_to_keep, first_col_to_keep:last_col_to_keep]
# This does require us to know exactly which column numbers we want.
df <- data[, 1:111]  # Columns 1 to 111, and all rows
chem <- data[, c(1, 112:117)]  # Columns 1 and 112 to 117, and all rows


# The method above uses what's known as `base R`. This means that no additional 
# libraries are needed for this code to run. Another way of subsetting is using 
# the `select` function from the `dplyr` package, which is part of the tidyverse.
# The `select` function allows us to choose specific columns that we want to keep
# and/or exclude.

df2 <- data %>%   
  # This line of code means "create df2 by starting with data and then..."
  select(`Metabolite Name`:`Pool 10`) 
# ... selecting columns `Metabolite Name` up to `Pool 10` (inclusive)

# A shortcut for the %>% is Ctrl+Shift+M (Windows) or Cmd+Shift+M (Mac)

# In this case, df and df2 should be identical as we've asked R to do the same 
# thing in two different ways. We'll use tidyverse functions quite a lot in this
# course, but for now let's delete the df2 object.
rm("df2")

# Let's also get rid of unidentified metabolites. We can look at df and chem to 
# see where we need to cut off the data before we enter those values into the 
# two lines of code below.
df <- df[1:116, ]
chem <- chem[1:116, ]




# At this point, our df data frame is the wrong way round: the rows need to be 
# columns and the columns need to be rows. 

df <- df %>% # Make an object called df by starting with df and then...
  column_to_rownames(var = "Metabolite Name") %>%  #... move the column "Metabolite Names" to be the row names then...
  t() %>% # ... transpose the data (which turns it into a matrix) and then...
  data.frame() %>% # ... convert it back to a data frame and then...
  rownames_to_column(var = "ID") # move the row names to be a column called "ID"

# Our data is now in the right shape for analysis - this is called "tidy" data, 
# where each row is a single observation and each column is a single variable. 
# To be more specific, this is wide format (we'll be using long format later on).

# There's one more issue with our data at this point - R currently doesn't think 
# that our metabolite intensity values are numbers. We can look at this here:
str(df)

# Each row of the output has "chr" listed next to it, which means that R thinks 
# it's a character variable, which is a string of text. We can easily change this
# using the `mutate` function. We want to change "Factors" to be a factor, and 
# all of our metabolites to be numeric.

# The mutate function can be used to create a new variable, or overwrite an 
# existing one (or multiple). We'll be overwriting in this case.

# We use the same pipe notation as before
df <- df %>% 
  mutate(Factors = factor(Factors, levels = c("MECFS", "Control")),
         # The above line is telling R that Factors should be a factor type 
         # variable, and that MECFS should be the first group in it, otherwise 
         # it defaults to alphabetical order
         across(xylose:X1.5.anhydroglucitol, ~as.numeric(.x)))
         # across is what we use if we want to make the same changes to multiple
         # variables at once. In this case, the only change we are making is to 
         # tell R that the data in these columns is numeric. Make sure to include
         # the ~ (tilda sign) and .x (not just x) when using the across function.


# Our data is now almost 100% ready to use - except we have 10 QC samples that we
# don't want, so we can remove them by using the `filter` function. This is like
# the `select` function, except we're now removing rows that we don't want. In 
# this case, we don't want rows where the Factors column contains "QC".
df <- df %>% 
  filter(Factors != "QC")

# The != means "does not equal". If we want something that does equal, we use ==
# i.e., two equals signs, not one. We can see that our filtering has worked 
# because when we run the table function below, it shows that we only have MECFS
# and Controls in the Factors column, and not QC (or anything else):
table(df$Factors)


# Two small last things: first, there's a row in the chem object that says 
# "Factors" and nothing else. We also want to filter that row out.

chem <- chem %>% 
  filter(`Metabolite Name` != "Factors")

# Second, look at the column names for df:
colnames(df)

# It looks like spaces have been replaced with dots, and metabolites that start
# with a number have an X put in front of them now. We can change that by 
# specifying column names for these columns. HINT: Think about the column numbers
# we want to replace, and where we can get our replacement column names from
colnames(df)[3:117] <- chem$`Metabolite Name`

# Now, after a lot of preparation, we're ready to actually work on our analyses!



# UNIVARIATE ANALYSES


# This first part here is optional - the code below will perform a log 
# transformation on all of the metabolite intensities to make the data normally 
# distributed, and we will also perform a # half-minimum missing value imputation. 
# Given this is GC-MS data, there are no missing values, so this is moot. However, 
# the code is interesting to look at and understand how it works. Plus, if you try 
# this on LC-MS data, it's very likely there will be some missing data.

# Like we did earlier when we changed our data to be numeric, we want to 
# `mutate` `across` multiple columns. This is because we want to make these 
# changes to every metabolite.

# This time, we're not going to # overwrite the df object because we don't want 
# to use log-transformed data when calculating fold changes later on. We'll call 
# the new object df_transform. To run a half-minimum missing value imputation, 
# we need to calculate the minimum value of each metabolite AFTER removing the 
# missing values, then multiply that value by 0.5:

df_transform <- df %>% 
  mutate(across(xylose:`1,5-anhydroglucitol`, ~replace_na(.x, 0.5 * min(.x, na.rm = TRUE))))

# The above code can be explained in a few steps, as it's actually four functions
# being run inside of each other. Starting from the left: 
# 1. mutate means we want to create/overwrite a variable
# 2. across means we want to work on multiple variables (which are xylose to 
# `1,5-anhydroglucitol`, inclusive) 
# 3. replace_na is going to replace the missing values with a value of our choice.
# When R has a missing value, it's written in to the data as NA.
# 4. min calculates the minimum value of the variable we're looking at, and the 
# na.rm = TRUE means "calculate the minimum after removing any missing data", 
# otherwise the minimum value will be a missing value.
# Once the minimum value is calculated, it's multiplied by 0.5 and this value 
# replaces any missing values for that particular variable. The function then 
# moves on to the next variable in the across function.


# Next, the log transformation. 

df_transform <- df_transform %>% 
  mutate(across(xylose:`1,5-anhydroglucitol`, ~log10(.x))) # Again, don't forget the ~ and .x





# Now, we can start by conducting independent-samples t-tests on the data. We're
# going to run a series of steps one at a time to show how each of them work, and 
# then re-run it all at once to show how it does it from scratch. 

# First, we'll create a new object called univariate_results that is based on 
# df_transform:
univariate_results <- df_transform

# If you look over in the environment pane, you'll see that univariate_results
# just appeared there, and it's the same size as df_transform. 

# Next, we're going to change our data from wide format to long format. I 
# previously explained that wide format is one row per observation, and one column
# per variable. Long format is one row per data point where there will be several
# columns saying what the data point is, and then the value itself. This uses
# a function called pivot_longer. We need to tell it which columns are going to be 
# used, what we should call the "name" column, and what we should call the "value"
# column". I'm calling them "Metabolite" and "Intensity", respectively.

univariate_results <- univariate_results %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity")

# If you look at the environment pane again, you'll see that the shape of 
# univariate_results has changed a lot - there are now many more rows and far 
# fewer columns. We have ID, Factors, Metabolite, and Intensity as our four columns.

# This next step doesn't visually do anything, but it tells R that we want to 
# split up what we do after it. This is important because the next thing we do is 
# run the t-test on each metabolite, and we don't want to run it on all the data 
# at once - it needs to be run on each metabolite separately.

univariate_results <- univariate_results %>% 
  group_by(Metabolite)

# Now, we can finally find out which metabolites are significantly altered in our
# data! We run a t-test. To do this, we specify the formula, which is 
# dependent_variable ~ independent_variable (so Intensity ~ Factors in our study).
# Note that I'm also explicitly stating that I don't want my p values to be 
# adjusted for multiple comparisons here. I'm going to do that later.

univariate_results <- univariate_results %>% 
  t_test(Intensity ~ Factors, p.adjust.method = "none")

# Now, let's add another column to the univariate_results that has the FDR-corrected
# p values too (by using the mutate function)
univariate_results <- univariate_results %>% 
  mutate(FDR = p.adjust(p, method = "fdr"))


# We don't really need all of these columns. In fact, you could argue that we only
# need the metabolite name, the raw p value, and the FDR-corrected p value. Let's
# select those columns to keep, and get rid of the others.
univariate_results <- univariate_results %>% 
  select(Metabolite, p, FDR)

# We now have the statistical significance of our univariate results done!
# Let's re-run that as one big operation but save it as univariate_results2 so we
# can double check that everything comes out the same

univariate_results2 <- df_transform %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity") %>% 
  group_by(Metabolite) %>% 
  t_test(Intensity ~ Factors, p.adjust.method = "none") %>% 
  mutate(FDR = p.adjust(p, method = "fdr")) %>% 
  select(Metabolite, p, FDR)

# Pretty straightforward once you know what you're doing!


# Next, let's calculate the fold change between the two groups for each metabolite.
# We'll follow the same process of going through each individual step, then re-run
# everything in one chunk.

# To start, we want to create an object called fold_changes by copying df. We do
# NOT want to use df_transform, as that has log-transformed data.
fold_changes <- df

# As with the t-test results, we pivot our data into long format to look at the 
# metabolite intensity
fold_changes <- fold_changes %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity")

# We use group_by again (which visually doesn't change anything), but this time 
# we group by both metabolite and Factors, as we're going to perform a calculation
# on every metabolite FOR BOTH groups

fold_changes <- fold_changes %>% 
  group_by(Factors, Metabolite)

# That calculation is calculating the median of each metabolite for each group.
# Our fold_changes data frame currently has three columns: Factors, Metabolite, 
# and M. Each row is the median intensity for a metabolite for one group.
fold_changes <- fold_changes %>% 
  summarise(M = median(Intensity))

# The next step is to pivot the data again, but to make it wide instead of long.
# The function needs to know our id column, which column the names are coming from,
# and which column the values are coming from. We want to have a data frame where
# each row is a metabolite, and there's one column for median MECFS intensity, and 
# another for median control intensity.
fold_changes <- fold_changes %>% 
  pivot_wider(id_cols = Metabolite, names_from = Factors, values_from = M)

# Next, let's create two new variables. The first will be called fc, where we 
# calculate the fold change (group X / group Y). We'll use the control group as 
# the reference group, which is the denominator in the ratio. The second variable
# can be called log2fc, and will be the log2-fold change. We'll use this for the 
# volcano plot. HINT: You can use a variable you've just made in mutate in the 
# same mutate function.
fold_changes <- fold_changes %>% 
  mutate(fc = MECFS/Control,
         log2fc = log2(fc))

# Lastly, we only need the metabolite name, fold change, and log2-fold change
# from this analysis, so let's select those three columns only
fold_changes <- fold_changes %>% 
  select(Metabolite, fc, log2fc)


# Let's put that all together in a single code chunk now and call the object 
# fold_changes2 so we can double check it worked the same
fold_changes2 <- df %>% 
  pivot_longer(cols = xylose:`1,5-anhydroglucitol`, names_to = "Metabolite", values_to = "Intensity") %>% 
  group_by(Factors, Metabolite) %>% 
  summarise(M = median(Intensity)) %>% 
  pivot_wider(id_cols = Metabolite, names_from = Factors, values_from = M) %>% 
  mutate(fc = MECFS/Control,
         log2fc = log2(fc)) %>% 
  select(Metabolite, fc, log2fc)


# We now have run all of the univariate analyses! But, the t-test results, the 
# fold changes, and metabolite metadata are all in separate data frames. We want 
# to put them all into one so we can save that as a single csv file. We can do 
# this with a join. A join adds columns from data frame y to data frame x, matching
# rows based on the key or keys.

# There are two main issues to consider with joins: What type of join? and What is
# my key?

# The four types of join we can easily run here are inner, left, right, and full. 
# An inner join only keeps observations from x that have a matching key in y. 
# The most important property of an inner join is that unmatched rows in either 
# input are not included in the end result. A left join keeps all observations in 
# x. A right join keeps all observations in y. A full join keeps all observations 
# in x and y.

# Because we've prepared ahead of time, we can actually use any of these joins 
# and it won't make any difference IN THIS SPECIFIC CASE. This is because in 
# univariate_results, fold_changes, and chem, all the data frames have the exact 
# same metabolite names and no unique ones. This is why I spend so much time at 
# the preparation stage before we get into the actual analysis.

# Second question is "what is my key?". This means "which column/s am I using to 
# match my rows across my data frames?". In our case, univariate_results and 
# fold_changes both have a column called "Metabolite" and chem has a column called
# "Metabolite Name" that contains the metabolites we've analyzed. This will be our
# key.

# Let's put everything into the univariate_results data frame. Note how the order
# I enter data frames into the join function will affect the order of the columns
# at the end.
univariate_results <- inner_join(univariate_results, fold_changes, by = "Metabolite") %>% 
  inner_join(chem, ., by = c("Metabolite Name" = "Metabolite"))

# In the first line of that code chunk, our key column has the same name in both 
# data frames. In the second, we want to put the columns from chem to the left of 
# our joined univariate_results/fold_changes (denoted by the dot), and then we 
# specify which column in x matches to which column in y.

# We now have a single data frame that contains all of our univariate results - 
# let's save it as a csv file
write.csv(univariate_results, file = "UnivariateResults.csv", row.names = FALSE)


# As a small aside, we can also save a couple of files that will be ready for 
# entry into MetaboAnalyst's pathway analysis module. This needs two files - 
# one of the statistically significant metabolites (as identified by KEGG IDs), 
# and another of every metabolite that was detected in the data (as identified by
# KEGG IDs). We can select columns and filter rows from univariate_results to do
# this:
all_keggs <- univariate_results %>% 
  filter(!is.na(KEGG)) %>% 
  select(KEGG)

sig_keggs <- univariate_results %>% 
  filter(!is.na(KEGG) & p < 0.05) %>% 
  select(KEGG)

write.table(all_keggs, file = "all_keggs.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sig_keggs, file = "sig_keggs.csv", row.names = FALSE, col.names = FALSE, sep = ",")



# Univariate analyses are done and saved, as well as files ready for pathway analysis!





# MULTIVARIATE ANALYSIS

# The typical multivariate analysis run on metabolomics data is PLSDA. This is 
# run through the mixOmics package - they have tutorials on their website at
# http://mixomics.org


plsda_res <- plsda(X = df[, 3:117], Y = df$Factors, ncomp = 10)


# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(plsda_res, 
          comp = 1:2, 
          group = df$Factors, 
          pch = df$Factors,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, 
          title = "PLSDA with confidence ellipses")


# The number of components to use is a crucial decision and is dictated by the 
# performance of the PLSDA model, i.e., its ability to correctly classify novel 
# samples. The perf() function is used for this. This is done with repeated 
# cross-validation. Based on the output of this function, the optimal number of 
# components to use can be identified.

# A five-fold, 10 repeat cross-validation procedure is utilized here. Generally, 
# for datasets with numerous samples, at least 10 folds is recommended. 3 or 5 
# folds is appropriate for smaller datasets and those with minimal samples should 
# use Leave-One-Out (LOO) validation. Consider using 50-100 repeats to reduce the 
# impact of the randomly allocated folds during each repeat.

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

# Kind of looks like there's some fatty acid alterations occurring?

auc_plsda <- auroc(final_plsda, roc.comp = 2, print = FALSE)
# ROC and AUC criteria may not be particularly insightful, or be in agreement 
# with the PLSDA performance, as the prediction threshold in PLS-DA is based on 
# specified distance. AUROC curves use a cutoff that maximises specificity and 
# sensitivity rather than this distance and hence should be used a merely a 
# complementary tool. 

# Lastly, let's calculate VIP scores for each metabolite and component, then join 
# them to the univariate_results and save the whole analysis as "results.csv"

vip_scores <- vip(final_plsda) %>% 
  rowname_to_column(var = "Metabolite Name")

results <- inner_join(univariate_results, vip_scores, by = "Metabolite Name")

write.csv(results, "results.csv", row.names = FALSE



# At this point, we've now completed our multivariate analysis, and found that
# our PLSDA model is bad at predicting MECFS vs. healthy controls. 



# VISUALIZATIONS


# Let's go back to the univariate analyses. We're going to make two kinds of plots
# here. First, we'll make a volcano plot, then we'll make boxplots. For the 
# boxplots, we can get a bit more advanced by telling R to only make it if the 
# univariate analysis was statistically significant.

# Volcano plot first though. A volcano plot shows the log2 fold change on the 
# x-axis and the -log10(p value) on the y-axis. We can use EnhancedVolcano to make
# and customize these plots extremely easily

# We'll save this plot as an object called p1. I'm going to put each argument on
# it's own line so I can briefly describe how it can be customized.
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

# There are plenty of other features that can be altered e.g., label size, 
# type ?EnhancedVolcano to get help on how to change things
ggsave("VolcanoPlot.png", p1, width = 6, height = 5)

# Let's move on to box plots. Before we create boxplots for all of them, let's 
# learn how to make just one in ggplot2. ggplot2 is the data visualization 
# package of the tidyverse. As a general rule, it follows a pattern of:
# 1. Defining the data
# 2. Specifying the type of graph you want
# 3. Making changes you want to it

# For example, we can make a simple box plot for xylose as follows:
ggplot(df, aes(x = Factors, y = xylose)) + 
  geom_boxplot()

# This has worked, but it's very plain. We can change pretty much anything we want
# though - let's start with adding colors to the boxes
ggplot(df, aes(x = Factors, y = xylose, fill = Factors)) + 
  geom_boxplot()

# Better already! It's important to note how the ggplot function works though:
# If something in your graph is going to change, such as values on either axis or
# colors between groups, you put that inside the aes parentheses. If it's going to 
# be consistent, it goes in the geom function below. For example:
ggplot(df, aes(x = Factors, y = xylose)) + 
  geom_boxplot(fill = "darkgreen")

# Let's go back to the previous version and change some axis labels and add a title
ggplot(df, aes(x = Factors, y = xylose, fill = Factors)) + 
  geom_boxplot() +
  labs(title = "Xylose", x = "Group", y = "Intensity")

# Now, let's change the overall theme (the gray background isn't great)
ggplot(df, aes(x = Factors, y = xylose, fill = Factors)) + 
  geom_boxplot() +
  labs(title = "Xylose", x = "Group", y = "Intensity") +
  theme_classic()

# We also don't need the legend at all, and all the text could be bigger
ggplot(df, aes(x = Factors, y = xylose, fill = Factors)) + 
  geom_boxplot() +
  labs(title = "Xylose", x = "Group", y = "Intensity") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16))

# I'm happy with that - let's save it as an object and then save it as a png file
p2 <- ggplot(df, aes(x = Factors, y = xylose, fill = Factors)) + 
  geom_boxplot() +
  labs(title = "Xylose", x = "Group", y = "Intensity") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("boxplot_xylose.png", p2, width = 6, height = 5)

# This process would be very time consuming for all 115 metabolites in our data.
# We're going to create a box plot for each of them, but given this, it's neater
# to create their own folder for them to go in
dir.create(paste0(getwd(), "/Boxplots/"))

# dir.create will create a directory
# paste0 joins two strings together, which are the current working directory 
# getwd() and Boxplots, which is the new folder we want to create.

# Next we're going to want to loop through each column of the data to make a 
# boxplot for each of them. There are definitely other ways to do this as well,
# but this is just one way. Plus, you get to learn loops.

for (i in 3:ncol(df)) { # This is saying "start at column 3 and continue until the last column in df
  
  # Create a two column data frame with Factors as one column, and the current metabolite as the other
  temp_df <- df %>% 
    select(2, all_of(i)) # you don't have to worry about what "all_of" means or does
  
  # Save the current metabolite name, and rename the second column of temp_df as "Intensity" so we can 
  # use it consistently over the loop
  current_metabolite <- colnames(temp_df)[2]
  colnames(temp_df)[2] <- "Intensity"
  
  # Use temp_df as the data frame in the plot, with "Intensity" as the y axis value
  temp_p <- ggplot(temp_df, aes(x = Factors, y = Intensity, fill = Factors)) +
    geom_boxplot() + 
    labs(title = current_metabolite, x = "Group", y = "Intensity") + # Put saved "current metabolite" as title
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 16))
  # Save the plot in the boxplots folder with the current metabolite as the file name
  ggsave(paste0(getwd(), "/Boxplots/", current_metabolite, ".png"), temp_p, width = 6, height = 5)
}


# That just created boxplots for every identified metabolite in our data set. If
# we have a data set with hundreds/thousands of metabolites (e.g., LC-MS data),
# we might want to refine our approach to only including those that are statistically
# significant. We can use the exact same loop as before, but we remove columns 
# from df that are non-significant metabolites.

# First, let's find which compounds were significant in univariate_results by 
# using filter
sig_metabolites <- univariate_results %>% 
  filter(p < 0.05)

# Then, let's only keep the metabolite names and no other columns
sig_metabolites <- sig_metabolites %>% 
  pull(`Metabolite Name`)

# That leaves us with three names which we can put into a select function to keep
significant_df <- df %>% 
  select(ID, Factors, all_of(sig_metabolites))



# Let's also make a new directory called "SignificantBoxplots" just to show this 
# works:
dir.create(paste0(getwd(), "/SignificantBoxplots/"))



for (i in 3:ncol(significant_df)) { # This is saying "start at column 3 and continue until the last column in significant_df
  
  # Create a two column data frame with Factors as one column, and the current metabolite as the other
  temp_df <- significant_df %>% 
    select(2, all_of(i)) # you don't have to worry about what "all_of" means or does
  
  # Save the current metabolite name, and rename the second column of temp_df as "Intensity" so we can 
  # use it consistently over the loop
  current_metabolite <- colnames(temp_df)[2]
  colnames(temp_df)[2] <- "Intensity"
  
  # Use temp_df as the data frame in the plot, with "Intensity" as the y axis value
  temp_p <- ggplot(temp_df, aes(x = Factors, y = Intensity, fill = Factors)) +
    geom_boxplot() + 
    labs(title = current_metabolite, x = "Group", y = "Intensity") + # Put saved "current metabolite" as title
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 16))
  # Save the plot in the boxplots folder with the current metabolite as the file name
  ggsave(paste0(getwd(), "/SignificantBoxplots/", current_metabolite, ".png"), temp_p, width = 6, height = 5)
}

# Now we have boxplots saved!


# PCA plots and PLSDA plots come from the mixOmics package, which is quite
# different to ggplot2. One key difference is that when you want to save the 
# plot itself, you need to put $graph at the end of the function. The ways to 
# change details of the plots are also all within plotIndiv

pca_res <- pca(df[, 3:117], center = TRUE, scale = TRUE)
p3 <- plotIndiv(pca_res, group = df$Factors, legend = TRUE, ellipse = TRUE, pch = df$Factors)$graph
ggsave("PCA_plot.png", p3, width = 6, height = 5)



p4 <- plotIndiv(final_plsda, 
                comp = 1:2, 
                group = df$Factors, 
                pch = df$Factors, # colour by class label
                ellipse = TRUE, # include 95% confidence ellipse 
                legend = TRUE, 
                title = "PLS-DA on MECFS, comp 1 & 2")$graph
ggsave("PLSDA_plotIndiv.png", p4, width = 6, height = 5)



# The plotVar for PLSDA is entirely different again - although it produces a plot,
# the object itself is actually data that needs to be entered into ggplot2 manually!

# Advanced challenge - make this look as close to the original as possible!


p5 <- plotVar(final_plsda, 
              comp = 1:2, 
              cutoff = 0.5, 
              cex = 3) %>% 
  ggplot(., aes(x = x, y = y, label = names), col = col) +
  geom_point() +
  geom_text_repel() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  theme_minimal() +
  labs(x = "Component 1", y = "Component 2")

ggsave("PLSDA_plotVar.png", p5, width = 6, height = 5)

 
auc_plsda <- auroc(final_plsda, roc.comp = 2, print = FALSE)$graph.Comp2
ggsave("PLSDA_AUC.png", auc_plsda, width = 6, height = 5)