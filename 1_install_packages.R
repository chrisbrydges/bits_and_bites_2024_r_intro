# This short chunk of code is to install the libraries that we're going to use in 
# this course. I've made this a separate file because installing libraries is 
# something you should only have to do once (assuming that you use R again on 
# the same computer).

# Libraries can generally (but not always) be found in one of three places:
# CRAN (the main R repository for libraries), Bioconductor (a large repository
# that specializes more in bioinformatics), and github (where developers often 
# share the latest version of their libraries before they go to CRAN/Bioconductor)

# We're going to install some libraries from CRAN first:
install.packages(c("readxl", "janitor", "tidyverse", "rstatix", "ggpubr", "BiocManager"))

# With BiocManager installed, we can now reach the Bioconductor repository and
# install libraries from there.

BiocManager::install(c("EnhancedVolcano", "mixOmics"))

