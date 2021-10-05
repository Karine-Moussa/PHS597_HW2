# Lasso - Cyclic Coordinate Descent comparison
# Author: Karine Moussa
# Due: October 6th, 2021

# Since this code takes a while, I've saved the major data objects
# You may load them using these commands. 
X <- readRDS("X.rds")
Y <- readRDS("Y.rds")
X.scale <- readRDS("X_scale.rds")
Y.scale <- readRDS("Y_scale.rds")
X.scale.subset <- readRDS("X_scale_subset.rds")
Y.scale.subset <- readRDS("Y_scale_subset.rds")
coef_lasso_fit <- readRDS("coef_lasso_fit.rds")
print(coef_lasso_fit)
intercept <- readRDS("intercept.rds")
beta_hat <- readRDS("beta_hat.rds")
print(beta_hat)

# Note, I unfortunately was unable to get the output of two approaches to match. 
# But here are some areas I think could be the cause.
#
# 1. Scaling the data. I chose to scale X and Y prior to inputting
# them in glmnet, as instructed by the SLS text book. 
# However I am unsure whether or not I should've scaled them.
# 2. glmnet. I had to choose a very low lambda in order to see any 
# values for the coefficients This low lambda may have thrown things off in the 
# cyclic coordinate descent functions
# 3. Cyclic coordinate descent. I emulated the equations that are in the text 
# (I've cited page numbers/equation numbers)
# however the final beta hat vector came out as an ordered vector, 
# so I likely am performing one of the operations incorrectly.

### SET WORKING DIRECTORY
setwd("~/Documents/PennState/FALL_PSU_2021/PHS 597 Data Mining/HW/LassoCyclic")

### PACKAGES
library(glmnet)

### GET DATA
## Expression:
expr <- read.delim("GEUVADIS_normalized_expression_chr20")
## Genotypes:
gt <- read.delim("GEUVADIS_chr20_processed.traw")

### SET UP Y and X FOR ONE GENE
## GET Y
gene_num <- 1 # which column gene we're looking at
Y <- t(expr[gene_num,-c(1,2,3,4)])
colnames(Y) <- expr[gene_num,1]
## GET X (500,000 bp around gene)
pos.min <- expr[1,3] - 500000
pos.max <- expr[1,3] + 500000
gt.subset <- gt[gt$POS >= pos.min & gt$POS <= pos.max,] # subset genotypes
X <- t(gt.subset[,7:ncol(gt.subset)]) # create X matrix
rownames(X) <- sub("_.*", "", rownames(X))
colnames(X) <- rep(paste0("X", 1:ncol(X)))
saveRDS(Y, "Y.rds")
saveRDS(X, "X.rds")

# Scale the data (using information at bottom of page 14)
X.scale <- scale(X, center = FALSE, scale = TRUE) # stdev = 1
Y.scale <- scale(Y, center = FALSE, scale = FALSE)
set.seed(1) # simplify the amount of data, only use 1/3
sub_rows <- sample(1:nrow(Y.scale), .33*nrow(Y.scale))
X.scale.subset <- X.scale[sub_rows,]
Y.scale.subset <- Y.scale[sub_rows]
saveRDS(X.scale, "X_scale.rds")
saveRDS(X.scale.subset, "X_scale_subset.rds")
saveRDS(Y.scale, "Y_scale.rds")
saveRDS(Y.scale.subset, "Y_scale_subset.rds")

### FIT USING GLMNET
lasso.fit <- glmnet(X.scale.subset, Y.scale.subset, 
                    family = "gaussian",
                    alpha = 1, # make sure we're performing LASSO
                    lambda = 0.05,
                    intercept = T, # this may not be correct
                    standardize = FALSE) 
coef(lasso.fit)
saveRDS(coef(lasso.fit), "coef_lasso_fit.rds")

### FIT USING CYCLIC COORDINATE DESCENT
# Get beta_hat j using equation 2.10
get_beta_hat_j <- function(gamma, lambda){
    if ( gamma > lambda ){
        beta_hat <- gamma - lambda
    }
    if ( abs(gamma) <= lambda ){
        beta_hat <- 0
    }
    if ( gamma < abs(lambda))  {
        beta_hat <- gamma + lambda
    }
    return(beta_hat)
}

# Set up initial variables
beta_hat <- as.matrix(rep(1, ncol(X.scale.subset))) # initial beta hat vector
lambda = 0.05 # use same lambda as glmnet
N = ncol(X.scale.subset)
M = nrow(X.scale.subset)
iterations = 100 
intercept = sum(Y.scale.subset - X.scale.subset%*%beta_hat) # set up initial intercept

for( i in 1:iterations ){
    for ( j in 1:N ){
        # Get beta j 
        # Let gamma be the partial_residual/N from 2.14 in text
        partial_res <- sum(Y.scale.subset - X.scale.subset%*%beta_hat) # page 16 of text 
        gamma <- partial_res/N
        beta_hat[j] <- get_beta_hat_j(gamma, lambda)
        # Update the intercept
        intercept <- sum(Y.scale.subset - X.scale.subset%*%beta_hat)
    }  
}
saveRDS(beta_hat, "beta_hat.rds")
saveRDS(intercept, "intercept.rds")
