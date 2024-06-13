# Clear the global environment
remove(list=ls())
cat("\f")

#----------------------------------------------------------
# 1: Preparation of the session
#----------------------------------------------------------

#----------------------------------------------------------
# Install/load packages and import/read data 
#----------------------------------------------------------
#install.packages("wbstats", dependencies = TRUE) 

# install.packages("plm")
# install.packages("sandwich")
# install.packages("lmtest")

library(plm)

library(sandwich)

library(lmtest)

library(wbstats)

library(stargazer)

library(tidyverse)

library(plm)

library(plyr)

library(dplyr)

library(car)


#----------------------------------------------------------
# Define folder structure
#----------------------------------------------------------

dir <- "~/Desktop/BAM/Thesis/Panel_data/"

dirData <- paste0(dir, "Data/")
dirProg <- paste0(dir, "Programs/")
dirRslt <- paste0(dir, "Results/")

#----------------------------------------------------------
# Import/read data
#----------------------------------------------------------

dfPanelData <- read.csv(file=paste0(dirData,"Panel_data.csv"),
                            stringsAsFactors=FALSE)

# Remove Country.iso codes (not used)
dfPanelData$B2 <- NULL
dfPanelData$B3 <- NULL
dfPanelData$B4 <- NULL
dfPanelData$B5 <- NULL
dfPanelData$B6 <- NULL
dfPanelData$B7 <- NULL
dfPanelData$B8 <- NULL
dfPanelData$B11 <- NULL
dfPanelData$B12 <- NULL

#-------------------------------------------------------------------------
# Create dummies
#-------------------------------------------------------------------------

# Replace "double_cover" with "double" in the Cover_Crop column
dfPanelData$Cover_crop <- gsub("double_cover", "double", dfPanelData$Cover_crop)

# Print the unique values in the Cover_Crop column to verify the replacement
unique(dfPanelData$Cover_crop)

# Create dummy variables for Cover_crop and class using one-hot encoding
cover_crop_dummies <- model.matrix(~ Cover_crop - 1, data = dfPanelData)  # Exclude intercept
class_dummies <- model.matrix(~ class - 1, data = dfPanelData)  # Exclude intercept

# Rename the dummy variables to remove the intercept notation
colnames(cover_crop_dummies) <- gsub("Cover_crop", "", colnames(cover_crop_dummies))
colnames(class_dummies) <- gsub("class", "", colnames(class_dummies))

# Combine dummy variables with numerical variables
dfPanelData <- cbind(dfPanelData, cover_crop_dummies, class_dummies)

# First, combine Year and Month into a single variable
dfPanelData$date <- as.Date(paste(dfPanelData$Year, dfPanelData$Month, "01", sep = "-"))

# Then, convert the combined variable to a factor to get unique indices
dfPanelData$date <- as.factor(dfPanelData$date)

# Now, convert the factor levels to integer indices
dfPanelData$date <- as.integer(dfPanelData$date)

# Print the head of the dataframe to check the new variable
head(dfPanelData$date)

#----------------------------------------------------------
# Data engineering to account for seasonal changes
#----------------------------------------------------------

dfPanelData <- dfPanelData %>%
  group_by(id) %>%
  dplyr::mutate(monthsStack = row_number()) %>%
  ungroup()

#----------------------------------------------------------
# Standardize remote sensing variables
#----------------------------------------------------------

# Calculate mean and standard deviation for each variable
mean_BSI <- mean(dfPanelData$BSI)
sd_BSI <- sd(dfPanelData$BSI)

mean_NDVI <- mean(dfPanelData$NDVI)
sd_NDVI <- sd(dfPanelData$NDVI)

mean_NDMI <- mean(dfPanelData$NDMI)
sd_NDMI <- sd(dfPanelData$NDMI)

mean_NDTI <- mean(dfPanelData$NDTI)
sd_NDTI <- sd(dfPanelData$NDTI)

mean_SI <- mean(dfPanelData$SI)
sd_SI <- sd(dfPanelData$SI)

mean_SSI <- mean(dfPanelData$SSI)
sd_SSI <- sd(dfPanelData$SSI)

mean_EVI <- mean(dfPanelData$EVI)
sd_EVI <- sd(dfPanelData$EVI)

# Standardize the variables
dfPanelData$BSI <- (dfPanelData$BSI - mean_BSI) / sd_BSI
dfPanelData$NDVI <- (dfPanelData$NDVI - mean_NDVI) / sd_NDVI
dfPanelData$NDMI <- (dfPanelData$NDMI - mean_NDMI) / sd_NDMI
dfPanelData$NDTI <- (dfPanelData$NDTI - mean_NDTI) / sd_NDTI
dfPanelData$SI <- (dfPanelData$SI - mean_SI) / sd_SI
dfPanelData$SSI <- (dfPanelData$SSI - mean_SSI) / sd_SSI
dfPanelData$EVI <- (dfPanelData$EVI - mean_EVI) / sd_EVI

# Calculate barren months (1 if BSI >= 0.5, 0 otherwise)
dfPanelData <- dfPanelData %>%
  mutate(barren_month = ifelse(BSI >= 0.5, 1, 0))

str(dfPanelData[c("longitude", "latitude")])
# Convert longitude and latitude to numeric
dfPanelData$longitude <- as.numeric(dfPanelData$longitude)
dfPanelData$latitude <- as.numeric(dfPanelData$latitude)

# Check the data types after conversion
str(dfPanelData[c("longitude", "latitude")])

# Find rows with missing values in longitude
longitude_missing_rows <- which(is.na(dfPanelData$longitude))

# Find rows with missing values in latitude
latitude_missing_rows <- which(is.na(dfPanelData$latitude))

# Print the row numbers where missing values occur
print("Rows with missing values in longitude:")
print(longitude_missing_rows)

print("Rows with missing values in latitude:")
print(latitude_missing_rows)

#----------------------------------------------------------
# Summary stats
#----------------------------------------------------------
model <- data.frame(dfPanelData$NDVI, dfPanelData$Cover_crop, dfPanelData$class, dfPanelData$BSI, dfPanelData$NDMI, dfPanelData$NDTI, dfPanelData$SI, dfPanelData$SSI)
stargazer(model, type = "text")

# Filter out only the numeric columns
numeric_columns <- dfPanelData[, sapply(dfPanelData, is.numeric)]

# Calculate the correlation matrix
cor_matrix <- cor(numeric_columns)
# View the correlation matrix
cor_matrix

#----------------------------------------------------------
# Regression models (NDVI)
#----------------------------------------------------------
mdlA <- NDVI ~  BSI + NDMI + NDTI + SI + SSI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton

# Pooled Model

pooled_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet +
                      monthsStack + constant_contrast + longitude + latitude,
                    data = dfPanelData,
                    index = c("id", "date"),
                    model = "pooling")
summary(pooled_model)

# Fixed Effects Model
fixed_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                     monthsStack + constant_contrast + longitude + latitude,
                                   data = dfPanelData,
                                   index = c("id", "date"),
                                   model = "within")
summary(fixed_model)

# Random Effects Model
random_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                      monthsStack + constant_contrast + longitude + latitude,
                                   data = dfPanelData,
                                   index = c("id", "date"),
                                   model = "random")
summary(random_model)

# Between Effects Model
between_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                       monthsStack + constant_contrast + longitude + latitude,
                     data = dfPanelData,
                     index = c("id", "date"),
                     model = "between")
summary(between_model)

stargazer(pooled_model, between_model, fixed_model, random_model,
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

#----------------------------------------------------------
# Regression models (NDVI) Robust standard errors
#----------------------------------------------------------

# Pooled Model with Clustered Standard Errors
pooled_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                      monthsStack + constant_contrast + longitude + latitude,
                    data = dfPanelData,
                    index = c("id", "date"),
                    model = "pooling")

# Clustered standard errors by plot (id)
pooled_se <- vcovHC(pooled_model, type = "HC1", cluster = "group")
pooled_test <- coeftest(pooled_model, vcov = pooled_se)

# Fixed Effects Model with Clustered Standard Errors
fixed_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                     monthsStack + constant_contrast + longitude + latitude,
                   data = dfPanelData,
                   index = c("id", "date"),
                   model = "within")

# Clustered standard errors by plot (id)
fixed_se <- vcovHC(fixed_model, type = "HC1", cluster = "group")
fixed_test <- coeftest(fixed_model, vcov = fixed_se)

# Random Effects Model with Clustered Standard Errors
random_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                      monthsStack + constant_contrast + longitude + latitude,
                    data = dfPanelData,
                    index = c("id", "date"),
                    model = "random")

# Clustered standard errors by plot (id)
random_se <- vcovHC(random_model, type = "HC1", cluster = "group")
random_test <- coeftest(random_model, vcov = random_se)

# Between Effects Model with Clustered Standard Errors
#between_model <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
#                       monthsStack + constant_contrast + longitude + latitude,
#                     data = dfPanelData,
#                     index = c("id", "date"),
#                     model = "between")
#
# Clustered standard errors by plot (id)
#between_se <- vcovHC(between_model, type = "HC1", cluster = "group")
# between_test <- coeftest(between_model, vcov = between_se)

stargazer(pooled_model, fixed_model, random_model,
          se = list(sqrt(diag(pooled_se)),
                    sqrt(diag(fixed_se)),
                    sqrt(diag(random_se))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")


#-------------------------------------------------------------------------
# Fixed effect regression model 
#-------------------------------------------------------------------------

# Fixed Effects Model
fixed_model_simple <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton,
                   data = dfPanelData,
                   index = c("id", "date"),
                   model = "within")

fixed_se1 <- vcovHC(fixed_model_simple, type = "HC1", cluster = "group")
fixed_test1 <- coeftest(fixed_model_simple, vcov = fixed_se)


fixed_model_control <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                        monthsStack + constant_contrast + longitude + latitude,
                   data = dfPanelData,
                   index = c("id", "date"),
                   model = "within")

fixed_se2 <- vcovHC(fixed_model_control, type = "HC1", cluster = "group")
fixed_test2 <- coeftest(fixed_model_control, vcov = fixed_se2)

fixed_model_interaction <- plm(NDVI ~ NDMI + SI + SSI + NDTI + 
                                 NDTI * cover + NDTI * double + NDTI * soy_corn + NDTI * soy_millet + NDTI * soy_cotton +
                                 barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                                 monthsStack + constant_contrast + longitude + latitude,
                               data = dfPanelData,
                               index = c("id", "date"),
                               model = "within")

fixed_se3 <- vcovHC(fixed_model_interaction, type = "HC1", cluster = "group")
fixed_test3 <- coeftest(fixed_model_interaction, vcov = fixed_se3)

# without robust standard errors
stargazer(fixed_model_simple, fixed_model_control, fixed_model_interaction,
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

# with robust standard errors
stargazer(fixed_model_simple, fixed_model_control, fixed_model_interaction,
          se = list(sqrt(diag(fixed_se1)),
                    sqrt(diag(fixed_se2)),
                    sqrt(diag(fixed_se3))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")

#-------------------------------------------------------------------------
# Random effect regression model 
#-------------------------------------------------------------------------

# Random Effects Model
random_model_simple <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton,
                          data = dfPanelData,
                          index = c("id", "date"),
                          model = "random")

random_se1 <- vcovHC(random_model_simple, type = "HC1", cluster = "group")
random_test1 <- coeftest(random_model_simple, vcov = random_se1)

random_model_control <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                              monthsStack + constant_contrast + longitude + latitude,
                           data = dfPanelData,
                           index = c("id", "date"),
                           model = "random")

random_se2 <- vcovHC(random_model_control, type = "HC1", cluster = "group")
random_test2 <- coeftest(random_model_control, vcov = random_se2)

random_model_interaction <- plm(NDVI ~ NDMI + SI + SSI + NDTI + 
                                 NDTI * cover + NDTI * double + NDTI * soy_corn + NDTI * soy_millet + NDTI * soy_cotton +
                                 barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                                  monthsStack + constant_contrast + longitude + latitude,
                               data = dfPanelData,
                               index = c("id", "date"),
                               model = "random")

random_se3 <- vcovHC(random_model_interaction, type = "HC1", cluster = "group")
random_test3 <- coeftest(random_model_interaction, vcov = random_se3)

# Without robust standard errors
stargazer(random_model_simple, random_model_control, random_model_interaction,
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

# With robust standard errors
stargazer(random_model_simple, random_model_control, random_model_interaction,
          se = list(sqrt(diag(random_se1)),
                    sqrt(diag(random_se2)),
                    sqrt(diag(random_se3))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")

#-------------------------------------------------------------------------
# Between effect regression model 
#-------------------------------------------------------------------------

# Random Effects Model
between_model_simple <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton,
                           data = dfPanelData,
                           index = c("id", "date"),
                           model = "between")
summary(between_model_simple)

between_model_control <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                               monthsStack + constant_contrast + longitude + latitude,
                            data = dfPanelData,
                            index = c("id", "date"),
                            model = "between")
summary(between_model_control)

between_model_interaction <- plm(NDVI ~ NDMI + SI + SSI + NDTI + 
                                  NDTI * cover + NDTI * double + NDTI * soy_corn + NDTI * soy_millet + NDTI * soy_cotton +
                                  barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                                   monthsStack + constant_contrast + longitude + latitude,
                                data = dfPanelData,
                                index = c("id", "date"),
                                model = "between")
summary(between_model_interaction)

stargazer(between_model_simple, between_model_control, between_model_interaction,
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

#-------------------------------------------------------------------------
# Pooled effect regression model 
#-------------------------------------------------------------------------

# Pooled Model

pooled_model_simple <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton,
                           data = dfPanelData,
                           index = c("id", "date"),
                           model = "pooling")

pooled_se1 <- vcovHC(pooled_model_simple, type = "HC1", cluster = "group")
pooled_test1 <- coeftest(pooled_model_simple, vcov = pooled_se1)

pooled_model_control <- plm(NDVI ~ NDMI + SI + SSI + NDTI + barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                              monthsStack + constant_contrast + longitude + latitude,
                            data = dfPanelData,
                            index = c("id", "date"),
                            model = "pooling")

pooled_se2 <- vcovHC(pooled_model_control, type = "HC1", cluster = "group")
pooled_test2 <- coeftest(pooled_model_control, vcov = pooled_se2)

pooled_model_interaction <- plm(NDVI ~ NDMI + SI + SSI + NDTI + 
                                  NDTI * cover + NDTI * double + NDTI * soy_corn + NDTI * soy_millet + NDTI * soy_cotton +
                                  barren_month + cover + double + soy_corn + soy_millet + soy_cotton +
                                  monthsStack + constant_contrast + longitude + latitude,
                                data = dfPanelData,
                                index = c("id", "date"),
                                model = "pooling")

pooled_se3 <- vcovHC(pooled_model_interaction, type = "HC1", cluster = "group")
pooled_test3 <- coeftest(pooled_model_interaction, vcov = pooled_se3)

# Without robust standard errors
stargazer(pooled_model_simple, pooled_model_control, pooled_model_interaction,
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

# With robust standard errors
stargazer(pooled_model_simple, pooled_model_control, pooled_model_interaction,
          se = list(sqrt(diag(pooled_se1)),
                    sqrt(diag(pooled_se2)),
                    sqrt(diag(pooled_se3))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")


# Evaluate the fixed effects model versus the pooled 
# regression model
pFtest(fixed_model_simple, pooled_model_simple)
pooltest(pooled_model_simple, fixed_model_simple)

pFtest(fixed_model_control, pooled_model_control)
pooltest(pooled_model_control, fixed_model_control)

pFtest(fixed_model_interaction, pooled_model_interaction)
pooltest(pooled_model_interaction, fixed_model_interaction)

## insignificant test: both models are consistent. 
# Rejection in favour of the fixed effect model so the fixed effect model 
# should be preffered over the pooled effects model

#----------------------------------------------------------
# Regression models (NDMI)
#----------------------------------------------------------

# Pooled Model

pooled_model_NDMI <- plm(NDMI ~ SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                           monthsStack + constant_contrast + longitude + latitude,
                         data = dfPanelData,
                         index = c("id", "date"),
                         model = "pooling")
pooled_se <- vcovHC(pooled_model_NDMI, type = "HC1", cluster = "group")
pooled_test <- coeftest(pooled_model_NDMI, vcov = pooled_se)

# Fixed Effects Model
fixed_model_NDMI <- plm(NDMI ~ SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                          monthsStack + constant_contrast + longitude + latitude,
                        data = dfPanelData,
                        index = c("id", "date"),
                        model = "within")
fixed_se <- vcovHC(fixed_model_NDMI, type = "HC1", cluster = "group")
fixed_test <- coeftest(fixed_model_NDMI, vcov = fixed_se)

# Random Effects Model
random_model_NDMI <- plm(NDMI ~ SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                           monthsStack + constant_contrast + longitude + latitude,
                         data = dfPanelData,
                         index = c("id", "date"),
                         model = "random")
random_se <- vcovHC(random_model_NDMI, type = "HC1", cluster = "group")
random_test <- coeftest(random_model_NDMI, vcov = random_se)

# Between Effects Model
#between_model_NDMI <- plm(NDMI ~ SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
#                            monthsStack + constant_contrast + longitude + latitude,
#                          data = dfPanelData,
#                          index = c("id", "date"),
#                          model = "between")
#between_se <- vcovHC(between_model_NDMI, type = "HC1", cluster = "group")
#between_test <- coeftest(between_model_NDMI, vcov = between_se)


# With robust standard errors
stargazer(pooled_model_NDMI, random_model_NDMI, fixed_model_NDMI,
          se = list(sqrt(diag(pooled_se)),
                    sqrt(diag(random_se)),
                    sqrt(diag(fixed_se))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")

#----------------------------------------------------------
# Regression models (BSI)
#----------------------------------------------------------

# Pooled Model

pooled_model_BSI <- plm(BSI ~ SI + SSI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                          monthsStack + constant_contrast + longitude + latitude,
                         data = dfPanelData,
                         index = c("id", "date"),
                         model = "pooling")
pooled_se <- vcovHC(pooled_model_BSI, type = "HC1", cluster = "group")
pooled_test <- coeftest(pooled_model_BSI, vcov = pooled_se)

# Fixed Effects Model
fixed_model_BSI <- plm(BSI ~ SI + SSI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                         monthsStack + constant_contrast + longitude + latitude,
                        data = dfPanelData,
                        index = c("id", "date"),
                        model = "within")
fixed_se <- vcovHC(fixed_model_BSI, type = "HC1", cluster = "group")
fixed_test <- coeftest(fixed_model_BSI, vcov = fixed_se)

# Random Effects Model
random_model_BSI <- plm(BSI ~ SI + SSI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                          monthsStack + constant_contrast + longitude + latitude,
                         data = dfPanelData,
                         index = c("id", "date"),
                         model = "random")
random_se <- vcovHC(random_model_BSI, type = "HC1", cluster = "group")
random_test <- coeftest(random_model_BSI, vcov = random_se)

# Between Effects Model
between_model_BSI <- plm(BSI ~ SI + SSI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                           monthsStack + constant_contrast + longitude + latitude,
                          data = dfPanelData,
                          index = c("id", "date"),
                          model = "between")
summary(between_model_BSI)

stargazer(pooled_model_BSI, between_model_BSI, random_model_BSI, fixed_model_BSI, 
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

# With robust standard errors
stargazer(pooled_model_BSI, random_model_BSI, fixed_model_BSI,
          se = list(sqrt(diag(pooled_se)),
                    sqrt(diag(random_se)),
                    sqrt(diag(fixed_se))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")


#----------------------------------------------------------
# Regression models (SSI)
#----------------------------------------------------------

# Pooled Model

pooled_model_SSI <- plm(SSI ~ BSI + SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                          monthsStack + constant_contrast + longitude + latitude,
                        data = dfPanelData,
                        index = c("id", "date"),
                        model = "pooling")
pooled_se <- vcovHC(pooled_model_SSI, type = "HC1", cluster = "group")
pooled_test <- coeftest(pooled_model_SSI, vcov = pooled_se)

# Fixed Effects Model
fixed_model_SSI <- plm(SSI ~ BSI + SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                         monthsStack + constant_contrast + longitude + latitude,
                       data = dfPanelData,
                       index = c("id", "date"),
                       model = "within")
fixed_se <- vcovHC(fixed_model_SSI, type = "HC1", cluster = "group")
fixed_test <- coeftest(fixed_model_SSI, vcov = fixed_se)

# Random Effects Model
random_model_SSI <- plm(SSI ~ BSI + SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                          monthsStack + constant_contrast + longitude + latitude,
                        data = dfPanelData,
                        index = c("id", "date"),
                        model = "random")
random_se <- vcovHC(random_model_SSI, type = "HC1", cluster = "group")
random_test <- coeftest(random_model_SSI, vcov = random_se)

# Between Effects Model
between_model_SSI <- plm(SSI ~ BSI + SI + NDTI + cover + double + soy_corn + soy_millet + soy_cotton +
                           monthsStack + constant_contrast + longitude + latitude,
                         data = dfPanelData,
                         index = c("id", "date"),
                         model = "between")
summary(between_model_SSI)

stargazer(pooled_model_SSI, between_model_SSI, random_model_SSI, fixed_model_SSI, 
          align=TRUE, no.space=TRUE, intercept.bottom = FALSE, type="text")

# With robust standard errors
stargazer(pooled_model_SSI, random_model_SSI, fixed_model_SSI,
          se = list(sqrt(diag(pooled_se)),
                    sqrt(diag(random_se)),
                    sqrt(diag(fixed_se))),
          align = TRUE, no.space = TRUE, intercept.bottom = FALSE, type = "text")
