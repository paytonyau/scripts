# Dependencies
# R. We used 3.3.1.
#install.packages('ggplot2')
library('ggplot2')
library('AICcmodavg')

# Theme for Birdy Study, and set for current session
birdtheme <- theme_bw() + theme(axis.title = element_text(size=9),
axis.text = element_text(size=8),
strip.background = element_rect(colour="transparent", fill="transparent"),
strip.text.x = element_text(size=10, face="bold"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_rect(fill="transparent")
) +
theme(axis.line.x = element_line(color="black", size = 0.5),
axis.line.y = element_line(color="black", size = 0.5))

theme_set(birdtheme)

# Read in the data and look at it
klepdat <- read.csv('GLMM_input_data.csv')
head(klepdat)

# Recode factor variables
klepdat$SiteFactor[klepdat$Site==1] <- "Billingsgate"
klepdat$SiteFactor[klepdat$Site==2] <- "Brancaster"
klepdat$SiteFactor <- as.factor(klepdat$SiteFactor)

klepdat$SpeciesFactor[klepdat$Species==13] <- "HG" # Herring Gull
klepdat$SpeciesFactor[klepdat$Species==43] <- "CG" # Common Gull
klepdat$SpeciesFactor[klepdat$Species==213] <- "BHG" # Black Headed Gull
klepdat$SpeciesFactor[klepdat$Species==322] <- "GBB" # Great Black Backed Gull
klepdat$SpeciesFactor <- as.factor(klepdat$SpeciesFactor)

klepdat$SeasonFactor[klepdat$Season==1] <- "breeding"
klepdat$SeasonFactor[klepdat$Season==2] <- "nonbreeding"
klepdat$SeasonFactor <- as.factor(klepdat$SeasonFactor)

# Remove all observations that include NAs
NAs <- klepdat == "999"
klepdat[NAs] <- NA
writeLines("Raw Data")
print(summary(klepdat)) # confirms 15 NAs in PreySize only

cleanklepdat <- na.omit(klepdat)
writeLines("\n\nCleaned Data")
print(summary(cleanklepdat))

# Create transformed and standardized versions of continuous variables
cleanklepdat$LogKlepRate3 <- log( 60 * (cleanklepdat$KleptoNum/cleanklepdat$PatchDur))
cleanklepdat$StdLogKlepRate3 <- scale(cleanklepdat$LogKlepRate3, center=T, scale=T)
cleanklepdat$StdPreySize <- scale(cleanklepdat$PreySize, center=T, scale=T)
cleanklepdat$StdPopDens <- scale(cleanklepdat$PopDens, center=T, scale=T)
writeLines("\n\nStandardized Data")
print(summary(cleanklepdat))

# Create and Display Betas for Main Effects Model for all potential predictors
MainEffectsModel <- glm(StdLogKlepRate3 ~ SiteFactor + SpeciesFactor + SeasonFactor + StdPreySize + StdPopDens, family=gaussian(link="identity"), data=cleanklepdat)
writeLines("\n\nMain Effects Model")
print(MainEffectsModel)
writeLines("Exact AIC for Main Effects Model")
print(AIC(MainEffectsModel))
writeLines("Exact Corrected AIC for Main Effects Model")
print(AICc(MainEffectsModel))

# Stepwise backward delta-AIC approach
writeLines("\n\nStepwise Subtraction from Main Effects Model")
step(MainEffectsModel, direction="backward")

# Create and Display Betas for Minimal Model with Main Effects only
MinimalModel <- glm(StdLogKlepRate3 ~ SiteFactor + StdPreySize + StdPopDens, family=gaussian(link="identity"), data = cleanklepdat)
writeLines("\n\nMinimal Model")
print(MinimalModel)
writeLines("Exact AIC for Minimal Model")
print(AIC(MinimalModel))
writeLines("Exact Corrected AIC for Minimal Model")
print(AICc(MinimalModel))

# make predictions based on model and display
# fitted is exactly the same as predict(MinimalModel,type="response")
cleanklepdat$PredKlep <- fitted(MinimalModel)
cleanklepdat$ResKlep <- residuals(MinimalModel)
writeLines("\n\nData with Predicted Values from Minimal Model")
print(summary(cleanklepdat))

# Write the output file (optional re-read data)
write.csv(cleanklepdat, file="GLMM_output_data.csv", quote=F, row.names=F)
#cleanklepdat <- read.csv('GLMM_output_data.csv')

# plotting
writeLines("\n\nPlotting now...")

# predicted against residuals plots (deprecated)
#tiff("ResPlot%03d.tiff", width=8, height=8, res=1200, units="in")
#plot(MinimalModel, pch=16)
#dev.off()

# plot the effect of PreySize on Standardized Klepto Rate
pdf("PreySize.pdf", width=6.77, height=3.39)
ps <- ggplot(cleanklepdat, aes(x=PreySize, y=StdLogKlepRate3)) + geom_point() + geom_smooth(colour="black", method=lm, formula=y~x)
ps <- ps + facet_grid(~ SiteFactor, scales="free_x") + xlab("Prey Size (bill lengths)") + ylab("Standardized Ln Kleptoparasitism Rate")
print(ps)
dev.off()
ggsave("PreySize.tiff", width=172, height=86, unit="mm", dpi=600)
#ppi <- 600
#tiff("PreySize.tiff", width=6.77*ppi, height=3.39*ppi, res=ppi)
#print(ps)
#dev.off()

# plot the effect of StdPopDens on Standardized Klepto Rate
pdf("PopDens.pdf", width=6.77, height=3.39) #res=600, units="in"
pd <- ggplot(cleanklepdat, aes(x=PopDens, y=StdLogKlepRate3)) + geom_point() + geom_smooth(colour="black", method=lm, formula=y~x)
pd <- pd + facet_grid(~ SiteFactor, scales="free_x") + xlab(expression(paste("Population Density (birds/", km^2, ")", sep=""))) + ylab("Standardized Ln Kleptoparasitism Rate")
print(pd)
dev.off()
ggsave("PopDens.tiff", width=172, height=86, unit="mm", dpi=600)
#tiff("PopDens.tiff", width=6.77*ppi, height=3.39*ppi)
#print(pd)
#dev.off()

# Create a Model adding all pairwise interactions to Minimal Model
InteractionsModel <- glm(StdLogKlepRate3 ~ (SiteFactor + StdPreySize + StdPopDens)^2, family=gaussian(link="identity"), data = cleanklepdat)
writeLines("\n\nInteractions Model")
print(InteractionsModel)
writeLines("Exact AIC for Interactions Model")
print(AIC(InteractionsModel))
writeLines("Exact Corrected AIC for Interactions Model")
print(AICc(InteractionsModel))

# Backward delta-AIC again
writeLines("\n\nStepwise Subtraction from Interactions Model")
print(drop1(InteractionsModel, test="F"))

# Create a Model with weak interactions only and report AIC
PlusModel <- glm(formula = StdLogKlepRate3 ~ SiteFactor + StdPreySize + StdPopDens + SiteFactor:StdPopDens + StdPreySize:StdPopDens, family=gaussian(link="identity"), data = cleanklepdat)
writeLines("\n\nExact AIC for Model with one interaction removed.")
print(AIC(PlusModel))
writeLines("Exact Corrected AIC for Model with one interaction removed.")
print(AICc(PlusModel))

# Site differences for (non-standardized) prey size and population density
writeLines("\n\nSite Differences")
print(tapply(cleanklepdat$PopDens, cleanklepdat$SiteFactor, length))
writeLines("\nPopulation Density")
print(tapply(cleanklepdat$PopDens, cleanklepdat$SiteFactor, summary))
writeLines("\nPrey Size")
print(tapply(cleanklepdat$PreySize, cleanklepdat$SiteFactor, summary))

# Mann-Whitney tests for population density and prey size between sites
writeLines("\n\nHypothesis test\nPopulation Density")
print(wilcox.test(cleanklepdat$PopDens ~ cleanklepdat$SiteFactor, paired = FALSE))
writeLines("\nPrey Size")
print(wilcox.test(cleanklepdat$PreySize ~ cleanklepdat$SiteFactor, paired = FALSE))
