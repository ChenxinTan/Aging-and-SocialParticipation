######################################################################################################
##  Status in Quo of Social Participation and Associated Sociodemographic Characteristics in China  ##
##                    Author: Chenxin Tan, Applied Quantitative Research at NYU                     ##
##                                      (chenxin.tan@nyu.edu)                                       ##
##        Data: China Health and Retirement Longitudinal Study (CHARLS) 2015 Cross-sectional        ##
##                             (http://charls.pku.edu.cn/index/en.html)                             ##
##                                   Last Modified: July 06, 2020                                   ##
######################################################################################################

# library packages
library(here)
library(haven)
library(dplyr)
library(poLCA)
library(reshape2)
library(ggplot2)
library(psych)
library(DescTools)
library(writexl)
library(gmodels)
library(ggpubr)
library(fastDummies)

## set working directory
setwd("~/Desktop/New York University/AQA I/finalpaper")

## load CHARLS 2015 health status and functioning data
healthsts <- read_dta("CHARLS2015/Health_Status_and_Functioning.dta")

# social participation items
healthsts <- healthsts %>%
  mutate(friends = as_factor(da057_1_), playcards = as_factor(da057_2_), helpothers = as_factor(da057_3_),
         gotoclub = as_factor(da057_4_), organization = as_factor(da057_5_), volunteer = as_factor(da057_6_),
         careothers = as.factor(da057_7_), takecourse = as.factor(da057_8_), otheract = as_factor(da057_11_))

# recode factor levels (R function)
lvls <- c("1 Almost Daily", "2 Almost Every Week", "3 Not Regularly", "4 Never")

levels(healthsts$friends)[3] <- "3 Not Regularly"
healthsts$friends <- factor(healthsts$friends, levels = lvls)
healthsts$friends[is.na(healthsts$friends)] <- "4 Never"

levels(healthsts$playcards)[3] <- "3 Not Regularly"
healthsts$playcards <- factor(healthsts$playcards, levels = lvls)
healthsts$playcards[is.na(healthsts$playcards)] <- "4 Never"

levels(healthsts$helpothers)[3] <- "3 Not Regularly"
healthsts$helpothers <- factor(healthsts$helpothers, levels = lvls)
healthsts$helpothers[is.na(healthsts$helpothers)] <- "4 Never"

levels(healthsts$gotoclub)[3] <- "3 Not Regularly"
healthsts$gotoclub <- factor(healthsts$gotoclub, levels = lvls)
healthsts$gotoclub[is.na(healthsts$gotoclub)] <- "4 Never"

levels(healthsts$organization)[3] <- "3 Not Regularly"
healthsts$organization <- factor(healthsts$organization, levels = lvls)
healthsts$organization[is.na(healthsts$organization)] <- "4 Never"

levels(healthsts$volunteer)[3] <- "3 Not Regularly"
healthsts$volunteer <- factor(healthsts$volunteer, levels = lvls)
healthsts$volunteer[is.na(healthsts$volunteer)] <- "4 Never"

levels(healthsts$careothers)[3] <- "3 Not Regularly"
healthsts$careothers <- factor(healthsts$careothers, levels = lvls)
healthsts$careothers[is.na(healthsts$careothers)] <- "4 Never"

levels(healthsts$takecourse)[3] <- "3 Not Regularly"
healthsts$takecourse <- factor(healthsts$takecourse, levels = lvls)
healthsts$takecourse[is.na(healthsts$takecourse)] <- "4 Never"

levels(healthsts$otheract)[3] <- "3 Not Regularly"
healthsts$otheract <- factor(healthsts$otheract, levels = lvls)
healthsts$otheract[is.na(healthsts$otheract)] <- "4 Never"

### latent class analysis (poLCA package) ###
set.seed(2048)

# store the goodness of fitness results as a datafram
lca.fitness <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(lca.fitness) <- c("N Class","G2", "BIC","Entropy")

# fit the LCA model looping through 1- to 6-class
spitems <- with(healthsts, cbind(friends, playcards, helpothers, gotoclub, organization, volunteer, careothers, takecourse, otheract))

for (i in 1:6) {
  lca.charls2015.ncls <- poLCA(spitems ~ 1, nclass = i, data = healthsts, nrep = 5, na.rm = TRUE, graphs = FALSE, maxiter = 100000)
  lca.fitness[i, "N Class"] = i
  lca.fitness[i, "G2"] = lca.charls2015.ncls$Gsq
  lca.fitness[i, "BIC"] = lca.charls2015.ncls$bic
  lca.fitness[i, "Entropy"] = poLCA.entropy(lca.charls2015.ncls)
}

lca.fitness #The 4-class model have the smallest BIC among all models

# predicted class from 4-class model
lca.charls2015.ncls4 <- poLCA(spitems ~ 1, nclass = 4, data = healthsts, nrep = 5, na.rm = TRUE, graphs = FALSE, maxiter = 100000)

plot(lca.charls2015.ncls4)

# reorder the latent classes
probs.start.new <- poLCA.reorder(lca.charls2015.ncls4$probs.start, order(c(4,1,3,2)))

lca.reorder <- poLCA(spitems ~ 1, nclass = 4, data = healthsts, nrep = 5, na.rm = TRUE, graphs = FALSE, maxiter = 100000, probs.start = probs.start.new)

plot(lca.reorder)

healthsts$spclass <- lca.reorder$predclass

# name the yielded latent classes
healthsts$spclass <- factor(healthsts$spclass, levels = c(1, 2, 3, 4), labels = c("1 society-oriented", "2 recreational", "3 occasional", "4 non-participated"))
healthsts$spclass <- relevel(healthsts$spclass, ref = 4)
table(healthsts$spclass, useNA = "ifany")

# plot the LCA results
lca.probs <- melt(lca.reorder$probs, level = 2)

lca.probs$Var1 <- recode(lca.probs$Var1, "class 1: " = "Society-oriented", "class 2: " = "Recreational", "class 3: " = "Occasional", "class 4: " = "Non-particiapted")
lca.probs$Var2 <- recode(lca.probs$Var2, "Pr(1)" = "Almost Daily", "Pr(2)" = "Almost Weekly", "Pr(3)" = "Not Regularly", "Pr(4)" = "Never")

zp1 <- ggplot(lca.probs, aes(x = L2, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Var1) +
  scale_fill_brewer(type = "seq", palette="Greys", direction = -1) +
  labs(x = "Social Particiaption Item", y = "Estimated Probablity", fill = "Item Response", caption = "Source: CHARLS 2015") +
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
        plot.caption = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 16),
        strip.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, angle = 45, vjust = 0.6),
        axis.text.y = element_text(size = 14))

print(zp1)

ggsave("table_and_figure/lcaplot1.png", plot = zp1, width = 13, height = 8) 

zp2 <- ggplot(lca.probs, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ L2) +
  scale_x_discrete("Latent Class", expand = c(0, 0)) +
  scale_y_continuous("Response Proportion", expand = c(0, 0)) +
  scale_fill_brewer(type = "seq", palette = "Greys", direction = -1, "Item Response") +
  labs(caption = "Source: CHARLS 2015") +
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = 1),
        plot.margin = margin(t = 15, r = 25, b = 15, l = 25, unit = "pt"),
        plot.caption = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 16),
        strip.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, angle = 15, vjust = 0.6),
        axis.text.y = element_text(size = 14))

print(zp2)

ggsave("table_and_figure/lcaplot2.png", plot = zp2, width = 13, height = 8)

# dummy variable of participation
healthsts$spitem.response.total <- rowSums(spitems)

healthsts$participation <- 0
for (col in colnames(spitems)){
  healthsts$participation <- replace(healthsts$participation, healthsts[[col]] == "1 Almost Daily" | healthsts[[col]] == "2 Almost Every Week", 1)
}
healthsts$participation <- factor(healthsts$participation, levels = c(0,1), labels = c("0 No", "1 Yes"), ordered = TRUE)

table(healthsts$participation)

## load CHARLS 2015 demography information data
demography <- read_dta("CHARLS2015/Demographic_Background.dta")

# year of birth
demography$birthyr <- ifelse(demography$ba002 == 2, demography$ba002_1, demography$ba004_w3_1)
demography$birthyr <- ifelse(is.na(demography$birthyr) == TRUE & is.na(demography$ba004_w3_1) == FALSE, demography$ba004_w3_1, demography$birthyr)

# age (2015 - year of age)
demography$age = 2015 - demography$birthyr

# hukou status (non-agricultural/urban hukou == 1) and highest educational level
# <! NOTE !> Since the CHARLS 2015 didn't provide the preload variables and only the transition of Hukou status has been asked (released) in the 2015 dataset,
#            we need to retrive information from (and merge) CHARLS 2013 where the preload variables were included

demo2015 <- demography %>%
  dplyr::select(ID, bc001_w3_1, bc002_w3_1, bc001_w3_2, bd001_w2_4)

demo2013 <- read_dta("CHARLS2013/Demographic_Background.dta")

demo2013 <- demo2013 %>%
  mutate(hukou2011 = zbc001, edu2011 = zbd001) %>%
  dplyr::select(ID, householdID, communityID, hukou2011, bc001, edu2011, bd001_w2_1, bd001_w2_2, bd001_w2_3, bd001_w2_4, bd001) %>%
  mutate(hukou2013 = ifelse(is.na(bc001), hukou2011, bc001)) %>%
  mutate(hukou2013 = ifelse(hukou2013 == 1, 1, ifelse(hukou2013 == 2 | hukou2013 == 3, 0, NA)))

demo2013$hukou2013 <- factor(demo2013$hukou2013, levels = c(1,0), labels = c("Agricultural", "Non-agricultural"))

demo2013$edu2013 <- ifelse(demo2013$bd001_w2_1 == 1, demo2013$edu2011, demo2013$bd001_w2_3)
demo2013$edu2013 <- ifelse(is.na(demo2013$edu2013 & demo2013$bd001_w2_4 != 12), demo2013$bd001_w2_4, demo2013$edu2013)
demo2013$edu2013 <- ifelse(is.na(demo2013$edu2013), demo2013$bd001, demo2013$edu2013)

demo2013 <- dplyr::select(demo2013, ID, hukou2013, edu2013)

demo2015_merge <- merge(demo2015, demo2013, by = "ID", all.x = TRUE)

demo2015_merge <- demo2015_merge %>%
  mutate(hukou2015 = ifelse(!is.na(bc002_w3_1), bc002_w3_1, hukou2013),
         edu2015 = ifelse(!is.na(bd001_w2_4) & bd001_w2_4 != 12, bd001_w2_4, edu2013)) %>%
  mutate(hukou2015 = ifelse(is.na(hukou2015) & !is.na(bc001_w3_2), bc001_w3_2, hukou2015),
         edu2015 = ifelse(is.na(edu2015) & !is.na(edu2013), edu2013, edu2015)) %>%
  mutate(hukou2015 = ifelse(hukou2015 == 1, 1, ifelse(hukou2015 == 2 | hukou2015 == 3, 0, NA)))

demo2015_merge$hukou2015 <- factor(demo2015_merge$hukou2015, levels = c(1,0), labels = c("Agricultural", "Non-agricultural"))

demo2015_merge <- demo2015_merge %>%
  mutate(edulv = case_when(edu2015 == 1 ~ 1,
                           edu2015 >= 2 & edu2015 <= 4 ~ 2,
                           edu2015 == 5 ~ 3,
                           edu2015 >= 6 & edu2015 <= 11 ~ 4)) %>%
  mutate(edulv = factor(edulv, levels = c(1,2,3,4), labels = c("1 No Education", "2 Elementary", "3 Middle", "4 HS or Higher"))) %>%
  mutate(edulv = relevel(edulv, ref = "1 No Education"))

mean(is.na(demo2015_merge$hukou2015)) * 100 #missing value inspection
mean(is.na(demo2015_merge$edulv)) * 100 #missing value inspection

table(demo2015_merge$edulv)

# marital status
demography$married <- demography$be001

demography$married[demography$married >= 1 & demography$married <= 3] <- "1 Married"
demography$married[demography$married >= 4 & demography$married <= 7] <- "0 Unmarried (inc. Widowed)"
demography$married <- factor(demography$married)

## Health and Function Wrangling

# gender (male == 1)
healthsts$male <- factor(ifelse(healthsts$xrgender == 1, 1, 0), levels = c(1, 0), labels = c("1 male", "0 female"))
healthsts$male <- relevel(healthsts$male, ref = "0 female")

# self-rated -5c- health status
healthsts$srhealth <- ifelse(!is.na(healthsts$da001), healthsts$da001, ifelse(!is.na(healthsts$da002), healthsts$da002, NA))
healthsts$srhealth <- factor(healthsts$srhealth, levels = c(1, 2, 3, 4, 5), labels = c("1 Excellent", "2 Very Good", "3 Good", "4 Fair", "5 Poor"), ordered = FALSE)
healthsts$srhealth <- relevel(healthsts$srhealth, ref = "5 Poor")

# physical functioning hardship
phyhard.alpha <- psych::alpha(dplyr::select(healthsts, db001, db004, db005, db006, db007, db008, db009), na.rm = TRUE) #Cronbach's alpha = 0.78

healthsts$phyfunction <- phyhard.alpha$scores

# ADL (activities of daily living)
healthsts <- healthsts %>%
  mutate(dressing   = case_when(db010 %in% c(2,3,4) ~ 1, db010 == 1 ~ 0),
         bathing    = case_when(db011 %in% c(2,3,4) ~ 1, db011 == 1 ~ 0),
         feeding    = case_when(db012 %in% c(2,3,4) ~ 1, db012 == 1 ~ 0),
         transfer   = case_when(db013 %in% c(2,3,4) ~ 1, db013 == 1 ~ 0),
         toileting  = case_when(db014 %in% c(2,3,4) ~ 1, db014 == 1 ~ 0),
         continence = case_when(db015 %in% c(2,3,4) ~ 1, db015 == 1 ~ 0)) %>%
  mutate(ADLscore = rowSums(cbind(dressing, bathing, feeding, transfer, toileting, continence), na.rm = TRUE),
         miss.ADL = rowSums(is.na(cbind(dressing, bathing, feeding, transfer, toileting, continence)))) %>%
  mutate(ADLscore = ifelse(miss.ADL >= 4, NA, ADLscore))

# IADL (instrumental activities of daily living)
healthsts <- healthsts %>%
  mutate(housekeeping = case_when(db016 %in% c(2,3,4) ~ 1, db016 == 1 ~ 0),
         cooking      = case_when(db017 %in% c(2,3,4) ~ 1, db017 == 1 ~ 0),
         shopping     = case_when(db018 %in% c(2,3,4) ~ 1, db018 == 1 ~ 0),
         telephone    = case_when(db035 %in% c(2,3,4) ~ 1, db035 == 1 ~ 0),
         medication   = case_when(db020 %in% c(2,3,4) ~ 1, db020 == 1 ~ 0),
         financemng   = case_when(db019 %in% c(2,3,4) ~ 1, db019 == 1 ~ 0)) %>%
  mutate(IADLscore = rowSums(cbind(housekeeping, cooking, shopping, telephone, medication, financemng), na.rm = TRUE),
         miss.IADL = rowSums(is.na(cbind(housekeeping, cooking, shopping, telephone, medication, financemng)))) %>%
  mutate(IADLscore = ifelse(miss.IADL >= 4, NA, IADLscore))

# MMSE (mini-mental state examination in CHARLS)

#1 time orientation
healthsts <- healthsts %>%
  mutate(crtyear   = ifelse(dc001s1 == 1, 1, 0),
         crtmonth  = ifelse(dc001s2 == 2, 1, 0),
         crtday    = ifelse(dc001s3 == 1, 1, 0),
         crtwkday  = case_when(dc002 == 1 ~ 1, dc002 == 2 ~ 0),
         crtseason = case_when(dc003 == 1 ~ 1, dc002 == 2 ~ 0))

#immediate word recall test
for (i in 1:10) {
  name <- paste("dc006s", i, sep = "")
  newvar <- paste("wordtest", i, sep = "")
  healthsts[[newvar]] <- ifelse(is.na(healthsts[[name]]), 1, 0)
}

healthsts <- healthsts %>%
  mutate(wordscore = healthsts %>%
           dplyr::select(starts_with("wordtest")) %>% rowSums())

healthsts$wordscore[healthsts$dc006s11 == 11] <- 0
healthsts$wordscore[healthsts$wordscore == 0 & is.na(healthsts$dc006s11) == TRUE] <- NA

#calculation test
answer <- 100 - 7
for (i in 19:23) {
  name <- paste("dc0", i, sep = "")
  newvar <- paste("calcutest", i - 18, sep = "")
  print(answer)
  healthsts[[newvar]] <- ifelse(healthsts[[name]] == answer, 1, ifelse(!is.na(healthsts[[name]]), 0, NA))
  answer <- answer - 7
}

healthsts <- healthsts %>%
  mutate(calcuscore = healthsts %>%
           dplyr::select(starts_with("calcutest")) %>% rowSums())

#drawing test
healthsts$drawing <- ifelse(healthsts$dc025 == 1, 1, ifelse(healthsts$dc025 == 2, 0, NA))

#delayed word recall test
for (i in 1:10) {
  name <- paste("dc027s", i, sep = "")
  newvar <- paste("delaywordtest", i, sep = "")
  healthsts[[newvar]] <- ifelse(is.na(healthsts[[name]]), 1, 0)
}

healthsts <- healthsts %>%
  mutate(delaywordscore = healthsts %>%
           dplyr::select(starts_with("delaywordtest")) %>% rowSums())

healthsts$delaywordscore[healthsts$dc027s11 == 11] <- 0
healthsts$delaywordscore[healthsts$delaywordscore == 0 & is.na(healthsts$dc027s11) == TRUE] <- NA

healthsts <- healthsts %>%
  mutate(MMSEscore = rowSums(cbind(crtyear, crtmonth, crtday, crtwkday, crtseason, wordscore, calcuscore, drawing, delaywordscore), na.rm = TRUE),
         miss.MMSE = rowSums(is.na(cbind(crtyear, crtmonth, crtday, crtwkday, crtseason, wordscore, calcuscore, drawing, delaywordscore))))

# Depression (measured by CES-D scale)
healthsts$depress1 <- healthsts$dc009 - 1
healthsts$miss.depress <- ifelse(is.na(healthsts$dc009), 1, 0)

for (i in 10:18) {
  name <- paste("dc0", i, sep = "")
  newvar <- paste("depress", i-8, sep = "")
  if (i %in% c(13, 16)) {
    healthsts[[newvar]] <- 4 - healthsts[[name]]
  } else {
    healthsts[[newvar]] <- healthsts[[name]] - 1
  }
  healthsts[["miss.depress"]] <- ifelse(is.na(healthsts[[name]]), healthsts[["miss.depress"]]+1, healthsts[["miss.depress"]])
}

healthsts <- healthsts %>%
  mutate(CESDscore = healthsts %>%
           dplyr::select(starts_with("depress")) %>% rowSums(., na.rm = TRUE)) %>%
  mutate(CESDscore = ifelse(miss.depress == 10, NA, CESDscore))

summary(healthsts$CESDscore)

# Life Satisfication
healthsts$lifesatisfy <- as_factor(healthsts$dc028)
healthsts$lifesatisfy <- relevel(healthsts$lifesatisfy, ref = "5 Not at All Satisfied")

## load Work_Retirement_and_Pension.dta
worknretire <- read_dta("CHARLS2015/Work_Retirement_and_Pension.dta")

# wages for employee
# variable         type           unit   (note)
# ff002_1          salary         year
# ff003braket_min  salary         year
# ff003braket_max  salary         year
# ff004_1          salary         month
# ff005braket_min  salary         month
# ff005braket_max  salary         month
# ff006            salary         week
# ff008            salary         day
# ff009braket_min  salary         day
# ff009braket_max  salary         day
# ff010            salary         hour
# ff011braket_min  salary         hour
# ff011braket_max  salary         hour
# ff012_1          salary         month   (contract/performance based)
# ff013braket_min  salary         month   (contract/performance based)
# ff013braket_max  salary         month   (contract/performance based)
# ff014            bonus          year
# ff015braket_min  bonus          year
# ff015braket_max  bonus          year
# fh010            self-employed  year
# fh011braket_min  self_employed  year
# fh011braket_max  self_employed  year
# fh020            family-busins  year
# fh020braket_min  family-busins  year
# fh020braket_max  family-busins  year
# fj003            side-job       month
# fj004braket_min  side-job       month
# fj004braket_max  side_job       month
# fl009            lastjob        month
# fl010braket_min  lastjob        month
# fl010braket_max  lastjob        month
# fl011            lastjob        year
# fl011braket_min  lastjob        year
# fl011braket_max  lastjob        year
# fn005_w2_1_      pension        month
# fn005_w2_2_      pension        month
# fn005_w2_3_      pension        month

worknretire <- worknretire %>%
  dplyr::select(ff002_1, ff003bracket_min, ff003bracket_max, ff004_1, ff005bracket_min, ff005bracket_max, ff006,
                ff008, ff009bracket_min, ff009bracket_max, ff010, ff011bracket_min, ff011bracket_max,
                ff012_1, ff013bracket_min, ff013bracket_max, ff014, ff015bracket_min, ff015bracket_max,
                fh010, fh011bracket_min, fh011bracket_max, fh020, fh020_bracket_min, fh020_bracket_max,
                fj003, fj004bracket_min, fj004bracket_max, fl009, fl010bracket_min, fl010bracket_max,
                fl011, fl011_bracket_min, fl011_bracket_max, fn005_w2_1_, fn005_w2_2_, fn005_w2_3_,
                ID, householdID, communityID) %>%
  mutate(ff003bracket = (ff003bracket_min + ff003bracket_max)/2,
         ff005bracket = (ff005bracket_min + ff005bracket_max)/2,
         ff009bracket = (ff009bracket_min + ff009bracket_max)/2,
         ff011bracket = (ff011bracket_min + ff011bracket_max)/2,
         ff013bracket = (ff013bracket_min + ff013bracket_max)/2,
         ff015bracket = (ff015bracket_min + ff015bracket_max)/2,
         fh011bracket = (fh011bracket_min + fh011bracket_max)/2,
         fh020bracket = (fh020_bracket_min + fh020_bracket_max)/2,
         fj004bracket = (fj004bracket_min + fj004bracket_max)/2,
         fl010bracket = (fl010bracket_min + fl010bracket_max)/2,
         fl011bracket = (fl011_bracket_min + fl011_bracket_max)/2) %>%
  mutate(miss.wage = rowSums(is.na(cbind(ff002_1, ff003bracket, ff004_1, ff005bracket, ff006, ff008, ff009bracket, ff010, ff011bracket,
                                         ff012_1, ff013bracket, ff014, ff015bracket, fh010, fh011bracket, fh020, fh020bracket,
                                         fj003, fj004bracket, fl009, fl010bracket, fl011, fl011bracket, fn005_w2_1_, fn005_w2_2_, fn005_w2_3_)))) %>%
  mutate(indiwage = rowSums(cbind(ff002_1, ff003bracket, ff004_1 * 12, ff005bracket * 12, ff006 * 48, ff008 * 5*48, ff009bracket * 5*48,
                                  ff010 * 8*5*48, ff011bracket * 8*5*48, ff012_1 * 12, ff013bracket * 12, ff014, ff015bracket,
                                  fh010, fh011bracket, fh020, fh020bracket, fj003 * 12, fj004bracket * 12, fl009 * 12, fl010bracket * 12,
                                  fl011, fl011bracket, fn005_w2_1_ * 12, fn005_w2_2_ * 12, fn005_w2_3_ * 12), na.rm = TRUE)) %>%
  mutate(indiwage = ifelse(miss.wage == 26, NA, indiwage)) %>%
  mutate(indiwage = ifelse(indiwage < 0, 0 , indiwage))

summary(worknretire$indiwage)

## load Individual_Income.dta
indincome <- read_dta("CHARLS2015/Individual_Income.dta")

indincome <- indincome %>%
  mutate(wage = case_when(!is.na(ga002) ~ ga002,
                          is.na(ga002) & !is.na(ga002_bracket_min) & !is.na(ga002_bracket_max) ~ (ga002_bracket_min + ga002_bracket_max)/2,
                          is.na(ga002) & is.na(ga002_bracket_min) & !is.na(ga002_bracket_max) ~ ga002_bracket_max,
                          is.na(ga002) & !is.na(ga002_bracket_min) & is.na(ga002_bracket_max) ~ ga002_bracket_min),
         compensation = rowSums(cbind(ga004_1_, ga004_2_, ga004_3_, ga004_4_, ga004_5_, ga004_6_, ga004_7_, ga004_8_, ga004_9_), na.rm = TRUE),
         miss.compensation = rowSums(is.na(cbind(ga004_1_, ga004_2_, ga004_3_, ga004_4_, ga004_5_, ga004_6_, ga004_7_, ga004_8_, ga004_9_)))) %>%
  mutate(compensation = ifelse(miss.compensation == 9, NA, compensation))

indincome$indinc <- case_when(!is.na(indincome$wage) & !is.na(indincome$compensation) ~ indincome$wage + indincome$compensation,
                              !is.na(indincome$wage) & is.na(indincome$compensation) ~ indincome$wage,
                              is.na(indincome$wage) & !is.na(indincome$compensation) ~ indincome$compensation)

summary(indincome$indinc)

## load Household_Income.dta
hhincome <- read_dta("CHARLS2015/Household_Income.dta")

# HH wage
for (i in 1:10) {
  yearinc <- paste("ga006_1_", i, "_", sep = "")
  monthinc <- paste("ga006_2_", i, "_", sep = "")
  newvar <- paste("hhmemberwage", i, sep = "")
  hhincome[[newvar]] <- hhincome[[yearinc]]
  hhincome[[newvar]] <- ifelse(is.na(hhincome[[newvar]]), hhincome[[monthinc]] * 12, hhincome[[newvar]])
  print(i)
}

for (i in 1:6) {
  bracketmin <- paste("ga006_bracket_", i, "_min", sep = "")
  bracketmax <- paste("ga006_bracket_", i, "_max", sep = "")
  newvar <- paste("hhmemberwage.bracket", i, sep = "")
  hhincome[[newvar]] <- case_when(!is.na(hhincome[[bracketmin]]) & !is.na(hhincome[[bracketmax]]) ~ (hhincome[[bracketmin]]+hhincome[[bracketmax]])/2,
                                  !is.na(hhincome[[bracketmin]]) &  is.na(hhincome[[bracketmax]]) ~ hhincome[[bracketmin]],
                                  is.na(hhincome[[bracketmin]]) & !is.na(hhincome[[bracketmax]]) ~ hhincome[[bracketmax]])
  print(i)
}

hhincome$hhwage.total <- rowSums(hhincome[, grep("hhmemberwage", names(hhincome))], na.rm = TRUE)
hhincome$miss.hhwage.total <- rowSums(is.na(hhincome[, grep("hhmemberwage", names(hhincome))]))
hhincome$hhwage.total[hhincome$miss.hhwage.total == 16] <- NA

summary(hhincome$hhwage.total)

# HH compensation
for (k in 1:9) {
  yearinc <- paste("hhmember.compensation.year.type", k, sep = "")
  monthinc <- paste("hhmember.compensation.month.type", k, sep = "")
  yeardf <- hhincome[ , grep(paste("ga008_", k, "b_", sep = ""), names(hhincome))]
  monthdf <- hhincome[ , grep(paste("ga008_", k, "c_", sep = ""), names(hhincome))]
  # by year for each ID
  hhincome[["misscount"]] <- rowSums(is.na(yeardf))
  hhincome[[yearinc]][hhincome$misscount != ncol(yeardf)] <- rowSums(yeardf, na.rm = TRUE)
  # by month for each ID
  hhincome[["misscount"]] <- rowSums(is.na(monthdf))
  hhincome[[monthinc]][hhincome$misscount != ncol(monthdf)] <- rowSums(monthdf, na.rm = TRUE) * 12
  print(k)
  
  hhincome[["misscount"]] <- rowSums(is.na(hhincome[ , grep("hhmember.compensation", names(hhincome))]))
  hhincome[["hh.compensation.total"]][hhincome$misscount != ncol(hhincome[ , grep("hhmember.compensation", names(hhincome))])] <-
    rowSums(hhincome[ , grep("hhmember.compensation.year", names(hhincome))], na.rm = TRUE) + rowSums(hhincome[ , grep("hhmember.compensation.month", names(hhincome))], na.rm = TRUE) * 12
}

summary(hhincome$hh.compensation.total)

# HH agricultural income
# crobs etc
hhincome$crobs.total <- hhincome$gb005_1 - ifelse(!is.na(hhincome$gb005_2), hhincome$gb005_2, 0)
hhincome$crobs.total <- ifelse(!is.na(hhincome$gb005_1) & is.na(hhincome$gb005_2) & !is.na(hhincome$gb005_3), hhincome$gb005_1*hhincome$gb005_3/100, hhincome$crobs.total)

hhincome$crobs.bracket <- case_when(!is.na(hhincome$gb005_bracket_min) & !is.na(hhincome$gb005_bracket_max) ~ (hhincome$gb005_bracket_min + hhincome$gb005_bracket_max)/2,
                                    !is.na(hhincome$gb005_bracket_min) &  is.na(hhincome$gb005_bracket_max) ~ hhincome$gb005_bracket_min,
                                    is.na(hhincome$gb005_bracket_min) & !is.na(hhincome$gb005_bracket_max) ~ hhincome$gb005_bracket_max)
hhincome$crobs.consume.bracket <- case_when(!is.na(hhincome$gb005_w2_bracket_min) & !is.na(hhincome$gb005_w2_bracket_max) ~ (hhincome$gb005_w2_bracket_min + hhincome$gb005_w2_bracket_max)/2,
                                            !is.na(hhincome$gb005_w2_bracket_min) &  is.na(hhincome$gb005_w2_bracket_max) ~ hhincome$gb005_bracket_min,
                                            is.na(hhincome$gb005_w2_bracket_min) & !is.na(hhincome$gb005_w2_bracket_max) ~ hhincome$gb005_bracket_max)
hhincome$crobs.bracket.dropconsume <- hhincome$crobs.bracket - ifelse(!is.na(hhincome$crobs.consume.bracket), hhincome$crobs.consume.bracket, 0)

hhincome$crobs.total <- ifelse(is.na(hhincome$crobs.total), hhincome$crobs.bracket.dropconsume, hhincome$crobs.total)

summary(hhincome$crobs.total)
ggplot(data = hhincome, aes(x = ifelse(crobs.total == 0, 0, log(crobs.total)))) + geom_histogram() + labs(x = "log household agricultural income")

#livestocks
hhincome$livestock.total <- hhincome$gb011_1 - ifelse(!is.na(hhincome$gb011_2), hhincome$gb011_2, 0)
hhincome$livestock.total <- ifelse(!is.na(hhincome$gb011_1) & is.na(hhincome$gb011_2) & !is.na(hhincome$gb011_3), hhincome$gb011_1*hhincome$gb011_3/100, hhincome$livestock.total)

hhincome$livestock.bracket <- case_when(!is.na(hhincome$gb011_bracket_min) & !is.na(hhincome$gb011_bracket_max) ~ (hhincome$gb011_bracket_min + hhincome$gb011_bracket_max)/2,
                                        !is.na(hhincome$gb011_bracket_min) &  is.na(hhincome$gb011_bracket_max) ~ hhincome$gb011_bracket_min,
                                        is.na(hhincome$gb011_bracket_min) & !is.na(hhincome$gb011_bracket_max) ~ hhincome$gb011_bracket_max)

hhincome$livestock.total <- ifelse(is.na(hhincome$livestock.total), hhincome$livestock.bracket, hhincome$livestock.total)

summary(hhincome$livestock.total)

#products of livestocks
hhincome <- hhincome %>%
  mutate(liveproduct.total = gb012_1 - ifelse(!is.na(gb012_2), gb012_2, 0)) %>%
  mutate(liveproduct.total = ifelse(!is.na(gb012_1) & !is.na(gb012_2) & !is.na(gb012_3), gb012_1*gb011_3/100, liveproduct.total),
         liveproduct.bracket = case_when(!is.na(gb012_bracket_min) & !is.na(gb012_bracket_max) ~ (gb012_bracket_min+gb012_bracket_max)/2,
                                         !is.na(gb012_bracket_min) &  is.na(gb012_bracket_max) ~ gb012_bracket_min/2,
                                         is.na(gb012_bracket_min) & !is.na(gb012_bracket_max) ~ gb012_bracket_max)/2) %>%
  mutate(liveproduct.total = ifelse(is.na(liveproduct.total), liveproduct.bracket, liveproduct.total))

summary(hhincome$liveproduct.total)

#HH self-employed
hhincome <- hhincome %>%
  mutate(selfemp.total = hhincome[ , c("gc005_1_", "gc005_2_", "gc005_3_")] %>% rowSums(., na.rm = TRUE)) %>%
  mutate(selfemp.total = ifelse(is.na(gc005_1_) & is.na(gc005_2_) & is.na(gc005_3_), NA, selfemp.total))

for (i in 1:3) {
  bracketmin <- paste("gc005_bracket_", i, "_min", sep = "")
  bracketmin <- paste("gc005_bracket_", i, "_max", sep = "")
  newvar <- paste0("selfemp.bracket", i)
  hhincome[[newvar]] = case_when(!is.na(hhincome[[bracketmin]]) & !is.na(hhincome[[bracketmax]]) ~ (hhincome[[bracketmin]]+hhincome[[bracketmax]])/2,
                                 !is.na(hhincome[[bracketmin]]) &  is.na(hhincome[[bracketmax]]) ~ hhincome[[bracketmin]],
                                 is.na(hhincome[[bracketmin]]) & !is.na(hhincome[[bracketmax]]) ~ hhincome[[bracketmax]])
}

hhincome$miss.selfemp <- rowSums(is.na(hhincome[ , grep("selfemp.bracket", names(hhincome))]))
hhincome$selfemp.total <- ifelse(hhincome$miss.selfemp == 3, NA, hhincome$selfemp.total + rowSums(hhincome[ , grep("selfemp.bracket", names(hhincome))], na.rm = TRUE))

summary(hhincome$selfemp.total)

#HH Transfer Income
hhincome <- hhincome %>%
  mutate(miss.transfer = rowSums(is.na(hhincome %>%
                                         dplyr::select(gd001, contains("gd002_"), contains("gd003_")))),
         transfer.total = hhincome %>%
           dplyr::select(gd001, contains("gd002_"), contains("gd003_")) %>% rowSums(na.rm = TRUE)) %>%
  mutate(transfer.total = replace(transfer.total, miss.transfer == 11, NA))

summary(hhincome$transfer.total)

# aggregate income from multiple resources (at household level)
hh.worknretire <- worknretire %>%
  group_by(householdID) %>%
  mutate(primary.wage.total = sum(indiwage)) %>%
  dplyr::select(householdID, primary.wage.total) %>%
  unique()

# another way for aggregating at household level
# aggregate(indiwage ~ householdID, data = worknretire, sum)

hh.indincome <- indincome %>%
  group_by(householdID) %>%
  mutate(hhinc.total = sum(indinc)) %>%
  dplyr::select(householdID, hhinc.total) %>%
  unique()

income.merge <- hh.worknretire %>%
  merge(hh.indincome, by = "householdID", all = TRUE) %>%
  merge(hhincome[ , c("householdID", "hhwage.total", "hh.compensation.total", "crobs.total", "livestock.total", "liveproduct.total", "selfemp.total", "transfer.total")],
        by = "householdID", all = TRUE)

income.merge <- income.merge %>%
  mutate(miss.hhincome.all = rowSums(is.na(income.merge %>% dplyr::select(ends_with(".total")))),
         hhincome.total = rowSums(income.merge %>% dplyr::select(ends_with(".total")), na.rm = TRUE)) %>%
  mutate(hhincome.total = replace(hhincome.total, miss.hhincome.all == 9, NA)) %>%
  mutate(hhincome.total = replace(hhincome.total, hhincome.total < 0, 0))

apply(income.merge[ , -1], 2, summary)

summary(income.merge$hhincome.total)
ggplot(data = income.merge, aes(x = ifelse(hhincome.total == 0, 0 ,log(hhincome.total)))) + geom_histogram() + labs(x = "Log Household Income Total")

## load Household_Member.dta
hhmember <- read_dta("CHARLS2015/Household_Member.dta")

hhmember <- hhmember %>%
  group_by(householdID) %>%
  mutate(hh.member = n()) %>%
  dplyr::select(householdID, hh.member) %>%
  unique()

## load Weights.dta
sample.weights <- read_dta("CHARLS2015/Weights.dta")

# merge (or join) all related datasets
charls2015 <- healthsts %>%
  dplyr::select(ID, householdID, communityID,
                friends, playcards, helpothers, gotoclub, organization, volunteer, careothers, takecourse, otheract, spclass, participation,
                male, srhealth, phyfunction, ADLscore, IADLscore, MMSEscore, CESDscore, lifesatisfy)%>%
  merge(dplyr::select(demography, ID, age, married), by = "ID", all = TRUE) %>%
  merge(dplyr::select(demo2015_merge, ID, hukou2015, edulv), by = "ID", all = TRUE) %>%
  left_join(dplyr::select(hhmember, householdID, hh.member), by = "householdID") %>%
  left_join(dplyr::select(income.merge, householdID, hhincome.total), by = "householdID") %>%
  merge(dplyr::select(sample.weights, ID, INDV_weight_ad2), by = "ID", all = TRUE) %>%
  mutate(income.per.hhmember = hhincome.total/hh.member,
         agesq = age*age,
         loghh.perincome = ifelse(income.per.hhmember == 0, 0, log(income.per.hhmember))) %>%
  filter(age >= 60)

# missing value inspection and listwise
table(rowSums(is.na(charls2015[ , -c(1,2,3)])))
apply(charls2015, 2, function(x) paste0(round(mean(is.na(x))*100, digits = 2), "%"))

complete.charls2015 <- charls2015 %>%
  dplyr::select(-hhincome.total) %>%
  na.omit()

# social participation
table(complete.charls2015$spclass)
MultinomCI(table(complete.charls2015$spclass), conf.level = 0.95)*100

table(complete.charls2015$participation)
BinomCI(length(complete.charls2015$participation[complete.charls2015$participation == "1 Yes"]), nrow(complete.charls2015), conf.level = 0.95) * 100

table(complete.charls2015$spclass, complete.charls2015$participation)
summary(table(complete.charls2015$spclass, complete.charls2015$participation))

spresponse <- melt(spitems) %>%
  group_by(Var2, value) %>%
  mutate(count = n()) %>%
  dplyr::select(-Var1) %>%
  unique() %>%
  arrange(Var2, value)

# aggregated item response
item.response <- ggplot(data = spresponse, aes(x = Var2, y = count, fill = as.factor(value))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Greys", direction = -1, labels = c("Almost Daily", "Almost Every Week", "Not Regularly", "Never")) + 
  labs(x = "Social Particiaption Item", y = "Response Proportion", fill = "Item Response", caption = "Source: CHARLS 2015") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.caption = element_text(size = 12, face = "italic"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, angle = 30, vjust = 0.6),
        axis.text.y = element_text(size = 14))

print(item.response)

ggsave("table_and_figure/item_response.png", plot = item.response, width = 10, height = 8)

# participation classification
plot.spclass <- ggplot(data = complete.charls2015, aes(x = factor(spclass))) +
  geom_bar(stat = "count", width = 0.6, fill = "white", colour = "black") + 
  geom_text(stat='count', aes(label=..count..), vjust = -0.6, size = 5) +
  ylim(0,5250) +
  labs(x = "Latent Class", y = "Count", caption = " ") +
  scale_x_discrete(labels = c("Society-oriented", "Recreational", "Occasional", "Non-participated")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        plot.caption = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14))

print(plot.spclass)

plot.participation <- ggplot(data = complete.charls2015, aes(x = factor(participation))) +
  geom_bar(stat = "count", width = 0.5, fill = "white", colour = "black") + 
  geom_text(stat='count', aes(label=..count..), vjust = -0.6, size = 5) +
  ylim(0,5250) +
  labs(x = "Participator", y = "Count", caption = "Source: CHARLS 2015") +
  scale_x_discrete(labels = c("Participator", "Non-participator")) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        plot.caption = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14))

print(plot.participation)

plot.spmembership <- ggarrange(plot.spclass, plot.participation, nrow = 1, widths = c(1.5,1))

print(plot.spmembership)

ggsave("table_and_figure/spmembership.png", plot = plot.spmembership, width = 13, height = 8)

## Descriptive Statistics

# the dichotomous indicator
desc.charls2015 <- complete.charls2015 %>%
  mutate(male = as.numeric(male) - 1, urban = as.numeric(hukou2015) - 1, married = as.numeric(married) - 1) %>%
  dplyr::select(spclass, participation, male, age, urban, married, hh.member, srhealth, phyfunction, ADLscore, IADLscore, MMSEscore, CESDscore, lifesatisfy, edulv, income.per.hhmember) %>%
  dummy_cols(select_columns = c("srhealth", "lifesatisfy", "edulv"), remove_first_dummy = FALSE, remove_selected_columns = TRUE)

desc.table.full <- desc.charls2015 %>%
  stat.desc(basic = TRUE) %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::select(mean, std.dev)

desc.table.participator <- desc.charls2015 %>%
  filter(participation == "1 Yes") %>%
  stat.desc(basic = TRUE) %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::select(mean, std.dev)

desc.table.nonparticipator <- desc.charls2015 %>%
  filter(participation == "0 No") %>%
  stat.desc(basic = TRUE) %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::select(mean, std.dev)

desc.table.dichotomous <- bind_cols(desc.table.full, desc.table.participator, desc.table.nonparticipator)
rownames(desc.table.dichotomous) = rownames(desc.table.full)

desc.table.dichotomous

# the categorical indicator
desc.table.categorical <- data.frame(matrix(nrow = 27))
for (level in levels(complete.charls2015$spclass)) {
  assign("desc.table.categorical.item", desc.charls2015 %>%
           filter(spclass == level) %>%
           stat.desc(basic = TRUE) %>%
           t(.) %>%
           as.data.frame(.) %>%
           dplyr::select(mean, std.dev))
  print(level)
  desc.table.categorical <- bind_cols(desc.table.categorical, desc.table.categorical.item)
}
rownames(desc.table.categorical) = rownames(desc.table.full)

desc.table.categorical

# export dataframes into .xlsx files
desc.table.dichotomous$variable = rownames(desc.table.dichotomous)
desc.table.categorical$varaible = rownames(desc.table.categorical)

write_xlsx(desc.table.dichotomous, path = "table_and_figure/descriptives_dichotomous.xlsx", col_names = TRUE)
write_xlsx(desc.table.categorical, path = "table_and_figure/descriptives_categorical.xlsx", col_names = TRUE)

# % of participator by age group and gender
age.specific <- complete.charls2015 %>%
  mutate(agegp = case_when(age >= 60 & age < 70 ~ 1,
                           age >= 70 & age < 80 ~ 2,
                           age >= 80 & age < 90 ~ 3,
                           age >= 90 ~ 4)) %>%
  mutate(agegp = factor(agegp, levels = c(1,2,3,4), labels = c("60 - 70", "70 - 80", "80 - 90", "90 and older"))) %>%
  group_by(male, agegp) %>%
  summarise(participate = mean(as.numeric(participation)-1))

plot.age.gender <- ggplot(data = age.specific, aes(x = agegp, y = participate, group = male)) +
  geom_point() + 
  geom_line(aes(linetype = male)) +
  scale_y_continuous(labels = scales::percent, limits = c(0.1,0.4)) + 
  labs(x = "Age Group", y = "% Participator", caption = "Source: CHARLS 2015") +
  scale_linetype_discrete(name = "Gender", labels = c("Female", "Male")) +
  theme(legend.position = c(0.5, 0.1),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.text = element_text(size = 14),
        panel.background = element_rect(fill = "white", color = "black"),
        plot.caption = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14))

print(plot.age.gender)

ggsave("table_and_figure/participator_by_age_gender.png", width = 10, height = 6)

age.gender.specific <- complete.charls2015 %>%
  mutate(agegp = case_when(age >= 60 & age < 70 ~ 1,
                           age >= 70 & age < 80 ~ 2,
                           age >= 80 & age < 90 ~ 3,
                           age >= 90 ~ 4)) %>%
  mutate(agegp = factor(agegp, levels = c(1,2,3,4), labels = c("60 - 70", "70 - 80", "80 - 90", "90 and older"))) %>%
  group_by(male, agegp)

spclass.age.gender <- ggplot(data = age.gender.specific, aes(x = agegp, fill = spclass)) +
  geom_bar(stat = "count", position = "fill", color = "black") +
  facet_wrap(~ male, labeller = labeller(male = c("0 female" = "Female", "1 male" = "Male"))) +
  scale_fill_brewer(palette = "Greys", direction = 1, labels = c("Non-participated", "Occasional", "Recreational", " Society-oriented")) + 
  labs(x = "Age Group", y = "Class Proportion", fill = "Latent Class", caption = "Source: CHARLS 2015") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.caption = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14))

print(spclass.age.gender)

ggsave("table_and_figure/spclass_by_age_gender.png", spclass.age.gender, width = 13, height = 10)

## Regressions

# <! NOTE > Since it is not straightforward to conduct regressions with sampling weights and clustered standard errors in R,
#           the following part will also be conducted in Stata (and the results are the same).

# Convert the data into .dta (Stata) data
write.dta(complete.charls2015, "CHARLS2015/complete_charls2015.dta", convert.factors = "labels")

# Logistic Model (with clustered SE but no sampling weights)
formula.logit.basic <- participation ~ male + age + agesq + hukou2015 + married + hh.member
model.logit.basic <- glm.cluster(formula.logit.basic, data = complete.charls2015, family = "binomial"(link = "logit"), cluster = "communityID")
summary(model.logit.basic)

formula.logit.physical <- participation ~ male + age + agesq + hukou2015 + married + hh.member + srhealth + phyfunction + ADLscore + IADLscore + MMSEscore + CESDscore + lifesatisfy
model.logit.physical <- glm.cluster(formula.logit.physical, data = complete.charls2015, family = "binomial"(link = "logit"), cluster = "communityID")
summary(model.logit.physical)

formula.logit.sociodemographic <- as.formula(participation ~ male + age + age^2 + hukou2015 + married + hh.member + srhealth + phyfunction + ADLscore + IADLscore + MMSEscore + CESDscore + lifesatisfy + edulv + loghh.perincome)
model.logit.sociodemographic <- glm.cluster(formula.logit.sociodemographic, data = complete.charls2015, family = "binomial"(link = "logit"), cluster = "communityID")
summary(model.logit.sociodemographic)

# Multinomial Model (no clustered SEs and sampling weights)
mlogitdf <- mlogit.data(complete.charls2015, shape = "wide", choice = "spclass")

formula.mlogit.basic <- mFormula(spclass ~ 1 | male + age + agesq + hukou2015 + married + hh.member)
model.mlogit.basic <- mlogit(formula = formula.mlogit.basic, data = mlogitdf, reflevel = "4 non-participated")
summary(model.mlogit.basic)

formula.mlogit.physical <- mFormula(spclass ~ 1 | male + age + agesq + hukou2015 + married + hh.member + srhealth + phyfunction + ADLscore + IADLscore + MMSEscore + CESDscore + lifesatisfy)
model.mlogit.physical <- mlogit(formula.mlogit.physical, data = mlogitdf, reflevel = "4 non-participated")
summary(model.mlogit.physical)

formula.mlogit.sociodemographic <- mFormula(spclass ~ 1 | male + age + age^2 + hukou2015 + married + hh.member + srhealth + phyfunction + ADLscore + IADLscore + MMSEscore + CESDscore + lifesatisfy + edulv + loghh.perincome)
model.mlogit.sociodemographic <- mlogit(formula.mlogit.sociodemographic, data = mlogitdf, reflevel = "4 non-participated")
summary(model.mlogit.sociodemographic)

##################################################################
##                      Author: Chenxin Tan                     ##
##                 Last Modified: July 06, 2020                 ##
##################################################################