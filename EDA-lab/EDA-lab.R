## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::latex()

## ----options, include=FALSE----------------------------------------------
options(digits=3, width=90, stringsAsFactors = FALSE)
opts_chunk$set(echo=TRUE,tidy=FALSE,include=TRUE,
               dev='png', fig.width = 8, fig.height = 6, comment = '#>  ',
               dpi =300,	cache = T, background="grey93" )

## ----required packages and data, echo = TRUE, cache=FALSE, warning=FALSE, message=FALSE----
load(url("http://www-huber.embl.de/users/klaus/BasicR/seqZyx.rda"))
library("TeachingDemos")
data(golub, package = "multtest")
library(biomaRt)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)
#library(xlsx)
library(vioplot)
library(beeswarm)

## ----ggplot-ex1----------------------------------------------------------
    summary(iris)
    p <- ggplot(iris, aes(Sepal.Length, Sepal.Width) )

## ----ggplot-ex2----------------------------------------------------------
 p + geom_point() 

## ----ggplot-geom-char-col------------------------------------------------
 p + geom_point(shape = 5, color = "coral3") 

## ----ggplot-ex-qplot, fig.width=8, fig.height=6--------------------------
qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)


## ----qplot-set-ex, fig.width=8, fig.height=6-----------------------------
qplot(Sepal.Length, Sepal.Width, data = iris, shape = Species,
      color = I("coral3"))

## ----ggplot-ex-smoother--------------------------------------------------
ggsmooth <- (qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)
+ stat_smooth(method = "lm"))
ggsmooth 

## ----ggplot-ex-smoother-only---------------------------------------------
ggsmoothPure <- (ggplot(iris, aes(x = Sepal.Length,  y = Sepal.Width, color = Species) )
+ stat_smooth(method = "lm"))
ggsmoothPure 

## ----ggplot-ex-transformation, results ='hide', fig.keep = 'none'--------
transformed.data <- as.list(print(ggsmooth))$data[[2]]

## ----smoothInh-----------------------------------------------------------
ggsmooth2 <- (qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)
+ stat_smooth(method = "lm", aes(group = 1))) 
ggsmooth2 

## ----smoothInh2----------------------------------------------------------
ggsmooth3 <- (qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)
+ stat_smooth(method = "lm", aes(group = 1))) + geom_line()
ggsmooth3 

## ----color_examples------------------------------------------------------
ggsmooth + scale_color_brewer(palette = "Dark2")
ggsmooth + scale_color_manual(values=c(setosa = "tan3", 
                                       versicolor = "chartreuse2", 
                                       virginica = "hotpink"))

## ----load-protein-data---------------------------------------------------
proteins <- read.csv("http://www-huber.embl.de/users/klaus/BasicR/proteins.csv")[,-1]
sample_n(proteins, 4)
proteins_pMek <- subset(proteins, proteins$Target == "pMEK")
proteins_pMek_sub <- subset(proteins_pMek, proteins_pMek$Condition == "10ng/mL HGF")

## ----plot-protein-data---------------------------------------------------
proteins_pMek_sub
qplot(min, Signal, data = proteins_pMek_sub, geom = "line")
# or
#ggplot(proteins_pMek_sub, aes(min, Signal)) + geom_line()

## ----plot-protein-data-2-------------------------------------------------
qplot(min, Signal, data = proteins_pMek, geom = "line", color = Condition)
# or
#ggplot(proteins_pMek, aes(min, Signal, color = Condition) ) + geom_line()

## ----limitng vs zooming--------------------------------------------------
## limit scale
(qplot(min, Signal, data = proteins_pMek, geom = "line", color = Condition)
+ ylim(c(0,1e9)))
## retrain scale
(qplot(min, Signal, data = proteins_pMek, geom = "line", color = Condition)
+ coord_cartesian(ylim = c(0, 1e9)))

## ----iris-panels---------------------------------------------------------
ggsmooth + facet_wrap( ~ Species)

## ----iris-panels qplot---------------------------------------------------
(qplot(Sepal.Length, Sepal.Width, data = iris, 
color = Species, facets = . ~  Species))

## ----spread_example_tidyr,   echo = TRUE, eval = TRUE--------------------
proteins_spread <- proteins %>%
  select(-Sigma) %>%
  spread(key = min, value = Signal)

 sample_n(proteins_spread, 4)

## ----gahter_example,   echo = TRUE, eval = TRUE--------------------------
proteins_gathered <- proteins_spread %>% 
  gather(key = min, value = Signal, -Target, -Condition) 

sample_n(proteins_gathered, 4)

## ----ensembl-query, eval=FALSE-------------------------------------------
## ensembl <-  useMart("ensembl",dataset="hsapiens_gene_ensembl")
## seqZyx <- getSequence(id = "ENSG00000159840", type = "ensembl_gene_id",
## mart  <- ensembl, seqType = "gene_exon")[1,]
## seqZyx <- strsplit(as.character(seqZyx[1,1]), split = character(0))
## seqZyx <- seqZyx[[1]]
## save(seqZyx, file = "seqZyx.rda")

## ----Zyx-table-----------------------------------------------------------
table(seqZyx) ## table
prop.table(table(seqZyx)) ## frequency table
barplot(table(seqZyx))

## ----Zyx-table-qplot-----------------------------------------------------
dataGG <- data.frame(bases = factor(seqZyx))
ggBarplot <- qplot(bases, data = dataGG)
ggBarplot


## ----Zyx-table-ggplot, eval=FALSE----------------------------------------
## ggBarplot <- ggplot(data = dataGG,  aes(x = bases))
## ggBarplot <- ggBarplot + xlab("") + geom_bar()
## ggBarplot

## ----Zyx-pie-------------------------------------------------------------
pie(table(seqZyx))

## ----Zyx-pie-ggplot1-----------------------------------------------------
dataGG <- data.frame(bases = factor(seqZyx))

ggPie <- (qplot(x = factor(1), fill = bases, data = dataGG, xlab = NULL)
+ scale_x_discrete(breaks=NULL))
ggPie

## ----Zyx-pie-ggplot2-----------------------------------------------------
ggPie + coord_polar(theta = "y")

## ----plotCCND3, echo = TRUE----------------------------------------------
plot(golub[1042,])

## ----stripchartCCND3, echo = TRUE----------------------------------------
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub[1042,] ~ gol.fac, method="jitter", vertical = TRUE)

## ----beeswarmCCND3, echo = TRUE------------------------------------------
beeswarm(golub[1042,] ~ gol.fac, method = "center")

## ----ggplot2scatter, echo = TRUE-----------------------------------------
CCND3 <- data.frame(expCCND3 = golub[1042,], patients = gol.fac)

ggScatter <- qplot(seq_along(expCCND3), expCCND3, data = CCND3, geom ="point") 
ggScatter + xlab("Index") 

ggStrip <-  ggplot(CCND3, aes(x=patients, y=expCCND3))
ggStrip + geom_point(position =  position_jitter(w = .3, h = .0))   

## ----beeswarmCCND3_ggplot, echo = TRUE-----------------------------------
bee <- beeswarm(golub[1042,] ~ gol.fac, method = "center")

qplot(x, y, data = bee) + scale_x_continuous(breaks = 1:2,
labels = levels(gol.fac), expand = c(0, 0.5))


## ----histCCND3, echo = TRUE----------------------------------------------
hist(golub[1042, gol.fac=="ALL"], col = "lavender")

## ----ggplot2hist, echo = TRUE--------------------------------------------
CCND3.ALL <- data.frame(CCND3.ALL = golub[1042, gol.fac=="ALL"])
ggplotAbsFreq <- qplot(golub[1042, gol.fac=="ALL"], binwidth = 0.5, 
                       fill = I("darkslategrey"))

ggplotHist <- ggplot(CCND3.ALL, aes(x=CCND3.ALL), guide=F)
ggplotHist <- ggplotHist + geom_histogram(aes(y = ..density..), 
                                          fill = I("darkslategrey"),
                                          binwidth = 0.5 )


ggplotAbsFreq
ggplotHist

## ----ggplot2density, echo = TRUE-----------------------------------------
CCND3.ALL <- data.frame(CCND3.ALL = golub[1042, gol.fac=="ALL"])
 qplot(golub[1042, gol.fac=="ALL"], geom = "density")

## ----boxplotCCND3, echo = TRUE-------------------------------------------
boxplot(golub[1042,] ~ gol.fac)

## ----ggBox---------------------------------------------------------------
ggBox <-  ggplot(CCND3, aes(x=patients, y=expCCND3, color = patients))
ggBox + geom_boxplot() 

## ----quantilesCCND3,  echo = TRUE----------------------------------------
quantile(golub[1042,])
## custom quantiles
quantile(golub[1042,], probs = c(0, 0.1, 0.5))

## ----boxplot.stats,  echo = TRUE-----------------------------------------
boxplot.stats(golub[1042, gol.fac == "ALL"])
### outliers are in the "out" element of the list 

## ----Violin plot example-------------------------------------------------
vioplot(split(golub[1042,],  gol.fac)[[1]], split(golub[1042,],  gol.fac)[[2]],
        names = c("ALL", "AML"))
ggBox + geom_violin() 

## ----ecdfCCND3, fig.height=5, echo = TRUE--------------------------------
plot(ecdf(golub[1042, ]) , main = "ecdf of CCND3")

## ----ggecdfCCND3, fig.height=5, echo = TRUE------------------------------
ggECDF <-  ggplot(CCND3, aes(x=expCCND3, color = patients))
ggECDF + stat_ecdf(geom = "step", size = 2) + ylab("Frequency")

## ----QQCCND3, fig.height=5, echo = TRUE----------------------------------
qqnorm(golub[1042, gol.fac=="ALL"])
qqline(golub[1042, gol.fac=="ALL"])

## ----QQCCND3vs, fig.height=5, echo = TRUE--------------------------------
qqplot(golub[1042, gol.fac=="ALL"], golub[1042, gol.fac=="AML"])

## ----ggQQplotCCND3, fig.height=5, echo = TRUE----------------------------
ggQQ <-  ggplot(CCND3.ALL, aes(sample=CCND3.ALL)) +  stat_qq()
#ggQQ 
temp <- cbind(CCND3.ALL, 
           as.data.frame((print(ggQQ)$data)[[1]][, c( "sample", "theoretical")]))
ggQQ + stat_smooth(mapping = aes(x = theoretical, y = sample), 
                   data = temp, method = "lm", se = F )


## ----meanMedianCCND3,   echo = TRUE--------------------------------------
mean(golub[1042, gol.fac=="ALL"])
median(golub[1042, gol.fac=="ALL"])

## ----spreadCCND3,   echo = TRUE------------------------------------------
sd(golub[1042, gol.fac=="ALL"])
IQR(golub[1042, gol.fac=="ALL"]) / 1.349
mad(golub[1042, gol.fac=="ALL"])

## ----sol-simple-ggplot,   echo = TRUE,  eval = TRUE, results='hide'------
#a
##########################################################
plot_1 <- ggplot(aes( x = min, y = Signal ), data = proteins_pMek_sub)
plot_1 <- plot_1 + geom_point()
plot_1 <- plot_1  + xlab("Time [min]") + ylab("pMEK Signal")

#b
##########################################################
plot_1 <- plot_1 + geom_line()
 
#c
##########################################################
plot_1 <- plot_1 + geom_errorbar(aes(ymax = Signal+2*Sigma, 
ymin = Signal-2*Sigma))

## ----sol-facets-ggplot,   echo = TRUE,  eval = TRUE, results='hide', fig.keep='none'----
(qplot(min, Signal, data = proteins_pMek, geom = "line", color = Condition) 
+ facet_wrap( ~ Condition) + geom_errorbar(aes(ymax = Signal+2*Sigma, 
ymin = Signal-2*Sigma)))

## ----sol-complex-ggplot,   echo = TRUE,  eval = TRUE, results='hide', fig.keep='none'----
(qplot(min, Signal, data = proteins, geom = "line", color = Condition) 
+ facet_wrap( ~ Target, ncol = 2) + geom_errorbar(aes(ymax = Signal+2*Sigma, 
                                                  ymin = Signal-2*Sigma)))

## ----sol-mean-sd,   echo = TRUE,  eval = TRUE, results='hide'------------
#Use 
x<- c(1,1.5,2,2.5,3) 

#and
mean(x) 

sd(x) 
#to obtain
#that the mean is 2 and the standard deviation is 0.79

#(b) Now the mean is 7.4 and  the standard deviation dramatically increased 
#(c) The outlier increased the mean as well as the standard deviation.

## ----sol-CCND3-plot,   echo = TRUE,  eval = TRUE, results='hide', fig.keep='none'----
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub[1042,] ~ gol.fac,method="jitter")
stripchart(golub[1042,] ~ gol.fac,method="jitter",vertical = TRUE)
stripchart(golub[1042,] ~ gol.fac,method="jitter",col=c("red", "blue"),
vertical = TRUE)
stripchart(golub[1042,] ~ gol.fac,method="jitter",col=c("red", "blue"), pch="*"
,vertical = TRUE)
title("CCND3 Cyclin D3 expression value for ALL and AMl patients")

## ----sol-GeneExprssion-ggplot2-plot,   echo = TRUE, eval = TRUE, results='hide'----
#a
###############################################################
golub.df <- as.data.frame(cbind(t(golub), gol.fac))
golub.df$gol.fac <- as.factor(golub.df$gol.fac)

#b
###############################################################
rand.sample <- c(sample(dim(golub)[1],6),3052)

#c
###############################################################
# with tidyr - gather and chaining
dataForPlot <- golub.df %>% 
  select(rand.sample) %>%
  gather(key = gene, value = exp, -gol.fac) 

#e
###############################################################
pdf("boxplots.pdf", width = 14, height = 10)
p <- qplot(gol.fac, exp,  data = dataForPlot) + 
  geom_boxplot(aes(fill = gol.fac)) +
  facet_wrap(~ gene) +
  geom_point(colour = 'black', alpha = 0.5) 
p
dev.off()

## ----sol-compare,   echo = TRUE, eval = FALSE----------------------------
## #Comparing two genes
## #(a) Use
## boxplot(golub[66,]~gol.fac)
## dev.new()
## boxplot(golub[790,]~gol.fac)
## #to observe that 790 has three
## #outliers and 66 has no outlier.
## # (dev.new() opens a new graphical window)
## 
## #(b) Use
## qqnorm(golub[66,gol.fac=="ALL"])
## qqline(golub[66,gol.fac=="ALL"])
## dev.new()
## qqnorm(golub[790,gol.fac=="ALL"])
## qqline(golub[790,gol.fac=="ALL"])
## #to observe that nearly all values of 66 are on the line, where as for
## #790 the three outliers are way of the normality line. Hypothesis:
## #The expression values of 66 are normally distributed, but those of
## #row 790 are not.
## 
## #(c) Use
## mean(golub[66,gol.fac=="ALL"])
## median(golub[790,gol.fac=="ALL"])
## #and#
## mean(golub[790,gol.fac=="ALL"])
## median(golub[790,gol.fac=="ALL"])
## #The mean (-1.174024) is larger than the median (-1.28137) due to
## #outliers on the right hand side. For the gene in row 66 the mean is
## #1.182503 and the median 1.23023. The differences are smaller.

## ----sessionInfo, cache=FALSE, results='asis'----------------------------
toLatex(sessionInfo())

