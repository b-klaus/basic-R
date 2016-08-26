## ----style, echo=FALSE, results='asis'-----------------------------------
BiocStyle::latex()

## ----options, include=FALSE----------------------------------------------
options(digits=3, width=85, stringsAsFactors = TRUE )
opts_chunk$set(echo=TRUE,tidy=FALSE,include=TRUE,
               fig.path='Rlab-',dev='png', 
		 fig.width = 12, fig.height = 9, comment = '#>', dpi = 300,
		cache = TRUE, lazy.load = FALSE, background="grey93", message = FALSE )

## ----required packages and data, echo = TRUE, cache = F------------------
library(TeachingDemos)
library(xlsx)
library(multtest)
library(Biobase)
library(plyr)
library(dplyr)
library(ggplot2)

## ----ex-O-1, echo = TRUE-------------------------------------------------
4 + 6

## ----ex-O-2, echo = TRUE-------------------------------------------------
x <- 6
y <- 4
z <- x+y
z

## ----ex-O-3, echo = TRUE-------------------------------------------------
ls()

## ----ex-O-4, echo = TRUE-------------------------------------------------
sqrt(16)

## ----ex-O-5, echo = TRUE-------------------------------------------------
z <- c(5,9,1,0)

## ----ex-O-5b, echo = TRUE------------------------------------------------
x <- c(5,9)
y <- c(1,0)
z <- c(x,y)

## ----ex-O-6, echo = TRUE-------------------------------------------------
seq(1,9,by=2)
seq(8,20,length=6)

## ----ex-O-7, echo = TRUE-------------------------------------------------
x <- seq(1,10)

## ----ex-O-8, echo = TRUE-------------------------------------------------
rep(1:3,6)
## repeat each element six times
rep(1:3,c(6,6,6))
## simplified
rep(1:3,rep(6,3))

## ----vec-----------------------------------------------------------------
x <- 1:5; y <- 5:1
x
y
x + y
x^2
# another example
x <- c(6,8,9)
y <- c(1,2,4)
x + y
x * y

## ----vec-2---------------------------------------------------------------
x <- c(6,8,9)
x + 2

## ----calcEx, eval = FALSE------------------------------------------------
## x + cos(pi/y)

## ----summary-1,   echo = TRUE--------------------------------------------
x <- c(7.5,8.2,3.1,5.6,8.2,9.3,6.5,7.0,9.3,1.2,14.5,6.2)
mean(x)
var(x)
summary(x)

## ----subscr,   echo = TRUE-----------------------------------------------
x[1:6]
x[7:12]
summary(x[1:6])
summary(x[7:12])

## ----subscr-2,   echo = TRUE---------------------------------------------
x[c(2,4,9)]
x[-(1:6)]
# compare to 
x[7:12]

## ----sort-rank,   echo = TRUE--------------------------------------------
x <- c(1.3,3.5,2.7,6.3,6.3)
sort(x)
order(x)
x[order(x)]
rank(x)

## ----object-examples,   echo = TRUE--------------------------------------
#assign value "9" to an object
a <- 9
# is a a string?
is.character(a) 
# is a a number?
is.numeric(a) 
# What's its type?
typeof(a)
# now turn it into a factor
a <- as.factor(a)
# Is it a factor?
is.factor(a)
# assign an string to a: 
a <- "NAME"
# what's a?
class(a)
str(a) 

## ----cbind-ex,   echo = TRUE---------------------------------------------
x <- c(5,7,9)
y <- c(6,3,4)
z <- cbind(x,y)
z
## dimensions: 3 rows and 2 columns
dim(z)
### matrix constructor
z <- matrix(c(5,7,9,6,3,4),nrow=3)

## ----Matrix-ex,   echo = TRUE--------------------------------------------
z <- matrix(c(5,7,9,6,3,4),nr=3,byrow=T)
z

## ----Matrix-op,   echo = TRUE--------------------------------------------
y <- matrix(c(1,3,0,9,5,-1),nrow=3,byrow=T)
y
y + z
y * z

## ----Matrix-op-2,   echo = TRUE------------------------------------------
x <- matrix(c(3,4,-2,6),nrow=2,byrow=T)
x
y %*% x

## ----Matrix-op-3,   echo = TRUE------------------------------------------
t(z)
solve(x)

## ----Matrix-op-4,   echo = TRUE------------------------------------------
z[1,1]
z[,2]
z[1:2,]
z[-1,]
z[-c(1,2),]

## ----load-Patients,   echo = TRUE----------------------------------------
pat <- read.csv("http://www-huber.embl.de/users/klaus/BasicR/Patients.csv")
pat
str(pat)

colnames(pat)
### equivalent
names(pat)

rownames(pat)

### a factor
pat$Gender

## ------------------------------------------------------------------------
pat <- read.csv("http://www-huber.embl.de/users/klaus/BasicR/Patients.csv",
              stringsAsFactors = FALSE)
### not a factor
pat$Gender


## ----load-Patients-excel,   echo = TRUE----------------------------------
pat.xls<-read.xlsx("Patients.xls", sheetIndex=1)
pat.xls
str(pat.xls)

## ----which-patients,   echo = TRUE---------------------------------------
## which patients are  less than 1.5 tall?
which(pat$Height<1.5) 
## How tall are they?
pat[which(pat$Height<1.5),]
### recommended alternative that is less error prone
pat[pat$Height<1.5,]
### access by rownames
low.height.rows <- subset(rownames(pat), pat$Height<1.5)
low.height.rows
pat[low.height.rows,]
### access via subset
subset(pat, Height<1.5)

## ----list-example,   echo = TRUE-----------------------------------------
L = list(one=1, two=c(1,2), five=seq(1, 4, length=5))
L
names(L)
 L$five + 10
## equivalent
L[[3]] + 10

## ----list-example-2,   echo = TRUE---------------------------------------
pat$Height
#equivalent
pat[[1]]

## ----apply-example,   echo = TRUE----------------------------------------
# Calculate mean for each of the first two columns
sapply(X = pat[,1:2], FUN = mean, na.rm = TRUE)
# Mean height separately for each gender
tapply(X = pat$Height, FUN = mean, INDEX = pat$Ge) 


## ----plot-example,   echo = TRUE-----------------------------------------
#pdf(file="plot-example.pdf", width=12, height=6)
x <- seq(-3,3, by = 0.1); y <- cos(x); z <- sin(x)
plot(x,y, type="l", col="darkviolet", main="Cos and Sin")
points(x, rep(0, length(x)), pch=3)
lines(x,z, type="l", col="magenta")
legend("topleft", c("cos","sin"), col=c("darkviolet",
"magenta"), lty=1)
#dev.off()

## ----ggplot-ex1----------------------------------------------------------
    summary(iris)
    p <- ggplot(iris, aes(Sepal.Length, Sepal.Width) )

## ----ggplot-ex2----------------------------------------------------------
 p + geom_point() 

## ----ggplot-ex-qplot, fig.width=8, fig.height=6--------------------------
qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)


## ----ggplot-ex-smoother--------------------------------------------------
ggsmooth <- (qplot(Sepal.Length, Sepal.Width, data = iris, color = Species)
+ stat_smooth(method = "lm"))
ggsmooth 

## ----ggplot-ex-transformation, results ='hide', fig.keep = 'none'--------
transformed.data <- as.list(print(ggsmooth))$data[[2]]

## ----curr-conv,   echo = TRUE, eval = TRUE-------------------------------
euro.calc<-function(x, currency="US") {
  ## currency has a default argrument "US"
  if(currency=="US") return(x*1.33)
  if(currency=="Pounds") return(x*0.85)
}
euro.calc(100) ## how many dollars are 100 Euros?

## ----if-example,   echo = TRUE, eval = TRUE------------------------------
w= 3
	if( w < 5 ){
	d=2
	}else{
	d=10
	}
d

## ----for-example,   echo = TRUE, eval = TRUE-----------------------------
h <- seq(from = 1, to = 8)
s <- integer() # create empty integer vector
	for(i in 2:10)
	{
	s[i] <- h[i] * 10
	}
s

## ----ifelse-example,   echo = TRUE, eval = TRUE--------------------------
s <- seq(from = 1, to = 10)
binary.s <- ifelse(s > 5, "high", "low")
binary.s

## ----sol-1,   echo = TRUE, results='hide'--------------------------------
rep(6,6)
rep(c(5,8),4)
c(rep(5,4), rep(8,4))

## ----calcExerc, results='hide'-------------------------------------------
x + cos(pi/y)

## ----sol-2,   echo = TRUE, results='hide'--------------------------------
x <-  -0.25
y <- 2
x + cos(pi/y)

## ----sol-3,   echo = TRUE, results = 'hide'------------------------------
y <- c(33,44,29,16,25,45,33,19,54,22,21,49,11,24,56)
# day of the week summary, example: Tuesday
Tuesday <- y[ (1:3) + 3 ]
summary(Tuesday)
## Shop 2 summary 
Shop2 <- y[ seq(2, length(y), by = 3 )  ]
summary(Shop2)

## ----sol-Matrix,   echo = TRUE, eval = FALSE-----------------------------
## #b
## ########################################################################
## a <- 3
## b <- 4.5
## 
## #c
## ########################################################################
## 
## is.numeric(a)
## is.character(b)
## 
## #d
## ########################################################################
## 
## A <- matrix(seq(1,9), nrow = 3, byrow = TRUE )
## A[3,3] <- 10
## #A<-matrix(c(1,2,3,4,5,6,7,8,10),nrow=3)
## 
## B <- matrix(seq(1,9), nrow = 3, byrow = FALSE)
## #B<-matrix(c(1,2,3,4,5,6,7,8,10),nrow=3, byrow=TRUE)
## 
## 
## y<-matrix(c(1,2,3), nrow=3)
## 
## #e
## ########################################################################
## a^2 + 1 /b
## a*A
## A %*% B
## 
## det(A)
## #if the determinant of a matrix is zero, it cannot be inverted
## solve(A)
## t(A)
## 
## #f
## ########################################################################
## 
## A[2,3]
## B[3,2]
## 
## #g
## ########################################################################
## 
## # element-wise
## #A[1,]*B[,2]
## 
## # a^T * b
## A[1,]%*%B[,2]
## 
## # b^T*a
## fr.A <- as.matrix(A[1,])
## sc.B <- as.matrix(t(B[,2]))
## 
## fr.A %*% sc.B
## 

## ----apply-test,   echo = TRUE, results = 'hide'-------------------------
#a
########################################################################

pat <- read.csv("http://www-huber.embl.de/users/klaus/BasicR/Patients.csv")
pat

#b
########################################################################

is.data.frame(pat)
summary(pat)
str(pat)

#c
########################################################################

head(pat)
pat
str(pat)
#d
########################################################################

summary(pat$Weight)
#There is a NA value, which is easy to spot, since the data set 
# is really small! 
#Otherwise, just use the which()-Funktion ...
which(is.na(pat$Weight))

## other NA methods
#?na.omit

## remove  NAs from Weight
na.omit(pat$Weight)

#e
########################################################################

#Pay attention to access the data set directly ...
pat[2,2] <- mean(pat$Weight, na.rm=TRUE)
pat


#f
########################################################################

BMI <- pat[,2] / pat[,1]^2

### alternatively

BMI = pat$Weight / pat$Height^2

pat$Weight[2] = mean(pat$Weight, na.rm=TRUE)

print(BMI)

### attach BMI to the data frame
pat <- cbind(pat, BMI)
pat

## ----sol-plot,   echo = TRUE, eval = FALSE-------------------------------
## x <- seq(from=-2, to=2, by=0.2)
## length(x)
## x
## stand.normal <- dnorm(x, mean=0, sd=1)
## # or: stand.normal<-dnorm(x)
## length(stand.normal)
## stand.normal
## 
## #visualize it
## #
## plot(x,stand.normal, type="l")
## 
## plot(x,stand.normal, type="b")
## plot(x,stand.normal, type="h", col = "darkgreen")
## plot(x,stand.normal, type="h", col = "darkgreen", main = "Standard Normal Density")
## 
## # use qplot
## qplot(x, stand.normal, color = I("darkgreen"))

## ----sol-read,   echo = TRUE, results = 'hide'---------------------------
#a
########################################################################

source("http://www-huber.embl.de/users/klaus/BasicR/readError.R")

test <- readError(1000)
## number of errors
sum(test == "error")
##  error probability
sum(test == "error") / 1000

prop.table(table(test))


#b
########################################################################
readError2 <- function(noBases){

  positions <- integer(noBases) ## initialize vector
		for(i in 1:noBases){	
		positions[i] <- rbinom(n=1, size = 1, prob = 0.05)
		}
  return(ifelse(positions, "correct", "error"))
	}



### equivalent function
rbinom(n=1000, size =1, prob = 0.05)

## ----sessionInfo, cache=FALSE, results='asis'----------------------------
toLatex(sessionInfo())

