#######################################################################
#
# QDA Analysis with 10 fold cross-validation
#
#######################################################################

# Set the working directory
# setwd( "/Volumes/GoogleDrive/My Drive/Research/Seashols/")
setwd( "/Users/ciararhodes/RStudio/miRNA in DNA/")

# Clear everything from memory
rm( list = ls() )

# Spread for the new "Other" data
Spread1 <- 10

# Load the Rio Package to read in the data
library( rio )
library( rpart )
library( rpart.plot )
library( stringr )
library( norm )
library( condMVNorm )
library( MASS )

# Read in the data
data1 <- import("miRNA in DNA extracts.xlsx")

data2a <- data.frame( data1, stringsAsFactors = TRUE )
data2a$Sample1 <- ifelse( str_detect( data1$Sample, "Feces"), "Feces", data2a$Sample )
data2a$Sample1 <- ifelse( str_detect( data1$Sample, "Urine"), "Urine", data2a$Sample1 )
data2a$Sample1 <- ifelse( str_detect( data1$Sample, "MB"), "Menstrual Blood", data2a$Sample1 )
data2a$Sample1 <- ifelse( str_detect( data1$Sample, "VF"), "Vaginal Fluid", data2a$Sample1 )
data2a$Color <- rep( "red", nrow( data2a) ) 
data2a$Color[data2a$Sample1 == "Feces" ] <- "brown"
data2a$Color[data2a$Sample1 == "Urine" ] <- "seagreen"
data2a$Color[data2a$Sample1 == "Menstrual Blood" ] <- "magenta"
data2a$Color[data2a$Sample1 == "Semen" ] <- "cyan"
data2a$Color[data2a$Sample1 == "Saliva" ] <- "blue"
data2a$Color[data2a$Sample1 == "Vaginal Fluid" ] <- "black"

# Get Descriptives
desc1 <- function( Var1, Fluid1, data2 ){
  data3 <- data2[ data2$Sample1 == Fluid1 , Var1 ]
  data4 <- data3[ !is.na( data3 ) ]
  n1 <- length( data4 )
  m1 <- length( data3 ) - n1
  mean1 <- mean( data4 )
  sd1 <- sd( data4 )
  min1 <- min( data4 )
  max1 <- max( data4 )
  res1 <- data.frame( Fluid1, Var1, n1, mean1, sd1, min1, max1, m1)
  return( res1 )
}

# Only want 320c 10b 891a 200b 141 412 205

# Run Blood
out1 <- desc1( "miR200b", "Blood", data2a )
out1 <- rbind( out1, desc1( "miR320c", "Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Blood", data2a ) )

# Saliva
out1 <- rbind( out1, desc1( "miR200b", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Saliva", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Saliva", data2a ) )

# Vaginal Fluid
out1 <- rbind( out1, desc1( "miR200b", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Vaginal Fluid", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Vaginal Fluid", data2a ) )

# Feces
out1 <- rbind( out1, desc1( "miR200b", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Feces", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Feces", data2a ) )

# Urine
out1 <- rbind( out1, desc1( "miR200b", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Urine", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Urine", data2a ) )

# Menstrual Blood
out1 <- rbind( out1, desc1( "miR200b", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Menstrual Blood", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Menstrual Blood", data2a ) )

# Semen
out1 <- rbind( out1, desc1( "miR200b", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR320c", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR10b", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR891a", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR141", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR412", "Semen", data2a ) )
out1 <- rbind( out1, desc1( "miR205", "Semen", data2a ) )

# Write the descriptives out to a table
write.csv( out1, "DescOutJan4other_12202021.csv", row.names = FALSE )

TotVar1 <- c( "miR200b", 
              "miR320c", 
              "miR10b", 
              "miR891a",
              "miR141",
              "miR412",
              "miR205")

# Variable to remove
filename1 <- "ConfusionJan4other_12202021.csv"

keepVar1 <- TotVar1#[ TotVar1 != RemVar1 ]
nVar1 <- length( keepVar1 )

# Build the formula
formula1 <- "Sample ~"
for( i in 1:nVar1 ){
  if( i < nVar1 ){
    formula1 <- paste( formula1, keepVar1[i], "+")
  }else{
    formula1 <- paste( formula1, keepVar1[i] )
  }
}

data2 <- data2a[,c("Sample1",keepVar1,"Color")]
data3 <- na.omit( data2 )

# Get the basic information so missing values can be plugged in.
data3f <- data3[ ,keepVar1]
data3c <- data3$Color
pairs( data3f )
pairs( data3f, col = data3c) 
plot( data3f$miR200b, data3$miR10b, col = as.factor( data3$Sample1 ), )
legend( -6, 15, unique(data3$Sample1), pch = c(1,1,1,1,1,1), col = as.factor(unique(data3$Sample1)),
        horiz = FALSE )
mean1 <- apply( data3f, 2, mean )
sigma1 <- var( data3f )

# Fill in the missing data with the conditional mean
n1 <- nrow( data2 )
outfullmean1 <- matrix( 0, nrow = n1, ncol = nVar1)
for( i in 1:n1 ){
  datatest1 <- data2[ i, keepVar1 ]
  datacheck <- sum( ifelse( is.na( datatest1 ),1,0) )
  if( datacheck > 0 & datacheck < nVar1 ){
    datatest2 <- datatest1
    Miss.ind1 <- (1:nVar1)[is.na( datatest2 )]
    Give.ind1 <- (1:nVar1)[!is.na( datatest2 ) ]
    Xgive1 <- as.vector(datatest2[ Give.ind1 ])
    CDG1 <- sigma1[Give.ind1,Give.ind1]
    CDM1 <- sigma1[Miss.ind1,Give.ind1]
    mean2 <- mean1[Miss.ind1] + CDM1%*%solve(CDG1)%*%t(Xgive1 - mean1[Give.ind1])
    out1 <- rep( 0, length( datatest2 ) )
    out1[ Miss.ind1 ] <- mean2
    out1[ Give.ind1 ] <- Xgive1 
    out2 <- as.matrix(unlist( out1 ), nrow = 1 )
  }else{
    out2 <- as.matrix( datatest1, nrow = 1 )
  }
  outfullmean1[ i , ] <- out2 
}
colnames( outfullmean1 ) <- colnames( datatest1 )


data4 <- data.frame( Sample = data2$Sample1[!is.na(outfullmean1[,1])],
                     outfullmean1[!is.na(outfullmean1[,1]),], Color = data2$Color[!is.na(outfullmean1[,1])],
                     stringsAsFactors = FALSE)
# pairs( data4[,2:6], col = data4$Color)

# Work on each of the body fluids
type1 <- unique( data4$Color)
cand1 <- cbind( runif( 10000, min( data4[,2]), max( data4[,2])),
                runif( 10000, min( data4[,3]), max( data4[,3])),
                runif( 10000, min( data4[,4]), max( data4[,4])),
                runif( 10000, min( data4[,5]), max( data4[,5])),
                runif( 10000, min( data4[,6]), max( data4[,6])),
                runif( 10000, min( data4[,7]), max( data4[,7])),
                runif( 10000, min( data4[,8]), max( data4[,8]))
)
accept1h <- matrix( 0, nrow = 10000, ncol = length( type1 ) )

for( i  in 1:length( type1 ) ){
  data4h <- data4[ data4$Color == type1[i], ]
  #head(outfullmean1)
  mean4Other <- apply( data4h[,2:8], 2, mean, na.rm = TRUE )
  Sig4Other <- var( data4h[,2:8] )
  accept1h[,i] <- dmvnorm( cand1, mean4Other, Spread1^2*Sig4Other, log = TRUE )
}
accept1max <- apply( accept1h, 2, max )
accept2a1 <- exp(accept1h[,1]-accept1max[1])
accept2a2 <- exp(accept1h[,2]-accept1max[2])
accept2a3 <- exp(accept1h[,3]-accept1max[3])
accept2a4 <- exp(accept1h[,4]-accept1max[4])
accept2 <- cbind( 1 - accept2a1, 1 - accept2a2, 1 - accept2a3, 1 - accept2a4 )
accept3 <- apply( accept2, 1, prod )
ind1 <- sample( 1:10000, 500, prob = accept3 )
samp2 <- cand1[ ind1, ]

# Add the new observations to the data
data4a <- data.frame( Sample = rep("Other", 500),
                      miR200b = samp2[,1],
                      miR320c = samp2[,2],
                      miR10b = samp2[,3],
                      miR891a = samp2[,4],
                      miR141 = samp2[,5],
                      miR412 = samp2[,6],
                      miR205 = samp2[,7],
                      Color = rep("lightblue", 500),
                      stringsAsFactors = FALSE)

data4 <- data.frame( Sample = data2$Sample1[!is.na(outfullmean1[,1])],
                     outfullmean1[!is.na(outfullmean1[,1]),], Color = data2$Color[!is.na(outfullmean1[,1])],
                     stringsAsFactors = FALSE)
#
data4 <- rbind( data4, data4a )
pairs( data4[,2:8], col = data4$Color)


###############################################################################
#
# Predictive performance
#
###############################################################################

# Get each of the classes into their own datasets for subsampling
Blood1 <- data4[ data4$Sample == "Blood", ]
Feces1 <- data4[ data4$Sample == "Feces", ]
MB1 <- data4[ data4$Sample == "Menstrual Blood", ]
Saliva1 <- data4[ data4$Sample == "Saliva", ]
Semen1 <- data4[ data4$Sample == "Semen", ]
Urine1 <- data4[ data4$Sample == "Urine", ]
VF1 <- data4[ data4$Sample == "Vaginal Fluid", ]
Other1 <- data4[ data4$Sample == "Other",]
# Get the sample sizes for each
n1.Blood1 <- nrow( Blood1 )
n1.Feces1 <- nrow( Feces1 )
n1.MB1 <- nrow( MB1 )
n1.Saliva1 <- nrow( Saliva1 )
n1.Semen1 <- nrow( Semen1 )
n1.Urine1 <- nrow( Urine1 )
n1.VF1 <- nrow( VF1 )
n1.Other1 <- nrow( Other1 )

ind1 <- rep( 1:10, 60 )

# Scramble the data in each type and assign the group
Blood2 <- cbind( Blood1[ sample( 1:n1.Blood1, n1.Blood1, replace = FALSE ),],
                 Group = ind1[1:n1.Blood1])
Feces2 <- cbind( Feces1[ sample( 1:n1.Feces1, n1.Feces1, replace = FALSE ),],
                 Group = ind1[1:n1.Feces1])
MB2 <- cbind( MB1[ sample( 1:n1.MB1, n1.MB1, replace = FALSE ),],
              Group = ind1[1:n1.MB1])
Saliva2 <- cbind( Saliva1[ sample( 1:n1.Saliva1, n1.Saliva1, replace = FALSE ),],
                  Group = ind1[1:n1.Saliva1])
Semen2 <- cbind( Semen1[ sample( 1:n1.Semen1, n1.Semen1, replace = FALSE ),],
                 Group = ind1[1:n1.Semen1])
Urine2 <- cbind( Urine1[ sample( 1:n1.Urine1, n1.Urine1, replace = FALSE ),],
                 Group = ind1[1:n1.Urine1])
VF2 <- cbind( VF1[ sample( 1:n1.VF1, n1.VF1, replace = FALSE ),],
              Group = ind1[1:n1.VF1])
Other2 <- cbind( Other1[ sample( 1:n1.Other1, n1.Other1, replace = FALSE ),],
                 Group = ind1[1:n1.Other1])

# Put all the data together
All.data1 <- rbind( Blood2, Feces2, MB2, Saliva2, Semen2, Urine2, VF2, Other2 )

# Create containers to hold the results
qda.out1 <- list()

PredictValid1 <- All.data1$Sample
# Number of groups here is 10 for 10fold cross validation.
for( j in 1:10 ){
  Train1 <- All.data1[ All.data1$Group != j, ]
  ValidIndex1 <- (1:length(All.data1$Group))[All.data1$Group == j ]
  Valid1 <- All.data1[ All.data1$Group == j, ]
  
  data4.qda <- qda( x = Train1[,keepVar1], 
                    grouping = Train1$Sample )
  pred1.qda <- predict( data4.qda, newdata = Valid1[keepVar1], method = "predictive" )$class
  PredictValid1[ ValidIndex1 ] <- as.character(pred1.qda)
  table( Valid1$Sample, pred1.qda )
  fit.table1.qda <- table( Valid1$Sample, pred1.qda ) 
  sum( diag( fit.table1.qda))/sum( fit.table1.qda )
  qda.out1[[j]] <- fit.table1.qda
  
}

# Get the summaries
qda.out2 <- qda.out1[[1]]
for( i in 2:10 ){
  qda.out2 <- qda.out2 + qda.out1[[i]]
}

sum( diag( qda.out2))/sum( qda.out2 )

qda.out2

write.csv( qda.out2, filename1)

PredictOut1 <- data.frame( All.data1[,1:7], Predicted = PredictValid1 )
write.csv( PredictOut1, "DataWPredictionsJan4other_12202021.csv", row.names = FALSE )