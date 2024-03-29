#################################################
## Date: 09/18/2014
## Title: Multiple linear regression block bootstrapping
#################################################

install.packages("hydroGOF") ## all GOF functions
install.packages("lubridate") ## date package
library(hydroGOF)
library(lubridate)

#All data file location.  Enter block size and change header to TRUE or FALSE
file.location = "G:/Statistics/R/mydata/modelEvaluation/silverRainfallpost1933.csv"
# file.location = "G:/Statistics/R/mydata/modelEvaluation/ECFT_Output.csv"

#Path for output files.  This directory must be existing
out.path <- "G:/Statistics/R/output/modelEvaluation/"
block_size <- 48
no_bootstrap_realiz <- 500
bootstrap = "Y" # (Y/N) bootstrap the results or just go through different consecutive periods "N"
dep.var <-  "Elev_M.0048"  # "Elev_M.0048", "Discharge_Silver.River.German"  y-variable in the text file
dep.type = "Water Level (NAVD88, ft)" #  "Water Level (NAVD88, ft)", "Discharge (cfs)"

calib.start = as.Date("1933-01-01")
calib.end = as.Date("1970-12-31")
resid.start = as.Date("1990-01-01")
resid.end = as.Date("1999-12-31")

alldata <- read.csv(file.location, header = TRUE)
alldata$YearMonth = as.Date(alldata$YearMonth, format = "%m/%d/%Y")
alldata$Month = as.factor(alldata$Month)
calib.data = subset(alldata, subset = alldata$YearMonth < calib.end & alldata$YearMonth> calib.start)
n_row = nrow(calib.data)
n_col<- ncol(alldata)

#convert all data to numeric so that it can be compared to pred.all
mat.alldata <- data.matrix(alldata[dep.var])
mat.calib.data <- data.matrix(calib.data[dep.var])
alldata.y <- as.numeric(mat.alldata)
calib.data.y <- as.numeric(mat.calib.data)
n_blocks <- as.integer(n_row/block_size)

#Make a data frame for gof stats so that all realizations can be easily compared
df.gof <- data.frame(matrix(nrow=21, ncol=no_bootstrap_realiz)) 
sum.stat.table = data.frame(date.start = as.Date(character()), 
                            date.end = as.Date(character()), 
                            R2 = numeric(), 
                            DF = numeric(), 
                            intercept = numeric(),
                            slope = numeric(),
                            RMSE = numeric())

#Loop for the number of bootstrap realizations
for (r in 1:no_bootstrap_realiz) {
  
  #Make a data frame to store the bootstrap sample in
  boot = data.frame(matrix(vector(), 0, length(names(calib.data)),
                           dimnames=list(c(), names(calib.data))), stringsAsFactors=F)
  df_boot = data.frame(matrix(vector(), 0, length(names(calib.data)),
                              dimnames=list(c(), names(calib.data))), stringsAsFactors=F)
  index = block_size
  
  if(bootstrap == "Y"){
    for (j in 1:n_blocks){
      block <- as.integer(runif(1, 0, n_blocks))
      boot = calib.data[seq((block*block_size),(block*block_size+block_size)),]
      df_boot = rbind(df_boot, boot)
    }
  }else{
    if ((r-1) > n_blocks){break}
    boot = calib.data[seq(((r-1)*block_size),((r-1)*block_size+block_size)),]
    df_boot = boot
  }
  #write the output bootstrap data frame to a text file
  #   out.file <- paste(out.path, "BlockbootRealiz", r, ".csv", sep= "")
  #   write.csv(df_boot, file=out.file)
  
  ####################Additional user input
  #Note: this equation will have to be changed based on the column names of the dependent variable and regressors
  fit_boot <- lm(Elev_M.0048 ~ 
                   # Precip24_Month,  # for simple linear regression 24 month accumlated time series
                   # Month +  # can include dummy variables but doesn't work
                   Precip_Ocala +
                   Precip_Lag1_Ocala + 
                   Precip_Lag2_Ocala + 
                   Precip_Lag3_Ocala + 
                   Precip_Lag4_Ocala + 
                   Precip_Lag5_Ocala + 
                   Precip_Lag6_Ocala + 
                   Precip_Lag7_Ocala + 
                   Precip_Lag8_Ocala + 
                   Precip_Lag9_Ocala + 
                   Precip_Lag10_Ocala + 
                   Precip_Lag11_Ocala + 
                   Precip_Lag12_Ocala + 
                   Precip_Lag13_Ocala + 
                   Precip_Lag14_Ocala + 
                   Precip_Lag15_Ocala + 
                   Precip_Lag16_Ocala + 
                   Precip_Lag17_Ocala + 
                   Precip_Lag18_Ocala + 
                   Precip_Lag19_Ocala + 
                   Precip_Lag20_Ocala + 
                   Precip_Lag21_Ocala + 
                   Precip_Lag22_Ocala + 
                   Precip_Lag23_Ocala + 
                   Precip_Lag24_Ocala + 
                   Precip_Lag25_Ocala + 
                   Precip_Lag26_Ocala + 
                   Precip_Lag27_Ocala + 
                   Precip_Lag28_Ocala + 
                   Precip_Lag29_Ocala + 
                   Precip_Lag30_Ocala + 
                   Precip_Lag31_Ocala + 
                   Precip_Lag32_Ocala + 
                   Precip_Lag33_Ocala + 
                   Precip_Lag34_Ocala + 
                   Precip_Lag35_Ocala,
                 na.action = na.omit, data=df_boot)
  
  ####################End additional user input
  
  # Export summary of regression
  summary.file <- paste(out.path, "RegSummary_", r, ".csv", sep= "")
  out<-capture.output(summary(fit_boot))
  cat(out,file=summary.file,sep="\n",append=TRUE)
  
  sum.stat = data.frame(date.start = min(df_boot$YearMonth), 
                        date.end = max(df_boot$YearMonth), 
                        R2 = format(reg.sum$adj.r.squared, digits = 3), 
                        DF = format(reg.sum$df[2]), 
                        intercept = fit_boot$coefficients[1], 
                        slope = fit_boot$coefficients[2], 
                        RMSE = format(sqrt(mean(fit_boot$residuals^2)), digits = 3))
  sum.stat.table = rbind(sum.stat.table, sum.stat)
  
  ## Export the predictions on the entire data set
  pred.all <- predict.lm(fit_boot, alldata) #this is the prediction on the entire data set
  pred.calib <- predict.lm(fit_boot,calib.data) #this is the prediction on the calibration data set
  #   pred.all.file <- paste(out.path,"AllDataPred_", r, ".csv", sep="")
  #   write.csv(pred.all, file=pred.file)
  
  #### FIGURES ###
  
  ## FIGURE - MLR Coefficient graphic
  jpeg.out.name <- paste(out.path,"Coefficient_", r, ".jpg", sep= "")
  jpeg(jpeg.out.name, width=8, height=8, units= 'in', res= 300)
  par(mfrow =c(2,1))  
  reg.sum = summary(fit_boot)
  reg.r2 = reg.sum$adj.r.squared
  reg.sum.p = 1-reg.sum$coefficients[,4]
  reg.sum.sum = reg.sum$coefficients[,1]
  barplot(height = as.vector(reg.sum.sum)[2:37], names.arg = c(0,seq(2:36)), 
          xlab = "Lag", ylab = expression(beta), main = substitute(paste("Coefficient Estimate (", beta, ")")))
  barplot(height = as.vector(reg.sum.p)[2:36], names.arg = c(0,seq(2:35)), xlab = "Lag", 
          main = substitute(paste("Significance (1-", italic(p), ")", sep = "")), 
          ylab = substitute(paste("Significance (1-", italic(p), ")", sep = "")))
  abline(h = 0.95, col = "red")
  dev.off()
  
  ## FIGURE - Observed vs Simulated Times series
  jpeg.out.name <- paste(out.path,"ObsVSimCalib_", r, ".jpg", sep= "")
  jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
  simCalib = zoo(as.vector(pred.calib[1:n_row]), calib.data$YearMonth[1:n_row])
  obsCalib = zoo(calib.data.y, calib.data$YearMonth[1:n_row])
  calib.resid = (simCalib-obsCalib)
  plot(calib.resid)  
  plot(obsCalib, calib.resid)
  ggof(simCalib,obsCalib, na.rm=TRUE, 
       ylab = dep.type, 
       xlab = "Date", tick.tstep = "years", lab.tstep = "years", lab.fmt = "%Y"
       #        ,cal.ini = c(min(df_boot$YearMonth, na.rm = TRUE), max(df_boot$YearMonth, na.rm = TRUE))
  )
  dev.off()
    
  ## FIGURE - Residuals time series with average resdual between resid.start and resid.end
  jpeg.out.name <- paste(out.path,"Residual_",r,"_", year(calib.start),".", year(calib.end), ".jpg", sep= "")
  jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
  sim = zoo(as.vector(pred.all), alldata$YearMonth)
  obs = zoo(alldata.y, alldata$YearMonth)
  resid = (sim - obs)
  avg.resid = round(median(resid,  na.rm = TRUE), digits = 2)
  resid.dates = seq(resid.start, resid.end, by = "day")  
  resid1 = resid[resid.dates]
  avg.resid1 = round(median(resid1,  na.rm = TRUE), digits = 2)
    plot(resid, ylab = "Residual(Sim - Obs)", xlab = "Year")
  legend("topleft", paste("Avg Resid (", year(resid.start), "-", year(resid.end),") =", avg.resid1, sep = ""), 
         bty = "n",col = c("black"), cex = 0.75, merge = F, border=F)
  abline(h = c(0), col = c("black"), lty = c(1))
  lines(x = c(resid.start, resid.end), y = c(avg.resid1, avg.resid1), type = "l", col = c("red"), lty = c(2))
  abline(v = c(resid.start,resid.end), col = c("red", "red"), lty = c(1,1))
  dev.off()
  
  ## FIGURE - CHECK RESDIDUALS FOR NON-LINEARITY
  #    plot(obsCalib,calib.resid, main = "Residuals vs. Observed") # shows non-linear relationship
  
  ## FIGURE - Observed vs simulated over whole period of record
  jpeg.out.name <- paste(out.path,"ObsVSim_", r, ".jpg", sep= "")
  jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
  ggof(sim,obs, na.rm=TRUE, ylab=dep.type, xlab = "Date", 
       tick.tstep = "years", lab.tstep = "years", lab.fmt = "%Y", 
       cal.ini = c(min(df_boot$YearMonth, na.rm = TRUE), 
                   max(df_boot$YearMonth, na.rm = TRUE)))  # vert red line calibration period
  dev.off()
  
  gof.stats = gof(simCalib, obsCalib, na.rm=TRUE) # gof stats for individual model
  df.gof[r]<- rbind(gof.stats, avg.resid) # creates data batle of all gof.stats for each model
  #   write.csv(gof.stats, file=gof.file)
  #   gof.file <- paste(out.path,"AllDataGOF_", r, ".csv", sep="")
}

## FIGURE - R2 model metric. The R2 can be that of only the bootstrap data "as.numeric(as.vector(sum.stat.table$R2))" or over the calibration period "R2.num"
mat.R2 <- data.matrix(df.gof[17, 1:no_bootstrap_realiz])
R2.num <- as.numeric(mat.R2)
jpeg.out.name <- paste(out.path,"Stats_R2.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
# hist(as.numeric(as.vector(sum.stat.table$R2)), 
#      xlab = substitute(paste("Coefficient of Determination ",R^2)), 
#      main = substitute(paste("Distribution ", R^2)))
hist(x = R2.num, xlab = substitute(paste("Coefficient of Determination ",R^2)), 
     main = substitute(paste("Distribution ", R^2)))
dev.off()

## FIGURE - Distirubtion of average residual between  resid.start and resid.end
jpeg.out.name <- paste(out.path,"Stats_avg.resid.jpg", sep= "")
jpeg(jpeg.out.name, width=12, height=4, units= 'in', res= 300)
par(mfrow = c(1,2))
t.df.gof = t(df.gof) # transpose the df.gof for next step
resid.good = as.numeric(subset(t.df.gof, subset = t.df.gof[,7] > 0.7)[,21]) # average residuals of models with R2 greater than 0.7
avg.resid.num <- as.numeric(data.matrix(df.gof[21, 1:no_bootstrap_realiz])) # residual of all models
hist(x = avg.resid.num, xlab = "Average Residual", 
     main = paste("Distribution Residual (", year(resid.start),"-", year(resid.end), ")", sep = ""))
boxplot(x = avg.resid.num, ylab = "Average Residual", 
        main = paste("Distribution Residual (", year(resid.start),"-", year(resid.end), ")", sep = ""))
dev.off()

# quantile(resid.good, c(1, 0.95, .75,.5,.25,0.05,0))
quantile(avg.resid.num, c(1, 0.95, .75,.5,.25,0.05,0)) # quantile of residuals

## FIGURES OF METRICS
mat.PBIAS <- data.matrix(df.gof[6, 1:no_bootstrap_realiz])
PBIAS.num <- as.numeric(mat.PBIAS)
jpeg.out.name <- paste(out.path,"Stats_PBIAS.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
hist(x = PBIAS.num,xlab = "Percent Bias (PBIAS)", main = "Distribution PBIAS")
dev.off()

mat.MAE <- data.matrix(df.gof[2, 1:no_bootstrap_realiz])
MAE.num <- as.numeric(mat.MAE)
jpeg.out.name <- paste(out.path,"Stats_MAE.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
hist(x = MAE.num,xlab = "Mean Absolute Error (MAE)", main = "Distribution MAE")
dev.off()

mat.MSE <- data.matrix(df.gof[3, 1:no_bootstrap_realiz])
MSE.num <- as.numeric(mat.MSE)
jpeg.out.name <- paste(out.path,"Stats_MSE.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
hist(x = MSE.num,xlab = "Mean Squared Error (MSE)", main = "Distribution MSE")
dev.off()

mat.RMSE <- data.matrix(df.gof[4, 1:no_bootstrap_realiz])
RMSE.num <- as.numeric(mat.RMSE)
jpeg.out.name <- paste(out.path,"Stats_RMSE.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
hist(x = RMSE.num,xlab = "Root Mean Squared Error (RMSE)", main = "Distribution RMSE")
dev.off()

mat.NSE <- data.matrix(df.gof[9, 1:no_bootstrap_realiz])
NSE.num <- as.numeric(mat.NSE)
jpeg.out.name <- paste(out.path,"Stats_NSE.jpg", sep= "")
jpeg(jpeg.out.name, width=8, height=4, units= 'in', res= 300)
hist(x = NSE.num,xlab="Nash Sutcliffe Efficiency (NSE)", main = "Distribution NSE")
dev.off()

gof.file <- paste(out.path,"AllDataGOF_Summary.csv", sep="")
write.csv(df.gof, file=gof.file)

## FIGURE - regression plot of one variable
# jpeg.out.name <- paste(out.path,"Regression.jpg", sep= "")
# jpeg(jpeg.out.name, width=7, height=7, units= 'in', res= 300)
# plot(calib.data[["Elev_M.0048"]] ~ calib.data[["Precip24_Month"]], xlab = "Precip 24 Month Accumulation", ylab = "Elev M-0048")
# abline(fit_boot)
# final.reg = fit_boot
# legend(x = "topleft", bty= "n",
#          legend = c(paste("R2 =", format(summary(final.reg)$adj.r.squared, digits = 3)), 
#                     paste("DF =", format(summary(final.reg)$df[2])),
#                     paste("RMSE =", format(summary(final.reg)$sigma, digits = 2))
#                     )
#   )
# dev.off()