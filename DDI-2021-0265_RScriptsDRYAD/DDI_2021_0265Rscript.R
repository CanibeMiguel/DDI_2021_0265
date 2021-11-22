library(raster)
library(rgdal)
library(RStoolbox)
###Read images, radiometric calibration, atmospheric correction
# mask, crop and save###

###June 2000###
metaL7Jun2000<- readMeta("D:/TFMData/Landsat/LE07_L1TP_204031_20000624_20170211_01_T1/LE07_L1TP_204031_20000624_20170211_01_T1_MTL.txt")
RadCor1<-stackMeta(metaL7Jun2000)
plotRGB(RadCor1, r=3, g=2, b=1, stretch="lin")
RadCor2<-radCor(RadCor1, metaL7Jun2000, method="apref")
RadCor3<-radCor(RadCor1, metaL7Jun2000, method="costz")
RBGX<-readOGR("D:/TFMData/Shapefiles/Limite", "RB-GeresXures")
RadCor1<-crop(RadCor1, RBGX)
RadCor1<-mask(RadCor1, RBGX)
RadCor2<-crop(RadCor2, RBGX)
RadCor2<-mask(RadCor2, RBGX)
RadCor3<-crop(RadCor3, RBGX)
RadCor3<-mask(RadCor3, RBGX)
#writeRaster (RadCor1, "D:/TFMData/RasterOutputs/2000/RadCor1_Jun2000.tiff")
#writeRaster (RadCor2, "D:/TFMData/RasterOutputs/2000/RadCor2_Jun2000.tiff")
#writeRaster (RadCor3, "D:/TFMData/RasterOutputs/2000/RadCor3_Jun2000.tiff")

###September 2000###
metaL5Sep2000<- readMeta("D:/TFMData/Landsat/LT05_L1TP_204031_20000904_20171211_01_T1/LT05_L1TP_204031_20000904_20171211_01_T1_MTL.txt")
RadCor1<-stackMeta(metaL5Sep2000)
plotRGB(RadCor1, r=3, g=2, b=1, stretch="lin")
RadCor2<-radCor(RadCor1, metaL5Sep2000, method="apref")
RadCor3<-radCor(RadCor1, metaL5Sep2000, method="costz")
RadCor1<-crop(RadCor1, RBGX)
RadCor1<-mask(RadCor1, RBGX)
RadCor2<-crop(RadCor2, RBGX)
RadCor2<-mask(RadCor2, RBGX)
RadCor3<-crop(RadCor3, RBGX)
RadCor3<-mask(RadCor3, RBGX)
#writeRaster (RadCor1, "D:/TFMData/RasterOutputs/2000/RadCor1_Sep2000.tiff")
#writeRaster (RadCor2, "D:/TFMData/RasterOutputs/2000/RadCor2_Sep2000.tiff")
#writeRaster (RadCor3, "D:/TFMData/RasterOutputs/2000/RadCor3_Sep2000.tiff")

###July 2010###
metaL5Jul2010<- readMeta("D:/TFMData/Landsat/LT05_L1TP_204031_20100730_20161014_01_T1/LT05_L1TP_204031_20100730_20161014_01_T1_MTL.txt")
RadCor1<-stackMeta(metaL5Jul2010)
plotRGB(RadCor1, r=3, g=2, b=1, stretch="lin")
RadCor2<-radCor(RadCor1, metaL5Jul2010, method="apref")
RadCor3<-radCor(RadCor1, metaL5Jul2010, method="costz")
RadCor1<-crop(RadCor1, RBGX)
RadCor1<-mask(RadCor1, RBGX)
RadCor2<-crop(RadCor2, RBGX)
RadCor2<-mask(RadCor2, RBGX)
RadCor3<-crop(RadCor3, RBGX)
RadCor3<-mask(RadCor3, RBGX)
#writeRaster (RadCor1, "D:/TFMData/RasterOutputs/2000/RadCor1_Jul2010.tiff")
#writeRaster (RadCor2, "D:/TFMData/RasterOutputs/2000/RadCor2_Jul2010.tiff")
#writeRaster (RadCor3, "D:/TFMData/RasterOutputs/2000/RadCor3_Jul2010.tiff")

###October 2010###
metaL5Oct2010<- readMeta("D:/TFMData/Landsat/LT05_L1TP_204031_20101018_20180128_01_T1/LT05_L1TP_204031_20101018_20180128_01_T1_MTL.txt")
RadCor1<-stackMeta(metaL5Oct2010)
plotRGB(RadCor1, r=3, g=2, b=1, stretch="lin")
RadCor2<-radCor(RadCor1, metaL5Oct2010, method="apref")
RadCor3<-radCor(RadCor1, metaL5Oct2010, method="costz")
RadCor1<-crop(RadCor1, RBGX)
RadCor1<-mask(RadCor1, RBGX)
RadCor2<-crop(RadCor2, RBGX)
RadCor2<-mask(RadCor2, RBGX)
RadCor3<-crop(RadCor3, RBGX)
RadCor3<-mask(RadCor3, RBGX)
#writeRaster (RadCor1, "D:/TFMData/RasterOutputs/2010/RadCor1_Oct2010.tiff")
#writeRaster (RadCor2, "D:/TFMData/RasterOutputs/2010/RadCor2_Oct2010.tiff")
#writeRaster (RadCor3, "D:/TFMData/RasterOutputs/2010/RadCor3_Oct2010.tiff")

library(raster)
#leyenda
#Agua = 1
#Bosque Caducifolio = 2
#Bosque Perene = 3
#Matorral = 4
#Prados y cultivos = 5
#Suelo sin vegetacion = 6

setwd("D:/TFMData/Layers/LCMaps/")
#for (iyear in c("2000", "2010")){
#  dir.create(iyear)
#}

#for(iyear in c("2000", "2010")){
#  for(imod in c("BR", "NN", "RF", "SVM")){
#    setwd(paste("D:/TFMData/Layers/LCMaps/", iyear, sep=""))
#    dir.create(imod)
#  }
#}

#for (iyear in c("2000","2010")){
#  for (imod in c("BR", "NN", "RF", "SVM")){
#    for (irad in c("RC1", "RC2", "RC3")){
#      setwd(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, sep=""))
#      dir.create(irad)
#    }
#  }
#}


###
library(rgdal)
St_Area<-readOGR("D:/TFMData/Shapefiles/St_Area","St_Area")
my.db<-read.csv("D:/TFMData/BirdData/41598_2019_40766_MOESM5_ESM.csv", header=TRUE, sep=";", dec=".")
my.data<-subset(my.db, type == "External_Temporal_Transferability")
my.points<-my.data[,c("POINT_X", "POINT_Y")]
coordinates(my.points)= ~ POINT_X + POINT_Y
###Raster reclassification: Year 2010, I don't think I can include both years
#in the same loop
for(iyear in c("2000","2010")){
 for (imod in c("BR", "NN", "RF", "SVM")){
  for (irad in c("RC1", "RC2", "RC3")){
      r1<-raster(paste("D:/TFMData/RasterOutputs/RC", iyear, "/", imod,"/", irad,"/", irad,"_", iyear, "_", imod,".tif", sep="" ))
      LCC_1<- reclassify (r1, c(0,1,100, 1,Inf,0))
      names(LCC_1)<-"Wa"
      LCC_2<- reclassify (r1, c(0,1,0, 1,2,100, 2,Inf,0))
      names(LCC_2)<- "DeFo"
      LCC_3<- reclassify (r1, c(0,2,0, 2,3,100, 3,Inf,0))
      names(LCC_3)<-"EvFo"
      LCC_4<- reclassify (r1, c(0,3,0, 3,4,100, 4,Inf,0))
      names(LCC_4)<-"Sh"
      LCC_5<- reclassify (r1, c(0,4,0, 4,5,100, 5,Inf,0))
      names(LCC_5)<-"CrGr"
      LCC_6<- reclassify (r1, c(0,5,0, 5,6,100))
      names(LCC_6)<-"SpVe"
    
      Pred30m<- stack(LCC_1, LCC_2, LCC_3, LCC_4, LCC_5, LCC_6)
      ref<- raster("D:/TFMData/Layers/DEM/Res200/Ref200.tif")
      Pred200m<-resample(Pred30m, ref, method="bilinear")
      Pred200m<-crop(Pred200m, my.points)
      #Pred200m<-mask(Pred200m, St_Area, method="bilinear")
      writeRaster(Pred200m, paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/.tif", sep=""), bylayer=TRUE, suffix="names", overwrite=TRUE)
      
    }
  }
}


library(raster)
setwd("D:/TFMData/SDM/pred/")

#leyenda
#Agua = 1
#Bosque Caducifolio = 2
#Bosque Perene = 3
#Matorral = 4
#Prados y cultivos = 5
#Suelo sin vegetacion = 6

###First, create directories###
#for (iyear in c("2000", "2010")){
#  dir.create(iyear)
#}

#for(iyear in c("2000", "2010")){
#  for(imod in c("BR", "NN", "RF", "SVM")){
#    setwd(paste("D:/TFMData/SDM/pred/", iyear, sep=""))
#    dir.create(imod)
#  }
#}

#for (iyear in c("2000","2010")){
#  for (imod in c("BR", "NN", "RF", "SVM")){
#    for (irad in c("RC1", "RC2", "RC3")){
#      setwd(paste("D:/TFMData/SDM/pred/", iyear, "/", imod, sep=""))
#      dir.create(irad)
#    }
#  }
#}


#Now let's get terrain data from the DEM
altitude<- raster("D:/TFMData/Layers/DEM/Res200/Ref200.tif")
St_Area<-readOGR("D:/TFMData/Shapefiles/St_Area","St_Area")
altitude<-crop(altitude,my.points)
names(altitude)<-"altitude"


slope<-terrain(altitude, opt="slope", unit="degrees", neighbors = 8)
aspect<-terrain(altitude, opt="aspect", unit="degrees", neighbors = 8)

plot(altitude)
points(my.points)
plot(slope)
plot(aspect)
points(my.points)

for (iyear in c("2000", "2010")){
  for (imod in c("BR", "NN", "RF", "SVM")){
    for (irad in c("RC1", "RC2", "RC3")){
      lcc1<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_Wa.tif", sep=""))
      names(lcc1)<-"Wa"
      lcc2<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_DeFo.tif", sep=""))
      names(lcc2)<-"DeFo"
      lcc3<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_EvFo.tif", sep=""))
      names(lcc3)<-"EvFo"
      lcc4<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_Sh.tif", sep=""))
      names(lcc4)<-"Sh"
      lcc5<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_CrGr.tif", sep=""))
      names(lcc5)<-"CrGr"
      lcc6<-raster(paste("D:/TFMData/Layers/LCMaps/", iyear, "/", imod, "/", irad, "/_SpVe.tif", sep=""))
      names(lcc6)<-"SpVe"
      
      predStack<-stack(lcc1, lcc2, lcc3, lcc4, lcc5, lcc6, altitude, slope, aspect)
      
      writeRaster(predStack, paste("D:/TFMData/SDM/pred/", iyear, "/", imod, "/", irad, "/pred.tif", sep=""), overwrite=TRUE)  
    }
  }
}




########################################
#########PREDICTORS LEGEND##############
########################################
#1=Water
#2=Deciduous Forest
#3=Evergreen Forest
#4=Shrublands
#5=Croplands&Grasslands
#6=Sparse Vegetation
#7=Altitude
#8=Slope
#9=Aspect
########################################
#Step 0. Create directories to keep everything organised and get the loop to work
#for(imod in c("BR", "NN", "RF", "SVM")){
#  setwd(paste("D:/TFMData/SDM/models/", sep=""))
#  dir.create(imod)
#}

#for (imod in c("BR", "NN", "RF", "SVM")){
#   for (irad in c("RC1", "RC2", "RC3")){
#      setwd(paste("D:/TFMData/SDM/models/",imod, "/", sep=""))
#      dir.create(irad)
#    }
#  }
#Later when we model, directories with the name of each species will be created,
#after that step we will create directories for plots and outputs

#Load required packages
library(biomod2)
library(raster)

#Read species Data, create subsets for each dataset. Read predictors.
my.db<-read.csv("D:/TFMData/BirdData/41598_2019_40766_MOESM5_ESM.csv", header=TRUE, sep=";", dec=".")
my.data<-subset(my.db, type == "External_Temporal_Transferability")
View(my.data)
#Prepare the other data subsets for evaluation (part of this will be in the
#loop since we need it to change)
my.ITT.data<-subset(my.db, type=="Internal_Temporal_Transferability")
my.ITT.XY<-my.ITT.data[,c("POINT_X", "POINT_Y")]
coordinates(my.ITT.XY)= ~ POINT_X + POINT_Y

my.past.data<-subset(my.db, type=="Calibration")
my.past.XY<-my.past.data[,c("POINT_X", "POINT_Y")]
coordinates(my.past.XY)= ~ POINT_X + POINT_Y

#Finally, load the DataSplitTable
load("D:/TFMData/SDM/myDataSplitTable")

#This is the modelling loop, to compute species one by one or in pairs
#(or all at once) just modify species code in line 54.
for (imod in c("BR", "RF","SVM", "NN")){
  for (irad in c("RC1", "RC2", "RC3")){
    for (isp in c("PIBE")){
setwd(paste("D:/TFMData/SDM/models/",imod, "/", irad, "/", sep=""))
myResp<-as.numeric(paste(my.data[,isp]))
myRespName<-paste(isp)
myRespXY<-my.data[,c("POINT_X", "POINT_Y")]
myExpl<-stack(paste("D:/TFMData/SDM/pred/2010/", imod, "/", irad, "/pred.tif", sep=""))

#Evaluation Datasets
Stack.ITT.values<-extract(myExpl, my.ITT.XY)
my.ITT.data<-cbind(my.ITT.data, Stack.ITT.values)
ITT.eval.data<-my.ITT.data[,c(paste(isp), "pred.1", "pred.2", "pred.3", "pred.4",
                                   "pred.5", "pred.6", "pred.7", "pred.8", "pred.9")]

myExplPast<-stack(paste("D:/TFMData/SDM/pred/2000/", imod, "/", irad, "/pred.tif", sep=""))
Stack.past.values<-extract(myExplPast, my.past.XY)
my.past.data<-cbind(my.past.data, Stack.past.values)
past.eval.data<-my.past.data[,c(paste(isp), "pred.1", "pred.2", "pred.3", "pred.4",
                              "pred.5", "pred.6", "pred.7", "pred.8", "pred.9")]

#Step 1. Format Data
myBiomodData<- BIOMOD_FormatingData(resp.var= myResp,
                                    expl.var=myExpl,
                                    resp.xy=myRespXY,
                                    resp.name=myRespName)

myBiomodData


## Step 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

##Step 3. Computing the Models
myBiomodModelOut<- BIOMOD_Modeling(myBiomodData,
                                   models=c("ANN","GBM","RF", "GLM"),
                                   models.options=myBiomodOption,
                                   VarImport=3,
                                   models.eval.meth=c("ROC", "TSS", "KAPPA"),
                                   SaveObj=TRUE,
                                   DataSplitTable=myDataSplitTable,
                                   modeling.id=paste(myRespName,"FirstModeling", sep=""))

myBiomodModelOut
#Now we can create the directories where plots and outputs will be saved
dir.create(paste(isp,"/plots", sep=""))
dir.create(paste(isp,"/outputs", sep=""))
dir.create(paste(isp, "/maps", sep=""))
# get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
capture.output(get_evaluations(myBiomodModelOut),
               file=paste(isp, "/outputs/ModelEval.txt", sep=""))
# print variable importances
get_variables_importance(myBiomodModelOut)
capture.output(get_variables_importance(myBiomodModelOut),
               file=paste(isp, "/outputs/VariableImportance.txt",sep=""))

#Evaluate model
evaluate(myBiomodModelOut, ITT.eval.data, stat=c("ROC", "TSS", "KAPPA"))
capture.output(evaluate(myBiomodModelOut, ITT.eval.data, stat=c("ROC", "TSS", "KAPPA")),
               file=paste(isp, "/outputs/ITT_Test.txt",sep=""))

evaluate(myBiomodModelOut, past.eval.data, stat=c("ROC", "TSS", "KAPPA"))
capture.output(evaluate(myBiomodModelOut, past.eval.data, stat=c("ROC", "TSS", "KAPPA")),
               file=paste(isp, "/outputs/Past_Test.txt",sep=""))

##Step4. Ensemble Modeling
myBiomodEM <- BIOMOD_EnsembleModeling (modeling.output = myBiomodModelOut,
                                       chosen.models = 'all',
                                       em.by = 'all',
                                       eval.metric = c('ROC'),
                                       eval.metric.quality.threshold = c(0.65),
                                       models.eval.meth = c('ROC','TSS', 'KAPPA'),
                                       prob.mean = T,
                                       prob.cv = F,
                                       prob.ci = F,
                                       prob.ci.alpha = 0.05,
                                       prob.median = T,
                                       committee.averaging = F,
                                       prob.mean.weight = T,
                                       prob.mean.weight.decay = 'proportional')

# print summary
myBiomodEM
# get evaluation scores
capture.output(get_evaluations(myBiomodEM),
               file=paste(isp, "/outputs/EMevaluation.txt", sep=""))

#Evaluate ensemble model
evaluate(myBiomodEM, ITT.eval.data, stat=c("ROC", "TSS", "KAPPA"))
capture.output(evaluate(myBiomodEM, ITT.eval.data, stat=c("ROC", "TSS", "KAPPA")),
               file=paste(isp, "/outputs/EM_ITT_Test.txt",sep=""))

evaluate(myBiomodEM, past.eval.data, stat=c("ROC", "TSS", "KAPPA"))
capture.output(evaluate(myBiomodEM, past.eval.data, stat=c("ROC", "TSS", "KAPPA")),
               file=paste(isp, "/outputs/EM_Past_Test.txt",sep=""))

##Step 5. Projections
#First let's over the whole study area under current conditions
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = myExpl,
                                  proj.name ='current',
                                  selected.models ='all',
                                  binary.meth =c("ROC","TSS"),
                                  compress ='xz',
                                  clamping.mask = FALSE,
                                  output.format ='.grd')
myBiomodProj

#Let's make a projection for Year 2000
myBiomodProjPast<-BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env = myExplPast,
                                    proj.name ='past',
                                    selected.models ='all',
                                    binary.meth =c("ROC","TSS"),
                                    compress ='xz',
                                    clamping.mask = FALSE,
                                    output.format ='.grd')
?BIOMOD_Projection
myComputedModels<-myBiomodModelOut@models.computed
myBiomodModelOut


##Step6. Let's now use the ensemble model made earlier for this same projections
myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj)


myBiomodEF
#We are going to save the maps both as plots and as raster files
pdf(file=paste(isp,"/plots/current_EF.pdf", sep=""), onefile=T)
plot(myBiomodEF)
dev.off()

myBiomodEF.raster<-get_predictions(myBiomodEF)
names(myBiomodEF.raster)<-c("EM_mean", "EM_median", "EM_wmean")
writeRaster(myBiomodEF.raster,
            paste("D:/TFMData/SDM/models/", imod, "/", irad, "/", isp,"/maps/Present.tif", sep=""),
            bylayer=TRUE, suffix="names", overwrite=TRUE)

myBiomodEFPast <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                             projection.output = myBiomodProjPast)

pdf(file=paste(isp, "/plots/past_EF.pdf", sep=""), onefile=T)
plot(myBiomodEFPast)
dev.off()

myBiomodEFPast.raster<-get_predictions(myBiomodEFPast)
names(myBiomodEFPast.raster)<-c("EM_mean","EM_median","EM_wmean")
writeRaster(myBiomodEFPast.raster,
            paste("D:/TFMData/SDM/models/", imod, "/", irad, "/", isp,"/maps/Past.tif", sep=""),
            bylayer=TRUE, suffix="names", overwrite=TRUE)
    }
  }
}



#Set the working directory for the data we want to use, either dropbox or local,
#remember to coment the one you are not going to use
setwd("C:/Users/migue/Dropbox/MiguelMaster/Traballo fin de máster/SDMs/models/")
#setwd("D:/TFMData/SDM/models/")

#Create vectors with the elements we are going to need for the loops
RadCor<- c("RC1","RC2", "RC3")
ClassAlg<-c("BR", "NN", "RF", "SVM")
sp.names<-c("CPAL",  "STUR",  "CCAN",  "PVIR",  "AARV",  "TTRO",  "PMOD",  "ERUB",  "STOR",  "TMER",  "SUND", 
            "SCOM",  "SATR",  "RIGN",  "PCRI",  "PATE",  "PMAJ",  "CBRA",  "OORI",  "LCOL",  "GGLA",  "FCOE", 
            "SSER",  "CCHL",  "CCANN", "ECIA",  "PIBE" )

#Create the final data frame, we are going to be merging it with the dfs
#we generate on each run of the loop
final.df <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("AUC", "Sensibilidad", "Especificidad", "Spp", "ModAlg", "Run", "Rad", "SCAlg", "Dataset")
colnames(final.df) <- x

#Now the loop
for (alg in ClassAlg){
  for (rad in RadCor){
    for (sp in sp.names){
      
      
AccRaw <- read.table(file = paste(alg, "/", rad, "/", sp,"/outputs/ModelEval.txt", sep=""), 
                      sep = "" , header = T , nrows = 10000, na.strings ="NA",  
                      fill = TRUE, blank.lines.skip = TRUE, stringsAsFactors= F)
      
Acc <- AccRaw[AccRaw$X. %in% "ROC", c(2,4,5)]
colnames(Acc) <- c("AUC", "Sensibilidad", "Especificidad")
Acc$Spp<-sp
Acc$ModAlg<-rep(c("ANN", "GBM", "RF", "GLM")) #Algoritmos dos modelos
Acc$Run <- rep(paste("RUN",1:20, sep=""), each=length(c('ANN','GBM', 'RF', 'GLM')))##Numero de RUN
Acc$Rad <- rad
Acc$SCAlg<- alg
Acc$DataSet<-"Crossval"
      
      
SpInd_Raw<-read.table(file = paste(alg, "/", rad, "/", sp,"/outputs/ITT_Test.txt", sep=""), 
                      sep = "" , header = F , nrows = 10000, na.strings ="NA",
                      fill = TRUE, blank.lines.skip = TRUE, stringsAsFactors= F)
      
SpInd <- SpInd_Raw[SpInd_Raw$V1 %in% "ROC",c(2,4,5)]
colnames(SpInd)<-c("AUC", "Sensibilidad", "Especificidad")
SpInd$Spp<- sp
SpInd$ModAlg<-rep(c("ANN", "GBM", "RF", "GLM"))
SpInd$Run<-rep(paste("RUN",1:20, sep=""), each=length(c('ANN','GBM', 'RF', 'GLM')))
SpInd$Rad<- rad
SpInd$SCAlg<- alg
SpInd$DataSet<-"SpInd"
      
TempInd_Raw<-read.table(file = paste(alg, "/", rad, "/", sp,"/outputs/Past_Test.txt", sep=""), 
                        sep = "" , header = F , nrows = 10000, na.strings ="NA",
                        fill = TRUE, blank.lines.skip = TRUE, stringsAsFactors= F)
      
TempInd <- TempInd_Raw[TempInd_Raw$V1 %in% "ROC",c(2,4,5)]
colnames(TempInd)<-c("AUC", "Sensibilidad", "Especificidad")
TempInd$Spp<- sp
TempInd$ModAlg<-rep(c("ANN", "GBM", "RF", "GLM"))
TempInd$Run<-rep(paste("RUN",1:20, sep=""), each=length(c('ANN','GBM', 'RF', 'GLM')))
TempInd$Rad<- rad
TempInd$SCAlg<- alg
TempInd$DataSet<-"TempInd"
      
      
final.df<-rbind(final.df, Acc, SpInd, TempInd)      
    }
  }
}

View(final.df)

write.table(final.df, "C:/Users/migue/Dropbox/MiguelMaster/Traballo fin de máster/SDMs/results/AUCtable.csv", dec=".", sep=";")



### read AUC values
AUC.table<-read.table("C:/Users/miguel/Dropbox/MiguelMaster/Traballo fin de máster/SDMs/results/AUCtable.csv", sep=";")

AUC.table <- na.omit(AUC.table)
### A Generelised Linear Mixed Model with using Template Model Builder (TMB)
library(glmmTMB)
library(performance)
library(effects)



fit2 <- glmmTMB(AUC ~ Rad + 
                    (1|DataSet:Spp:Run), data=AUC.table, family = gaussian())


summary(fit2)
check_model(fit2)
y<-model_performance(fit2)
fit2

plot(allEffects(fit2))

### A Generelised Linear Mixed Model at species level
setwd("F:/TFMData/SDM/glmm.sp.effects/")
sp.names <- levels(AUC.table$Spp)
for (sp.n in sp.names){
  fit <- glmmTMB(AUC ~ ModAlg + Rad + SCAlg + 
            (1|DataSet:Run), data=AUC.table[validmodels$Spp %in% sp.n,], family = gaussian())
  print(summary(fit))
  print(model_performance(fit))
  pdf(file=paste(sp.n,"_effect.pdf", sep=""), onefile=T)
  plot(allEffects(fit))
  dev.off()
}

Spmean<-aggregate(AUC.table[, 1], list(AUC.table$Spp), mean)
Spsd<-aggregate(AUC.table[, 1], list(AUC.table$Spp), sd)
Spmin<-aggregate(AUC.table[, 1], list(AUC.table$Spp), min)
Spmax<-aggregate(AUC.table[, 1], list(AUC.table$Spp), max)
Sp.stats<-cbind(Spmean,Spsd,Spmin, Spmax)
Sp.stats<-Sp.stats[,c(1, 2, 4, 6, 8)]
colnames(Sp.stats)<-c("Spp", "Mean", "SD", "Min", "Max")          
is.num<-sapply(Sp.stats, is.numeric)
Sp.stats[is.num]<-lapply(Sp.stats[is.num], round, 3)

#Compute model performance for each species
cv<-subset(AUC.table, DataSet=="Crossval")
Spmean<-aggregate(cv[, 1], list(cv$Spp), mean)
Spsd<-aggregate(cv[, 1], list(cv$Spp), sd)
Spmin<-aggregate(cv[, 1], list(cv$Spp), min)
Spmax<-aggregate(cv[, 1], list(cv$Spp), max)
cv.stats<-cbind(Spmean,Spsd,Spmin, Spmax)
cv.stats<-cv.stats[,c(1, 2, 4, 6, 8)]
colnames(cv.stats)<-c("Spp", "Mean", "SD", "Min", "Max")          
is.num<-sapply(cv.stats, is.numeric)
cv.stats[is.num]<-lapply(cv.stats[is.num], round, 3)


#Commonality analysis
df <- data.frame(matrix(ncol=1, nrow=0))
names <- c("R2_marginal")
colnames(df) <- names
sp.names<- levels(as.factor(AUC.table$Spp))
for (sp in sp.names){
  mod <- glmmTMB (AUC  ~ModAlg + (1|DataSet:Spp:Run), data=AUC.table[AUC.table$Spp %in% sp,], family=gaussian())
  ma <- model_performance(mod)
  ma <- subset(ma, select=R2_marginal)
  rad <- glmmTMB(AUC ~Rad + (1|DataSet:Spp:Run), data=AUC.table[AUC.table$Spp %in% sp,], family = gaussian())
  r <- model_performance(rad)
  r<- subset(r, select=R2_marginal)
  sc <- glmmTMB ( AUC ~ SCAlg + (1|DataSet:Spp:Run), data=AUC.table[AUC.table$Spp %in% sp,], family = gaussian())
  sa <- model_performance(sc)
  sa <- subset(sa, select=R2_marginal)
  full <- glmmTMB (AUC ~ ModAlg + Rad + SCAlg + (1|DataSet:Spp:Run), data=AUC.table[AUC.table$Spp %in% sp,], family = gaussian())
  f <- model_performance(full)
  f <- subset(f, select=R2_marginal)
  df <- rbind(df, ma, r, sa, f)
}

data <- df
data$Spp <- rep(sp.names, each=4)
data$model <- rep(c("ModAlg", "Rad", "SCAlg", "full"))   
library(reshape2)

data<- melt(data)
d2 <- dcast(data, Spp~model)
d2 <- d2 %>% mutate_if(is.numeric, round, digits=4)
write.table(d2, "F:/TFMData/dataframes/AUC_R2.csv", dec=".", sep=";")




library(biomod2)
library(raster)
library(dplyr)
AUC.table<-read.table("C:/Users/PC/Dropbox/MiguelMaster/Traballo fin de máster/SDMs/results/AUCtable.csv", sep=";")
Sp.names<-levels(as.factor(AUC.table$Spp))
#Create empty vector to store all values
values<-vector()
setwd("D:/MCanibe/TFMData/SDM/models/BR/RC1/")

for (iSC in  c("BR", "NN", "RF", "SVM")){
  for (iRad in c("RC1", "RC2", "RC3")){
    for (iSp in Sp.names){
      setwd(paste("D:/MCanibe/TFMData/SDM/models/", iSC, "/", iRad, "/", sep=""))
      myBinProjCurrent<-brick(paste(iSp, "/proj_current/proj_current_", iSp, "_ROCbin", sep=""))
      myBinProjPast<-brick(paste(iSp, "/proj_past/proj_past_", iSp, "_ROCbin", sep=""))
      a1<-cellStats(myBinProjCurrent, "sum")
      a2<-cellStats(myBinProjPast, "sum")
      b<-((a1-a2)/(a2))*100
      values<-c(values, b)
    }
  }
}

Suitable<-subset(AUC.table, AUC.table$DataSet == "Crossval")
Suitable<-na.omit(Suitable)
Suitable<-cbind(Suitable, values)
names(Suitable)[10]<-"Suitable_Diff"
Suitable <- Suitable[, -c(1,2,3,9)]
Suitable<-Suitable %>% mutate_if(is.numeric, round, digits=3)
write.table(Suitable, "D:/MCanibe/TFMData/dataframes/Suitable.csv", dec=".", sep=";")

library(dplyr)
Suitable<-Suitable %>% mutate_if(is.numeric, round, digits=3)

Suitable <- read.table("F:/TFMData/dataframes/Suitable.csv", sep=";", header = T)
library(ggplot2)
View(Suitable)
Suit_AARV <- subset(Suitable, Spp=="AARV")
#levels(AUC.table$Rad)[levels(AUC.table$Rad)=="RC1"] <- "DN"
Suitable$Rad <- as.factor(Suitable$Rad)
levels(Suitable$Rad)[levels(Suitable$Rad)=="RC1"] <- "DN"
levels(Suitable$Rad)[levels(Suitable$Rad)=="RC2"] <- "TOA"
levels(Suitable$Rad)[levels(Suitable$Rad)=="RC3"] <- "SRef"

png("F:/TFMData/Eng_Paper_Plots/Suit_AARV.png", units="px", width=2560, height=1600, res=200)
myplot <- ggplot(Suit_AARV, aes(x=Rad, y=Suitable_Diff)) +
  geom_boxplot(aes(fill=SCAlg)) +
  facet_wrap (~ModAlg)+
  scale_fill_brewer(palette="Reds") + theme_light()+
  theme(axis.text=element_text(size=16, face="bold"),
        axis.title=element_text(size=16,face="bold"))+
  theme(strip.text.x = element_text(size=16, color="black",
                                    face="bold.italic"),
        strip.text.y = element_text(size=16, color="black",
                                    face="bold.italic"))+
  ggtitle(isp)+
  xlab("Preprocessing") + 
  ylab("Relative change in suitable area")+
  theme(legend.text=element_text(size=14))

myplot
dev.off()


Suitable<-read.table("C:/Users/PC/Dropbox/MiguelMaster/Traballo fin de máster/Result_Dataframes/Suitable.csv", sep=";", header = T)
### A Generelised Linear Mixed Model with using Template Model Builder (TMB)
library(glmmTMB)
library(performance)
library(effects)
library(see)
library(gridExtra)



fit2 <- glmmTMB(Suitable_Diff ~ SCAlg + Rad + ModAlg +
                  (1|Spp) + (1|Run), data=Suitable, family = gaussian())

hist(Suitable$Suitable_Diff)
summary(fit2)
check_model(fit2)
z<-model_performance(fit2)
z
fit2


setwd("E:/TFM_WD/GLMM/RCSA_bySpp/")
for (model in c("full", "rad", "SC", "mod")) {
  for (sp in sp.names) {
    setwd(paste("E:/TFM_WD/GLMM/RCSA_bySpp/", model, "/", sep = ""))
    dir.create(sp)
  }
}
setwd("E:/TFM_WD/GLMM/RCSA_bySpp/")
for (sp.n in sp.names){
  full_model <- glmmTMB(Suitable_Diff ~ ModAlg + Rad + SCAlg +
                          (1|Run), data=Suitable[Suitable$Spp %in% sp.n,], family = gaussian(), REML = FALSE)
  null_rad <- glmmTMB(Suitable_Diff ~ ModAlg + SCAlg +
                        (1|Run), data=Suitable[Suitable$Spp %in% sp.n,], family = gaussian(), REML=FALSE)
  null_SC <- glmmTMB(Suitable_Diff ~ ModAlg + Rad +
                       (1|Run), data=Suitable[Suitable$Spp %in% sp.n,], family = gaussian(), REML=FALSE)
  null_mod <- glmmTMB(Suitable_Diff ~ Rad + SCAlg +
                        (1|Run), data=Suitable[Suitable$Spp %in% sp.n,], family = gaussian(), REML=FALSE)
  
  sum_full <- summary(full_model)
  coef_full <- as.data.frame(sum_full$coefficients$cond)
  write.table(coef_full, file = paste("full/", sp.n, "/Coefs.csv", sep = ""), dec = ".", sep = ";")
  perf_full <- model_performance(full_model)
  write.table(perf_full, file = paste("full/", sp.n, "/Perf.csv", sep = ""), dec = ".", sep = ";")
  
  sum_rad <- summary(null_rad)
  coef_rad <- as.data.frame(sum_rad$coefficients$cond)
  write.table(coef_rad, file = paste("rad/", sp.n, "/Coefs.csv", sep = ""), dec = ".", sep = ";")
  perf_rad <- model_performance(null_rad)
  write.table(perf_rad, file = paste("rad/", sp.n, "/Perf.csv", sep = ""), dec = ".", sep = ";")
  
  sum_SC <- summary(null_SC)
  coef_SC <- as.data.frame(sum_SC$coefficients$cond)
  write.table(coef_SC, file = paste("SC/", sp.n, "/Coefs.csv", sep = ""), dec = ".", sep = ";")
  perf_SC <- model_performance(null_SC)
  write.table(perf_SC, file = paste("SC/", sp.n, "/Perf.csv", sep = ""), dec = ".", sep = ";")
  
  sum_mod <- summary(null_mod)
  coef_mod <- as.data.frame(sum_mod$coefficients$cond)
  write.table(coef_mod, file = paste("mod/", sp.n, "/Coefs.csv", sep = ""), dec = ".", sep = ";")
  perf_mod <- model_performance(null_mod)
  write.table(perf_mod, file = paste("mod/", sp.n, "/Perf.csv", sep = ""), dec = ".", sep = ";")
  
  lrt_rad <- anova(full_model, null_rad)
  write.table(lrt_rad, file = paste("rad/", sp.n, "/LRT.csv", sep = ""), dec = ".", sep = ";")
  lrt_SC <- anova(full_model, null_SC)
  write.table(lrt_SC, file = paste("SC/", sp.n, "/LRT.csv", sep = ""), dec = ".", sep = ";")
  lrt_mod <- anova(full_model, null_mod)
  write.table(lrt_mod, file = paste("mod/", sp.n, "/LRT.csv", sep = ""), dec = ".", sep = ";")
  
}

#Now gather significance on LRT analysis
names <- colnames(lrt_mod)
df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(df) <- names
setwd("E:/TFM_WD/GLMM/")
for (model in c("mod", "rad","SC")) {
  for (sp in sp.names) {
    LRT <- read.table(paste("RCSA_bySpp/", model, "/", sp, "/LRT.csv", sep = ""), sep=";")
    df <- rbind(df,LRT)
  }
}
df <- na.omit(df)
df$Spp <- rep(sp.names)
df$factor <- rep(c("mod", "rad", "SC"), each=27)
df$Pvalue <- cut(df$Pr..Chisq., breaks = c(0, 0.001, 0.01, 0.05, 1),
                 labels = c("***", "**", "*", "."))
write.table(df, file = "RCSA_bySpp/LRT_results.csv", dec=".", sep=";")










      
     



