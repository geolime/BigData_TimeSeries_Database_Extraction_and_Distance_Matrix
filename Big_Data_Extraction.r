# Created by Leif Holmquist
# Date: 2019-10-01
# This R script is designed to extract millions of records over multiple time series with differing variables standards.'
# This script additionally calculates the distance between the nearest train stations of an input file with station coordinates.
# This was designed for analysing the commuting distances and demographics of workers who commuted to Malmö, SE and analysing 
# commuting distances and demographics of workers who lived in Malmö and commuted within Malmö or elsewhere in Sweden.


# Set working directory
setwd("")

# Load libraries
library(RODBC);
library(tidyverse);
library(dplyr);
library(taRifx);
library(bazar);
library(ggplot2);
library(raster);
library(foreign);

# Open and set the database connection to DBMS
conn <- odbcConnect("");

# Load the coordinate file of the Malmö stations
StationCoords <- read.csv(file = paste("StationCoords.csv", sep=""))
############################################################################################################
# Employees who live in Malmö from 1990-1997
for (year in 1990:1997){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who live in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE Municipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Malmö
  EmployeeM <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for all businesses in Sweden
  WorkplaceS <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceS)[names(WorkplaceS) == "Grid"] <- "wGrid"
  names(WorkplaceS)[names(WorkplaceS) == "Municipality"] <- "wMunicipality"
  names(WorkplaceS)[names(WorkplaceS) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "City"] <- "wCity"
  names(WorkplaceS)[names(WorkplaceS) == "Namn"] <- "wNamn"
  names(EmployeeM)[names(EmployeeM) == "XKOORD"] <- "eYKOORD"
  names(EmployeeM)[names(EmployeeM) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceS <- subset.data.frame(WorkplaceS, Organization_Nr > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wGrid > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wYKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wXKOORD > 0)
  WorkplaceS <- unique(WorkplaceS)
  EmployeeM <- subset.data.frame(EmployeeM, Resident_Nr > 0)
  EmployeeM <- subset.data.frame(EmployeeM, Grid > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eYKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eXKOORD > 0)
  EmployeeM <- unique(EmployeeM)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceS, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeM, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between Employees, their workplaces and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i,"eYKOORD"])
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i, "wYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(EmployeeCoords, WorkplaceCoords, lonlat=FALSE), pointDistance(EmployeeCoords, Centralen, lonlat=FALSE), pointDistance(EmployeeCoords, Triangeln, lonlat=FALSE), pointDistance(EmployeeCoords, Hyllie, lonlat=FALSE), pointDistance(EmployeeCoords, Svagertorp, lonlat=FALSE), pointDistance(EmployeeCoords, Rosengard, lonlat=FALSE), pointDistance(EmployeeCoords, Persborg, lonlat=FALSE), pointDistance(EmployeeCoords, Ostra, lonlat=FALSE), pointDistance(EmployeeCoords, Oxie, lonlat=FALSE))
  } 
  
  # Merge Euclidean distance to main variable
  colnames(EuclideanDist) <- c("Resident_Nr", "WorkDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Resident_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Resident_Distances",year,".dta", sep=""));
}
############################################################################################################
# Employees who live in Malmö from 1998-2003
for (year in 1998:2003){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who live in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE Municipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Malmö
  EmployeeM <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for all businesses in Sweden
  WorkplaceS <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceS)[names(WorkplaceS) == "Grid"] <- "wGrid"
  names(WorkplaceS)[names(WorkplaceS) == "Municipality"] <- "wMunicipality"
  names(WorkplaceS)[names(WorkplaceS) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "City"] <- "wCity"
  names(WorkplaceS)[names(WorkplaceS) == "Namn"] <- "wNamn"
  names(EmployeeM)[names(EmployeeM) == "XKOORD"] <- "eYKOORD"
  names(EmployeeM)[names(EmployeeM) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceS <- subset.data.frame(WorkplaceS, Organization_Nr > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wGrid > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wYKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wXKOORD > 0)
  WorkplaceS <- unique(WorkplaceS)
  EmployeeM <- subset.data.frame(EmployeeM, Resident_Nr > 0)
  EmployeeM <- subset.data.frame(EmployeeM, Grid > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eYKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eXKOORD > 0)
  EmployeeM <- unique(EmployeeM)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceS, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeM, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between Employees, their workplaces and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i,"eYKOORD"])
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i, "wYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(EmployeeCoords, WorkplaceCoords, lonlat=FALSE), pointDistance(EmployeeCoords, Centralen, lonlat=FALSE), pointDistance(EmployeeCoords, Triangeln, lonlat=FALSE), pointDistance(EmployeeCoords, Hyllie, lonlat=FALSE), pointDistance(EmployeeCoords, Svagertorp, lonlat=FALSE), pointDistance(EmployeeCoords, Rosengard, lonlat=FALSE), pointDistance(EmployeeCoords, Persborg, lonlat=FALSE), pointDistance(EmployeeCoords, Ostra, lonlat=FALSE), pointDistance(EmployeeCoords, Oxie, lonlat=FALSE))
  } 
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "WorkDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Resident_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Resident_Distances",year,".dta", sep=""));
}
############################################################################################################
# Employees who live in Malmö from 2004
for (year in 2004){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who live in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE Municipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Malmö
  EmployeeM <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for all businesses in Sweden
  WorkplaceS <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceS)[names(WorkplaceS) == "Grid"] <- "wGrid"
  names(WorkplaceS)[names(WorkplaceS) == "Municipality"] <- "wMunicipality"
  names(WorkplaceS)[names(WorkplaceS) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "City"] <- "wCity"
  names(WorkplaceS)[names(WorkplaceS) == "Namn"] <- "wNamn"
  names(EmployeeM)[names(EmployeeM) == "XKOORD"] <- "eYKOORD"
  names(EmployeeM)[names(EmployeeM) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceS <- subset.data.frame(WorkplaceS, Organization_Nr > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wGrid > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wYKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wXKOORD > 0)
  WorkplaceS <- unique(WorkplaceS)
  EmployeeM <- subset.data.frame(EmployeeM, Resident_Nr > 0)
  EmployeeM <- subset.data.frame(EmployeeM, Grid > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eYKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eXKOORD > 0)
  EmployeeM <- unique(EmployeeM)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceS, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeM, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between Employees, their workplaces and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i,"eYKOORD"])
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i, "wYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(EmployeeCoords, WorkplaceCoords, lonlat=FALSE), pointDistance(EmployeeCoords, Centralen, lonlat=FALSE), pointDistance(EmployeeCoords, Triangeln, lonlat=FALSE), pointDistance(EmployeeCoords, Hyllie, lonlat=FALSE), pointDistance(EmployeeCoords, Svagertorp, lonlat=FALSE), pointDistance(EmployeeCoords, Rosengard, lonlat=FALSE), pointDistance(EmployeeCoords, Persborg, lonlat=FALSE), pointDistance(EmployeeCoords, Ostra, lonlat=FALSE), pointDistance(EmployeeCoords, Oxie, lonlat=FALSE))
  } 
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "WorkDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Resident_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Resident_Distances",year,".dta", sep=""));
}
############################################################################################################
# Employees who live in Malmö from 2005-2014
for (year in 2005:2014){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who live in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE Municipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Malmö
  EmployeeM <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for all businesses in Sweden
  WorkplaceS <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceS)[names(WorkplaceS) == "Grid"] <- "wGrid"
  names(WorkplaceS)[names(WorkplaceS) == "Municipality"] <- "wMunicipality"
  names(WorkplaceS)[names(WorkplaceS) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "City"] <- "wCity"
  names(WorkplaceS)[names(WorkplaceS) == "Namn"] <- "wNamn"
  names(EmployeeM)[names(EmployeeM) == "XKOORD"] <- "eYKOORD"
  names(EmployeeM)[names(EmployeeM) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceS <- subset.data.frame(WorkplaceS, Organization_Nr > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wGrid > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wYKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wXKOORD > 0)
  WorkplaceS <- unique(WorkplaceS)
  EmployeeM <- subset.data.frame(EmployeeM, Resident_Nr > 0)
  EmployeeM <- subset.data.frame(EmployeeM, Grid > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eYKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eXKOORD > 0)
  EmployeeM <- unique(EmployeeM)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceS, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeM, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between Employees, their workplaces and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i,"eYKOORD"])
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i, "wYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(EmployeeCoords, WorkplaceCoords, lonlat=FALSE), pointDistance(EmployeeCoords, Centralen, lonlat=FALSE), pointDistance(EmployeeCoords, Triangeln, lonlat=FALSE), pointDistance(EmployeeCoords, Hyllie, lonlat=FALSE), pointDistance(EmployeeCoords, Svagertorp, lonlat=FALSE), pointDistance(EmployeeCoords, Rosengard, lonlat=FALSE), pointDistance(EmployeeCoords, Persborg, lonlat=FALSE), pointDistance(EmployeeCoords, Ostra, lonlat=FALSE), pointDistance(EmployeeCoords, Oxie, lonlat=FALSE))
  } 
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "WorkDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income_K04 <- as.integer(DistMerge$Disposable_Income_K04)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Resident_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Resident_Distances",year,".dta", sep=""));
}
############################################################################################################
# Employees who live in Malmö from 2015-2017
for (year in 2015:2017){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who live in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE Municipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Malmö
  EmployeeM <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for all businesses in Sweden
  WorkplaceS <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceS)[names(WorkplaceS) == "P0794_Organization_Nr"] <- "Organization_Nr"
  names(WorkplaceS)[names(WorkplaceS) == "Grid"] <- "wGrid"
  names(WorkplaceS)[names(WorkplaceS) == "WorkMunicipality"] <- "wMunicipalityKOD"
  names(WorkplaceS)[names(WorkplaceS) == "Ruta_XKOORD"] <- "wYKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "Ruta_YKOORD"] <- "wXKOORD"
  names(WorkplaceS)[names(WorkplaceS) == "AstDeso"] <- "wDeso"
  names(WorkplaceS)[names(WorkplaceS) == "AstCity"] <- "wCity"
  names(WorkplaceS)[names(WorkplaceS) == "astCity_2015_Namn"] <- "wNamn"
  names(EmployeeM)[names(EmployeeM) == "P0794_Resident_Nr"] <- "Resident_Nr"
  names(EmployeeM)[names(EmployeeM) == "Municipality"] <- "MunicipalityKOD"
  names(EmployeeM)[names(EmployeeM) == "XKOORD"] <- "eYKOORD"
  names(EmployeeM)[names(EmployeeM) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceS <- subset.data.frame(WorkplaceS, Organization_Nr > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wGrid > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wYKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wXKOORD > 0)
  WorkplaceS <- subset.data.frame(WorkplaceS, wMunicipalityKOD != "Rest")
  WorkplaceS <- unique(WorkplaceS)
  EmployeeM <- subset.data.frame(EmployeeM, Resident_Nr > 0)
  EmployeeM <- subset.data.frame(EmployeeM, Grid > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eYKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, eXKOORD > 0)
  EmployeeM <- subset.data.frame(EmployeeM, MunicipalityKOD != "Rest")
  EmployeeM <- unique(EmployeeM)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceS, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeM, WorkplaceMerge, by="Resident_Nr")
  
  # Fix Municipality code
  EmployeeMerge$Municipality <- paste(EmployeeMerge$Municipality, EmployeeMerge$MunicipalityKOD, sep="")
  EmployeeMerge$wMunicipality <- paste(EmployeeMerge$WorkMunicipality, EmployeeMerge$wMunicipalityKOD, sep="")
  EmployeeMerge <- subset.data.frame(EmployeeMerge, select = -c(wDeso,wMunicipalityKOD,MunicipalityKOD))
  
  # Calculate the Euclidean distance between Employees, their workplaces and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i,"eYKOORD"])
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i, "wYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(EmployeeCoords, WorkplaceCoords, lonlat=FALSE), pointDistance(EmployeeCoords, Centralen, lonlat=FALSE), pointDistance(EmployeeCoords, Triangeln, lonlat=FALSE), pointDistance(EmployeeCoords, Hyllie, lonlat=FALSE), pointDistance(EmployeeCoords, Svagertorp, lonlat=FALSE), pointDistance(EmployeeCoords, Rosengard, lonlat=FALSE), pointDistance(EmployeeCoords, Persborg, lonlat=FALSE), pointDistance(EmployeeCoords, Ostra, lonlat=FALSE), pointDistance(EmployeeCoords, Oxie, lonlat=FALSE))
  } 
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "WorkDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income_K04 <- as.integer(DistMerge$Disposable_Income_K04)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  DistMerge$Capital_Income <- as.integer(DistMerge$Capital_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Resident_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Resident_Distances",year,".dta", sep=""));
}
############################################################################################################
# Malmö workplaces for employees in Sweden from 1990-1997
for (year in 1990:1997){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who work in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE WorkMunicipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Sweden
  EmployeeS <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for businesses in Malmö
  WorkplaceM <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceM)[names(WorkplaceM) == "Grid"] <- "wGrid"
  names(WorkplaceM)[names(WorkplaceM) == "Municipality"] <- "wMunicipality"
  names(WorkplaceM)[names(WorkplaceM) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "City"] <- "wCity"
  names(WorkplaceM)[names(WorkplaceM) == "Namn"] <- "wNamn"
  names(EmployeeS)[names(EmployeeS) == "XKOORD"] <- "eYKOORD"
  names(EmployeeS)[names(EmployeeS) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceM <- subset.data.frame(WorkplaceM, Organization_Nr > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wGrid > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wYKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wXKOORD > 0)
  WorkplaceM <- unique(WorkplaceM)
  EmployeeS <- subset.data.frame(EmployeeS, Resident_Nr > 0)
  EmployeeS <- subset.data.frame(EmployeeS, Grid > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eYKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eXKOORD > 0)
  EmployeeS <- unique(EmployeeS)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceM, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeS, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between the Workplaces, their employees and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i,"wYKOORD"])
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i, "eYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(WorkplaceCoords, EmployeeCoords, lonlat=FALSE), pointDistance(WorkplaceCoords, Centralen, lonlat=FALSE), pointDistance(WorkplaceCoords, Triangeln, lonlat=FALSE), pointDistance(WorkplaceCoords, Hyllie, lonlat=FALSE), pointDistance(WorkplaceCoords, Svagertorp, lonlat=FALSE), pointDistance(WorkplaceCoords, Rosengard, lonlat=FALSE), pointDistance(WorkplaceCoords, Persborg, lonlat=FALSE), pointDistance(WorkplaceCoords, Ostra, lonlat=FALSE), pointDistance(WorkplaceCoords, Oxie, lonlat=FALSE))
  }
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "EmployeeDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Workplace_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Workplace_Distances",year,".dta", sep=""));
}
############################################################################################################
# Malmö workplaces for employees in Sweden from 1998-2003
for (year in 1998:2003){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who work in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE WorkMunicipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Sweden
  EmployeeS <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for businesses in Malmö
  WorkplaceM <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceM)[names(WorkplaceM) == "Grid"] <- "wGrid"
  names(WorkplaceM)[names(WorkplaceM) == "Municipality"] <- "wMunicipality"
  names(WorkplaceM)[names(WorkplaceM) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "City"] <- "wCity"
  names(WorkplaceM)[names(WorkplaceM) == "Namn"] <- "wNamn"
  names(EmployeeS)[names(EmployeeS) == "XKOORD"] <- "eYKOORD"
  names(EmployeeS)[names(EmployeeS) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceM <- subset.data.frame(WorkplaceM, Organization_Nr > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wGrid > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wYKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wXKOORD > 0)
  WorkplaceM <- unique(WorkplaceM)
  EmployeeS <- subset.data.frame(EmployeeS, Resident_Nr > 0)
  EmployeeS <- subset.data.frame(EmployeeS, Grid > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eYKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eXKOORD > 0)
  EmployeeS <- unique(EmployeeS)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceM, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeS, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between the Workplaces, their employees and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i,"wYKOORD"])
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i, "eYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(WorkplaceCoords, EmployeeCoords, lonlat=FALSE), pointDistance(WorkplaceCoords, Centralen, lonlat=FALSE), pointDistance(WorkplaceCoords, Triangeln, lonlat=FALSE), pointDistance(WorkplaceCoords, Hyllie, lonlat=FALSE), pointDistance(WorkplaceCoords, Svagertorp, lonlat=FALSE), pointDistance(WorkplaceCoords, Rosengard, lonlat=FALSE), pointDistance(WorkplaceCoords, Persborg, lonlat=FALSE), pointDistance(WorkplaceCoords, Ostra, lonlat=FALSE), pointDistance(WorkplaceCoords, Oxie, lonlat=FALSE))
  }
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "EmployeeDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Workplace_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Workplace_Distances",year,".dta", sep=""));
}
############################################################################################################
# Malmö workplaces for employees in Sweden from 2004
for (year in 2004){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who work in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE WorkMunicipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Sweden
  EmployeeS <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for businesses in Malmö
  WorkplaceM <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceM)[names(WorkplaceM) == "Grid"] <- "wGrid"
  names(WorkplaceM)[names(WorkplaceM) == "Municipality"] <- "wMunicipality"
  names(WorkplaceM)[names(WorkplaceM) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "City"] <- "wCity"
  names(WorkplaceM)[names(WorkplaceM) == "Namn"] <- "wNamn"
  names(EmployeeS)[names(EmployeeS) == "XKOORD"] <- "eYKOORD"
  names(EmployeeS)[names(EmployeeS) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceM <- subset.data.frame(WorkplaceM, Organization_Nr > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wGrid > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wYKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wXKOORD > 0)
  WorkplaceM <- unique(WorkplaceM)
  EmployeeS <- subset.data.frame(EmployeeS, Resident_Nr > 0)
  EmployeeS <- subset.data.frame(EmployeeS, Grid > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eYKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eXKOORD > 0)
  EmployeeS <- unique(EmployeeS)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceM, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeS, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between the Workplaces, their employees and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i,"wYKOORD"])
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i, "eYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(WorkplaceCoords, EmployeeCoords, lonlat=FALSE), pointDistance(WorkplaceCoords, Centralen, lonlat=FALSE), pointDistance(WorkplaceCoords, Triangeln, lonlat=FALSE), pointDistance(WorkplaceCoords, Hyllie, lonlat=FALSE), pointDistance(WorkplaceCoords, Svagertorp, lonlat=FALSE), pointDistance(WorkplaceCoords, Rosengard, lonlat=FALSE), pointDistance(WorkplaceCoords, Persborg, lonlat=FALSE), pointDistance(WorkplaceCoords, Ostra, lonlat=FALSE), pointDistance(WorkplaceCoords, Oxie, lonlat=FALSE))
  }
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "EmployeeDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt <- as.numeric(DistMerge$Gender_Wt)
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id <- as.integer(DistMerge$Disposable_Income_per_Id)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income <- as.integer(DistMerge$Disposable_Income)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Workplace_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Workplace_Distances",year,".dta", sep=""));
}
############################################################################################################
# Malmö workplaces for employees in Sweden from 2005-2014
for (year in 2005:2014){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who work in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE WorkMunicipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Sweden
  EmployeeS <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for businesses in Malmö
  WorkplaceM <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year," WHERE Municipality LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceM)[names(WorkplaceM) == "Grid"] <- "wGrid"
  names(WorkplaceM)[names(WorkplaceM) == "Municipality"] <- "wMunicipality"
  names(WorkplaceM)[names(WorkplaceM) == "XKOORD"] <- "wYKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "YKOORD"] <- "wXKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "City"] <- "wCity"
  names(WorkplaceM)[names(WorkplaceM) == "Namn"] <- "wNamn"
  names(EmployeeS)[names(EmployeeS) == "XKOORD"] <- "eYKOORD"
  names(EmployeeS)[names(EmployeeS) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceM <- subset.data.frame(WorkplaceM, Organization_Nr > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wGrid > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wYKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wXKOORD > 0)
  WorkplaceM <- unique(WorkplaceM)
  EmployeeS <- subset.data.frame(EmployeeS, Resident_Nr > 0)
  EmployeeS <- subset.data.frame(EmployeeS, Grid > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eYKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eXKOORD > 0)
  EmployeeS <- unique(EmployeeS)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceM, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeS, WorkplaceMerge, by="Resident_Nr")
  
  # Calculate the Euclidean distance between the Workplaces, their employees and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i,"wYKOORD"])
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i, "eYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(WorkplaceCoords, EmployeeCoords, lonlat=FALSE), pointDistance(WorkplaceCoords, Centralen, lonlat=FALSE), pointDistance(WorkplaceCoords, Triangeln, lonlat=FALSE), pointDistance(WorkplaceCoords, Hyllie, lonlat=FALSE), pointDistance(WorkplaceCoords, Svagertorp, lonlat=FALSE), pointDistance(WorkplaceCoords, Rosengard, lonlat=FALSE), pointDistance(WorkplaceCoords, Persborg, lonlat=FALSE), pointDistance(WorkplaceCoords, Ostra, lonlat=FALSE), pointDistance(WorkplaceCoords, Oxie, lonlat=FALSE))
  }
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "EmployeeDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income_K04 <- as.integer(DistMerge$Disposable_Income_K04)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Workplace_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Workplace_Distances",year,".dta", sep=""));
}
############################################################################################################
# Workplaces for individuals in Sweden from 2015-2017
for (year in 2015:2017){
  # Query statistics databasse table for employment status and save to employed variable
  EmployedStat <- sqlQuery(conn, paste("SELECT TOP (1) c.name AS column_name FROM sys.views AS t INNER JOIN sys.columns c ON t.OBJECT_ID = c.OBJECT_ID WHERE t.name = 'StateDB",year,"' AND c.name LIKE 'EmployedStat%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query statistics databse table for residents who work in Malmö municipality (code = 1280) AND Employment status is TRUE
  StateDB <- sqlQuery(conn, paste("SELECT * FROM dbo.StateDB",year," WHERE WorkMunicipality = 1280 AND ",EmployedStat," = 1;", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query housing geodatabase for residents in Sweden
  EmployeeS <- sqlQuery(conn, paste("SELECT * FROM dbo.housing_geodata",year,";", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  # Query business geodatabase for businesses in Malmö
  WorkplaceM <- sqlQuery(conn, paste("SELECT * FROM dbo.business_geodata",year," WHERE AstDeso LIKE '1280%';", sep=""), stringsAsFactors = FALSE, as.is = TRUE);
  
  # Rename fields and Fix X and Y Coordinates
  names(WorkplaceM)[names(WorkplaceM) == "P0794_Organization_Nr"] <- "Organization_Nr"
  names(WorkplaceM)[names(WorkplaceM) == "Grid"] <- "wGrid"
  names(WorkplaceM)[names(WorkplaceM) == "WorkMunicipality"] <- "wMunicipalityKOD"
  names(WorkplaceM)[names(WorkplaceM) == "Ruta_XKOORD"] <- "wYKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "Ruta_YKOORD"] <- "wXKOORD"
  names(WorkplaceM)[names(WorkplaceM) == "AstDeso"] <- "wDeso"
  names(WorkplaceM)[names(WorkplaceM) == "AstCity"] <- "wCity"
  names(WorkplaceM)[names(WorkplaceM) == "astCity_2015_Namn"] <- "wNamn"
  names(EmployeeS)[names(EmployeeS) == "P0794_Resident_Nr"] <- "Resident_Nr"
  names(EmployeeS)[names(EmployeeS) == "Municipality"] <- "MunicipalityKOD"
  names(EmployeeS)[names(EmployeeS) == "XKOORD"] <- "eYKOORD"
  names(EmployeeS)[names(EmployeeS) == "YKOORD"] <- "eXKOORD"
  
  # Remove erroneous data
  WorkplaceM <- subset.data.frame(WorkplaceM, Organization_Nr > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wGrid > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wYKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wXKOORD > 0)
  WorkplaceM <- subset.data.frame(WorkplaceM, wMunicipalityKOD != "Rest")
  WorkplaceM <- unique(WorkplaceM)
  EmployeeS <- subset.data.frame(EmployeeS, Resident_Nr > 0)
  EmployeeS <- subset.data.frame(EmployeeS, Grid > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eYKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, eXKOORD > 0)
  EmployeeS <- subset.data.frame(EmployeeS, MunicipalityKOD != "Rest")
  #EmployeeS <- unique(EmployeeS)
  
  # Merge employee and workplace information
  WorkplaceMerge <- merge.data.frame(WorkplaceM, StateDB, by="Organization_Nr")
  EmployeeMerge <- merge.data.frame(EmployeeS, WorkplaceMerge, by="Resident_Nr")
  
  # Fix Municipality code for Workplace and Employees
  EmployeeMerge$Municipality <- paste(EmployeeMerge$Municipality, EmployeeMerge$MunicipalityKOD, sep="")
  EmployeeMerge$wMunicipality <- paste(EmployeeMerge$WorkMunicipality, EmployeeMerge$wMunicipalityKOD, sep="")
  EmployeeMerge <- subset.data.frame(EmployeeMerge, select = -c(wDeso,wMunicipalityKOD,MunicipalityKOD))
  
  # Calculate the Euclidean distance between the Workplaces, their employees and the 7 stations in Malmö Municipality
  EuclideanDist <- matrix(ncol=10, nrow=nrow(EmployeeMerge))
  for (i in 1: nrow(EmployeeMerge)){
    WorkplaceCoords <- c(EmployeeMerge[i,"wXKOORD"], EmployeeMerge[i,"wYKOORD"])
    EmployeeCoords <- c(EmployeeMerge[i,"eXKOORD"], EmployeeMerge[i, "eYKOORD"])
    Centralen <- c(StationCoords[8,3], StationCoords[8,4])
    Triangeln <- c(StationCoords[7,3], StationCoords[7,4])
    Hyllie <- c(StationCoords[6,3], StationCoords[6,4])
    Svagertorp <- c(StationCoords[5,3], StationCoords[5,4])
    Persborg <- c(StationCoords[4,3], StationCoords[4,4])
    Rosengard <- c(StationCoords[3,3], StationCoords[3,4])
    Ostra <- c(StationCoords[2,3], StationCoords[2,4])
    Oxie <- c(StationCoords[1,3], StationCoords[1,4])
    
    EuclideanDist[i,] <- c(EmployeeMerge$Resident_Nr[i], pointDistance(WorkplaceCoords, EmployeeCoords, lonlat=FALSE), pointDistance(WorkplaceCoords, Centralen, lonlat=FALSE), pointDistance(WorkplaceCoords, Triangeln, lonlat=FALSE), pointDistance(WorkplaceCoords, Hyllie, lonlat=FALSE), pointDistance(WorkplaceCoords, Svagertorp, lonlat=FALSE), pointDistance(WorkplaceCoords, Rosengard, lonlat=FALSE), pointDistance(WorkplaceCoords, Persborg, lonlat=FALSE), pointDistance(WorkplaceCoords, Ostra, lonlat=FALSE), pointDistance(WorkplaceCoords, Oxie, lonlat=FALSE))
  }
  
  # Merge Euclidean distance to main file
  colnames(EuclideanDist) <- c("Resident_Nr", "EmployeeDistance", "MalmoC", "Triangeln", "Hyllie", "Svagertorp", "Rosengard", "Persborg", "Ostra", "Oxie")
  DistMerge <- merge.data.frame(EuclideanDist, EmployeeMerge, by="Resident_Nr")
  
  # Set ID's as Character type (String)
  DistMerge$Resident_Nr <- as.character(DistMerge$Resident_Nr)
  DistMerge$Sp_Organization_Nr <- as.character(DistMerge$Sp_Organization_Nr)
  DistMerge$Id_Nr <- as.character(DistMerge$Id_Nr)
  DistMerge$Organization_Nr <- as.character(DistMerge$Organization_Nr)
  
  # Make sure certain fields are integer/numeric
  DistMerge$Gender_Wt04 <- as.numeric(DistMerge$Gender_Wt04)
  DistMerge$Social_Id <- as.integer(DistMerge$Social_Id)
  DistMerge$Income <- as.integer(DistMerge$Income)
  DistMerge$Disposable_Income_per_Id04 <- as.integer(DistMerge$Disposable_Income_per_Id04)
  DistMerge$Disposable_Income_K <- as.integer(DistMerge$Disposable_Income_K)
  DistMerge$Disposable_Income_K04 <- as.integer(DistMerge$Disposable_Income_K04)
  DistMerge$Disposable_Income04 <- as.integer(DistMerge$Disposable_Income04)
  DistMerge$Capital_Income <- as.integer(DistMerge$Capital_Income)
  
  # Save to CSV file
  write.csv(DistMerge, file = paste("Malmo_Workplace_Distances",year,".csv", sep=""), row.names = FALSE);
  
  #Save to DTA STATA file
  write.dta(DistMerge, file = paste("Malmo_Workplace_Distances",year,".dta", sep=""));
}
# Close the Database connection
close(conn)