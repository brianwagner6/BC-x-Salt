# Project: Biochar x Salt Greenhouse Project (2019)
# Location: Morton Arboretum Soils Lab, Lisle, Illinois
# Adviser & Corresponding Author: Dr. Meghan Midgley (mmidgley@mortonarb.org)
# By: Brian Wagner (brianandrewwagner@gmail.com)
# Created: September 18, 2022
# Updated: October 2, 2023
# Pre-print published August 5, 2023
# DOI: 10.1101/2023.08.03.551785

rm(list=ls(all=TRUE)) #Clears the Environment tab so all variables are unique to this Script

##### Load Packages and Files ####
{
  library(lmerTest)
  library(lsmeans)
  library(ggplot2)
  library(emmeans) 
  library(multcompView)
  library(tidyr)
  library(dplyr)
  library(pbkrtest)
  library(cowplot)
}

#### Import Master Data ####
{
  Master.df<-read.csv("C:\\Users\\mmidgley\\Downloads\\Master Data.csv", header = TRUE)
  # Make Species Block column that crosses the species and greenhouse groups (for t-tests)
  # (e.g., HON_3 represents one of the six HON in Greenhouse Group 3)
  Master.df$Species.Block <- paste(Master.df$Species, "_", Master.df$Greenhouse.Group)
  
  # Convert high sodium values (>200) to N/A becuase they are not correct 
  Master.df$Week.2.Leachate.NaCl..mg.[Master.df$Week.2.Leachate.NaCl..mg. > 200] <- NA
  
  # Convert high Above/Below Ground Biomass Ratio values above 10 becuase they are extraneous
  Master.df$New.Above.Below.Ground.Biomass.Ratio[Master.df$New.Above.Below.Ground.Biomass.Ratio > 10] <- NA
  
  # Create Roots and Shoots Biomass (g) column in Master.df
  # equal to the combined dry biomass of dry Fine Roots and Above Ground Biomass
  Master.df$Roots.and.Shoots.Biomass..g. <- Master.df$Fine.Root.Mass..g. + Master.df$New.Above.Ground.Biomass..g.
}

##### Create Average.SPAD column in Master.df ####
{
  Master.df$Average.SPAD <- rowMeans(Master.df[, c("Week.6.SPAD", "End.SPAD")], na.rm = TRUE)
  Master.df$Average.SPAD[is.nan(Master.df$Average.SPAD)] <- NA #Converts NaN (not a number) into NA
}

##### Adding LICOR (A.Net, g.sw) data to Master.df ####
{
  #Sugar Maple LICOR Data
  
  # Comments from Excel Workbook
  # Most leaves on Trees 120, 110, 109, 124, 85, 104, 93, 95, 118, 88 were 50% brown/crunchy
  # All leaves on Trees 100, 87, 115, 90, 106, 125, 97, 114,  were 50% brown/crunchy
  # On tree 119, all top leaves were dead. The leaves near the base were small but alive. The sample leave broke off during the test
  # Tree 122 had only one green leaf 
  # Trees 123, 94, 112, 103 were dead
  # Tree 117 was mostly dead
  
  # Import LICOR Data File
  LICOR_MAP.df <- read.csv("C:\\Users\\mmidgley\\Downloads\\LICOR MAP.csv", header = FALSE)
  
  #Extract the rows with data and labelled columns of the dataframe and override original variable 
  LICOR_MAP.df <- LICOR_MAP.df[15:nrow(LICOR_MAP.df),]
  
  #Dataframe with only data-containing rows and only columns 6 and 9 (for Tree ID and net assimilation rate [A_net = A], respectively)
  MAP_A_net.df <- LICOR_MAP.df[3:nrow(LICOR_MAP.df),c(6,9,14)]
  
  #Adds column titles
  colnames(MAP_A_net.df) <- c("TreeID", "A_Net", "g_sw")
  
  #Initiates the all_LICOR.df where all the LICOR data will be appended (using rbind())
  all_LICOR.df <- MAP_A_net.df
  
  #Honey Locust LICOR Data
  
  # Comments from Excel Workbook
  #All trees (Except ID 140) in this file need to be adjusted for leaf surface area (see pictures on Google Drive)
  
  # Import LICOR Data File
  LICOR_HON.df <- read.csv("C:\\Users\\mmidgley\\Downloads\\LICOR HON.csv", header = FALSE)
  
  #Extract the rows with data and labelled columns of the dataframe and override original variable 
  LICOR_HON.df <- LICOR_HON.df[14:nrow(LICOR_HON.df),] 
  
  #Dataframe with only data-containing rows and only columns 6 and 9 (for Tree ID and net assimilation rate [A_net = A], respectively)
  HON_A_net.df <- LICOR_HON.df[3:nrow(LICOR_HON.df),c(6,9,14)]
  
  #Adds column titles
  colnames(HON_A_net.df) <- c("TreeID", "A_Net", "g_sw")
  
  #Append data to all_LICOR.df
  all_LICOR.df <- rbind(all_LICOR.df,HON_A_net.df)
  
  #Red Oak LICOR Data
  
  # Comments from Excel Workbook
  # Trees 160, 137, 132, were all dead
  # Tree 153 had only 2 small leaves
  # Tree 163 had only 3 small leaves
  # Tree 154 had only 2 small leaves
  # Tree 156 had only 3 small leaves
  # The leaves on Trees 146 and 140 were less than 2 square cm
  # Trees 140 and 146 are in the file with the Honey Locusts 
  
  # Import LICOR Data File
  LICOR_OAK.df <- read.csv("C:\\Users\\mmidgley\\Downloads\\LICOR OAK.csv", header = FALSE)
  
  #Extract the rows with data and labelled columns of the dataframe and override original variable 
  LICOR_OAK.df <- LICOR_OAK.df[14:nrow(LICOR_OAK.df),] 
  
  #Dataframe with only data-containing rows and only columns 6 and 9 (for Tree ID and net assimilation rate [A_net = A], respectively)
  OAK_A_net.df <- LICOR_OAK.df[3:nrow(LICOR_OAK.df),c(6,9,14)]
  
  #Adds column titles
  colnames(OAK_A_net.df) <- c("TreeID", "A_Net", "g_sw")
  
  #Append data to all_LICOR.df
  all_LICOR.df <- rbind(all_LICOR.df,OAK_A_net.df)
  
  #Logdata 1 LICOR Data
  
  # Comments from Excel Workbook
  #Trees 12 and 37's leaves were curled and not the greatest
  
  # Import LICOR Data File
  LICOR_Data_1.df <- read.csv("C:\\Users\\mmidgley\\Downloads\\LICOR Logdata 1.csv", header = FALSE)
  
  #Extract the rows with data and labelled columns of the dataframe and override original variable 
  LICOR_Data_1.df <- LICOR_Data_1.df[14:nrow(LICOR_Data_1.df),] 
  
  #Dataframe with only data-containing rows and only columns 6 and 9 (for Tree ID and net assimilation rate [A_net = A], respectively)
  Data_1_A_net.df <- LICOR_Data_1.df[3:nrow(LICOR_Data_1.df),c(6,9,14)]
  
  #Adds column titles
  colnames(Data_1_A_net.df) <- c("TreeID", "A_Net", "g_sw")
  
  #Append data to all_LICOR.df
  all_LICOR.df <- rbind(all_LICOR.df,Data_1_A_net.df)
  
  #Logdata 2 LICOR Data
  
  # Comments from Excel Workbook
  #Tree 36 had bad data, will redo w/ different leaf 
  #The method that I used to merge the data in R overrides any existing data, so the good datapoint will be the one 
  #input into the data frame because the good datapoint is the second datapoint (since it was the redone value)
  
  # Import LICOR Data File
  LICOR_Data_2.df <- read.csv("C:\\Users\\mmidgley\\Downloads\\LICOR Logdata 2.csv", header = FALSE)
  
  #Extract the rows with data and labelled columns of the dataframe and override original variable 
  LICOR_Data_2.df <- LICOR_Data_2.df[14:nrow(LICOR_Data_2.df),] 
  
  #Dataframe with only data-containing rows and only columns 6 and 9 (for Tree ID and net assimilation rate [A_net = A], respectively)
  Data_2_A_net.df <- LICOR_Data_2.df[3:nrow(LICOR_Data_2.df),c(6,9,14)]
  
  #Adds column titles
  colnames(Data_2_A_net.df) <- c("TreeID", "A_Net", "g_sw")
  
  #Append data to all_LICOR.df
  all_LICOR.df <- rbind(all_LICOR.df,Data_2_A_net.df)
  
  #Merge LICOR data into Master.df and adjust it for Fraction of Licor covered 
  
  #Converts TreeID column into int rather than character
  all_LICOR.df$TreeID <- as.numeric(as.character(all_LICOR.df$TreeID))
  
  #Sorts numerically by TreeID
  all_LICOR.df <- arrange(all_LICOR.df,TreeID)
  
  #Add data to Master.df
  Master.df <- merge(Master.df, all_LICOR.df, by.x = "ID", by.y = "TreeID", all.x = TRUE)
  
  #Adjust A_Net based on the Fraction of Licor covered by Leaf
  
  #Create new column in dataframe
  Master.df$Adjusted.A.Net <- NA
  
  #Convert A_Net column to numeric 
  Master.df$A_Net <- as.numeric(as.character(Master.df$A_Net))
  
  #Convert g_sw column to numeric 
  Master.df$g_sw <- as.numeric(as.character(Master.df$g_sw))
  
  #Multiply the A_Net and Fraction columns, storing the products in the new Adjusted.A.Net column
  Master.df$Adjusted.A.Net <- Master.df$A_Net / Master.df$Fraction.of.LICOR.covered.by.Leaf
  
  # Convert high Adjusted.A.Net values to N/A becuase they are not correct
  # This only changes Adjusted.A.Net, not the original A.Net, but only Adjusted.A.Net values
  # are used for data analysis 
  Master.df$Adjusted.A.Net[Master.df$Adjusted.A.Net > 20] <- NA
  
  #Calculate WUE
  Master.df$WUE <- Master.df$A_Net / Master.df$g_sw
  
  #Create columns for squared and cubed Adusted.A.Net for Model 1 Statistics Table
  Master.df$Adjusted.A.Net.Squared <- (Master.df$Adjusted.A.Net)^2
  Master.df$Adjusted.A.Net.Cubed <- (Master.df$Adjusted.A.Net)^3
}

##### Create leachate.df ####
{
  # leachate.df is used for graphing leachate vs. time, such as in All Data Model 3 
  
  #Week Vector
  week.vector <- (rep(1,168)) #initiates the week vector to have one-hundred-sixty-eight 1's
  for(i in 2:8){ #completes the week vector so that each week is represented 168 times in a row
    week.vector <- c(week.vector,rep(i,168))
  }
  #Greenhouse Group vector 
  greenhouse.vector <- as.vector(Master.df$Greenhouse.Group)
  greenhouse.vector <- c(rep(greenhouse.vector,8))
  # Salt.Addition vector
  salt_addition.vector <- as.vector(Master.df$Salt.Addition)
  salt_addition.vector <- c(rep(salt_addition.vector,8))
  # Species vector
  species.vector <- as.vector(Master.df$Species)
  species.vector <- c(rep(species.vector,8))
  # Soil.Type vector
  soil_type.vector <- as.vector(Master.df$Soil.Type)
  soil_type.vector <- c(rep(soil_type.vector,8))
  #Leachate Vector (contains all the leachate data in ID order within a week, and all weeks are concatenated)
  NaCl_leachate.vector <- as.vector(as.matrix(Master.df[,c("Week.1.Leachate.NaCl..mg.", "Week.2.Leachate.NaCl..mg.", 
                                                           "Week.3.Leachate.NaCl..mg.", "Week.4.Leachate.NaCl..mg.", 
                                                           "Week.5.Leachate.NaCl..mg.", "Week.6.Leachate.NaCl..mg.",
                                                           "Week.7.Leachate.NaCl..mg.", "Week.8.Leachate.NaCl..mg.")]))
  # The following vector has each ID repeated 8 times (once for each week)
  leachate.df <- data.frame("ID" = rep(1:168,8), "Species" = species.vector, "Soil Type" = soil_type.vector,
                            "Salt Addition" = salt_addition.vector, "Greenhouse Group" = greenhouse.vector, 
                            "Week" = week.vector, "NaCl Leachate (mg)" = NaCl_leachate.vector)  
  #Converts Week column into character rather than integer
  leachate.df$Week <- as.character(as.numeric(leachate.df$Week))
  
  # Convert incorrect data point to N/A
  leachate.df$NaCl.Leachate..mg.[leachate.df$Week == 2 & leachate.df$NaCl.Leachate..mg. > 200] <- NA
}

#### Analyze functions ####
{
  # model_1_salt_addition_x_soil.function is a generic function that takes input for the y axis only because the dataframe is 
  # always master.df, etc. It can be used for any of the overall analyses in this file
  # This does not separate species but rather combines all trees into one "data pool"
  # species_Y_N is a string-form boolean for which Y will run the lmer incorporating "+ data.df$Species" to account for 
  # species differences
  # Example calling of the model_1_salt_addition_x_soil.function() : 
  # model_1_salt_addition_x_soil.function(Master.df, Master.df$New.Above.Ground.Biomass..g., "New.Above.Ground.Biomass..g.", "Y")
  model_1_salt_addition_x_soil.function <- function(data.df, y_axis, string_y, species_Y_N) {
    # Y-axis data vs. Species Graph
    print(ggplot(data.df, aes(Soil.Type, {{y_axis}}, fill = Salt.Addition, na.rm = TRUE)) + geom_boxplot() + 
            stat_boxplot(geom = "errorbar") + ggtitle(paste(string_y)))
    
    # lmer for Y-axis Data vs. Species 
    if(species_Y_N == "Y"){
      model_1_salt_addition_x_soil.lmer <- lmer(data.df[,string_y] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                                  data.df$Species + (1|Greenhouse.Group), na.action = na.exclude, 
                                                data = data.df)
    }
    else{
      model_1_salt_addition_x_soil.lmer <- lmer(data.df[,string_y] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                                  (1|Greenhouse.Group), na.action = na.exclude, data = data.df)
    }
    
    # ANOVA
    print(anova(model_1_salt_addition_x_soil.lmer)) 
    
    #Ask if want to run post-hoc test
    input <- readline("Run POST HOC test? Enter 'Y' or 'N': ")
    if(input == "Y") {
      input_2 <- readline("Which variable for the post hoc test? Enter 1, 2, or 3 (1: 'Soil.Type', 2: 'Species', or 3: 'Salt.Addition * Soil.Type'): ")
      if(input_2 == "1"){ #1: Soil.Type
        emmeans(model_1_salt_addition_x_soil.lmer, list(pairwise ~ Soil.Type), adjust = "tukey")
      } else if(input_2 == "2") {#2: Species
        emmeans(model_1_salt_addition_x_soil.lmer, list(pairwise ~ Species), adjust = "tukey")
      } else if(input_2 == "3") {#3: Soil.Type * Salt.Addition
        emmeans(model_1_salt_addition_x_soil.lmer, list(pairwise ~ Soil.Type * Salt.Addition), adjust = "tukey")
      } else {
        print("Error")
      }
    }
    else {
      print("No POST HOC test was run")
    }
  }
}

##### Create species-specific data frames ####
{
  CAT.df <- subset(Master.df, Species == "CAT")
  HON.df <- subset(Master.df, Species == "HON")
  MAP.df <- subset(Master.df, Species == "MAP")
  OAK.df <- subset(Master.df, Species == "OAK")
}

#### Add columns in data frames for testing for any necessary data transformations ####
{
  # Master data transformations
  {
    # Master Diameter Growth
    {
      # Inverse
      Master.df$Diameter.Growth.Inverse <- (Master.df$Diameter.Growth....)^(-1)
      
      # Squared
      Master.df$Diameter.Growth.Squared <- (Master.df$Diameter.Growth....)^2
      
      # Square Root
      # Smallest original value is -23.38251986
      # Add 23.5 to every value so that there are no negative values before square rooting
      Master.df$Diameter.Growth.Sqrt <- sqrt((Master.df$Diameter.Growth....)+23.5)
      
      # Log
      # Smallest original value is -23.38251986
      # Add 23.5 to every value so that there are no negative values before taking log 
      Master.df$Diameter.Growth.Log <- log((Master.df$Diameter.Growth....)+23.5)
    }
    
    # Master Stem Mass Growth
    {
      # Inverse
      Master.df$Stem.Mass.Growth.Inverse <- (Master.df$Stem.Mass.Growth....)^(-1)
      #Replace "Inf" with NA
      Master.df$Stem.Mass.Growth.Inverse[Master.df$Stem.Mass.Growth.Inverse == "Inf"] <- NA
      
      # Squared
      Master.df$Stem.Mass.Growth.Squared <- (Master.df$Stem.Mass.Growth....)^2
      
      # Square Root
      # Smallest original value is -60.9375
      # Add 61 to every value so that there are no negative values before square rooting
      Master.df$Stem.Mass.Growth.Sqrt <- sqrt((Master.df$Stem.Mass.Growth....)+61)
      
      # Log
      # Smallest original value is -60.9375
      # Add 61 to every value so that there are no negative values before taking the log
      Master.df$Stem.Mass.Growth.Log <- log((Master.df$Stem.Mass.Growth....)+61)
    }
    
    # Master Coarse Root Growth
    {
      # Inverse
      Master.df$Coarse.Root.Growth.Inverse <- (Master.df$Coarse.Root.Growth....)^(-1)
      #Replace "Inf" with NA
      Master.df$Coarse.Root.Growth.Inverse[Master.df$Coarse.Root.Growth.Inverse == "Inf"] <- NA
      
      # Squared
      Master.df$Coarse.Root.Growth.Squared <- (Master.df$Coarse.Root.Growth....)^2
      
      # Square Root
      # Smallest original value is -92.06349206
      # Add 92.5 to every value so that there are no negative values before square rooting
      Master.df$Coarse.Root.Growth.Sqrt <- sqrt((Master.df$Coarse.Root.Growth....)+92.5)
      
      # Log
      # Smallest original value is -92.06349206
      # Add 92.5 to every value so that there are no negative values before taking the log
      Master.df$Coarse.Root.Growth.Log <- log((Master.df$Coarse.Root.Growth....)+92.5)
    }
    
    # Master Pre Existing Biomass Growth
    {
      # Inverse
      Master.df$Pre.Existing.Biomass.Growth.Inverse <- (Master.df$Pre.Existing.Biomass.Growth....)^(-1)
      
      # Squared
      Master.df$Pre.Existing.Biomass.Growth.Squared <- (Master.df$Pre.Existing.Biomass.Growth....)^2
      
      # Square Root
      # Smallest original value is -69.48118006
      # Add 70 to every value so that there are no negative values before square rooting
      Master.df$Pre.Existing.Biomass.Growth.Sqrt <- sqrt((Master.df$Pre.Existing.Biomass.Growth....)+70)
      
      # Log
      # Smallest original value is -69.48118006
      # Add 70 to every value so that there are no negative values before taking the log
      Master.df$Pre.Existing.Biomass.Growth.Log <- log((Master.df$Pre.Existing.Biomass.Growth....)+70)
    }
    
    # Master New Above Ground Biomass
    {
      # Inverse 
      Master.df$New.Above.Ground.Biomass..g.Inverse <- (Master.df$New.Above.Ground.Biomass..g.)^(-1)
      #Replace "Inf" with NA
      Master.df$New.Above.Ground.Biomass..g.Inverse[Master.df$New.Above.Ground.Biomass..g.Inverse == "Inf"] <- NA
      
      # Squared
      Master.df$New.Above.Ground.Biomass..g.Squared <- (Master.df$New.Above.Ground.Biomass..g.)^2
      
      # Square Root
      Master.df$New.Above.Ground.Biomass..g.sqrt <- sqrt(Master.df$New.Above.Ground.Biomass..g.)
      
      # Log 
      # Smallest original value is 0
      # Add 0.5 to every value so that there are no negative values before taking the log
      Master.df$New.Above.Ground.Biomass..g.Log <- log((Master.df$New.Above.Ground.Biomass..g.)+0.5)
    }
    
    # Master Fine Root Mass
    {
      # Inverse 
      Master.df$Fine.Root.Mass..g.Inverse<- (Master.df$Fine.Root.Mass..g.)^(-1)
      #Replace "Inf" with NA
      Master.df$Fine.Root.Mass..g.Inverse[Master.df$Fine.Root.Mass..g.Inverse == "Inf"] <- NA
      
      # Squared 
      Master.df$Fine.Root.Mass..g.Squared <- (Master.df$Fine.Root.Mass..g.)^2
      
      # Square Root 
      Master.df$Fine.Root.Mass..g.Sqrt <- sqrt(Master.df$Fine.Root.Mass..g.)
      
      # Log
      # Smallest original value is 0
      # Add 0.5 to every value so that there are no negative values before square rooting
      Master.df$Fine.Root.Mass..g.Log <- log((Master.df$Fine.Root.Mass..g.)+0.5)
    }
  }
  
  # CAT data transformations
  {
    # CAT Diameter Growth
    {
      # Inverse
      CAT.df$Diameter.Growth.Inverse <- (CAT.df$Diameter.Growth....)^(-1)
      
      # Squared
      CAT.df$Diameter.Growth.Squared <- (CAT.df$Diameter.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -12.99734748
      # Add 13 to every value so that there are no negative values before square rooting
      CAT.df$Diameter.Growth.Sqrt <- sqrt((CAT.df$Diameter.Growth....)+13)
      
      # Log
      # Smallest original value is -12.99734748
      # Add 13 to every value so that there are no negative values before taking the log
      CAT.df$Diameter.Growth.Log <- log((CAT.df$Diameter.Growth....)+13)
    }
    
    # CAT Stem Mass Growth 
    {
      # No transformation needed
    }
    
    # CAT Coarse Root Growth
    {
      # Inverse
      CAT.df$Coarse.Root.Growth.Inverse <- (CAT.df$Coarse.Root.Growth....)^(-1)
      
      # Squared
      CAT.df$Coarse.Root.Growth.Squared <- (CAT.df$Coarse.Root.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -78.57142857
      # Add 79 to every value so that there are no negative values before square rooting
      CAT.df$Coarse.Root.Growth.Sqrt <- sqrt((CAT.df$Coarse.Root.Growth....)+79)
      
      # Log
      # Smallest original value is -78.57142857
      # Add 79 to every value so that there are no negative values before taking the log 
      CAT.df$Coarse.Root.Growth.Log <- log((CAT.df$Coarse.Root.Growth....)+79)
    }
    
    # CAT Pre-Existing Biomass Growth 
    {
      # Inverse
      CAT.df$Pre.Existing.Biomass.Growth.Inverse <- (CAT.df$Pre.Existing.Biomass.Growth....)^(-1)
      
      # Squared
      CAT.df$Pre.Existing.Biomass.Growth.Squared <- (CAT.df$Pre.Existing.Biomass.Growth....)^(2)
      
      # Square Root
      CAT.df$Pre.Existing.Biomass.Growth.Sqrt <- sqrt(CAT.df$Pre.Existing.Biomass.Growth....)
      
      # Log
      CAT.df$Pre.Existing.Biomass.Growth.Log <- log(CAT.df$Pre.Existing.Biomass.Growth....)
    }
    
    # CAT New Above Ground Biomass
    {
      # Inverse 
      CAT.df$New.Above.Ground.Biomass..g.Inverse <- (CAT.df$New.Above.Ground.Biomass..g.)^(-1)
      
      # Squared 
      CAT.df$New.Above.Ground.Biomass..g.Squared <- (CAT.df$New.Above.Ground.Biomass..g.)^(2)
      
      # Square Root
      CAT.df$New.Above.Ground.Biomass..g.Sqrt <- sqrt(CAT.df$New.Above.Ground.Biomass..g.)
      
      # Log
      CAT.df$New.Above.Ground.Biomass..g.Log <- log(CAT.df$New.Above.Ground.Biomass..g.)
    }
    
    # CAT Fine Root Mass (g)
    {
      # Inverse 
      CAT.df$Fine.Root.Mass..g.Inverse <- (CAT.df$Fine.Root.Mass..g.)^(-1)
      
      # Squared 
      CAT.df$Fine.Root.Mass..g.Squared <- (CAT.df$Fine.Root.Mass..g.)^(2)
      
      # Square Root
      CAT.df$Fine.Root.Mass..g.Sqrt <- sqrt(CAT.df$Fine.Root.Mass..g.)
      
      # Log
      CAT.df$Fine.Root.Mass..g.Log <- log(CAT.df$Fine.Root.Mass..g.)
    }
    
    # CAT Total Biomass Growth
    {
      # Inverse
      CAT.df$Total.Biomass.Growth.Inverse <- (CAT.df$Total.Biomass.Growth....)^(-1)
      
      # Squared
      CAT.df$Total.Biomass.Growth.Squared <- (CAT.df$Total.Biomass.Growth....)^2
      
      # Square Root
      CAT.df$Total.Biomass.Growth.Sqrt <- sqrt((CAT.df$Total.Biomass.Growth....))
      
      # Log
      CAT.df$Total.Biomass.Growth.Log <- log((CAT.df$Total.Biomass.Growth....))
    }
    
    # CAT Roots and Shoots Biomass (g)
    {
      # Inverse 
      CAT.df$Roots.and.Shoots.Biomass..g.Inverse <- (CAT.df$Roots.and.Shoots.Biomass..g.)^(-1)
      
      # Squared 
      CAT.df$Roots.and.Shoots.Biomass..g.Squared <- (CAT.df$Roots.and.Shoots.Biomass..g.)^(2)
      
      # Square Root
      CAT.df$Roots.and.Shoots.Biomass..g.Sqrt <- sqrt(CAT.df$Roots.and.Shoots.Biomass..g.)
      
      # Log
      CAT.df$Roots.and.Shoots.Biomass..g.Log <- log(CAT.df$Roots.and.Shoots.Biomass..g.)
    }
  }
  
  # MAP data transformations
  {
    # MAP Diameter Growth
    {
      # Inverse
      MAP.df$Diameter.Growth.Inverse <- (MAP.df$Diameter.Growth....)^(-1)
      
      # Squared
      MAP.df$Diameter.Growth.Squared <- (MAP.df$Diameter.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -13.46389229
      # Add 13.5 to every value so that there are no negative values before square rooting
      MAP.df$Diameter.Growth.Sqrt <- sqrt((MAP.df$Diameter.Growth....)+13.5)
      
      # Log 
      # Smallest original value is -13.46389229
      # Add 13.5 to every value so that there are no negative values before taking the log
      MAP.df$Diameter.Growth.Log <- log((MAP.df$Diameter.Growth....)+13.5)
    }
    
    # MAP Stem Mass Growth
    {
      # No transformation needed
    }
    
    # MAP Coarse Root Growth
    {
      # No transformation needed
    }
    
    # MAP Pre-Existing Biomass Growth
    {
      # No transformation needed
    }
    
    # MAP New Above Ground (g)
    {
      # No transformation needed
    }
    
    # MAP Fine Root (g)
    {
      # No transformation needed
    }
    
    # MAP Total Biomass Growth 
    {
      # No transformation needed
    }
  }
  
  # HON data transformations
  {
    # HON Diameter Growth 
    {
      # Inverse
      HON.df$Diameter.Growth.Inverse <- (HON.df$Diameter.Growth....)^(-1)
      
      # Squared
      HON.df$Diameter.Growth.Squared <- (HON.df$Diameter.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -9.143807149
      # Add 9.5 to every value so that there are no negative values before square rooting
      HON.df$Diameter.Growth.Sqrt <- sqrt((HON.df$Diameter.Growth....)+9.5)
      
      # Log
      # Smallest original value is -9.143807149
      # Add 9.5 to every value so that there are no negative values before taking the log
      HON.df$Diameter.Growth.Log <- log((HON.df$Diameter.Growth....)+9.5)
    }
    
    # HON Stem Mass Growth
    {
      # No data
    }
    
    # HON Coarse Root Growth
    {
      # Inverse
      HON.df$Coarse.Root.Growth.Inverse <- (HON.df$Coarse.Root.Growth....)^(-1)
      
      # Squared
      HON.df$Coarse.Root.Growth.Squared <- (HON.df$Coarse.Root.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -44.44444444
      # Add 44.5 to every value so that there are no negative values before square rooting
      HON.df$Coarse.Root.Growth.Sqrt <- sqrt((HON.df$Coarse.Root.Growth....)+44.5)
      
      # Log
      # Smallest original value is -44.44444444
      # Add 44.5 to every value so that there are no negative values before taking the log
      HON.df$Coarse.Root.Growth.Log <- log((HON.df$Coarse.Root.Growth....)+44.5)
    }
    
    # HON Pre-Existing Biomass Growth
    {
      # Inverse
      HON.df$Pre.Existing.Biomass.Growth.Inverse <- (HON.df$Pre.Existing.Biomass.Growth....)^(-1)
      
      # Squared
      HON.df$Pre.Existing.Biomass.Growth.Squared <- (HON.df$Pre.Existing.Biomass.Growth....)^(2)
      
      # Smallest original value is -38.96484375
      # Add 39 to every value so that there are no negative values before square rooting
      HON.df$Pre.Existing.Biomass.Growth.Sqrt <- sqrt((HON.df$Pre.Existing.Biomass.Growth....)+39)
      
      # Log
      # Smallest original value is -38.96484375
      # Add 39 to every value so that there are no negative values before taking the log
      HON.df$Pre.Existing.Biomass.Growth.Log <- log((HON.df$Pre.Existing.Biomass.Growth....)+39)
    }
    
    # HON New Above Ground Biomass (g)
    {
      # Squared
      HON.df$New.Above.Ground.Biomass..g.Squared <- (HON.df$New.Above.Ground.Biomass..g.)^(2)
      
      # Square Root
      HON.df$New.Above.Ground.Biomass..g.Sqrt <- sqrt(HON.df$New.Above.Ground.Biomass..g.)
      
      # Log
      # Smallest original value is 0
      # Add 0.5 to every value so that there are no zero values before taking the log
      HON.df$New.Above.Ground.Biomass..g.Log <- log((HON.df$New.Above.Ground.Biomass..g.)+0.5)
    }
    
    # HON Fine Root (g)
    {
      # Squared
      HON.df$Fine.Root.Mass..g.Squared <- (HON.df$Fine.Root.Mass..g.)^(2)
      
      # Square Root
      HON.df$Fine.Root.Mass..g.Sqrt <- sqrt(HON.df$Fine.Root.Mass..g.)
      
      # Log
      # Smallest original value is 0
      # Add 0.5 to every value so that there are no zero values before taking the log
      HON.df$Fine.Root.Mass..g.Log <- log((HON.df$Fine.Root.Mass..g.)+0.5)
    }
    
    # HON Total Biomass Growth 
    {
      # Inverse
      HON.df$Total.Biomass.Growth.Inverse <- (HON.df$Total.Biomass.Growth....)^(-1)
      
      # Squared
      HON.df$Total.Biomass.Growth.Squared <- (HON.df$Total.Biomass.Growth....)^2
      
      # Square Root
      # Smallest original value is -38.59667969
      # Add 39 to every value so that there are no negative values before square rooting
      HON.df$Total.Biomass.Growth.Sqrt <- sqrt((HON.df$Total.Biomass.Growth....)+39)
      
      # Log
      # Smallest original value is -38.59667969
      # Add 39 to every value so that there are no negative values before taking the log
      HON.df$Total.Biomass.Growth.Log <- log((HON.df$Total.Biomass.Growth....)+39)
    }
  }
  
  # OAK data transformations
  {
    # OAK Diameter Growth 
    {
      # Inverse
      OAK.df$Diameter.Growth.Inverse <- (OAK.df$Diameter.Growth....)^(-1)
      
      # Squared
      OAK.df$Diameter.Growth.Squared <- (OAK.df$Diameter.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -23.38251986
      # Add 23.5 to every value so that there are no negative values before square rooting
      OAK.df$Diameter.Growth.Sqrt <- sqrt((OAK.df$Diameter.Growth....)+23.5)
      
      # Log
      # Smallest original value is -23.38251986
      # Add 23.5 to every value so that there are no negative values before taking the log
      OAK.df$Diameter.Growth.Log <- log((OAK.df$Diameter.Growth....)+23.5)
    }
    
    # OAK Stem Mass Growth
    {
      # Inverse
      OAK.df$Stem.Mass.Growth.Inverse <- (OAK.df$Stem.Mass.Growth....)^(-1)
      # Remove 1/0 (undefined) values and convert to NA
      OAK.df$Stem.Mass.Growth.Inverse[OAK.df$Stem.Mass.Growth.Inverse == "Inf"] <- NA
      
      # Squared
      OAK.df$Stem.Mass.Growth.Squared <- (OAK.df$Stem.Mass.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -37.5
      # Add 38 to every value so that there are no negative values before square rooting
      OAK.df$Stem.Mass.Growth.Sqrt <- sqrt((OAK.df$Stem.Mass.Growth....)+38)
      
      # Log
      # Smallest original value is -37.5
      # Add 38 to every value so that there are no negative values before square rooting
      OAK.df$Stem.Mass.Growth.Log <- log((OAK.df$Stem.Mass.Growth....)+38)
    }
    
    # OAK Coarse Root Growth
    {
      # No transformation needed
    }
    
    # OAK Pre-Existing Biomass (g)
    {
      # Inverse
      OAK.df$Pre.Existing.Biomass.Growth.Inverse <- (OAK.df$Pre.Existing.Biomass.Growth....)^(-1)
      
      # Squared
      OAK.df$Pre.Existing.Biomass.Growth.Squared <- (OAK.df$Pre.Existing.Biomass.Growth....)^(2)
      
      # Square Root
      # Smallest original value is -42.04425711
      # Add 42.05 to every value so that there are no negative values before square rooting
      OAK.df$Pre.Existing.Biomass.Growth.Sqrt <- sqrt((OAK.df$Pre.Existing.Biomass.Growth....)+42.5)
      
      # Log
      # Smallest original value is -42.04425711
      # Add 42.5 to every value so that there are no negative values before taking the log
      OAK.df$Pre.Existing.Biomass.Growth.Log <- log((OAK.df$Pre.Existing.Biomass.Growth....)+42.5)
    }
    
    # OAK Fine Root Biomass (g)
    {
      # Inverse 
      OAK.df$Fine.Root.Mass..g.Inverse <- (OAK.df$Fine.Root.Mass..g.)^(-1)
      OAK.df$Fine.Root.Mass..g.Inverse[OAK.df$Fine.Root.Mass..g.Inverse == "Inf"] <-- NA
      
      # Squared 
      OAK.df$Fine.Root.Mass..g.Squared <- (OAK.df$Fine.Root.Mass..g.)^(2)
      
      # Square Root
      OAK.df$Fine.Root.Mass..g.Sqrt <- sqrt(OAK.df$Fine.Root.Mass..g.)
      
      # Log
      # Smallest original value is 0
      # Add 0.5 to every value so that there are no negative values before taking the log
      OAK.df$Fine.Root.Mass..g.Log <- log((OAK.df$Fine.Root.Mass..g.)+0.5)
    }
    
    # OAK Total Biomass Growth
    {
      # Inverse
      OAK.df$Total.Biomass.Growth.Inverse <- (OAK.df$Total.Biomass.Growth....)^(-1)
      
      # Squared
      OAK.df$Total.Biomass.Growth.Squared <- (OAK.df$Total.Biomass.Growth....)^2
      
      # Square Root
      # Smallest original value is -38.32903056
      # Add 39 to every value so that there are no negative values before square rooting
      OAK.df$Total.Biomass.Growth.Sqrt <- sqrt((OAK.df$Total.Biomass.Growth....)+39)
      
      # Log
      # Smallest original value is-38.32903056
      # Add 39 to every value so that there are no negative values before taking the log
      OAK.df$Total.Biomass.Growth.Log <- log((OAK.df$Total.Biomass.Growth....)+39)
    }
  }
}

#### Testing Data for Normalization, Transformation as needed ####
{
  qqtest.function <- function(data)
  {
    hist(data, main = deparse(substitute(data))) # adds title as data's object name
    qqnorm(data)
    qqline(data)
    shapiro.test(data)
  }
  
  # Tests residuals for the Master.df data set
  qqtest_overall_residuals.function <- function(data.df, string_y_axis, final_or_not) 
  {
    if(final_or_not == "Y"){
      input <- readline("Input initial data as string, e.g., 'Final.Diameter..mm.': ")
      initial.data <- data.df[,input]
      qqtest_overall_residuals.lmer <- lmer(data.df[,string_y_axis] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                              data.df$Species + initial.data + (1|Greenhouse.Group), 
                                            na.action = na.exclude, data = data.df)
    }
  
    else{
      qqtest_overall_residuals.lmer <- lmer(data.df[,string_y_axis] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                              data.df$Species + (1|Greenhouse.Group), na.action = na.exclude, 
                                            data = data.df)
    }
    
    # Test for normality
    qqnorm(residuals(qqtest_overall_residuals.lmer))
    qqline(residuals(qqtest_overall_residuals.lmer))
    print(shapiro.test(residuals(qqtest_overall_residuals.lmer)))
    
    # Test for homoscedasticity 
    plot(residuals(qqtest_overall_residuals.lmer) ~ fitted(qqtest_overall_residuals.lmer))
  }
  
  # Tests residuals for an input species data set
  qqtest_species_residuals.function <- function(data.df, string_y_axis, final_or_not) 
  {
    if(final_or_not == "Y"){
      input <- readline("Input initial data as string, e.g., 'Final.Diameter..mm.' ")
      initial.data <- data.df[,input]
      qqtest_overall_residuals.lmer <- lmer(data.df[,string_y_axis] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                            initial.data + (1|Greenhouse.Group), 
                                            na.action = na.exclude, data = data.df)
    }
    else{
      qqtest_overall_residuals.lmer <- lmer(data.df[,string_y_axis] ~ data.df$Salt.Addition * data.df$Soil.Type + 
                                              (1|Greenhouse.Group), na.action = na.exclude, 
                                            data = data.df)
    }
    
    # Test for normality
    qqnorm(residuals(qqtest_overall_residuals.lmer))
    qqline(residuals(qqtest_overall_residuals.lmer))
    print(shapiro.test(residuals(qqtest_overall_residuals.lmer)))
    
    # Test for homoscedasticity 
    plot(residuals(qqtest_overall_residuals.lmer) ~ fitted(qqtest_overall_residuals.lmer))
  }
  
  # Overall Data Transformation QQTESTS
  {
    #Adjusted.A.Net --> Per Meghan's playing around with the data: use A.Net.Cubed 
    
    #Diameter Growth --> Use Diameter.Growth.Log
    {
      # Original Data
      hist(Master.df$Diameter.Growth....)
      qqnorm(Master.df$Diameter.Growth....)
      qqline(Master.df$Diameter.Growth....)
      shapiro.test(Master.df$Diameter.Growth....) # p = 2.391e-15
      
      # Inverse
      hist(Master.df$Diameter.Growth.Inverse)
      qqnorm(Master.df$Diameter.Growth.Inverse)
      qqline(Master.df$Diameter.Growth.Inverse)
      shapiro.test(Master.df$Diameter.Growth.Inverse) # p < 2.2e-16
      
      # Squared 
      hist(Master.df$Diameter.Growth.Squared)
      qqnorm(Master.df$Diameter.Growth.Squared)
      qqline(Master.df$Diameter.Growth.Squared)
      shapiro.test(Master.df$Diameter.Growth.Squared) # p < 2.2e-16
      
      # Square Root
      hist(Master.df$Diameter.Growth.Sqrt)
      qqnorm(Master.df$Diameter.Growth.Sqrt)
      qqline(Master.df$Diameter.Growth.Sqrt)
      shapiro.test(Master.df$Diameter.Growth.Sqrt) # p = 5.833e-10 
      
      # Log 
      hist(Master.df$Diameter.Growth.Log) # Looks good 
      qqnorm(Master.df$Diameter.Growth.Log)
      qqline(Master.df$Diameter.Growth.Log)
      shapiro.test(Master.df$Diameter.Growth.Log) # p = 2.167e-14
      # Testing residuals 
      qqtest_overall_residuals.function(Master.df, "Diameter.Growth.Log", "N") # p < 2.2e-16
      
      #Conclusion: take the log of Diameter.Growth....
    }
    
    #Stem Mass Growth --> Use Original Data
    {
      # Original Data 
      hist(Master.df$Stem.Mass.Growth....)
      qqnorm(Master.df$Stem.Mass.Growth....)
      qqline(Master.df$Stem.Mass.Growth....)
      shapiro.test(Master.df$Stem.Mass.Growth....) # = 1.752e-10
      
      # Inverse
      hist(Master.df$Stem.Mass.Growth.Inverse)
      qqnorm(Master.df$Stem.Mass.Growth.Inverse)
      qqline(Master.df$Stem.Mass.Growth.Inverse)
      shapiro.test(Master.df$Stem.Mass.Growth.Inverse) # p = 2.198e-12
      
      # Squared
      qqtest.function(Master.df$Stem.Mass.Growth.Squared) # Shapiro p < 2.2e-16
      
      # Square Root
      hist(Master.df$Stem.Mass.Growth.Sqrt)
      qqnorm(Master.df$Stem.Mass.Growth.Sqrt)
      qqline(Master.df$Stem.Mass.Growth.Sqrt)
      shapiro.test(Master.df$Stem.Mass.Growth.Sqrt) # p = 1.648e-05
      
      # Log 
      qqtest.function(Master.df$Stem.Mass.Growth.Log) # p = 1.473e-13
      
      #Conclusion: Use original data
    }
    
    #Coarse Root Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(Master.df$Coarse.Root.Growth....) # p = 8.745e-13
      
      # Inverse
      qqtest.function(Master.df$Coarse.Root.Growth.Inverse) # p < 2.2e-16
      
      # Squared
      qqtest.function(Master.df$Coarse.Root.Growth.Squared) # Shapiro p < 2.2e-16
      
      # Square Root 
      qqtest.function(Master.df$Coarse.Root.Growth.Sqrt) # p = 0.0006588
      
      # Log
      qqtest.function(Master.df$Coarse.Root.Growth.Log) # p = 9.135e-11
      
      #Conclusion: Use original for Coarse Root 
    }
    
    #Pre-Existing Biomass Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(Master.df$Pre.Existing.Biomass.Growth....) # p = 1.145e-14
  
      # Inverse
      qqtest.function(Master.df$Pre.Existing.Biomass.Growth.Inverse) # Shapiro p < 2.2e-16
      
      # Squared
      qqtest.function(Master.df$Pre.Existing.Biomass.Growth.Squared) # Shapiro p < 2.2e-16
      
      # Square Root
      qqtest.function(Master.df$Pre.Existing.Biomass.Growth.Sqrt) # p = 4.984e-07
      
      # Log
      qqtest.function(Master.df$Pre.Existing.Biomass.Growth.Log)  # Shapiro p = = 2.214e-07
      
      #Conclusion: Use Original Data for Pre-Existing Biomass
    }
    
    #New Above Ground Biomass --> Use Original Data
    {
      # Original Data 
      qqtest.function(Master.df$New.Above.Ground.Biomass..g.) # p = 3.161e-14
      
      # Inverse 
      qqtest.function(Master.df$New.Above.Ground.Biomass..g.Inverse) # Shapiro p < 2.2e-16
      
      # Squared
      qqtest.function(Master.df$New.Above.Ground.Biomass..g.Squared) # Shapiro p < 2.2e-16
      
      # Square Root
      qqtest.function(Master.df$New.Above.Ground.Biomass..g.sqrt) # Shapiro p = 8.643e-06
      
      # Log 
      qqtest.function(Master.df$New.Above.Ground.Biomass..g.Log) # Shapiro p = 6.142e-06
      
      #Conclusion: Use Original Data for New Above Ground Biomass
      
    }
    
    #Fine Root Mass --> Use Original Data
    {
      # Original Data
      qqtest.function(Master.df$Fine.Root.Mass..g.) # Shapiro = 1.634e-12
      
      # Inverse 
      qqtest.function(Master.df$Fine.Root.Mass..g.Inverse) # Shapiro p < 2.2e-16
      
      # Squared 
      qqtest.function(Master.df$Fine.Root.Mass..g.Squared) # Shapiro p < 2.2e-16
      
      # Square Root 
      qqtest.function(Master.df$Fine.Root.Mass..g.Sqrt) # Shapiro p = 2.202e-05
      
      # Log
      qqtest.function(Master.df$Fine.Root.Mass..g.Log) # Shapiro p = 1.093e-07
      
      # Conclusion: use Original Data for Fine.Root
      
    }
    
    # Total Biomass Growth --> Use Original Data 
    {
      # Original Data
      qqtest.function(Master.df$Total.Biomass.Growth....) # p < 2.2e-16 
      
      # Inverse
      qqtest.function(Master.df$Total.Biomass.Growth.Inverse) # Shapiro p = 4.441e-10
      
      # Squared
      qqtest.function(Master.df$Total.Biomass.Growth.Squared) # Shapiro p < 2.2e-16
      
      # Square Root
      qqtest.function(Master.df$Total.Biomass.Growth.Sqrt) # Shapiro p = 4.346e-12
      
      # Log
      qqtest.function(Master.df$Total.Biomass.Growth.Log)  # Shapiro p = 2.11e-10
      
      #Conclusion: Use Original Data for Total Biomass
    }
  }
  
  # CAT Data Transformation QQTESTS
  {
    # CAT Diameter Growth --> Use Diameter.Growth.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$Diameter.Growth....) # Shapiro p = 0.0134
      
      # Inverse
      qqtest.function(CAT.df$Diameter.Growth.Inverse) # Shapiro p = 1.549e-07
      
      # Squared
      qqtest.function(CAT.df$Diameter.Growth.Squared) # Shapiro p = 7.756e-09
      
      # Square Root
      qqtest.function(CAT.df$Diameter.Growth.Sqrt) # Shapiro p = 0.4181
      # Test residuals
      qqtest_species_residuals.function(CAT.df, "Diameter.Growth.Sqrt", "N") # Shapiro p = 0.0977
                                        
      # Log
      qqtest.function(CAT.df$Diameter.Growth.Log) # Shapiro p = 2.563e-11
      
      # Conclusion: Use Sqrt for CAT Diameter
      
    }
    
    # CAT Stem Mass Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(CAT.df$Stem.Mass.Growth....) # Shapiro p = 0.1675
    }
    
    # CAT Coarse Root Growth --> Use Coarse.Root.Growth.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$Coarse.Root.Growth....) # Shapiro p = 0.002243
      
      # Inverse
      qqtest.function(CAT.df$Coarse.Root.Growth.Inverse) # Shapiro p = 0.0002483
      
      # Squared
      qqtest.function(CAT.df$Coarse.Root.Growth.Squared) # Shapiro p = 2.038e-10
      
      # Square Root
      qqtest.function(CAT.df$Coarse.Root.Growth.Sqrt) # Shapiro p = 0.1286
      # Test residuals
      qqtest_species_residuals.function(CAT.df, "Coarse.Root.Growth.Sqrt", "N") # Shapiro p = 0.1786
      
      # Log
      qqtest.function(CAT.df$Coarse.Root.Growth.Log) # Shapiro p = 2.817e-09
      
      # Conclusion: Use Sqrt for CAT Coarse Root
      
    }
    
    # CAT Pre-Existing Biomass Growth --> Use Pre.Existing.Biomass.Growth.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$Pre.Existing.Biomass.Growth....) # Shapiro p = 0.001144
      
      # Inverse
      qqtest.function(CAT.df$Pre.Existing.Biomass.Growth.Inverse) # Shapiro p = 3.084e-12
      
      # Squared
      qqtest.function(CAT.df$Pre.Existing.Biomass.Growth.Squared) # Shapiro p = 6.964e-09
      
      # Square Root
      qqtest.function(CAT.df$Pre.Existing.Biomass.Growth.Sqrt) # Shapiro p = 0.5409
      # Test residuals
      qqtest_species_residuals.function(CAT.df, "Pre.Existing.Biomass.Growth.Sqrt", "N") # Shapiro p = 0.7095
      
      # Log
      qqtest.function(CAT.df$Pre.Existing.Biomass.Growth.Log) # Shapiro p = 0.0008192
      
      # Conclusion: Use Sqrt for CAT Pre-Existing Biomass
    }
    
    # CAT New Above Ground Biomass --> Use New.Above.Ground.Biomass..g.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$New.Above.Ground.Biomass..g.) # Shapiro p = 0.006991
      
      # Inverse 
      qqtest.function(CAT.df$New.Above.Ground.Biomass..g.Inverse) # Shapiro p = 1.101e-05
      
      # Squared 
      qqtest.function(CAT.df$New.Above.Ground.Biomass..g.Squared) # Shapiro p = 3.842e-05
      
      # Square Root
      qqtest.function(CAT.df$New.Above.Ground.Biomass..g.Sqrt) # Shapiro p = 0.0358
      # Test residuals 
      qqtest_species_residuals.function(CAT.df, "New.Above.Ground.Biomass..g.Sqrt", "N") # Shaprio p = 0.4884
      
      # Log
      qqtest.function(CAT.df$New.Above.Ground.Biomass..g.Log) # Shapiro p = 0.0183
      
      # Conclusion: Use Sqrt for CAT New Above Ground Biomass 
      
    }
    
    # CAT Fine Root Biomass --> Use Fine.Root.Mass..g.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$Fine.Root.Mass..g.) # Shapiro p = 0.02063
      
      # Inverse 
      qqtest.function(CAT.df$Fine.Root.Mass..g.Inverse) # Shapiro p = 8.391e-06
      
      # Squared 
      qqtest.function(CAT.df$Fine.Root.Mass..g.Squared) # Shapiro p = 5.729e-07
      
      # Square Root
      qqtest.function(CAT.df$Fine.Root.Mass..g.Sqrt) # Shapiro p = 0.6145
      # Test residuals 
      qqtest_species_residuals.function(CAT.df, "Fine.Root.Mass..g.Sqrt", "N") # Shapiro p = 0.1855
      
      # Log
      qqtest.function(CAT.df$Fine.Root.Mass..g.Log) # Shapiro p = 0.2677
      # Test residuals
      qqtest_species_residuals.function(CAT.df, "Fine.Root.Mass..g.Log", "N") # Shapiro p = 0.4329
      
      # Conclusion: Use Sqrt for CAT Fine Root  
      
    }
    
    # CAT Total Biomass Growth --> Use Total.Biomass.Growth.Sqrt
    {
      # Original Data
      qqtest.function(CAT.df$Total.Biomass.Growth....) # Shaprio p = 0.0006936
      
      # Inverse
      qqtest.function(CAT.df$Total.Biomass.Growth.Inverse) # Shapiro p = 9.488e-10
      
      # Squared
      qqtest.function(CAT.df$Total.Biomass.Growth.Squared) # Shapiro p = 2.833e-09
      
      # Square Root
      qqtest.function(CAT.df$Total.Biomass.Growth.Sqrt) # Shapiro p = 0.6434
      # Test residuals 
      qqtest_species_residuals.function(CAT.df, "Total.Biomass.Growth.Sqrt", "N") # Shapiro p = 0.9913
      
      # Log
      qqtest.function(CAT.df$Total.Biomass.Growth.Log)  # Shapiro p  = 0.05382
      # Test residuals 
      qqtest_species_residuals.function(CAT.df, "Total.Biomass.Growth.Log", "N") # Shaprio p = 0.1562
      
      #Conclusion: Use Square Root for Total Biomass
    }
  }
  
  # MAP Data Transformation QQTESTS
  {
    # MAP Diameter Growth --> Use Original Data
    {
      # Original Data 
      qqtest.function(MAP.df$Diameter.Growth....) # p = 9.621e-09
      
      # Inverse
      qqtest.function(MAP.df$Diameter.Growth.Inverse) # Shapiro p = 5.538e-11
      
      # Squared
      qqtest.function(MAP.df$Diameter.Growth.Squared) # Shapiro p = 8.98e-14
      
      # Square Root
      qqtest.function(MAP.df$Diameter.Growth.Sqrt) # Shapiro p = 4.597e-05
      
      # Log 
      qqtest.function(MAP.df$Diameter.Growth.Log) # Shapiro p = 2.17e-10
      
      # Conclusion: Use Original Data
      
    }
    
    # MAP Stem Mass Growth --> Use Original Data
    {
      # Original Data 
      qqtest.function(MAP.df$Stem.Mass.Growth....) # p = 0.131
      
      # Conclusion: Use Original Data
      
    }
    
    # MAP Coarse Root Growth --> Use Original Data
    {
      # Original Data 
      qqtest.function(MAP.df$Coarse.Root.Growth....) # p = 4005
      
      # Conclusion: Use Original Data
      
    }
    
    # MAP Pre-Existing Biomass Growth --> Use Original Data
    {
      # Original Data 
      qqtest.function(MAP.df$Pre.Existing.Biomass.Growth....) # p = 0.05369
      
      # Conclusion: Use Original Data
      
    }
    
    # MAP New Above Ground --> Use Original Data
    {
      # Original Data
      qqtest.function(MAP.df$New.Above.Ground.Biomass..g.) # Shapiro p = 0.3255
      
      # Conclusion: Use Original Data
    }
    
    # MAP Fine Root --> Use Original Data
    {
      # Original Data
      qqtest.function(MAP.df$Fine.Root.Mass..g.) # Shapiro p = 0.3746
      
      # Conclusion: Use Original Data
    }
    
    # MAP Total Biomass Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(MAP.df$Total.Biomass.Growth....) # Shaprio p = 0.2838
      
      #Conclusion: Use Original Data for Total Biomass
    }
  }
  
  # HON Data Transformation QQTESTS
  {
    # HON Diameter Growth --> Use Diameter.Growth.Sqrt
    {
      # Original Data
      qqtest.function(HON.df$Diameter.Growth....) # Shapiro p = 0.0002828
      
      # Inverse
      qqtest.function(HON.df$Diameter.Growth.Inverse) # Shapiro p = 1.031e-06
      
      # Squared
      qqtest.function(HON.df$Diameter.Growth.Squared) # Shapiro p = 4.777e-10
      
      # Square Root
      qqtest.function(HON.df$Diameter.Growth.Sqrt) # Shapiro p = 0.2134
      # Test residuals
      qqtest_species_residuals.function(HON.df, "Diameter.Growth.Sqrt", "N") # Shapiro p = 0.3275
      
      # Log
      qqtest.function(HON.df$Diameter.Growth.Log) # Shapiro p = 0.001145
      
      # Conclusion: Use Sqrt for HON Diameter Growth
    }
    
    # HON Stem Mass Growth --> No Data
    
    # HON Coarse Root --> Use Original Data
    {
      # Original Data
      qqtest.function(HON.df$Coarse.Root.Growth....) # Shapiro p = 6.103e-08
      
      # Inverse
      qqtest.function(HON.df$Coarse.Root.Growth.Inverse) # Shapiro p = 9.886e-06
      
      # Squared
      qqtest.function(HON.df$Coarse.Root.Growth.Squared) # Shapiro p = 4.428e-13
      
      # Square Root
      qqtest.function(HON.df$Coarse.Root.Growth.Sqrt) # Shapiro p = 0.0008516
      
      # Log
      qqtest.function(HON.df$Coarse.Root.Growth.Log) # Shapiro p = 1.603e-08
      
      # Conclusion: Use Original Data for HON Coarse Root Growth 
    }
    
    # HON Pre-Existing Biomass Growth --> Use Pre.Existing.Biomass.Growth.Sqrt
    {
      # Original Data
      qqtest.function(HON.df$Pre.Existing.Biomass.Growth....) # Shapiro p = 0.001507
      
      # Inverse
      qqtest.function(HON.df$Pre.Existing.Biomass.Growth.Inverse) # Shapiro p = 8.313e-12
      
      # Squared
      qqtest.function(HON.df$Pre.Existing.Biomass.Growth.Squared) # Shapiro p = 8.23e-05
      
      # Square Root
      qqtest.function(HON.df$Pre.Existing.Biomass.Growth.Sqrt) # Shapiro p = 0.4706
      # Test residuals 
      qqtest_species_residuals.function(HON.df, "Pre.Existing.Biomass.Growth.Sqrt", "N") # Shaprio p = 0.6589
      
      # Log
      qqtest.function(HON.df$Pre.Existing.Biomass.Growth.Log) # Shapiro p = 1.502e-08
      
      # Conclusion: Use Log for HON Pre-Existing Biomass Growth  
    }
    
    # HON New Above Ground Biomass --> Use Original Data
    {
      # Original Data
      qqtest.function(HON.df$New.Above.Ground.Biomass..g.) # Shapiro p = 2.89e-08
      
      # Squared
      qqtest.function(HON.df$New.Above.Ground.Biomass..g.Squared) # Shapiro p = 1.721e-09
      
      # Square Root
      qqtest.function(HON.df$New.Above.Ground.Biomass..g.Sqrt) # Shapiro p = 4.579e-07
      
      # Log
      qqtest.function(HON.df$New.Above.Ground.Biomass..g.Log) # Shapiro p = 1.654e-07
      
      # Conclusion: Use Original Data for HON New Above Ground Biomass   
    }
    
    # HON Fine Root --> Use Original Data
    {
      # Original Data
      qqtest.function(HON.df$Fine.Root.Mass..g.) # Shapiro p = 6.914e-09
      
      # Squared
      qqtest.function(HON.df$Fine.Root.Mass..g.Squared) # Shapiro p = 5.216e-10
      
      # Square Root
      qqtest.function(HON.df$Fine.Root.Mass..g.Sqrt) # Shapiro p = 2.252e-07
      
      # Log
      qqtest.function(HON.df$Fine.Root.Mass..g.Log) # Shapiro p = 3.534e-08
      
      # Conclusion: Use Original Data for HON New Above Ground Biomass   
    }
    
    # HON Total Biomass Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(HON.df$Total.Biomass.Growth....) # Shaprio p = 1.602e-06
      
      # Inverse
      qqtest.function(HON.df$Total.Biomass.Growth.Inverse) # Shapiro p = 0.0002801
      
      # Squared
      qqtest.function(HON.df$Total.Biomass.Growth.Squared) # Shapiro p = 4.128e-11
      
      # Square Root
      qqtest.function(HON.df$Total.Biomass.Growth.Sqrt) # Shapiro p = 0.00501
      
      # Log
      qqtest.function(HON.df$Total.Biomass.Growth.Log)  # Shapiro p = 0.005192
      
      #Conclusion: Use Original Data for Total Biomass
    }
  }
  
  # OAK Data Transformation QQTESTS
  {
    # OAK Diameter Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(OAK.df$Diameter.Growth....) # Shapiro p = 0.006552
      
      # Inverse
      qqtest.function(OAK.df$Diameter.Growth.Inverse) # Shapiro p = 8.217e-09
      
      # Squared
      qqtest.function(OAK.df$Diameter.Growth.Squared) # Shapiro p = 4e-11
      
      # Square Root
      qqtest.function(OAK.df$Diameter.Growth.Sqrt) # Shapiro p = 0.003236
      
      # Log
      qqtest.function(OAK.df$Diameter.Growth.Log) # Shapiro p = 8.296e-10
      
      # Conclusion: Use Original Data
      
    }
    
    # OAK Stem Mass Growth --> Use Stem.Mass.Growth.Sqrt
    {
      # Original Data 
      qqtest.function(OAK.df$Stem.Mass.Growth....) # p = 4.877e-05
      
      # Inverse
      qqtest.function(OAK.df$Stem.Mass.Growth.Inverse) # Shapiro p = 0.0003653
      
      # Squared
      qqtest.function(OAK.df$Stem.Mass.Growth.Squared) # Shapiro p = 1.221e-11
      
      # Square Root
      qqtest.function(OAK.df$Stem.Mass.Growth.Sqrt) # Shapiro p = 0.2128
      # Test residuals 
      qqtest_species_residuals.function(OAK.df, "Stem.Mass.Growth.Sqrt", "N") # Shapiro p = 0.3324
      
      # Log
      qqtest.function(OAK.df$Stem.Mass.Growth.Log) # Shapiro p = 4.011e-06
      
      # Conclusion: Use OAK Square Root
      
    }
    
    # OAK Coarse Root Growth --> Use Original Data
    {
      # Original Data
      qqtest.function(OAK.df$Coarse.Root.Growth....) # Shapiro p = 0.3276
      
      # Conclusion: Use Original Data for OAK Coarse Root
    }
    
    # OAK Pre-Existing Biomass --> Use Pre.Existing.Biomass.Growth.Sqrt
    {
      # Original Data
      qqtest.function(OAK.df$Pre.Existing.Biomass.Growth....) # Shapiro p = 0.003739
      
      # Inverse
      qqtest.function(OAK.df$Pre.Existing.Biomass.Growth.Inverse) # Shapiro p = 2.391e-11
      
      # Squared
      qqtest.function(OAK.df$Pre.Existing.Biomass.Growth.Squared) # Shapiro p = 0.04935
      
      # Square Root
      qqtest.function(OAK.df$Pre.Existing.Biomass.Growth.Sqrt) # Shapiro p = 0.6
      # Test residuals 
      qqtest_species_residuals.function(OAK.df, "Pre.Existing.Biomass.Growth.Sqrt", "N") # Shapiro p = 0.9895
      
      # Log
      qqtest.function(OAK.df$Pre.Existing.Biomass.Growth.Log) # Shapiro p = 0.000515
      
      # Conclusion: Use Sqrt for OAK Pre-Existing Biomass
    }
    
    # OAK New Above Ground Biomass --> Use Original Data
    {
      # Original Data
      qqtest.function(OAK.df$New.Above.Ground.Biomass..g.) # Shapiro p = 0.0008717
      
      # Inverse 
      qqtest.function(OAK.df$New.Above.Ground.Biomass..g.Inverse) # Shapiro p = 8.37e-09
      
      # Squared 
      qqtest.function(OAK.df$New.Above.Ground.Biomass..g.Squared) # Shapiro p = 5.747e-08
      
      # Square Root
      qqtest.function(OAK.df$New.Above.Ground.Biomass..g.Sqrt) # Shapiro p = 0.02118
      
      # Log
      qqtest.function(OAK.df$New.Above.Ground.Biomass..g.Log) # Shapiro p = 0.01851
      
      # Conclusion: Use Original Data for OAK New Above Ground Biomass 
      
    }
    
    # OAK Fine Root Biomass --> Use Original Data
    {
      # Original Data
      qqtest.function(OAK.df$Fine.Root.Mass..g.) # Shapiro p = 8.056e-11
      
      # Inverse 
      qqtest.function(OAK.df$Fine.Root.Mass..g.Inverse) # Shapiro p = 5.855e-09
      
      # Squared 
      qqtest.function(OAK.df$Fine.Root.Mass..g.Squared) # Shapiro p = 4.894e-13
      
      # Square Root
      qqtest.function(OAK.df$Fine.Root.Mass..g.Sqrt) # Shapiro p = 1.437e-05
      
      # Log
      qqtest.function(OAK.df$Fine.Root.Mass..g.Log) # Shapiro p = 2.42e-07
      
      # Conclusion: Use Original Data for OAK Fine Root  
      
    }
    
    # OAK Total Biomass Growth --> Use Total.Biomass.Growth.Sqrt
    {
      # Original Data
      qqtest.function(OAK.df$Total.Biomass.Growth....) # Shaprio p = 0.0001928
      
      # Inverse
      qqtest.function(OAK.df$Total.Biomass.Growth.Inverse) # Shapiro p = 0.000159
      
      # Squared
      qqtest.function(OAK.df$Total.Biomass.Growth.Squared) # Shapiro p = 0.01434
      
      # Square Root
      qqtest.function(OAK.df$Total.Biomass.Growth.Sqrt) # Shapiro p = 0.1848
      # Test residuals 
      qqtest_species_residuals.function(OAK.df, "Total.Biomass.Growth.Sqrt", "N") # Shapiro p = 0.2081
      
      # Log
      qqtest.function(OAK.df$Total.Biomass.Growth.Log)  # Shapiro p = 0.02737
      
      #Conclusion: Use Sqrt for Total Biomass
    }
  }
}

#### All Data Models 1 and 3 ####
{
  # Leachate Model 1: Leachate ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Week 1 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.1.Leachate.NaCl..mg., "Week.1.Leachate.NaCl..mg.", "Y")
    # Week 2 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.2.Leachate.NaCl..mg., "Week.2.Leachate.NaCl..mg.", "Y")
    # Week 3 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.3.Leachate.NaCl..mg., "Week.3.Leachate.NaCl..mg.", "Y")
    # Week 4 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.4.Leachate.NaCl..mg., "Week.4.Leachate.NaCl..mg.", "Y")
    # Week 5 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.5.Leachate.NaCl..mg., "Week.5.Leachate.NaCl..mg.", "Y")
    # Week 6 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.6.Leachate.NaCl..mg., "Week.6.Leachate.NaCl..mg.", "Y")
    # Week 7 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.7.Leachate.NaCl..mg., "Week.7.Leachate.NaCl..mg.", "Y")
    # Week 8 Leachate vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Week.8.Leachate.NaCl..mg., "Week.8.Leachate.NaCl..mg.", "Y")
  }
  
  # Leachate Model 3: Leachate ~ Week + (1|Greenhouse.Group)
  {
    #lmer analysis for Overall Leachate by Week
    leachate_NaCl.lmer <- lmer(NaCl.Leachate..mg. ~ Week + (1|Greenhouse.Group), na.action = na.exclude, data = leachate.df)
    #Run ANOVA for leachate_NaCl.lmer, which returns a table in the Console
    anova(leachate_NaCl.lmer) #The result from this tells us that the weeks a significant factor in Na leachate 
    #Post hoc on Week
    emmeans(leachate_NaCl.lmer, list(pairwise ~ Week), adjust = "tukey")
    
    #lmer analysis for Overall Leachate by Week * Salt.Addition
    leachate_with_fill_NaCl.lmer <- lmer(NaCl.Leachate..mg. ~ Week * Salt.Addition + (1|Greenhouse.Group), na.action = na.exclude, data = leachate.df)
    #Run ANOVA for leachate_with_fill_NaCl.lmer, which returns a table in the Console
    anova(leachate_with_fill_NaCl.lmer) #The result from this tells us that the weeks a significant factor in Na leachate 
    #Post hoc on Week * Salt.Addition
    emmeans(leachate_with_fill_NaCl.lmer, list(pairwise ~ Week * Salt.Addition), adjust = "tukey")
  }
  
  # SPAD Model 1: SPAD ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Average SPAD vs. Soil.Type * Salt.Addition
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Average.SPAD, "Average.SPAD", "Y")
  }
  
  # Photosynthesis Model 1: Photosynthesis ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Adjusted A_Net vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Adjusted.A.Net, "Adjusted.A.Net", "Y")
    # Adjusted.A.Net.Cubed for Model 1 Statistics Table (DON'T USE GRAPH, just post hoc test)
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Adjusted.A.Net.Cubed, "Adjusted.A.Net.Cubed", "Y")
  }
  
  # Diameter.Growth.... Model 1: Diameter.Growth.... ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Untransformed Data
    #----
    # Diameter.Growth.... vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Diameter.Growth...., "Diameter.Growth....", "Y")
    #----
    # Transformed Data (Log)
    #----
    # Diameter.Growth.Log vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Diameter.Growth.Log, "Diameter.Growth.Log", "Y")
    #----
  }
  
  # Stem.Mass.Growth.... Model 1: Stem.Mass.Growth.... ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Stem.Mass.Growth.... vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Stem.Mass.Growth...., "Stem.Mass.Growth....", "Y")
  }
 
  # Coarse.Root.Growth.... Model 1: Coarse.Root.Growth.... ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Coarse.Root.Growth.... vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Coarse.Root.Growth...., "Coarse.Root.Growth....", "Y")
  }
 
  # Pre.Existing.Biomass.Growth.... Model 1: Pre.Existing.Biomass.Growth.... ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Pre.Existing.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Pre.Existing.Biomass.Growth...., "Pre.Existing.Biomass.Growth....", "Y")
  }
  
  # Total.Biomass.Growth.... Model 1: Total.Biomass.Growth.... ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Total.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Total.Biomass.Growth...., "Total.Biomass.Growth....", "Y")
  }
  
  # New.Above.Ground.Biomass..g. Model 1: New.Above.Ground.Biomass..g. ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # New.Above.Ground.Biomass..g. vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$New.Above.Ground.Biomass..g., "New.Above.Ground.Biomass..g.", "Y")
  }
  
  # Fine.Root.Mass..g. Model 1: Fine.Root.Mass..g. ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Fine.Root.Mass..g. vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Fine.Root.Mass..g., "Fine.Root.Mass..g.", "Y")
  }
  
  # Roots and Shoots Biomass Model 1: Roots.and.Shoots.Biomass..g. ~ Salt.Addition * Soil.Type + Species + (1|Greenhouse.Group)
  {
    # Roots and Shoots Biomass vs. Soil.Type * Salt.Addition 
    model_1_salt_addition_x_soil.function(Master.df, Master.df$Roots.and.Shoots.Biomass..g., "Roots.and.Shoots.Biomass..g.", "Y") 
  }

}

#### CAT Model 1 ####
{
  # CAT.df 
  
  # CAT SPAD Model 1: SPAD vs. Soil.Type * Salt.Addition 
  {
    # Average.SPAD
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Average.SPAD, "Average.SPAD", "N")
  }
  
  # CAT Adjusted.A.Net Model 1: Adjusted.A.Net vs. Soil.Type * Salt.Addition 
  {
    # Adjusted.A.Net.Squared for Model 1 Statistics Table (DO NOT USE GRAPH, just post hoc test)
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Adjusted.A.Net.Squared, "Adjusted.A.Net.Squared", "N")
  }

  # CAT Diameter.Growth.... Model 1: Diameter.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Diameter.Growth.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Diameter.Growth.Sqrt, "Diameter.Growth.Sqrt", "N")
  }
 
  # CAT Stem.Mass.Growth.... Model 1: Stem.Mass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Stem.Mass.Growth....
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Stem.Mass.Growth...., "Stem.Mass.Growth....", "N")
  }

  # CAT Coarse.Root.Growth.... Model 1: Coarse.Root.Growth.... vs. Soil.Type * Salt.Addition 
  {
  
    # Coarse.Root.Growth.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Coarse.Root.Growth.Sqrt, "Coarse.Root.Growth.Sqrt", "N")
  }
  
  # CAT Pre.Existing.Biomass.Growth.... Model 1: Pre.Existing.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Pre.Existing.Biomass.Growth.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Pre.Existing.Biomass.Growth.Sqrt, "Pre.Existing.Biomass.Growth.Sqrt", "N")
  }
  
  # CAT New.Above.Ground.Biomass..g. Model 1: New.Above.Ground.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # New.Above.Ground.Biomass..g.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$New.Above.Ground.Biomass..g.Sqrt, "New.Above.Ground.Biomass..g.Sqrt", "N")
  }
  
  # CAT Fine.Root.Mass..g. Model 1: Fine.Root.Mass..g. vs. Soil.Type * Salt.Addition 
  {
    # Fine.Root.Mass..g.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Fine.Root.Mass..g.Sqrt, "Fine.Root.Mass..g.Sqrt", "N")
  }
  
  # CAT Total.Biomass.Growth.... Model 1: Total.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Total.Biomass.Growth....Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Total.Biomass.Growth.Sqrt, "Total.Biomass.Growth.Sqrt", "N")
  }
  
  # CAT Roots.and.Shoots.Biomass..g. Model 1: Roots.and.Shoots.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # Roots.and.Shoots.Biomass..g.Sqrt
    model_1_salt_addition_x_soil.function(CAT.df, CAT.df$Roots.and.Shoots.Biomass..g.Sqrt, "Roots.and.Shoots.Biomass..g.Sqrt", "N")
  }

}

#### HON Model 1 ####
{
  # HON.df 
  
  # HON SPAD Model 1: SPAD vs. Soil.Type * Salt.Addition 
  {
   # Average.SPAD
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Average.SPAD, "Average.SPAD", "N")
  }
  
  # HON Adjusted.A.Net Model 1: Adjusted.A.Net vs. Soil.Type * Salt.Addition 
  {
    # Adjusted.A.Net
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Adjusted.A.Net, "Adjusted.A.Net", "N")
  }
  
  # HON Diameter.Growth.... Model 1: Diameter.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Diameter.Growth.Sqrt
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Diameter.Growth.Sqrt, "Diameter.Growth.Sqrt", "N")
  }
  
  # HON Stem.Mass.Growth.... Model 1: Stem.Mass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # No Data
  }
  
  # HON Coarse.Root.Growth.... Model 1: Coarse.Root.Growth.... vs. Soil.Type * Salt.Addition 
  {
  
    # Coarse.Root.Growth....
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Coarse.Root.Growth...., "Coarse.Root.Growth....", "N")
  }
  
  # HON Pre.Existing.Biomass.Growth.... Model 1: Pre.Existing.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Pre.Existing.Biomass.Growth.Sqrt
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Pre.Existing.Biomass.Growth.Sqrt, "Pre.Existing.Biomass.Growth.Sqrt", "N")
  }
  
  # HON New.Above.Ground.Biomass..g. Model 1: New.Above.Ground.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # New.Above.Ground.Biomass..g.
    model_1_salt_addition_x_soil.function(HON.df, HON.df$New.Above.Ground.Biomass..g., "New.Above.Ground.Biomass..g.", "N")
  }
  
  # HON Fine.Root.Mass..g. Model 1: Fine.Root.Mass..g. vs. Soil.Type * Salt.Addition 
  {
    # Fine.Root.Mass..g.
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Fine.Root.Mass..g., "Fine.Root.Mass..g.", "N")
  }
  
  # HON Total.Biomass.Growth.... Model 1: Total.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Total.Biomass.Growth....
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Total.Biomass.Growth...., "Total.Biomass.Growth....", "N")
  }
  
  # HON Roots.and.Shoots.Biomass..g. Model 1: Roots.and.Shoots.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # Untransformed Data
    #----
    # Roots.and.Shoots.Biomass..g.
    model_1_salt_addition_x_soil.function(HON.df, HON.df$Roots.and.Shoots.Biomass..g., "Roots.and.Shoots.Biomass..g.", "N")
    #----
    
    # Significant Variables' Graphs (Raw Data)
    #----
    # None
    #----
    
  }

}
 
#### MAP Model 1 ####
{
  # MAP.df 
  
  # MAP SPAD Model 1: SPAD vs. Soil.Type * Salt.Addition 
  {
    # Average.SPAD
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Average.SPAD, "Average.SPAD", "N")
  }
  
  # MAP Adjusted.A.Net Model 1: Adjusted.A.Net vs. Soil.Type * Salt.Addition 
  {
    # Adjusted.A.Net
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Adjusted.A.Net, "Adjusted.A.Net", "N")
  }
  
  # MAP Diameter.Growth.... Model 1: Diameter.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Diameter.Growth....
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Diameter.Growth...., "Diameter.Growth....", "N")
  }
  
  # MAP Stem.Mass.Growth.... Model 1: Stem.Mass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Stem.Mass.Growth....
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Stem.Mass.Growth...., "Stem.Mass.Growth....", "N")
  }
  
  # MAP Coarse.Root.Growth.... Model 1: Coarse.Root.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Coarse.Root.Growth....
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Coarse.Root.Growth...., "Coarse.Root.Growth....", "N")
  }
  
  # MAP Pre.Existing.Biomass.Growth.... Model 1: Pre.Existing.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Pre.Existing.Biomass.Growth....
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Pre.Existing.Biomass.Growth...., "Pre.Existing.Biomass.Growth....", "N")
  }
 
  # MAP New.Above.Ground.Biomass..g. Model 1: New.Above.Ground.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # New.Above.Ground.Biomass..g.
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$New.Above.Ground.Biomass..g., "New.Above.Ground.Biomass..g.", "N")
  }
 
  # MAP Fine.Root.Mass..g. Model 1: Fine.Root.Mass..g. vs. Soil.Type * Salt.Addition 
  {
    # Fine.Root.Mass..g.
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Fine.Root.Mass..g., "Fine.Root.Mass..g.", "N")
  }
 
  # MAP Total.Biomass.Growth.... Model 1: Total.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Total.Biomass.Growth....
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Total.Biomass.Growth...., "Total.Biomass.Growth....", "N")
  }
  
  # MAP Roots.and.Shoots.Biomass..g. Model 1: Roots.and.Shoots.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # Roots.and.Shoots.Biomass..g.
    model_1_salt_addition_x_soil.function(MAP.df, MAP.df$Roots.and.Shoots.Biomass..g., "Roots.and.Shoots.Biomass..g.", "N")
  
  }
  
}

#### OAK Model 1 ####
{
  # OAK.df 
  
  # OAK SPAD Model 1: SPAD vs. Soil.Type * Salt.Addition 
  {
    # Average.SPAD
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Average.SPAD, "Average.SPAD", "N")
  }
 
  # OAK Adjusted.A.Net Model 1: Adjusted.A.Net vs. Soil.Type * Salt.Addition 
  {
    # Adjusted.A.Net.Squared for Model 1 Statistics Table (DON'T USE THE GRAPH, just the post hoc test)
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Adjusted.A.Net.Squared, "Adjusted.A.Net.Squared", "N")
  }
 
  # OAK Diameter.Growth.... Model 1: Diameter.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Diameter.Growth....
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Diameter.Growth...., "Diameter.Growth....", "N")
  }
  
  # OAK Stem.Mass.Growth.... Model 1: Stem.Mass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Stem.Mass.Growth.Sqrt
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Stem.Mass.Growth.Sqrt, "Stem.Mass.Growth.Sqrt", "N")
  }
  
  # OAK Coarse.Root.Growth.... Model 1: Coarse.Root.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Coarse.Root.Growth....
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Coarse.Root.Growth...., "Coarse.Root.Growth....", "N")
  }
  
  # OAK Pre.Existing.Biomass.Growth.... Model 1: Pre.Existing.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    # Pre.Existing.Biomass.Growth.Sqrt
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Pre.Existing.Biomass.Growth.Sqrt, "Pre.Existing.Biomass.Growth.Sqrt", "N")
  }
  
  # OAK New.Above.Ground.Biomass..g. Model 1: New.Above.Ground.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    # New.Above.Ground.Biomass..g.
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$New.Above.Ground.Biomass..g., "New.Above.Ground.Biomass..g.", "N")
  }
  
  # OAK Fine.Root.Mass..g. Model 1: Fine.Root.Mass..g. vs. Soil.Type * Salt.Addition 
  {
    # Fine.Root.Mass..g.
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Fine.Root.Mass..g., "Fine.Root.Mass..g.", "N")
  }
  
  # OAK Total.Biomass.Growth.... Model 1: Total.Biomass.Growth.... vs. Soil.Type * Salt.Addition 
  {
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Total.Biomass.Growth.Sqrt, "Total.Biomass.Growth.Sqrt", "N")
  }
  
  # OAK Roots.and.Shoots.Biomass..g. Model 1: Roots.and.Shoots.Biomass..g. vs. Soil.Type * Salt.Addition 
  {
    model_1_salt_addition_x_soil.function(OAK.df, OAK.df$Roots.and.Shoots.Biomass..g.Sqrt, "Roots.and.Shoots.Biomass..g.Sqrt", "N")
  }
  
}

#### Beautified Graphs ####
{
  #Plots the Overall Leachate by Week with fill = Salt.Addition
  {
    # Beautified graph: Overall Leachate by Week with fill = Salt.Addition
    #----
    leachate_x_week_salt_fill.plot <- ggplot(leachate.df, aes(x = Week, y = NaCl.Leachate..mg., fill = Salt.Addition)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=16), 
            axis.text.x= element_text(margin = margin(t = 10), size=12),
            axis.text.y=element_text(margin = margin(r = 10), size=12), 
            axis.ticks.length=unit(-0.25, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            legend.title = element_text(colour="black", size=14),
            legend.text = element_text(size=12), 
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "Sodium Leaching (mg Na+)", x = "Week") +
      scale_x_discrete (limits = c("1", "2", "3", "4", "5", "6", "7", "8"), 
                        labels = c("1", "2", "3", "4", "5", "6", "7", "8"))+
      scale_fill_manual(name= "Salt Addition",
                        values=c("gray", "white"), 
                        limits = c("Y", "N"),
                        labels = c("Salt Added", "No Salt Added"))+
      geom_point(position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1),  
                 cex = .75) +
      theme(legend.position = c(0.8,0.8)) #Gets rid of the legend that specifies which color is which group
    
    leachate_x_week_salt_fill.plot
    
    ggsave("leachate_x_week_salt_fill.jpg", dpi = 600, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    
  }
  
  #Plots the Overall Leachate by Week with no fill
  {
    # Beautified graph: Overall Leachate by Week with no fill
    #----
    leachate_x_week_salt_no_fill.plot <- ggplot(leachate.df, aes(x = Week, y = NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=14),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=14), 
            axis.text.x= element_text(margin = margin(t = 10), size=12),
            axis.text.y=element_text(margin = margin(r = 10), size=10), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +  
            labs(x="Week", y=expression(Sodium~Leaching~(mg~Na^{"+"})))+
      scale_x_discrete (limits = c("1", "2", "3", "4", "5", "6", "7", "8"), 
                        labels = c("1", "2", "3", "4", "5", "6", "7", "8")) +
      
      geom_jitter(color="black", size=0.1, alpha=0.9, width = 0.1)
    
    leachate_x_week_salt_no_fill.plot
    
    ggsave("leachate_x_week_salt_no_fill.jpg", dpi = 1200, width = 9,
           height = 9,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    #----
  }
  
  # 9 Panel Graph
  {
    #Beautified 9-panel 1.1
    #----
    # Beautified graph "Week 1 Sodium Leachate vs. Soil Type"
    w1_leachate.df <- subset(leachate.df, Week == "1")
  
    panel_1.1_9_panel.plot <- ggplot(w1_leachate.df, aes(Soil.Type, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,350) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y=expression(Wk~1~Leaching~(mg~Na^{"+"})))+
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    
    panel_1.1_9_panel.plot
    
    ggsave("panel_1.1_9_panel.jpg", width = 6, height = 6, units = "cm", dpi = 1200, 
           path = "C:\\Users\\brian\\Downloads")
    
  
    #----
    #Beautified 9-panel 1.2
    #----
    # Beautified graph "Week 1 Sodium Leachate vs. Species"
    panel_1.2_9_panel.plot <- ggplot(w1_leachate.df, aes(Species, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,350) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_9_panel.plot
    
    ggsave("panel_1.2_9_panel.jpg", width = 6, height = 6, units = "cm", dpi = 1200, 
           path = "C:\\Users\\brian\\Downloads")
    #----
    #Beautified 9-panel 1.3
    #----
    # Beautified graph "Week 1 Sodium Leachate vs. Salt Addition"
    panel_1.3_9_panel.plot <- ggplot(w1_leachate.df, aes(Salt.Addition, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,350) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      scale_x_discrete (limits = c("N", "Y")) +
      
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.3_9_panel.plot
    
    ggsave("panel_1.3_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 2.1
    #----
    # Beautified graph "Week 4 Sodium Leachate vs. Soil Type"
    w4_leachate.df <- subset(leachate.df, Week == "4")
    panel_2.1_9_panel.plot <- ggplot(w4_leachate.df, aes(Soil.Type, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,80) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y=expression(Wk~4~Leaching~(mg~Na^{"+"})))+
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_9_panel.plot
    
    ggsave("panel_2.1_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 2.2
    #----
    # Beautified graph "Week 4 Sodium Leachate vs. Species"
    panel_2.2_9_panel.plot <- ggplot(w4_leachate.df, aes(Species, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,80) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_9_panel.plot
    
    ggsave("panel_2.2_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 2.3
    #----
    # Beautified graph "Week 4 Sodium Leachate vs. Salt Addition"
    panel_2.3_9_panel.plot <- ggplot(w4_leachate.df, aes(Salt.Addition, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,80) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      scale_x_discrete (limits = c("N", "Y")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.3_9_panel.plot
    
    ggsave("panel_2.3_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 3.1
    #----
    # Beautified graph "Week 8 Sodium Leachate vs. Soil Type"
    w8_leachate.df <- subset(leachate.df, Week == "8")
    panel_3.1_9_panel.plot <- ggplot(w8_leachate.df, aes(Soil.Type, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,35) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y=expression(Wk~8~Leaching~(mg~Na^{"+"})), x="Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_3.1_9_panel.plot
    
    ggsave("panel_3.1_9_panel.jpg",width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 3.2
    #----
    # Beautified graph "Week 8 Sodium Leachate vs. Species"
    panel_3.2_9_panel.plot <- ggplot(w8_leachate.df, aes(Species, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,35) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      # labs(y = "Sodium Leaching (mg Na+)", x = "Species") +
      labs(y = "", x = "Tree Species") + 
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_3.2_9_panel.plot
    
    ggsave("panel_3.2_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Beautified 9-panel 3.3
    #----
    # Beautified graph "Week 8 Sodium Leachate vs. Salt Addition"
    panel_3.3_9_panel.plot <- ggplot(w8_leachate.df, aes(Salt.Addition, NaCl.Leachate..mg.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,35) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            legend.title = element_text(colour="black", size=14),
            legend.text = element_text(size=14), 
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(x = "Salt Addition") + 
      scale_x_discrete (limits = c("N", "Y")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_3.3_9_panel.plot
    
    ggsave("panel_3.3_9_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 1,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #Full Beautified 9-panel 
    #----
    full_9_panel <- plot_grid(labels = "AUTO", panel_1.1_9_panel.plot, panel_1.2_9_panel.plot, panel_1.3_9_panel.plot,
                              panel_2.1_9_panel.plot, panel_2.2_9_panel.plot, panel_2.3_9_panel.plot,
                              panel_3.1_9_panel.plot, panel_3.2_9_panel.plot, panel_3.3_9_panel.plot, 
                              rel_heights = c(1,1,1.15),
                              rel_widths = c(1.15,1,1), 
                              label_size = 12, 
                              label_x = 0, label_y = 1.05
    )
    
    full_9_panel
    
    
    ggsave("full_9_panel.jpg", dpi = 1200, width = 18,
           height = 16,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    #----
  }
  
  # 4 Panel Graph
  {
    # Beautified 4-panel 1.1 
    #----
    # Beautified graph "Average SPAD vs. Species"
    panel_1.1_4_panel.plot <- ggplot(Master.df, aes(Species, Average.SPAD)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,47) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "SPAD") +
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +    
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_4_panel.plot
    
    ggsave("panel_1.1_4_panel.jpg", width = 2, height = 2.5, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 4-panel 1.2
    #----
    # Beautified graph "CAT Average SPAD vs. Soil Type"
    panel_1.2_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Average.SPAD)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,32) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "", x = "") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_4_panel.plot
    
    ggsave("panel_1.2_4_panel.jpg", width = 2, height = 2.5, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    # Beautified 4-panel 2.1
    #----
    # Beautified graph "Photosynthesis vs. Species"
    panel_2.1_4_panel.plot <- ggplot(Master.df, aes(Species, Adjusted.A.Net)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,8) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y=expression('Photosynthesis Rate\n'(mu*mol ~ C ~ m^{-2} ~ s^{-1})),
           x  = "Tree Species") + 
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +     
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_4_panel.plot
    
    ggsave("panel_2.1_4_panel.jpg", width = 2, height = 2.5, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 4-panel 2.2
    #----
    # Beautified graph "MAP Photosynthesis vs. Soil Type"
    panel_2.2_4_panel.plot <- ggplot(MAP.df, aes(Soil.Type, Adjusted.A.Net, fill = Salt.Addition)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,1.5) +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            legend.title = element_text(colour="black", size=7),
            legend.text = element_text(size=5),
            legend.box="horizontal",
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = " ", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      scale_fill_manual(name= "Salt Addition",
                        values=c("white", "gray"), 
                        limits = c("N", "Y"))+
      guides(color = guide_legend(nrow = 1), byrow=
               TRUE)+
      geom_point(size = 0.2, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
      theme(legend.position = c(0.5,0.82))
    
    panel_2.2_4_panel.plot
    
    ggsave("panel_2.2_4_panel.jpg", width = 2, height = 2.5, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Full Beautified 4-panel
    #----
    full_4_panel <- plot_grid(labels = "AUTO", panel_1.1_4_panel.plot, panel_1.2_4_panel.plot,
                              panel_2.1_4_panel.plot, panel_2.2_4_panel.plot,  
                              rel_heights = c(1,1.15),
                              rel_widths = c(1.10,1), 
                              label_size = 12,
                              align="V"
    )
    
    full_4_panel
    
    
    ggsave("full_4_panel.jpg", dpi = 1200, width = 12,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    #----
  }
  
  # 6-Panel Graph 
  {
    # Beautified 6-Panel Graph 1.1
    #----
    panel_1.1_6_panel.plot <- ggplot(Master.df, aes(Soil.Type, Pre.Existing.Biomass.Growth....)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,400) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "Stem and Coarse Root \nBiomass (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_6_panel.plot
    
    ggsave("panel_1.1_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 6-Panel Graph 1.2
    #----
    panel_1.2_6_panel.plot <- ggplot(Master.df, aes(Species, Pre.Existing.Biomass.Growth....)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,400) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +      
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_6_panel.plot
    
    ggsave("panel_1.2_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 6-Panel Graph 1.3
    #----
    panel_1.3_6_panel.plot <- ggplot(Master.df, aes(Salt.Addition, Pre.Existing.Biomass.Growth....)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,400) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=8), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      scale_x_discrete (limits = c("N", "Y")) +
      
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.3_6_panel.plot
    
    ggsave("panel_1.3_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 6-Panel Graph 2.1
    #----
    panel_2.1_6_panel.plot <- ggplot(Master.df, aes(Soil.Type, Roots.and.Shoots.Biomass..g.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,17) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "Shoot and Fine Root \nBiomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_6_panel.plot
    
    ggsave("panel_2.1_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 6-Panel Graph 2.2
    #----
    panel_2.2_6_panel.plot <- ggplot(Master.df, aes(Species, Roots.and.Shoots.Biomass..g.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,17) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(x = "Tree Species") +
      scale_x_discrete (limits = c("CAT", "MAP", "HON", "OAK"), labels = c("CASP", "ACSA", "GLTR", "QURU")) +      
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_6_panel.plot
    
    ggsave("panel_2.2_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified 6-Panel Graph 2.3  
    #----
    panel_2.3_6_panel.plot <- ggplot(Master.df, aes(Salt.Addition, Roots.and.Shoots.Biomass..g.)) + 
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,17) + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_blank(), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) + 
      labs(y = "", x = "Salt Addition") +
      scale_x_discrete (limits = c("N", "Y")) +
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.3_6_panel.plot
    
    ggsave("panel_2.3_6_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Full Beautified 6-Panel
    #----
    full_6_panel <- plot_grid(labels = "AUTO", panel_1.1_6_panel.plot, panel_1.2_6_panel.plot, panel_1.3_6_panel.plot,
                              panel_2.1_6_panel.plot, panel_2.2_6_panel.plot, panel_2.3_6_panel.plot,
                              rel_heights = c(1,1.15),
                              rel_widths = c(1.18,1,1), 
                              label_size = 12, 
                              label_x = 0, label_y = 1.05
    )
    
    full_6_panel
    
    ggsave("full_6_panel.jpg", dpi = 1200, width = 18,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    
    #----
    #----
  }
  
  # CAT 4-Panel Graph 
  {
    # Beautified CAT 4-Panel Graph 1.1
    #----
    panel_1.1_CAT_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Stem.Mass.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw()  +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Stem Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +    
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_CAT_4_panel.plot
    
    ggsave("panel_1.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified CAT 4-Panel Graph 1.2
    #----
    panel_1.2_CAT_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Coarse.Root.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Coarse Root Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_CAT_4_panel.plot
    
    ggsave("panel_1.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified CAT 4-Panel Graph 2.1
    #----
    panel_2.1_CAT_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, New.Above.Ground.Biomass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,13) +  
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Shoot Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_CAT_4_panel.plot
    
    ggsave("panel_2.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified CAT 4-Panel Graph 2.2
    #----
    panel_2.2_CAT_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Fine.Root.Mass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.25, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Fine Root Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_CAT_4_panel.plot
    
    ggsave("panel_2.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Full Beautified CAT 4-Panel
    #----
    full_CAT_4_panel <- plot_grid(labels = "AUTO", panel_1.1_CAT_4_panel.plot, panel_1.2_CAT_4_panel.plot,
                                  panel_2.1_CAT_4_panel.plot, panel_2.2_CAT_4_panel.plot, 
                                  rel_heights = c(1,1.15), 
                                  rel_widths = c(1,1.02), 
                                  label_size = 12, 
                                  align="v")
  
    full_CAT_4_panel
  
  
    ggsave("full_CAT_4_panel.jpg", dpi = 1200, width = 13.1,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    
    #----
  }
  
  # HON 4-Panel Graph 
  {
    # Original 4-Panel Graph
    #----
    # HON 4-Panel Graph Row 1
    # HON 4-Panel Graph 1.1: Stem Mass Growth (%) vs. Soil.Type
    # BECAUSE THERE IS NO INITIAL STEM DATA FOR HON, THIS GRAPH CANNOT SHOW ANY DATA
    ggplot(HON.df, aes(Soil.Type, Stem.Mass.Growth...., na.rm = TRUE)) + geom_boxplot() + 
      stat_boxplot(geom = "errorbar") + ggtitle("HON Stem Growth vs. Soil Type")
    
    # HON 4-Panel Graph 1.2: Coarse Root Growth (%) vs. Soil.Type
    ggplot(HON.df, aes(Soil.Type, Coarse.Root.Growth...., na.rm = TRUE)) + geom_boxplot() + 
      stat_boxplot(geom = "errorbar") + ggtitle("HON Coarse Root Growth vs. Soil Type")
    
    # HON 4-Panel Graph Row 2
    # HON 4-Panel Graph 2.1: New Above Ground Biomass (g) vs. Soil.Type
    ggplot(HON.df, aes(Soil.Type, New.Above.Ground.Biomass..g., na.rm = TRUE)) + geom_boxplot() + 
      stat_boxplot(geom = "errorbar") + ggtitle("HON New Above Ground Biomass vs. Soil Type")
    
    # HON 4-Panel Graph 2.2: Fine Root Mass (g) vs. Soil.Type
    ggplot(HON.df, aes(Soil.Type, Fine.Root.Mass..g., na.rm = TRUE)) + geom_boxplot() + 
      stat_boxplot(geom = "errorbar") + ggtitle("HON Fine Root Biomass vs. Soil Type")
    #----
    
    # Beautified HON 4-Panel Graph 1.1
    #----
    panel_1.1_HON_4_panel.plot <- ggplot(HON.df, aes(Soil.Type, Stem.Mass.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Stem Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +    
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_HON_4_panel.plot
    
    ggsave("panel_1.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified HON 4-Panel Graph 1.2
    #----
    panel_1.2_HON_4_panel.plot <- ggplot(HON.df, aes(Soil.Type, Coarse.Root.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Coarse Root Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_HON_4_panel.plot
    
    ggsave("panel_1.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified HON 4-Panel Graph 2.1
    #----
    panel_2.1_HON_4_panel.plot <- ggplot(HON.df, aes(Soil.Type, New.Above.Ground.Biomass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Shoot Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_HON_4_panel.plot
    
    ggsave("panel_2.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified HON 4-Panel Graph 2.2
    #----
    panel_2.2_HON_4_panel.plot <- ggplot(HON.df, aes(Soil.Type, Fine.Root.Mass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.25, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Fine Root Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_HON_4_panel.plot
    
    ggsave("panel_2.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Full Beautified HON 4-Panel
    #----
    full_HON_4_panel <- plot_grid(labels = "AUTO", panel_1.1_HON_4_panel.plot, panel_1.2_HON_4_panel.plot,
                                  panel_2.1_HON_4_panel.plot, panel_2.2_HON_4_panel.plot, 
                                  rel_heights = c(1,1.15), 
                                  rel_widths = c(1,1.02), 
                                  label_size = 12, 
                                  align="v")
    
    full_HON_4_panel
    
    
    ggsave("full_HON_4_panel.jpg", dpi = 1200, width = 13.1,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\mmidgley\\Downloads")
    #----
    # Beautified MAP 4-Panel Graph 1.2
    #----
    panel_1.1_MAP_4_panel.plot <- ggplot(MAP.df, aes(Soil.Type, Stem.Mass.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Stem Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +    
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_MAP_4_panel.plot
    
    ggsave("panel_1.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified MAP 4-Panel Graph 1.2
    #----
    panel_1.2_MAP_4_panel.plot <- ggplot(MAP.df, aes(Soil.Type, Coarse.Root.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Coarse Root Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_MAP_4_panel.plot
    
    ggsave("panel_1.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified MAP 4-Panel Graph 2.1
    #----
    panel_2.1_MAP_4_panel.plot <- ggplot(MAP.df, aes(Soil.Type, New.Above.Ground.Biomass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Shoot Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_MAP_4_panel.plot
    
    ggsave("panel_2.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Beautified MAP 4-Panel Graph 2.2
    #----
    panel_2.2_MAP_4_panel.plot <- ggplot(MAP.df, aes(Soil.Type, Fine.Root.Mass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.25, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Fine Root Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_MAP_4_panel.plot
    
    ggsave("panel_2.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Documents\\Personal\\Morton Arboretum\\SIR\\Data\\Images & Graphs from R\\Beautified Graphs Fall 2022\\CAT 4 Panel Graph")
    #----
    # Full Beautified HON 4-Panel
    #----
    full_MAP_4_panel <- plot_grid(labels = "AUTO", panel_1.1_MAP_4_panel.plot, panel_1.2_MAP_4_panel.plot,
                                  panel_2.1_MAP_4_panel.plot, panel_2.2_MAP_4_panel.plot, 
                                  rel_heights = c(1,1.15), 
                                  rel_widths = c(1,1.02), 
                                  label_size = 12, 
                                  align="v")
    
    full_MAP_4_panel
    
    
    ggsave("full_MAP_4_panel.jpg", dpi = 1200, width = 13.1,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\mmidgley\\Downloads")
    #----
  }
  
  # OAK 4-Panel Graph 
  {
    # Beautified OAK 4-Panel Graph 1.1
    #----
    panel_1.1_OAK_4_panel.plot <- ggplot(OAK.df, aes(Soil.Type, Stem.Mass.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() + 
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Stem Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +    
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.1_OAK_4_panel.plot
    
    ggsave("panel_1.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified OAK 4-Panel Graph 1.2
    #----
    panel_1.2_OAK_4_panel.plot <- ggplot(OAK.df, aes(Soil.Type, Coarse.Root.Growth....)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Coarse Root Growth (%)") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_1.2_OAK_4_panel.plot
    
    ggsave("panel_1.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified OAK 4-Panel Graph 2.1
    #----
    panel_2.1_OAK_4_panel.plot <- ggplot(OAK.df, aes(Soil.Type, New.Above.Ground.Biomass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.15, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Shoot Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.1_OAK_4_panel.plot
    
    ggsave("panel_2.1_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Beautified OAK 4-Panel Graph 2.2
    #----
    panel_2.2_OAK_4_panel.plot <- ggplot(OAK.df, aes(Soil.Type, Fine.Root.Mass..g.)) +
      stat_boxplot(geom = "errorbar") + theme_linedraw() +
      geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.border = element_rect(fill=NA, colour = "black", size=.7), 
            axis.title.x = element_text(margin = margin(t = 10, b=5), size=10),
            axis.title.y = element_text(margin = margin(l = 5, r=5), size=10), 
            axis.text.x= element_text(margin = margin(t = 10), size=9),
            axis.text.y=element_text(margin = margin(r = 10), size=9), 
            axis.ticks.length=unit(-0.25, "cm"),
            #axis.ticks.margin=unit(0.5, "cm"),
            axis.ticks = element_line(colour = "black", size = 0.4), 
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10)) +
      labs(y = "Fine Root Biomass (g)", x = "Soil Treatment") +
      scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) + 
      geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)
    
    panel_2.2_OAK_4_panel.plot
    
    ggsave("panel_2.2_CAT_4_panel.jpg", width = 2.5, height = 2, units = "in", dpi = 50, scale = 2,
           path = "C:\\Users\\brian\\Downloads")
    #----
    # Full Beautified OAK 4-Panel
    #----
    full_OAK_4_panel <- plot_grid(labels = "AUTO", panel_1.1_OAK_4_panel.plot, panel_1.2_OAK_4_panel.plot,
                                  panel_2.1_OAK_4_panel.plot, panel_2.2_OAK_4_panel.plot, 
                                  rel_heights = c(1,1.15), 
                                  rel_widths = c(1,1.02), 
                                  label_size = 12, 
                                  align="v")
    
    full_OAK_4_panel
    
    
    ggsave("full_OAK_4_panel.jpg", dpi = 1200, width = 13.1,
           height = 10.92,
           units = "cm",
           path = "C:\\Users\\brian\\Downloads")
    #----
  }

}

#WUE for CAT
model_1_salt_addition_x_soil.function(CAT.df, CAT.df$WUE, "WUE", "N")

model_1_salt_addition_x_soil.function(HON.df, HON.df$WUE, "WUE", "N")

model_1_salt_addition_x_soil.function(MAP.df, MAP.df$WUE, "WUE", "N")

model_1_salt_addition_x_soil.function(OAK.df, OAK.df$WUE, "WUE", "N")

model_1_salt_addition_x_soil.function(CAT.df, CAT.df$g_sw, "g_sw", "N")

model_1_salt_addition_x_soil.function(HON.df, HON.df$g_sw, "g_sw", "N")

model_1_salt_addition_x_soil.function(MAP.df, MAP.df$g_sw, "g_sw", "N")

model_1_salt_addition_x_soil.function(OAK.df, OAK.df$g_sw, "g_sw", "N")

model_1_salt_addition_x_soil.function(Master.df, Master.df$g_sw, "g_sw", "Y")
#----
# Beautified 4-panel 1.1 
#----
# Beautified graph "CAT Average SPAD vs. Soil Type"
panel_1.1_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Average.SPAD)) + 
  stat_boxplot(geom = "errorbar") + theme_linedraw() + ylim(NA,32) + 
  geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_blank(),
        axis.text.x= element_text(margin = margin(t = 10), size=9),
        axis.text.y=element_text(margin = margin(r = 10), size=9), 
        axis.ticks.length=unit(-0.15, "cm"),
        #axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) + 
  labs(y="SPAD",
       x  = "") +     scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
  geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)

panel_1.1_4_panel.plot

#----
# Beautified 4-panel 1.2
#----
# Beautified graph "CAT Photosynthesis vs. Soil Type"
panel_1.2_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, Adjusted.A.Net)) + 
  stat_boxplot(geom = "errorbar") + theme_linedraw()  + 
  geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.title.x = element_blank(),
        axis.text.x= element_text(margin = margin(t = 10), size=9),
        axis.text.y=element_text(margin = margin(r = 10), size=9), 
        axis.ticks.length=unit(-0.15, "cm"),
        #axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) + 
  labs(y=expression('e\n'(mu*mol ~ C ~ m^{-2} ~ s^{-1})),
       x  = "") + 
  scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
  geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)

panel_1.2_4_panel.plot

# Beautified graph "CAT Stomatal Conductance vs. Soil Type"
panel_2.1_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, g_sw)) + 
  stat_boxplot(geom = "errorbar") + theme_linedraw()  + ylim(NA,.45) +  
  geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.text.x= element_text(margin = margin(t = 10), size=9),
        axis.text.y=element_text(margin = margin(r = 10), size=9), 
        axis.ticks.length=unit(-0.15, "cm"),
        #axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) + 
  labs(y=expression('Stomatal Conductance\n'(mol ~ H[2]*O~ m^{-2} ~ s^{-1})),
       x  =  "Soil Treatment") + 
  scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
  geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)

panel_2.1_4_panel.plot

#----
# Beautified graph "CAT WUE vs. Soil Type"
panel_2.2_4_panel.plot <- ggplot(CAT.df, aes(Soil.Type, WUE)) + 
  stat_boxplot(geom = "errorbar") + theme_linedraw()  + ylim(NA,150) +  
  geom_boxplot(notch=F, lwd=.5, colour='black', stat="boxplot", outlier.shape=NA) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", size=.7), 
        axis.text.x= element_text(margin = margin(t = 10), size=9),
        axis.text.y=element_text(margin = margin(r = 10), size=9), 
        axis.ticks.length=unit(-0.15, "cm"),
        #axis.ticks.margin=unit(0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.4), 
        axis.ticks.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) + 
  labs(y=expression('WUE'~(mu*mol ~ C~mol ~ H[2]*O^{-1})),
       x  = "Soil Treatment") + 
  scale_x_discrete (limits=c("Control", "BCMI", "BCOT"), labels = c("Control", "BC-I", "BC-TD")) +
  geom_jitter(color="black", size=0.2, alpha=0.9, width = 0.1)

panel_2.2_4_panel.plot

full_4_panel <- plot_grid(labels = "AUTO", panel_1.1_4_panel.plot, panel_1.2_4_panel.plot,
                          panel_2.1_4_panel.plot, panel_2.2_4_panel.plot, 
                              rel_heights = c(1,1.15), 
                              rel_widths = c(1,1.02), 
                              label_size = 12, 
                              align="v")

full_4_panel


png(file="C:\\Users\\mmidgley\\Downloads\\full_4_panel3.png",
    height=10.92, width=13.1, units="cm", res=1200)
plot(full_4_panel)
dev.off()

ggsave("full_4_panel", dpi = 1200, width = 13.1,
       height = 10.92,
       units = "cm",
       path = "C:\\Users\\mmidgley\\Downloads")
#----
}
