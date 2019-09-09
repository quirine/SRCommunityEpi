# Load Data ---------------------------------------------------------------

Data05 <- read.csv("../Data/Longevity_05FAR.csv",skip=4, nrows=41)
Data075 <- read.csv("../Data/Longevity_075FAR.csv",skip=4, nrows=41)
Data1 <- read.csv("../Data/Longevity_1FAR.csv",skip=4, nrows=41)
Data125 <- read.csv("../Data/Longevity_125FAR.csv",skip=4, nrows=41)
Data15 <- read.csv("../Data/Longevity_15FAR.csv",skip=4, nrows=41)


# Prepare the Data for Analysis -------------------------------------------
Data05 <- Data05[,1:20]; Data05$X. <- NULL; Data05$X..1 <- NULL; Data05$X <- NULL; 
Numbers.released <- read.csv("../Data/Longevity_05FAR.csv",skip=49, nrows=1)
Numbers.released$X <- NULL; Numbers.released$X.1 <- NULL; Numbers.released$X.2 <- NULL; Numbers.released$X.3 <- NULL; Numbers.released$X.4 <- NULL
Numbers.released$X.5 <- NULL; Numbers.released$X.6 <- NULL
temp <- as.numeric(Numbers.released); newrow <- c(temp[1:7],mean(temp[1:6]),temp[8:14],mean(temp[8:13])); newrow <- c(-1,newrow)
Data05 <- rbind(newrow,Data05)  # adding an extra NULL row to capture the instant deaths at t=0

Data075 <- Data075[,1:20]; Data075$X. <- NULL; Data075$X..1 <- NULL; Data075$X <- NULL; 
Numbers.released <- read.csv("../Data/Longevity_075FAR.csv",skip=49, nrows=1)
Numbers.released$X <- NULL; Numbers.released$X.1 <- NULL; Numbers.released$X.2 <- NULL; Numbers.released$X.3 <- NULL; Numbers.released$X.4 <- NULL
Numbers.released$X.5 <- NULL; Numbers.released$X.6 <- NULL
temp <- as.numeric(Numbers.released); newrow <- c(temp[1:7],mean(temp[1:6]),temp[8:14],mean(temp[8:13])); newrow <- c(-1,newrow)
Data075 <- rbind(newrow,Data075)  # adding an extra NULL row to capture the instant deaths at t=0

Data1 <- Data1[,1:20]; Data1$X. <- NULL; Data1$X..1 <- NULL; Data1$X <- NULL; 
Numbers.released <- read.csv("../Data/Longevity_1FAR.csv",skip=49, nrows=1)
Numbers.released$X <- NULL; Numbers.released$X.1 <- NULL; Numbers.released$X.2 <- NULL; Numbers.released$X.3 <- NULL; Numbers.released$X.4 <- NULL
Numbers.released$X.5 <- NULL; Numbers.released$X.6 <- NULL
temp <- as.numeric(Numbers.released); newrow <- c(temp[1:7],mean(temp[1:6]),temp[8:14],mean(temp[8:13])); newrow <- c(-1,newrow)
Data1 <- rbind(newrow,Data1)  # adding an extra NULL row to capture the instant deaths at t=0

Data125 <- Data125[,1:20]; Data125$X. <- NULL; Data125$X..1 <- NULL; Data125$X <- NULL; 
Numbers.released <- read.csv("../Data/Longevity_125FAR.csv",skip=49, nrows=1)
Numbers.released$X <- NULL; Numbers.released$X.1 <- NULL; Numbers.released$X.2 <- NULL; Numbers.released$X.3 <- NULL; Numbers.released$X.4 <- NULL
Numbers.released$X.5 <- NULL; Numbers.released$X.6 <- NULL
temp <- as.numeric(Numbers.released); newrow <- c(temp[1:7],mean(temp[1:6]),temp[8:14],mean(temp[8:13])); newrow <- c(-1,newrow)
Data125 <- rbind(newrow,Data125)  # adding an extra NULL row to capture the instant deaths at t=0

Data15 <- Data15[,1:20]; Data15$X. <- NULL; Data15$X..1 <- NULL; Data15$X <- NULL; 
Numbers.released <- read.csv("../Data/Longevity_15FAR.csv",skip=49, nrows=1)
Numbers.released$X <- NULL; Numbers.released$X.1 <- NULL; Numbers.released$X.2 <- NULL; Numbers.released$X.3 <- NULL; Numbers.released$X.4 <- NULL
Numbers.released$X.5 <- NULL; Numbers.released$X.6 <- NULL
temp <- as.numeric(Numbers.released); newrow <- c(temp[1:7],mean(temp[1:6]),temp[8:14],mean(temp[8:13])); newrow <- c(-1,newrow)
Data15 <- rbind(newrow,Data15)  # adding an extra NULL row to capture the instant deaths at t=0

# Derive Time to Event Data -----------------------------------------------

SurvSR <- rbind(get.TTE(Data05,0.5), get.TTE(Data075,0.75), get.TTE(Data1,1.00), get.TTE(Data125,1.25), get.TTE(Data15,1.5))
I <- which(SurvSR$TTE > con.time)
SurvSR.init = SurvSR[-I,]
SurvSR.con = SurvSR[I,]

SurvSR.con <- upData(SurvSR.con, labels = c(ID = 'Mosquito id', TTE = 'Days at risk',Event = 'Died or survived',
                                    Treat = 'SR or control', Dose = 'Dosage relative to FAR',
                                    Group = 'lab batch'),
                 levels = list(Treat = list('control' = 0, 'SR' = 1),
                               Event = list('survived' = 0, 'died' = 1)))

SurvSR.init <- upData(SurvSR.init, labels = c(ID = 'Mosquito id', TTE = 'Days at risk',Event = 'Died or survived',
                                              Treat = 'SR or control', Dose = 'Dosage relative to FAR',
                                              Group = 'lab batch'),
                      levels = list(Treat = list('control' = 0, 'SR' = 1),
                                    Event = list('survived' = 0, 'died' = 1)))