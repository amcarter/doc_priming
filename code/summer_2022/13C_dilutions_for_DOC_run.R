# Calculate dilutions for running 13C DOC samples through picarro
#
# Target DOC concentrations: 2-10 mg/L
# Target d13C: 200 ppt

# delta to AF functions:####
deltoAF<- function (del) {
    Rst<-0.0112372
    R<-(del/1000+1)*Rst
    AF<-R/(1+R)

    return(AF)
}

AFtodel <- function (AF) {
    Rst<-0.0112372
    R<- AF/(1-AF)
    del<- (R/Rst-1)*1000

    return(del)
}


# Constants and target concentrations: ####
# assumed del values of each carbon source:

lake_del = -25
glu_del = -12   # appx del for corn
target_del = 200

lake_AF = deltoAF(lake_del)
glu_AF = deltoAF(glu_del)
target_AF = deltoAF(target_del)

# approximate concentrations (mg/L)
lake_conc = 1.2
glu_conc = 1
lechate_conc = 200

lechate_add = lechate_conc * 0.75/1000 # add 0.75 ml/L

# Calculations ####
# lake water + lechate treatment (LW_L)

# calculate total mg 13C/total C (AF)

LW_L_AF = (lechate_add + lake_conc * lake_AF)/
    (lechate_add + lake_conc)

LW_L_del = AFtodel(LW_L_AF)


# calculate glucose needed to dilute to a del of 200
# target_AF = (lechate_add + lake_conc * lake_AF + glu * glu_AF)/
#               (lechate_add + lake_conc + glu)

# target_AF*(lechate_add + lake_conc + glu) =
#           (lechate_add + lake_conc * lake_AF + glu * glu_AF)

# (target_AF - glu_AF) * glu = (lechate_add + lake_conc * lake_AF) - target_AF*(lechate_add + lake_conc)

# glu = (lechate_add * (1- target_AF) + lake_conc *(lake_AF - target_AF))/(target_AF - glu_AF)
glu = (lechate_add * (1- target_AF) + lake_conc *(lake_AF - target_AF))/(target_AF - glu_AF)


# expected del of sample:
glu = 60
exp_AF = (lechate_add + lake_conc * lake_AF + glu * glu_AF)/(lechate_add + lake_conc + glu)
exp_del = AFtodel(exp_AF)

# Dilution Calculations ####

# Glucose added to 4 ml of sample:

glu_add = glu * 4/1000

sample_conc = (glu_add + lake_conc * 4/1000 + lechate_add * 4/1000) /40 *1000

# stoc_glu (mgC /L)* 36ml * (1 L/1000 ml) = glu_add (mg C)

stock_glu = glu_add /36 *1000

# target 6.5 mg C/L for stock solution, made using a 750 mgC/L stock:

# ml_glu * 750 = 6.5 * 2000
ml_glu = 6.5 * 2000 / 750

# target phosphoric acid solution:

# 80 uL * 85% H2PO4 = 40 ml * conc PA
# 80 ul / 40 ml * (1 ml/1000 ul) * 85%
conc_PA = 80/40 * (1/1000) * 0.85

ml_85_PA = 2000/40 *80/1000

# Lechate dilution: ####
# make a 1 L solution of 6.5 mg/L glucose with lechate added to be
glu_lechate_dilution = 6.5 * 40/1000
# target_AF = (lechate + glu_lechate_dilution * glu_AF)/(lechate + glu_lechate_dilution)
# lechate + glu_lechate_dilution * glu_AF = target_AF *( lechate + glu_lechate_dilution)
# lechate * (1 - target_AF) = target_AF * glu_lechate_dilution - glu_lechate_dilution * glu_AF

lechate = (glu_lechate_dilution *(target_AF - glu_AF))/ (1 - target_AF)
lechate / 200

# 200 mg/L * ml_lechate = 6 mg/L * 40 ml
ml_lechate = 2 * 40 / 200
