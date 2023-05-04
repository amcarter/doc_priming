
#initial isotope calcs

deltoAF<- function (del) {
  Rst<-0.0112372
  R<-(del/1000+1)*Rst
  AF<-R/(1+R)
  AF
  
}


AFtodel <- function (AF) {
  Rst<-0.0112372
  R<- AF/(1-AF)
  del<- (R/Rst-1)*1000
  
  del
}

##added 0.5 mL leachate to 7 L water.  Assume 200ug/mL
doc_13c<- 0.5*200/7  #ugC/L
doc_13c/1000

#assume 2meq/L alk and therefore DOC

DIC<- 12* 2  ##mgC/L

## DIC went from -16 to say -10
deltoAF(-10)*DIC - deltoAF(-16) *DIC 
##1.5 ug has entered the DIC pool

#what is the del of the DOC pool?  Assume 2 mgC/L

AFtodel ( (2*deltoAF(-27) + 0.014)/2.0014 )
