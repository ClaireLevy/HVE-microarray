#making a design matrix based on limma vignette pg 42


#set up the targets frame, making sure the sample order
#is the same as in the data (3 controls, 3 treatments)
targets<-data.frame("CellLine"=c("HVE1","HVE2","HVE3",
                                 "HVE1","HVE2","HVE3"),
                    "Treatment"=c("control","control","control",
                                  "drug","drug","drug"))

#using the targets frame, make a design matrix that accounts for
#pairing within Cell lines (i.e. HVE1 dose=0 and HVE1 dose=500 
#go together). see limma guide pg 42

CellLine<-factor(targets$CellLine, levels=c("HVE1","HVE2","HVE3"))
Treatment<-factor(targets$Treatment, levels=c("control","drug"))
#make the design matrix
design<-model.matrix(~CellLine + Treatment)
