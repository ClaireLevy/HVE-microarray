#DAVID data has the GO id and term together as
# GO:xxxxxx~term and this function extracts the GO id part
# and puts it in a new column called "Pathway Id"
# It is called pathway id because that is what the GO id column
# is called in the innateDB data and this makes it easier to
#merge the data sets





getGO<-function(df){
  mutate(df,
         Pathway.Id= ifelse(str_detect(df$Term,"GO:")==TRUE,
                            substr(df$Term,1,10),NA))
}
