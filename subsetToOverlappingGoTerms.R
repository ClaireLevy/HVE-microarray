subsetToOverlappingGoTerms <- function(df, direction, withDay1) {
   if (withDay1) {
      toWrite <- df[df$direction == direction, ]
      toWrite <- toWrite %>%
         group_by(Term)%>%
         filter(n()==4)%>% # because there are 4 days
         arrange(Day,Term,PValue)%>%
         ungroup()
   } else {
      toWrite <- df[df$direction == direction & df$Day != 1, ]
      toWrite <- toWrite %>%
         group_by(Term)%>%
         filter(n()==3)%>% # because there are 3 days
         arrange(Day,Term,PValue)%>%
         ungroup()
   }
}