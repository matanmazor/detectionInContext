library('tidyverse')

read.csv('..\\experiments\\pilots\\hangman2_pilot\\data\\jatos_results_batch1.csv',na.strings=c(""," ","NA")) %>%
  group_by(PROLIFIC_PID) %>%
  summarise(a=1)%>%
  pivot_wider(values_from=a,names_from=PROLIFIC_PID) %>%
  write.table(file='..\\experiments\\Hangman\\data\\subj_list.csv', sep=",", row.names = FALSE, quote=FALSE)
