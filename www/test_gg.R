
tmp11 <- tmp1 %>% gather( 'CHF', 'Pneumonia', 'ILD_IPF', 'ILD_CTD', 'ILD_other', 'LungScars', key="Disease",value="Prob")
tmp111 <- tmp11 %>% dplyr::filter(ID == 'Init_prob' | ID =='Final_prob' )
ggplot(tmp111, aes(fill=ID, y=Prob, x=Disease)) + 
  geom_bar(position="dodge", stat="identity")

