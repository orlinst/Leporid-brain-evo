require(caper)
require(car)
require(plotly)
require(ggplot)
lepdata <-read.table("lepdata_trimmed2.csv", sep = ",", header = TRUE)
leptree <-read.nexus("Matthee.nex")
#or  
leptree <-read.tree("leptree.txt")

row.names(lepdata) <- lepdata$Name_phyl

leporid <- comparative.data(phy = leptree, data = lepdata, names.col = Name_phyl, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)


#models

model.Se<-pgls(log(ECV_OB)/log(ECV_T) ~ (log(SeT)+log(SeP))*log(AdultBodyMass), data = leporid, lambda='ML')
summary(model.Se) # oB but not ROB or T (both absolute and relative work here!)
plot(model.Se)



model.GR<-pgls(log(ECV_T) ~ (log(GR_Area)):log(AdultBodyMass), data = leporid, lambda='ML')
summary(model.GR) # T and ROB but not OB or relative OB (Home Range NS) - I think nothing from here






model.pgls.res.lep<-pgls(log(ECV_OB) ~ log(AdultBodyMass), data = leporid, lambda='ML')

model.pgls.res.lep<- residuals(model.pgls.res.lep, phylo = TRUE)


mod.l <- pgls.profile(model.pgls, 'lambda')

model.pgls<-pgls(log(ECV_T) ~ (log(ST)+log(SP))*log(AdultBodyMass), data = leporid, lambda='ML')
summary(model.pgls)


#plot
ggplot(lepdata, aes(x = log(ECV_OB)/log(ECV_T), y = log(SeP)*log(SeT))) + geom_point()


#or
p <- plot_ly(lepdata, x = ~log(ECV_OB)/log(ECV_T), y = ~log(SeT), z = ~log(AdultBodyMass), color = ~Activity, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'RelBrain'),
                      yaxis = list(title = 'SeT'),
                      zaxis = list(title = 'BodyMass')))




#Check for vifs

DataCA=data.frame(lepdata$ECV_OB, lepdata$SeT, lepdata$SeP, lepdata$AdultBodyMass)
Lm_CA=lm(lepdata.ECV_OB ~ (lepdata.SeT+lepdata.SeP)*lepdata.AdultBodyMass, data=DataCA)
#or disregarding body size
Lm_CA=lm(lepdata.ECV_OB ~ lepdata.SeT+lepdata.SeP, data=DataCA)
vif(Lm_CA) # low vifs


cor.test(lepdata$SeT, lepdata$SeP) #low cross correlation between SeT & SeP


