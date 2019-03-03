require(caper)
require(car)
require(plotly)
require(ggplot2)
lepdata <-read.table("lepdata_trimmed2.csv", sep = ",", header = TRUE)
leptree <-read.nexus("Matthee.nex")
#or  
leptree <-read.tree("leptree.txt")

row.names(lepdata) <- lepdata$Name_phyl

leporid <- comparative.data(phy = leptree, data = lepdata, names.col = Name_phyl, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#obtain residuals

#for total brain/body
model.pgls.res <- pgls(log(ECV_T) ~ log(AdultBodyMass), data = leporid, lambda='ML')
res <- model.pgls.res$res

#for OB rel to T
model.pgls.res.ob <- pgls(log(ECV_OB) ~ log(ECV_T), data = leporid, lambda='ML')
res.ob <- model.pgls.res.ob$res


#models

model.Se1<-pgls(log(ECV_OB) ~ log(SeT)*log(AdultBodyMass), data = leporid, lambda='ML')
summary(model.Se1) # oB but not ROB or T (both absolute and relative work here!)


model.Se2<-pgls(res.ob ~ log(SeP), data = leporid, lambda='ML')
summary(model.Se2)


model.GR<-pgls(log(ECV_T) ~ (log(GR_Area))*log(AdultBodyMass), data = leporid, lambda='ML')
summary(model.GR) # T and ROB but not OB or relative OB (Home Range NS) - I think nothing from here


#plot
ggplot(lepdata, aes(x = res.ob, y = log(SeP))) + geom_point()

       #or
p <- plot_ly(lepdata, x = ~log(ECV_OB)/log(ECV_T), y = ~log(SeT), z = ~log(AdultBodyMass), color = ~Activity, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'RelBrain'),
                      yaxis = list(title = 'SeT'),
                      zaxis = list(title = 'BodyMass')))




#Check for vifs

DataCA=data.frame(res.ob, lepdata$SeT, lepdata$SeP, lepdata$AdultBodyMass)
Lm_CA=lm(res.ob ~ lepdata.SeP*lepdata.AdultBodyMass, data=DataCA)
vif(Lm_CA) # low vifs


cor.test(lepdata$SeT, lepdata$SeP) #low cross correlation between SeT & SeP


