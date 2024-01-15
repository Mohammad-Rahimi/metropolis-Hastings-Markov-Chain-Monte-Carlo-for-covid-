
setwd('G:/My Drive/wvu/2nd semester/601/codes/homework11/project')
newjerseycases<-read.csv("./data_table_for_daily_case_trends__new_jersey.csv")
newdaily<-newjerseycases$New.Cases

date<-as.Date(newjerseycases$Date,format="%b %d %Y")
plot(date,newdaily,type="l", main="NewJersey new daily cases")
axis(1,date,format(date, "%b %d"), cex.axis=0.8)

