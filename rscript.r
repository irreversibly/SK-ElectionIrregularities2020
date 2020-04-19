setwd("C:/Users/user/Dropbox/ProduceKR")
library(eforensics)
set.seed(12345)
mcmc = list(burn.in=10000, n.adapt=10000, n.iter=10000, n.chains=10)

fnames <- c('Jongno2020','Jongno2016',
'JungguSeongdongguA2020','JungguSeongdongguA2016',
'JungguSeongdongguB2020','JungguSeongdongguB2016',
'Yongsan2020','Yongsan2016',
'GwangjinA2020','GwangjinA2016',
'GwangjinB2020','GwangjinB2016',
'SongpaB2020','SongpaB2016',
'GoyangsiA2020','GoyangsiA2016')

for(fname in fnames) {

	data <-read.csv(file=paste('Data/',fname,'.csv',sep=''))

	samples = eforensics(
		n.w ~ 1,
		n.a ~ 1,
		mu.iota.m ~ ind,
		mu.iota.s ~ ind,
		mu.chi.m  ~ ind,
		mu.chi.s  ~ ind,
		data=data,
		eligible.voters="n.rv",
		model="qbl",
		mcmc=mcmc,
		parameters = "all",
		parComp = TRUE,
		autoConv = TRUE,
		max.auto = 0,
		mcmc.conv.diagnostic = "MCMCSE",
		mcmc.conv.parameters = c("pi"),
		mcmcse.conv.precision = .05,
		mcmcse.combine = FALSE
	)

	coefs <- summary(samples, join.chains=T)

	num_z <- samples[[1]]$piZi*100
	for (i in 2:mcmc$n.chains) {
		num_z <- num_z + samples[[i]]$piZi*100
	}
	num_z <- num_z/mcmc$n.chains
	max_z <- apply(num_z,1,which.max)
	
	write.csv(coefs,paste('Outputs/',fname,'-coefs.csv',sep=''))
	write.csv(num_z,paste('Outputs/',fname,'-Zvalues.csv',sep=''))

}

