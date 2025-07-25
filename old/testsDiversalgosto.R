
k = 0.86;l=0.56;rho = 0.6
r = 30
contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r )

resNaif = SampleCovOnline(data$Z)

norm(resNaif$Sigma - Sigma0,"F")

table(data$labelsVrais,resNaif$outliers_labels)

resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)

norm(resUsOnline$Sigma[n,,] - Sigma0,"F")

table(data$labelsVrais,resUsOnline$outlier_labels)

resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))

norm(resUsStreaming$Sigma[n,,] - Sigma0,"F")

table(data$labelsVrais,resUsStreaming$outlier_labels)


k = 8.59;l=32;rho = 0.975
r = 39
contParam = ParmsF1(m1, k, l, rho1)

data = genererEchantillon(n,n,mu1 = mu0,mu2 = contParam$mu1,Sigma1 = Sigma0,Sigma2 = contParam$Sigma1,r )

resNaif = SampleCovOnline(data$Z)

norm(resNaif$Sigma - Sigma0,"F")

table(data$labelsVrais,resNaif$outliers_labels)

resUsOnline= StreamingOutlierDetection(data$Z,batch = 1)

norm(resUsOnline$Sigma[n,,] - Sigma0,"F")

table(data$labelsVrais,resUsOnline$outlier_labels)

resUsStreaming= StreamingOutlierDetection(data$Z,batch = ncol(data$Z))

norm(resUsStreaming$Sigma[n,,] - Sigma0,"F")

table(data$labelsVrais,resUsStreaming$outlier_labels)

eigen(resUsOnline$Sigma[n,,])$values

eigen(Sigma0)$values

