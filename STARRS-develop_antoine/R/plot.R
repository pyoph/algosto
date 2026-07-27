
#' Robust Mixture Model
#'
#' @param x Output of \link{multipleRobustMM}
#' @param showOutliers Logical. If FALSE (by default), outliers are removed from the plot
#' @param graph Character vector. The type of plot to be displayed. Possible values are:
#' \describe{
#' \item{"Two_Dim"}{Scatter plot of the projection of all the data points on the principal plan (along the two principal components) with color hue depending on clustering}
#' \item{"Two_Dim_Uncertainty"}{Scatter plot of the projection of all the data points on the principal plan (along the two principal components) with color hue depending on clustering and size depending on Uncertainty}
#' \item{"ICL"}{Plot of the ICL criterion as a function of the number of clusters}
#' \item{"BIC"}{Plot of the BIC criterion as a function of the number of clusters}
#' \item{"Profiles"}{Plot of the profiles of the data points with the median profile of each cluster}
#' \item{"Uncertainty"}{Plot of the Uncertainty of the data points}
#' }
#' @param K Integer. The number of clusters to be displayed. If NULL (by default), the best result is displayed
#' @param ... Additional parameters
#'
#' @return A plot for each item in the `graph` parameter
#'
#' @md
#'
#' @importFrom ggplot2 .data
#' @export
#' @method plot RMM
plot.RMM <- function(x,
                    showOutliers=FALSE,
                      graph=c('Two_Dim','Two_Dim_Uncertainty','ICL','BIC','Profiles','Uncertainty'),
                      K=NULL, ...)
{
  X <- x$data
  d <- ncol(X)
  checkmate::assert_int(d, lower=2)
  if (('BIC' %in% graph) | ('ICL' %in% graph)){
    # TODO dans ce cas, mettre un message explicite, comme quoi demander icl/bic pour un seul K n'a pas de sens
    checkmate::assert_int(length(x$nclust), lower=2)
  }

  # extract only the two first components (those with higher eigen values)
  # independant des bestresult ou K
  Gmed <- WeiszfeldRobustPCA(X,scores = 2)
  vec <- Gmed$vectors
  med0 <- Gmed$median

  med <- c(sum(med0*vec[,1]),sum(med0*vec[,2]))
  Xtrans <- matrix(0,nrow=nrow(X),ncol=2)
  for(i in 1:nrow(X))
  {
    Xtrans[i,1] <- sum(X[i,]*vec[,1])-med[1]
    Xtrans[i,2] <- sum(X[i,]*vec[,2])-med[2]
  }

  if (!is.null(K)){
    Ksel <- K
    Knclust <- which(x$nclust==Ksel)
    res <- x$allresults[[Knclust]]
  }
  else {
    Ksel <- x$Kopt
    res <- x$bestresult
  }

  cluster <- res$clusters

  tau <- res$tau
  tau <- tau*log(tau)
  Uncertainty <- -rowSums(tau)/max(abs(-rowSums(tau))) #TODO corriger: remplacer max(abs(-rowSums(tau))) par log(K) (a voir avec stephane)

  dataplot <- data.frame(x=Xtrans[,1],y=Xtrans[,2],K=as.character(cluster), Uncertainty=Uncertainty)

  # si on a reperé des outliers on les vire pourqu'ils ne pourissent pas le graphique
  # on pourrait les garder en les projetant sur les bords avec une autre couleur
  if (showOutliers==FALSE){
    if (length(res$outliers)>0){
      I <- res$outliers
      dataplot <- dataplot[-I,]
    }
  }


  #---------------------------------------------------------
  # TODO a decouper
  # scatter plot of the projection of all the data points on the principal plan (along the two principal components)
  # with color hue, depending on clustering
  if ('Two_Dim' %in% graph){
    plt <- ggplot(data=dataplot) +
      scale_x_continuous('First Principal Component',  limits = c(min(dataplot$x),max(dataplot$x))) +
      scale_y_continuous('Second Principal Component', limits = c(min(dataplot$y),max(dataplot$y))) +
      geom_point(aes(x = .data$x, y = .data$y,color = .data$K))
    print(plt)
  }

  #---------------------------------------------------------
  # affiche les points avec en plus Uncertainty calculé (qui change la taille des ppoints)
  if ('Two_Dim_Uncertainty' %in% graph){

    plt2 <- ggplot(data=dataplot) +
      scale_x_continuous('First Principal Component', limits = c(min(dataplot$x),max(dataplot$x))) +
      scale_y_continuous('Second Principal Component',limits = c(min(dataplot$y),max(dataplot$y))) +
      scale_size_continuous(range=c(1,4))+
      geom_point(aes(x = .data$x, y = .data$y,color = .data$K,size = .data$Uncertainty))
    print(plt2)
  }

  #------------------------------
  # TODO voir mclust
  if ('ICL' %in% graph){
    data_icl <- data.frame(K=x$nclust,ICL=x$ICL)
    plt_icl <- ggplot(data=data_icl)+
      geom_line(aes(x=.data$K,y=.data$ICL))
    print(plt_icl)
  }

  #------------------------------

  if ('BIC' %in% graph){
    data_bic <- data.frame(K=x$nclust,BIC=x$BIC)
    plt_bic <- ggplot(data=data_bic)+
      geom_line(aes(x=.data$K,y=.data$BIC))
    print(plt_bic)
  }

  # ------------------------------------
  # on fait un plot par classe avec tous les profils et le profil median
  if ('Profiles' %in% graph){
    Xt <- X
    cluster2 <- cluster

    if (showOutliers==FALSE)
    {
      if (length(res$outliers)>0){
        I <- res$outliers
        Xt <- Xt[-I,]
        cluster2 <- cluster2[-I]
      }
    }

    if (d>2){
      dataclust <- data.frame(x=1:d,y=t(Xt))
      datacen <- data.frame(xt=1:d,t(res$centers)) #TODO t per forza ?
      ind <- c()
      Kclu <- c()
      for (o in 1:nrow(Xt))
      {
        ind <- c(ind,rep(o,d))
        Kclu <-c(Kclu,rep(as.character(cluster2[o]),d))
      }
      Kcen <- c()
      indK <- c()
      col <- c()
      for (o in 1:Ksel)
      {
        I <- length(which(Kclu==as.character(o)))
        Kcen <- c(Kcen,rep(as.character(o),d))
        indK <- c(indK,rep((o+max(ind)),d))
      }
      data_new <- melt(dataclust, id = c("x"))
      data_cen <- melt(datacen, id = c("xt"))
      datanew <- data.frame(x=c(data_new$x,data_cen$xt),y=c(data_new$value,data_cen$value),
                          couleur=c(ind,indK),Kclu=c(Kclu,Kcen),
                          coul=c(rep('grey',length(data_new$x)),rep('red',length(data_cen$xt))))

      plt_profile <- ggplot(data=datanew,aes(x=.data$x , y=.data$y,color=.data$coul))+
        scale_color_manual(values=c('grey','red'))+
        scale_x_continuous('') +
        scale_y_continuous('Profiles',limits=c(min(Xt),max(Xt))) +
        geom_line(mapping=aes(x=.data$x , y=.data$y , group=.data$couleur),show.legend = F,alpha=0.4)+
        #  geom_line(aes(color=(coul)))+
        facet_wrap(~as.factor(Kclu))
      print(plt_profile)
    }
    if (d==2){
      dataclust <- data.frame(x=1:d,y=t(Xt))
      Kclu <- c()
      for (o in 1:nrow(Xt))
      {
        Kclu <- c(Kclu,as.character(cluster2[o]))
      }
      Kcen <- c()
      for (o in 1:Ksel)
      {
        Kcen <- c(Kcen,as.character(o))
      }
      data_clust <- data.frame(x=c(Xt[,1],(res$centers)[,1]),y=c(Xt[,2],(res$centers)[,2]),
                            K=c(Kclu,Kcen),coul=c(rep('grey',nrow(Xt)),rep('red',Ksel)))
      plt_profile <- ggplot(data=data_clust,aes(x=.data$x , y=.data$y,color=.data$coul))+
        scale_color_manual(values=c('grey','red'))+
        scale_x_continuous('') +
        scale_y_continuous('Profile',limits=c(min(Xt),max(Xt))) +
        geom_point(mapping=aes(x=.data$x , y=.data$y),show.legend = F,alpha=0.4)+
        #  geom_line(aes(color=(coul)))+
        facet_wrap(~as.factor(K))
      print(plt_profile)
    }
  }

  #----------------------------------------
  if ('Uncertainty' %in% graph)
  {
    ord <- order(Uncertainty)
    Uncertainty <- Uncertainty[ord]

    data_unc <- data.frame(X=1:nrow(X),Uncertainty=Uncertainty)
    plt_unc <- ggplot(data=data_unc,aes(x=.data$X,y=.data$Uncertainty))+
      geom_bar(stat='identity')
    print(plt_unc)
  }
}
