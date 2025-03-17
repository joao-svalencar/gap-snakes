#How to simply calculate Petchey and Gaston's Functional Diversity (FD)
#First, install the picante and ade4 packages

#install.packages("picante")
#install.packages("ade4")
library(picante)
library(ade4)
library(FD)
library(vegan)
library(ape)
library(SYNCSA)

comp<-read.csv("comm_sbv.csv", head=T, sep=";") #Load the 'community' file
traits<-read.csv("traits_sbv.csv", head=T, sep=";") # Load the 'traits' file
comp
traits
leaf_area <- rnorm(53,1.5,0.4)
leaf_area
traits$leaf_area2 <- leaf_area
sla <- rnorm(53,2.5,0.8)
traits$sla2 <- sla
str(traits)
traits <-traits[,c(2:6)]
head(traits)
traits
row.names(traits) <-traits$genus_sp

#Matriz de dist?ncia
dist.func<-vegdist(traits,method= "gower")
dist.func
?vegdist

#Use this distance matrix to construct a dendrogram that represents the similarity/dissimilarity among species according to their ecological traits
tree<-hclust(dist.func,"average") #Average ? o m?todo UPGMA
tree
plot(tree)

#Now we have to transform the dendrogram into a '.phylo' file
tree.p<-as.phylo(tree)
tree.p
par(mfrow=c(1,1), mai=c(0,1,0,1))
plot(tree.p, cex=1.0) #cex representa o tamanho da fonte
dev.off()

#Finally, calculate the functional diversity for each community (see 'help' to see how to include or not the root of the dendrogram; the default is to include it).
FD<-pd(comp, tree.p)
FD #retorna os valores de FD
plot(FD$SR, FD$PD)

#write.csv(FD,"fd_result.csv") #Salvar resultado em planilha Excel (.csv)

#You can see the values for each community on the column "PD"
#(this is correct, as Petchey and Gaston's FD is calculated with the same method as Faith's Phylogenetic Diversity - PD).
#You can also find the species richness ("SR") for each community.
## ?ndices de diversidade funcional baseados em dist?ncia (veja a ajuda das fun??es para maiores detalhes)

# ?ndice MPD (Webb 2000, Amer Nat; Webb et al. 2002 Ann Rev Ecol Syst):

MPD<-mpd(comp,cophenetic(tree.p))
MPD

#MPD.d <- mpd(comp,as.matrix(D)) #A diferen?a ? que a dist?ncia aqui n?o vem da ?rvore filogen?tica (cofen?tica), mas sim da matriz de dist?ncias (fun??o dist.ktab); o ideal ? utilizar as dist?ncias da matriz!
#MPD.d

#plot(MPD, MPD.d)
#abline(0, 1, lty=2)
#write.csv(MPD.d,"mpd_result.csv")

#plot(FD$SR, MPD.d) #MPD tende a n?o ter rela??o com riqueza

### ?ndice MNTD (Webb 2000, Amer Nat; Webb et al. 2002 Ann Rev Ecol Syst):

#MNTD<-mntd(comp,cophenetic(tree.p))
#MNTD.d<-mntd(comp,as.matrix(D)) #MNTD funcional pode ser interessante para casos em que h? invas?o por esp?cies e se quer saber se invasoras s?o similares ?s nativas (aumentam ou n?o espa?o funcional)
#MNTD.d

#plot(MNTD, MNTD.d)
#abline(0,1,lty=2)

#plot(FD$SR, MNTD.d) #Tende a ter uma rela??o negativa com a riqueza, pois bancos maiores de esp?cies possuem maior chance de adicionar uma esp?cie similar a outra j? ocorrente em uma comunidade rica

#MNTD.d <-data.frame(MNTD.d)
#colnames(MNTD.d ) <- "MNTD.d"

#write.csv(MNTD.d,"mntd_result.csv")

## A fun??o "dbFD" do pacote "FD" permite calcular v?rios ?ndices de diversidade funcional baseados em dist?ncia (e CWM) em uma s? linha de comando.
 # Antes de mais nada, ler "Details" na ajuda da fun??o.

?dbFD
fix(comp)
# FRic (Cornwell et al. 2006, Ecology; Vill?ger et al. 2008 Ecology)
# FEve (Vill?ger et al. 2008, Ecology)
# FDiv (Vill?ger et al. 2008, Ecology)
# FDis (Lalibert? & Legendre 2010, Ecology)
# Rao, Q (Botta-Duk?t 2005, J Veg Sci)

DBFD<-dbFD(D,comp, m=10, print.pco = T, corr = c("cailliez")) #A fun??o dbFD trabalha com matrizes de dist?ncias do tipo euclidiano. Como "dist.func" foi criada a partir do m?todo de Pavoine et al. 2008 ("Gower generalizado"), ent?o aqui foi necess?rio fazer uma corre??o da matriz de dist?ncia para transform?-la em euclidiana. Ver outras op??es de corre??o em "Details" da ajuda da fun??o.

DBFD

ls(DBFD)


DBFD$FRic
DBFD$FEve
DBFD$FDiv
DBFD$RaoQ
DBFD$FDis

#write.csv(DBFD$FRic,"FRic_result.csv")
#write.csv(DBFD$FEve,"FEve_result.csv")
#write.csv(DBFD$FDiv,"FDiv_result.csv")
#write.csv(DBFD$RaoQ,"RaoQ_result.csv")
#write.csv(DBFD$FDis,"FDis_result.csv")

DBFD$nbsp #vector listing the number of species in each community
DBFD$qual.FRic #quality of the reduced-space representation required to compute FRic and FDiv.
#DBFD$sing.sp #vector listing the number of functionally singular species in each community. If all species are functionally different, sing.sp will be identical to nbsp.
#DBFD$x.axes #PCoA axes. Only returned if print.pco is TRUE.
#DBFD$x.values #eigenvalues from the PCoA. Only returned if print.pco is TRUE.

plot(DBFD$nbsp, DBFD$FRic) #Volume funcional (FRic) satura com o aumento do n?mero de esp?cies (Mouchet et al. 2010, Funct. Ecol.)
plot(DBFD$nbsp, DBFD$FEve) #FEve em geral n?o ? relacionado com a riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)
plot(DBFD$nbsp, DBFD$FDiv) #FDiv em geral n?o ? relacionado com a riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)

plot(DBFD$nbsp, DBFD$RaoQ)  #RaoQ em geral ? pouco influenciado pela riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)porque seus valores s?o m?dias
plot(DBFD$nbsp, DBFD$FDis)

#At? que ponto ?ndices de diversidade funcional nos informam coisas distintas?

index_all <- data.frame(MPD.d, MNTD.d, DBFD$FRic, DBFD$FEve, DBFD$FDiv, DBFD$RaoQ,
           DBFD$FDis)

index_all
#plot(index_all)

#plot(MPD.d, DBFD$RaoQ)
#plot(DBFD$FDis, DBFD$RaoQ)
#plot(MPD.d, DBFD$FDis)
#plot(DBFD$RaoQ, DBFD$FDiv)

###SES ? uma estat?stica Z, tamb?m usada como transforma??o vetorial (normaliza??o), quando (obs - m?dia.esp) / sd.esp.
#Independentswap ? um modelo nulo mais restritivo do que os demais.
#Signific?ncia de z contra m?dia zero ? pobre. A maioria dos valores em geral n?o difere do nulo. 
#P-Fisher como s?ntese da tend?ncia dos Ps ? uma alternativa.
#Usar a tend?ncia dos valores de Z (ses) em an?lises de gradientes ? bem aceito.
#Tamanho de efeito forte (significativo) ? aquele menor do que -1.96 ou maior do 1.96 (retirado de uma pop te?rica com g.l. infinitos)
## Modelos nulos com FD, MPD and MNTD  #S?o melhores para filogen?tica
# Modelo nulo - FD

ses.FD<-ses.pd(comp,tree.p,null.model="independentswap",runs=1000,iterations=100) 
ses.FD
#ses.FD_lib<-ses.pd(comp,tree.p,null.model="taxa.labels",runs=100,iterations=10) 


#write.csv(ses.FD,"ses.FD_indsw_result.csv")
#write.csv(ses.FD_lib,"ses.FD_txlab_result.csv")

#write.csv(ses.FD,"ses.FD_indsw_result.csv")
#write.csv(ses.FD_lib,"ses.FD_txlab_result.csv")

#MPD null model:

#ses.MPD<-ses.mpd(comm,dist.func,null.model="independentswap",runs=1000,iterations=10) #d? para calcular com matriz de dist?ncias

#write.csv(ses.MPD,"ses.MPD_result.csv")

#MNTD null model:

#ses.MNTD<-ses.mntd(comm,dist.func,null.model="independentswap",runs=1000,iterations=10)

#write.csv(ses.MNTD,"mntd_nm_result.csv")

# CWM (Garnier et al. 2004, Ecology); Community-weigthed trait means n?o ? um ?ndice de diversidade funcional, mas sim uma m?dia do atributo no n?vel de comunidade ponderada pelas abund?ncias relativas das esp?cies
# A fun??o "matrix.t" do pacote SYNCSA gera uma matriz de CWM para v?rios atributos de uma s? vez (Pillar et al. 2009, J Veg Sci)

summary(colnames(comp)==rownames(traits))

resu_matrix.t <- matrix.t(comp, traits, scale = FALSE, notification = TRUE)
resu_matrix.t

ls(resu_matrix.t)

cwm_matrix <- as.data.frame(resu_matrix.t$matrix.T)
cwm_matrix

# Testando mudan?as m?dias do atributo ao longo do gradiente ("trait mean shifts")

#test.grad <- data.frame(rnorm(45, 50, 20), row.names = rownames(cwm_matrix)) #Inventa um gradiente ambiental
#colnames(test.grad) <- "test.grad"

#plot(test.grad$test.grad, cwm_matrix$BodyMass.Value)
#plot(test.grad$test.grad, cwm_matrix$Diet.Inv)
#plot(test.grad$test.grad, cwm_matrix$Diurnal)
