library(ade4)
library(ape)
library(picante)
library(vegan)
library(FD)
library(SYNCSA)

#abrindo planilhas método read.table
atributos<-read.table(here::here("functional-meet", "atributos.txt"), h=T)
atributos
comunidade<-read.table(here::here("functional-meet", "comunidade.txt"), h=T)
comunidade

#logartimizando dados de abundância
#?decostand
comunidadelog <- log(comunidade+1)  # quando os dados são quantitativos na mesma escala de mensuração
comunidadelog

#matriz de distância de Gower para dados mistos
traits <- FD::gowdis(atributos, asym.bin = NULL, ord = c("classic"))
traits
?gowdis

#Criando matriz de distância de Gower 'generalizado' (Pavoine et a. 2009)	

#atividade-AT	
#microhabitat-MH	
#vocalization-VO 	
#reproductive-R	
#Total- TL
#Mouth-M	
#Forearm-F	
#Femur-FE

#AT:
at.trait <- data.frame(AT=atributos[,1])
#MH
mh.trait <- data.frame(MH=atributos[,2])
#VO
vo.trait <- data.frame(VO=atributos[,3])
#R
r.trait <- data.frame(R=atributos[,4])
#TL
tl.trait <- data.frame(TL=atributos[,5])
#M
m.trait <- data.frame(M=atributos[,6])
#F
f.trait <- data.frame(F=atributos[,7])
#FE
fe.trait <- data.frame(FE=atributos[,8])


#Fazendo a lista de variáveis:
#Aqui adicionei o argumento rownames pra colocar o nome das especies
ktab <- ade4::ktab.list.df(list(at.trait,mh.trait,vo.trait,r.trait,tl.trait,m.trait,f.trait,fe.trait), rownames = row.names(atributos))
?ktab.list.df

#Making the distance matrix according Pavoine et al. (2009):
#matriz de distância feita com a medida de pavoine, "N" quer dizer variáveis nominais
p.dist <- ade4::dist.ktab(ktab, type=c("N","N","N","N","Q","Q","Q","Q"), option=c("scaledBYrange"))
?dist.ktab

#A função dbFD trabalha com matrizes de distâncias do tipo euclidiano. 
#Como "dist.func" foi criada a partir do método de Pavoine et al. 2008 ("Gower generalizado")
#Então aqui foi necessário fazer uma correção da matriz de distância para transformá-la em euclidiana.
#Ver outras opções de correção em "Details" da ajuda da função.

#Calcular índices de DF:
DBFD <- FD::dbFD(p.dist, comunidadelog, m=8, print.pco = T, calc.CWM= FALSE, corr = c("sqrt")) 

#Matriz antiga (errada): note o uso de "traits" ao inves de "p.dist"
#DBFD <- FD::dbFD(traits, comunidadelog, m=8, print.pco = T, calc.CWM= FALSE, corr = c("sqrt")) 

?dbFD

ls(DBFD)

DBFD$FRic
DBFD$FEve
DBFD$FDiv
DBFD$RaoQ
DBFD$FDis

write.csv(DBFD$FRic,"FRic_result.csv")
write.csv(DBFD$FEve,"FEve_result.csv")
write.csv(DBFD$FDiv,"FDiv_result.csv")
write.csv(DBFD$RaoQ,"RaoQ_result.csv")
write.csv(DBFD$FDis,"FDis_result.csv")

DBFD$nbsp #vector listing the number of species in each community
DBFD$qual.FRic #quality of the reduced-space representation required to compute FRic and FDiv.
DBFD$sing.sp #vector listing the number of functionally singular species in each community. If all species are functionally different, sing.sp will be identical to nbsp.
DBFD$x.axes #PCoA axes. Only returned if print.pco is TRUE.
DBFD$x.values #eigenvalues from the PCoA. Only returned if print.pco is TRUE.

plot(DBFD$nbsp, DBFD$FRic) #Volume funcional (FRic) satura com o aumento do n?mero de esp?cies (Mouchet et al. 2010, Funct. Ecol.)
plot(DBFD$nbsp, DBFD$FEve) #FEve em geral n?o ? relacionado com a riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)
plot(DBFD$nbsp, DBFD$FDiv) #FDiv em geral n?o ? relacionado com a riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)
plot(DBFD$nbsp, DBFD$RaoQ)  #RaoQ em geral ? pouco influenciado pela riqueza de esp?cies (Mouchet et al. 2010, Funct. Ecol.)porque seus valores s?o m?dias
plot(DBFD$nbsp, DBFD$FDis)
index_all <- data.frame(DBFD$FRic, DBFD$FEve, DBFD$FDiv, DBFD$RaoQ,
           DBFD$FDis)

index_all
plot(index_all)




