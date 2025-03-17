library(ade4)
library(ape)
library(picante)
library(vegan)
library(FD)
library(SYNCSA)

#abrindo planilhas método clipboard

atributos<-read.table('clipboard',header=T) #row.names=considera a primeira linha como label (legenda e não atributo)
atributos
comunidade<-read.table('clipboard',header=T)#riqueza/abundância de sp
comunidade

#abrindo planilhas método read.table

atributos<-read.table("atributos.txt", h=T)
atributos
comunidade<-read.table("comunidade.txt", h=T)
comunidade

#logartimizando dados de abundância
#?decostand
log(comunidade+1)->comunidadelog # quando os dados são quantitativos na mesma escala de mensuração
comunidadelog

#matriz de distância de Gower para dados mistos
traits<-gowdis(atributos, asym.bin = NULL, ord = c("classic"))
traits

#Criando matriz de distância de Gower 'generalizado' (Pavoine et a. 2009)	

#atividade-AT	
#microhabitat-MH	
#vocalization-VO 	
#reproductive-R	
#Total- TL
#Mouth-M	
#Forearm-F	
#Femur-FE


#a primeira coluna são as espécies, começar com a segunda, que é a primeira característica
#AT:
at.trait<-atributos[,1]#deve transformar todas as colunas que são fatores em data.frame pra dar certo na função da pavoine
head(at.trait)
class(at.trait)
at.trait<-data.frame(rownames(atributos),atributos[,1]) 
rownames(at.trait)<-at.trait[,1]#row.names- d? nome para as linhas
at.trait<-at.trait[,-1]
at.trait=data.frame(AT=at.trait)
rownames(at.trait)=atributos[,1]
class(at.trait)
fix(at.trait)#mostra o que modificou na planilha

#MH
mh.trait<-atributos[,2]
head(mh.trait)
class(mh.trait)
mh.trait<-data.frame(rownames(atributos),atributos[,2]) 
rownames(mh.trait)<-mh.trait[,1]
mh.trait<-mh.trait[,-1]
mh.trait=data.frame(MH=mh.trait)
rownames(mh.trait)=atributos[,2]
class(mh.trait)
fix(mh.trait)

#VO
vo.trait<-atributos[,3]
head(vo.trait)
class(vo.trait)
vo.trait<-data.frame(rownames(atributos),atributos[,3]) 
rownames(vo.trait)<-vo.trait[,1]
vo.trait<-vo.trait[,-1]
vo.trait=data.frame(VO=vo.trait)
rownames(vo.trait)=atributos[,3]
class(vo.trait)
fix(vo.trait)

#R
r.trait<-atributos[,4]
head(r.trait)
class(r.trait)
r.trait<-data.frame(rownames(atributos),atributos[,4]) 
rownames(r.trait)<-r.trait[,1]
r.trait<-r.trait[,-1]
r.trait=data.frame(R=r.trait)
rownames(r.trait)=atributos[,4]
class(r.trait)
fix(r.trait)

#TL
tl.trait<-atributos[,5]
head(tl.trait)
class(tl.trait)
tl.trait<-data.frame(rownames(atributos),atributos[,5]) 
rownames(tl.trait)<-tl.trait[,1]
tl.trait<-tl.trait[,-1]
tl.trait=data.frame(TL=tl.trait)
rownames(tl.trait)=atributos[,5]
class(tl.trait)
fix(tl.trait)

#M
m.trait<-atributos[,6]
head(m.trait)
class(m.trait)
m.trait<-data.frame(rownames(atributos),atributos[,6]) 
rownames(m.trait)<-m.trait[,1]
m.trait<-m.trait[,-1]
m.trait=data.frame(M=m.trait)
rownames(m.trait)=atributos[,6]
class(m.trait)
fix(m.trait)

#F
f.trait<-atributos[,7]
head(f.trait)
class(f.trait)
f.trait<-data.frame(rownames(atributos),atributos[,7])
rownames(f.trait)<-f.trait[,1]
f.trait<-f.trait[,-1]
f.trait=data.frame(F=f.trait)
rownames(f.trait)=atributos[,7]
class(f.trait)
fix(f.trait)


#FE
fe.trait<-atributos[,8]
head(fe.trait)
class(fe.trait)
fe.trait<-data.frame(rownames(atributos),atributos[,8]) 
rownames(fe.trait)<-fe.trait[,1]
fe.trait<-fe.trait[,-1]
fe.trait=data.frame(FE=fe.trait)
rownames(fe.trait)=atributos[,8]
class(fe.trait)
fix(fe.trait)


#Fazendo a lista de variáveis:

ktab<-ktab.list.df(list(at.trait,mh.trait,vo.trait,r.trait,tl.trait,m.trait,f.trait,fe.trait))

#Making the distance matrix according Pavoine et al. (2009):
p.dist<-dist.ktab(ktab,type=c("N","N","N","N","Q","Q","Q","Q"),option=c("scaledBYrange"))#matriz de distância feita com a medida de pavoine, "N" quer dizer variáveis nominais
p.dist

#Calcular índices de DF:
DBFD<-dbFD(traits,comunidadelog, m=8, print.pco = T, calc.CWM= FALSE, corr = c("sqrt")) #A função dbFD trabalha com matrizes de distâncias do tipo euclidiano. Como "dist.func" foi criada a partir do método de Pavoine et al. 2008 ("Gower generalizado"), então aqui foi necessário fazer uma correção da matriz de distância para transformá-la em euclidiana. Ver outras opções de correção em "Details" da ajuda da função.
DBFD


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




