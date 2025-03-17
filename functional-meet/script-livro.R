install.packages("FD") 
install.packages("ade4")
install.packages("ecodados")
devtools::install_github("paternogbc/ecodados")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("tidyverse")
install.packages("picante")
install.packages("vegan")
install.packages("SYNCSA")
install.packages("GGally")
install.packages("betapart")
install.packages("nlme")
install.packages("ape")
install.packages("TPD")
install.packages("cati")
install.packages("kableExtra")

## Pacotes 
library(FD)
library(ade4)
library(ecodados)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(picante)
library(vegan)
library(SYNCSA)
library(GGally)
library(betapart)
library(nlme)
library(ape)
library(TPD)
library(cati)
library(kableExtra)

## Dados e funções
comun_fren_dat <- ecodados::fundiv_frenette2012a_comu
ambie_fren_dat <- ecodados::fundiv_frenette2012a_amb
trait_fren_dat <- ecodados::fundiv_frenette2012a_trait
trait_dat      <- ecodados::fundiv_barbaro2009a_trait
comun_dat      <- ecodados::fundiv_barbaro2009a_comu
ambie_dat      <- ecodados::fundiv_barbaro2009a_amb
trait_baselga  <- ecodados::trait_baselga
comm_baselga   <- ecodados::comm_baselga
anuros_comm    <- ecodados::anuros_comm
traits         <- ecodados::traits
env            <- ecodados::env
# ecodados::wITV # funtion: wITV

## PCoA dos atributos contínuos

# 1. Padronização dos dados
trait_pad <- decostand(trait_fren_dat, "standardize")
?decostand

euclid_dis <- vegdist(trait_pad, "euclidean")
?vegdist

# 2. PCoA
# Resultados são idênticos aos resultados de uma PCA.
pcoa_traits_cont <- pcoa(euclid_dis, correction = "cailliez") 

# 3. Exportandos dados para gráfico
# Ao usar '[,1:2]' você irá selecionar os dois primeiros eixos.
eixos_cont <- as.data.frame(pcoa_traits_cont$vectors[, 1:2]) 

# 4. Gráfico de ordenação 
plot_trait_cont <- ggplot(eixos_cont, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(pch = 21, size = 4, color = "black", alpha = 0.7, fill = "red2") + 
  geom_text_repel(aes(Axis.1, Axis.2, label = rownames(eixos_cont))) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "PCO 1", y = "PCO 2", title = "Dados contínuos") + 
  tema_livro()
plot_trait_cont

## PCoA dos atributos categóricos

# 1. Selecionar somente os atributos categóricos
trait_cat <- trait_dat %>% 
  dplyr::select_if(is.character)

# 2. Calcular a distância de Gower
dist_categ <- gowdis(trait_cat)

# 3. PCoA da matriz de distância funcional (Gower)
pcoa_traits_cat <- pcoa(dist_categ, correction = "cailliez")
?pcoa
# 4. Exportar dados (escores) para ggplot
eixos_cat <- as.data.frame(pcoa_traits_cat$vectors[,1:2]) # Selecionar os dois primeiros eixos

# 5. Gráfico de ordenação
plot_trait_cat <- ggplot(eixos_cat, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(pch = 21, size = 4, alpha = 0.7, color = "black", fill = "cyan4") + 
  geom_text_repel(aes(Axis.1, Axis.2, label = rownames(eixos_cat))) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "PCO 1", y = "PCO 2", title = "Dados categóricos") + 
  tema_livro()
plot_trait_cat
