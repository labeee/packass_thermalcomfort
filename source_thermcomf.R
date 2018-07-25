#---Script que incorpora as fuções desenvolvidos no ambiente R.
#---O usuário deve clicar em "Source" para ativar as funções.
#---As funções ficarão disponibilizadas no "Global Environment" do R Studio e 
#---podem ser utilizadas através do "Console".

#---Set the working directory to this file path, to ensure the correct linking
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("f(x) comfortvectorfull.R")
source("f(x) comfortvector.R")
source("f(x) comfortindices.R")
source("f(x) pmvboth.R")
source("f(x) pmvelev.R")
source("f(x) pmvlow.R")
source("f(x) setashrae.R")
source("f(x) optemplimits.R")