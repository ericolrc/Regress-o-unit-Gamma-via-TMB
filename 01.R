# Carregando pacotes e dados
library(simplexreg)
data("sdac", package = "simplexreg")
head(sdac)

# Configurações do diretório de trabalho e compilação do modelo em C++
library(TMB)
setwd("C:\\Users\\Estudante\\Desktop\\Iniciação Científica\\Gama Unitária\\Aplicação final")

compile("gama_unitaria.cpp")
dyn.load(dynlib("gama_unitaria"))

# Carregando a matriz do modelo e preparando dados
matriz_modelo = read.table("C:\\Users\\Estudante\\Desktop\\Iniciação Científica\\Gama Unitária\\Aplicação final\\MatrizModeloGU.txt",header = F)
X = as.matrix(matriz_modelo[,-1])
Y = sdac$rcd
dados_gu = list(Y = Y, X = X)

# Definindo parâmetros iniciais para o modelo
parametros = list(beta = c(rep(0, 3)), logphi = 0)

# Ajuste do modelo com TMB
gu_TMB = MakeADFun(dados_gu, parametros, DLL = "gama_unitaria", hessian = TRUE, silent = TRUE)
opt = nlminb(start = gu_TMB$par, obj = gu_TMB$fn, gr = gu_TMB$gr)
rep = sdreport(gu_TMB)

# Exibindo resultados do ajuste
round(summary(rep, "fixed", p.value = TRUE), 4)
round(summary(rep, "report", p.value = TRUE), 4)

# Cálculo do BIC e AIC
BIC_gu = -2 * (-197.0202) + 3 * log(239) 
AIC_gu = -2 * (-197.0202) + 3 * 2

# Predições para o modelo GU com função logito
Xgama = 1.00263938 + sdac$ageadj * 0.0165 + sdac$chemo * 0.2441
gu_predict = exp(Xgama) / (1 + exp(Xgama))

# Cálculos intermediários baseados nas predições
yt = sdac$rcd
k = 3
n = length(yt)

phi_hat = exp(2.3683)
mu_hat = gu_predict
alpha_hat = (mu_hat^(1/phi_hat)) / (1 - mu_hat^(1/phi_hat))
mu_star = alpha_hat / mu_hat^((1/phi_hat) + 1)
y_star = (alpha_hat^2 * log(yt)) / (phi_hat * mu_hat^((1/phi_hat) + 1))
a_hat = alpha_hat / mu_hat^((1/phi_hat) + 1)
g_l = 1 / (mu_hat * (1 - mu_hat))
Wbb = diag(((a_hat^2) / phi_hat) * (1 / g_l)^2)
Wbb_12 = sqrt(Wbb)
s_hat = mu_star + y_star
H = Wbb_12 %*% X %*% solve(t(X) %*% Wbb %*% X) %*% t(X) %*% Wbb_12
ht = diag(H)

# Função de densidade e distribuição da GU
dGU = function(y, mu, phi, log = FALSE) {
  alpha = mu^(1/phi) / (1 - mu^(1/phi))
  fy = (alpha^(phi) / gamma(phi)) * y^(alpha - 1) * log(1/y)^(phi - 1)
  logfy = (phi * log(alpha)) - log(gamma(phi)) + (alpha - 1) * log(y) + (phi - 1) * log(log(1/y))
  if (log) { return(log(fy)) } else { return(fy) }
}

pGU = function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  value = numeric(length(q))
  for (i in 1:length(q)) {
    if (lower.tail) {
      value[i] = integrate(function(x) dGU(x, mu[i], phi), 0, q[i])$value
    } else {
      value[i] = 1 - integrate(function(x) dGU(x, mu[i], phi), 0, q[i])$value
    }
  }
  if (log.p) {
    return(log(value))
  } else {
    return(value)
  }
}

# Cálculo dos resíduos quantílicos
N = nrow(sdac)
resq = numeric(N)
for (i in 1:N) {
  resq[i] = qnorm(pGU(yt[i], mu_hat[i], 2.3683)) / sqrt(1 - ht)
}
dados_resq = data.frame(resq = resq)

# Gráficos dos resíduos e análise de influência
library(ggplot2)
library(patchwork)
library(car)
library(qqplotr)

# Gráfico quantil-quantil
g1 = ggplot(dados_resq, aes(sample = resq)) +
  stat_qq_band(bandType = "pointwise", fill = "gray", alpha = 0.8) +
  stat_qq_line(colour = "red") +
  stat_qq_point() +
  theme_minimal() +
  labs(x = "Quantis da Normal padrão", y = "Resíduo Quantílico", title = "")

# Resíduo ponderado padronizado
r_pond_pd = s_hat / sqrt(((mu_star^2) / phi_hat) * (1 - ht))

# Gráfico: Resíduo ponderado padronizado X Índice
g2 = ggplot() +
  geom_point(aes(y = r_pond_pd, x = 1:length(yt))) +
  labs(x = "Índice das observações", y = "Resíduo ponderado padronizado", title = "") +
  ylim(-4, 4) +
  geom_hline(aes(yintercept = 0), col = "black", linetype = 2) +
  geom_hline(aes(yintercept = 3), col = "red", linetype = 2) +
  geom_hline(aes(yintercept = -3), col = "red", linetype = 2) +
  theme_minimal()

# Gráfico: Resíduo ponderado padronizado X Preditor linear
g3 = ggplot() +
  geom_point(aes(x = mu_hat, y = r_pond_pd)) +
  labs(x = "Preditor linear", y = "Resíduo ponderado padronizado", title = "") +
  geom_hline(aes(yintercept = 0), col = "black", linetype = 2) +
  ylim(-4, 4) +
  geom_hline(aes(yintercept = 3), col = "red", linetype = 2) +
  geom_hline(aes(yintercept = -3), col = "red", linetype = 2) +
  theme_minimal()

# Gráfico: Distância de Cook
DC2 = (r_pond_pd^2) * (ht / (k * (1 - ht)))
g4 = ggplot() +
  geom_point(aes(y = DC2, x = 1:length(yt))) +
  labs(x = "Índice das observações", y = "Distância de Cook", title = "") +
  ylim(0, 0.30) +
  geom_hline(aes(yintercept = 2 * k / n), col = "red", linetype = 2) +
  theme_minimal()

# Gráfico: Alavanca
g5 = ggplot() +
  geom_point(aes(y = ht, x = 1:length(yt))) +
  labs(x = "Índice das observações", y = "htt", title = "") +
  geom_hline(aes(yintercept = 2 * k / n), col = "red", linetype = 2) +
  ylim(0, 0.30) +
  theme_minimal()

# Exibir gráficos combinados
(g4 + g2) / (g5 + g1)

################################################################################################
# Carregando as bibliotecas necessárias
library("simplexreg")
library(ggplot2)
library(xtable)
library(patchwork)

# Carregando os dados
data("sdac", package = "simplexreg")
head(sdac, 5)

## Análise descritiva
desc = function(x){
  round(rbind(min(x), mean(x), median(x),
              max(x), sd(x), var(x)), 4)
}

# Exibindo a descrição da variável 'rcd'
desc(sdac$rcd)
print(xtable(desc(sdac$rcd)), format.args = list(decimal.mark = ","))

# Gráfico 1: Histograma da taxa de recuperação de células CD34+
g1 = ggplot(sdac, aes(x = rcd, y = after_stat(density))) + 
  geom_histogram(colour = "black", fill = "gray", stat = "bin") +
  geom_vline(aes(xintercept = mean(rcd)), color = "red", linetype = "dashed", size = 1) +
  labs(x = "Taxa de recuperação de células CD34+", y = "Densidade") +
  theme_bw()

# Gráfico 2: Boxplot da taxa de recuperação de células CD34+
g2 = ggplot(sdac, aes(x = rcd)) + 
  geom_boxplot(color = "black", fill = "gray") +
  labs(x = "Taxa de recuperação de células CD34+", y = "") +
  theme_bw()

# Combinando os gráficos
g1/g2
