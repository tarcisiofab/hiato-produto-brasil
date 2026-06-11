# =============================================================================
# HIATO DO PRODUTO — BCB RI Jun/2024
# GRUPO I  — 7 Métodos Univariados (I a VI + VII Beveridge-Nelson KMW)
# GRUPO II — Multivariado I: Função de Produção (Combinação Simples)
# GRUPO II — Multivariado II: Areosa (2008)
# GRUPO II — Multivariado III: CBO (Shackleton, 2018)
# GRUPO II — Multivariado IV: Jarocinski & Lenza (2018)
# GRUPO II — Multivariado V: PCA
#
# Série PIB:       BCB/SGS 22109 (índice dessaz., trimestral)
# Série NUCI:      BCB/SGS 28561 (%, dessaz., trimestral)
# Série Desemp.:   BCB/SGS 24369 (PNADC %, dessaz., trimestral)
# Autor: Tarcísio | 2026
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite); library(tidyverse)
  library(mFilter);  library(KFAS)
  library(scales)
})

# ── FUNÇÃO AUXILIAR ───────────────────────────────────────────────────────────
get_sgs <- function(serie, inicio = "01/01/1996") {
  url <- sprintf(
    "https://api.bcb.gov.br/dados/serie/bcdata.sgs.%d/dados?formato=json&dataInicial=%s",
    serie, inicio)
  raw <- jsonlite::fromJSON(url)
  tibble(data  = as.Date(raw$data, "%d/%m/%Y"),
         valor = as.numeric(gsub(",", ".", raw$valor))) %>%
    filter(!is.na(valor))
}

# ════════════════════════════════════════════════════════════════════════════
# DADOS
# ════════════════════════════════════════════════════════════════════════════
cat("Baixando dados...\n")
pib_raw   <- get_sgs(22109)
nuci_raw  <- get_sgs(28561, "01/01/2001")
desemp_raw <- get_sgs(24369, "01/01/2012")
caged_raw <- get_sgs(28763, "01/01/2010")   # Novo Caged (MTE): admissões - demissões (desde 2010)
cat(sprintf("  PIB: %s–%s | NUCI: %s–%s | Desemp: %s–%s | Caged: %s–%s\n",
    min(pib_raw$data), max(pib_raw$data),
    min(nuci_raw$data), max(nuci_raw$data),
    min(desemp_raw$data), max(desemp_raw$data),
    min(caged_raw$data), max(caged_raw$data)))

# Série PIB para os univariados
pib <- pib_raw %>%
  arrange(data) %>%
  mutate(log_pib = log(valor), t = row_number())
n      <- nrow(pib)
y      <- pib$log_pib
pib_ts <- ts(y, start = c(1996, 1), frequency = 4)

# ════════════════════════════════════════════════════════════════════════════
# GRUPO I — MÉTODOS UNIVARIADOS
# ════════════════════════════════════════════════════════════════════════════

# ── I. Tendência Quadrática com Quebras ──────────────────────────────────────
cat("\nI.  Tendência Quadrática com Quebras (BCB)...\n")
quebras_datas <- as.Date(c("2000-01-01","2008-10-01","2013-07-01","2020-04-01"))
idx_q <- sapply(quebras_datas, function(d) which(pib$data >= d)[1])
k <- length(idx_q)
ini_s <- c(1, idx_q); fim_s <- c(idx_q - 1, n)
X <- matrix(0, n, 3*(k+1))
for (j in 1:(k+1)) {
  tl <- pib$t[ini_s[j]:fim_s[j]]
  X[ini_s[j]:fim_s[j], 3*(j-1)+1] <- 1
  X[ini_s[j]:fim_s[j], 3*(j-1)+2] <- tl
  X[ini_s[j]:fim_s[j], 3*(j-1)+3] <- tl^2
}
pib$hiato_TQ <- residuals(lm(y ~ X - 1)) * 100
cat(sprintf("    Hiato atual: %+.2f p.p.\n", tail(pib$hiato_TQ, 1)))

# ── II. Tendência Não-Paramétrica (LOESS) ────────────────────────────────────
cat("II.  Tendência Não-Paramétrica (LOESS, bwidth=0.4)...\n")
loess_bcb <- function(y, bwidth = 0.4) {
  T_ <- length(y); kk <- floor((T_*bwidth - 0.5)/2); tend <- numeric(T_)
  for (i in seq_len(T_)) {
    t_m <- max(1,i-kk); t_p <- min(T_,i+kk); idx <- t_m:t_p
    Delta <- max(t_p-i, i-t_m, 1)
    w <- (1 - pmin(abs(idx-i)/Delta, 1)^3)^3
    xl <- idx - i; Xl <- cbind(1, xl, xl^2); W <- diag(w)
    beta <- tryCatch(solve(t(Xl)%*%W%*%Xl, t(Xl)%*%W%*%y[idx]),
                     error = function(e) MASS::ginv(t(Xl)%*%W%*%Xl) %*% (t(Xl)%*%W%*%y[idx]))
    tend[i] <- beta[1]
  }
  tend
}
pib$hiato_NP <- (y - loess_bcb(y)) * 100
cat(sprintf("    Hiato atual: %+.2f p.p.\n", tail(pib$hiato_NP, 1)))

# ── III. Hodrick-Prescott (λ=1600) ───────────────────────────────────────────
cat("III. Hodrick-Prescott (λ=1600)...\n")
hp <- mFilter::hpfilter(pib_ts, freq = 1600, type = "lambda")
pib$hiato_HP <- as.numeric(hp$cycle) * 100
cat(sprintf("    Hiato atual: %+.2f p.p.\n", tail(pib$hiato_HP, 1)))

# ── IV. Tendência ℓ₁ (Kim et al., 2009) ─────────────────────────────────────
cat("IV.  Tendência ℓ₁ (Kim et al., 2009)...\n")
soft_thr <- function(x, k) sign(x) * pmax(abs(x) - k, 0)
l1_filter <- function(y, lambda, rho_mult = NULL, max_iter = 20000, tol = 1e-9) {
  T_ <- length(y)
  D  <- matrix(0, T_-2, T_)
  for (i in 1:(T_-2)) { D[i,i] <- 1; D[i,i+1] <- -2; D[i,i+2] <- 1 }
  if (is.null(rho_mult)) rho_mult <- lambda / (max(abs(D %*% y)) / 5)
  rho <- rho_mult * lambda
  Ach <- chol(diag(T_) + rho * t(D) %*% D)
  mu <- y; z <- rep(0, T_-2); u <- rep(0, T_-2)
  for (iter in 1:max_iter) {
    mu_old <- mu
    mu  <- backsolve(Ach, forwardsolve(t(Ach), y + rho * t(D) %*% (z - u)))
    Dmu <- D %*% mu; z <- soft_thr(Dmu + u, lambda/rho); u <- u + Dmu - z
    if (max(abs(mu - mu_old)) < tol) { cat(sprintf("(conv.%d) ", iter)); break }
  }
  mu
}
D_mat <- matrix(0, n-2, n)
for (i in 1:(n-2)) { D_mat[i,i] <- 1; D_mat[i,i+1] <- -2; D_mat[i,i+2] <- 1 }
lam_max <- max(abs(solve(D_mat %*% t(D_mat), D_mat %*% y)))
lam_l1  <- (0.5^5) * lam_max
rm_l1   <- lam_l1 / (max(abs(D_mat %*% y)) / 5)
pib$hiato_L1 <- (y - as.numeric(l1_filter(y, lam_l1, rho_mult = rm_l1))) * 100
cat(sprintf("    Hiato atual: %+.2f p.p.\n", tail(pib$hiato_L1, 1)))

# ── V. HP Modificada (Andrle/BCB AR(2)) ──────────────────────────────────────
cat("V.   HP Modificada (Andrle/BCB AR(2))...\n")
g_ss <- mean(diff(y))
build_hpm <- function(theta) {
  sg <- exp(theta[1]); sc <- exp(theta[2])
  rho <- plogis(theta[3]); r <- plogis(theta[4])*0.99; om <- plogis(theta[5])*pi
  phi1 <- 2*r*cos(om); phi2 <- -r^2; cd <- g_ss*(1-rho)
  Tt <- matrix(c(1,1,0,0,0, 0,rho,0,0,cd, 0,0,phi1,phi2,0,
                 0,0,1,0,0, 0,0,0,0,1), 5, 5, byrow = TRUE)
  Rt <- matrix(c(0,0,1,0,0,1,0,0,0,0), 5, 2, byrow = TRUE)
  Qt <- diag(c(sg^2, sc^2))
  Zt <- matrix(c(1,0,1,0,0), 1, 5)
  P1inf <- matrix(0,5,5); P1inf[1,1] <- 1
  P1 <- matrix(0,5,5)
  P1[2,2] <- sg^2 / max(1-rho^2, 1e-8)
  P1[3,3] <- sc^2 / max(1-phi1^2-phi2^2, 1e-8); P1[4,4] <- P1[3,3]
  SSModel(pib_ts ~ -1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, P1inf=P1inf, P1=P1),
          H = matrix(1e-8))
}
nll_hpm <- function(theta) tryCatch(-as.numeric(logLik(build_hpm(theta))), error=function(e) 1e10)
sd_dy <- sd(diff(y))
th0 <- c(log(sd_dy/20), log(sd_dy*2), qlogis(0.8), qlogis(0.9/0.99), qlogis(1/8))
cat("    Otimizando MLE...")
opt  <- optim(th0, nll_hpm, method="Nelder-Mead", control=list(maxit=10000,reltol=1e-12))
opt2 <- tryCatch(optim(opt$par, nll_hpm, method="BFGS",
                       control=list(maxit=3000,reltol=1e-12)), error=function(e) opt)
if (is.finite(opt2$value) && opt2$value < opt$value) opt <- opt2
ks_hpm <- KFS(build_hpm(opt$par), smoothing="state")
pib$hiato_HPM <- as.numeric(ks_hpm$alphahat[,3]) * 100
if (sd(pib$hiato_HPM, na.rm=TRUE) < 0.1)
  pib$hiato_HPM <- (y - as.numeric(ks_hpm$alphahat[,1])) * 100
cat(sprintf(" OK | Hiato atual: %+.2f p.p.\n", tail(pib$hiato_HPM, 1)))

# ── VI. Band-Pass Christiano-Fitzgerald (8-32 trimestres) ────────────────────
cat("VI.  Band-Pass Christiano-Fitzgerald (8-32 trimestres)...\n")
cf_filter <- function(y, pl=8, pu=32) {
  # Filtro CF ASSIMÉTRICO (Christiano & Fitzgerald, 2003) — versão random walk
  # Usa toda a amostra disponível em cada t (janela assimétrica). Coeficientes
  # de borda B~_k = -B_0/2 - Σ_{j=1}^{k-1} B_j garantem Σ pesos = 0 sem zerar
  # os extremos — o último ponto da amostra (hiato corrente) é estimado, não nulo.
  T_ <- length(y); wl <- 2*pi/pu; wu <- 2*pi/pl
  # Remoção de drift (CF 2003 recomendam): torna a anulação de tendência
  # linear EXATA também nas bordas da janela assimétrica
  g_drift <- (y[T_] - y[1]) / (T_ - 1)
  y <- y - g_drift * (0:(T_-1))
  B  <- numeric(T_); B[1] <- (wu - wl)/pi              # B_0
  for (k in 1:(T_-1)) B[k+1] <- (sin(wu*k) - sin(wl*k))/(pi*k)
  Btil <- function(k) {                                 # B~(0) = -B_0/2
    -B[1]/2 - (if (k >= 2) sum(B[2:k]) else 0)
  }
  cycle <- numeric(T_)
  for (t in 1:T_) {
    w    <- numeric(T_)
    w[t] <- B[1]
    # Leads: interiores j=1..(T-t-1); borda em y_T com B~(T-t)
    if (T_ - t >= 1) {
      if (T_ - t >= 2) for (j in 1:(T_-t-1)) w[t+j] <- w[t+j] + B[j+1]
      w[T_] <- w[T_] + Btil(T_ - t)
    } else {
      w[T_] <- w[T_] + Btil(0)      # t = T: meio peso no próprio ponto
    }
    # Lags: interiores j=1..(t-2); borda em y_1 com B~(t-1)
    if (t - 1 >= 1) {
      if (t - 1 >= 2) for (j in 1:(t-2)) w[t-j] <- w[t-j] + B[j+1]
      w[1] <- w[1] + Btil(t - 1)
    } else {
      w[1] <- w[1] + Btil(0)        # t = 1
    }
    cycle[t] <- sum(w * y)
  }
  cycle
}
pib$hiato_BP <- cf_filter(y) * 100
cat(sprintf("    Hiato atual: %+.2f p.p.\n", tail(pib$hiato_BP, 1)))

# ── VII. Beveridge-Nelson modificado (Kamber, Morley & Wong, 2018) ──────────
# Tendência BN = limite da esperança condicional do PIB em horizonte longo.
# AR(12) sobre o crescimento (constante = média), reparametrizado em
# rho = soma dos coef. AR. Razão sinal-ruído delta = sig2_dtau/sig2_dy
# imposta MAIS BAIXA (grid em delta <= 1; amplitude do hiato é crescente em
# delta, logo o ótimo fica na borda delta=1, que pina rho* por raiz única).
# Correção de heteroscedasticidade na pandemia: WLS em 2 etapas com
# variância relativa calibrada pelos resíduos (Morley et al., 2023).
cat("VII. Beveridge-Nelson modificado (Kamber et al., 2018)...\n")

bn_solve <- function(x, rho, p, w) {
  T_ <- length(x); idx <- (p+1):T_
  yreg <- x[idx] - rho * x[idx-1]
  Xreg <- sapply(1:(p-1), function(j) x[idx-j] - x[idx-j-1])
  fit  <- lm.wfit(Xreg, yreg, w)
  psi  <- fit$coefficients
  phi  <- numeric(p)
  phi[1] <- rho + psi[1]
  for (j in 2:(p-1)) phi[j] <- psi[j] - psi[j-1]
  phi[p] <- -psi[p-1]
  F_ <- rbind(phi, cbind(diag(p-1), 0))
  if (max(Mod(eigen(F_, only.values=TRUE)$values)) >= 0.999) return(NULL)
  M  <- F_ %*% solve(diag(p) - F_)
  cyc <- rep(NA_real_, T_)
  for (t in p:T_) cyc[t] <- -as.numeric(M[1,] %*% x[t:(t-p+1)])
  list(cycle = cyc, e = fit$residuals,
       sig2e = sum(w * fit$residuals^2) / sum(w))
}

bn_kmw <- function(yl, datas, p = 12) {
  dy <- diff(yl * 100); x <- dy - mean(dy); d_dy <- datas[-1]
  # Pesos COVID (2 etapas)
  idx_reg <- (p+1):length(x)
  w <- rep(1, length(idx_reg))
  covid_t <- d_dy[idx_reg] >= as.Date("2020-01-01") &
             d_dy[idx_reg] <= as.Date("2021-06-30")
  est0 <- bn_solve(x, 0.5, p, w)
  if (!is.null(est0) && any(covid_t)) {
    sig_pre <- sqrt(mean(est0$e[!covid_t]^2))
    kappa   <- pmax(1, (est0$e / sig_pre)^2)
    w[covid_t] <- 1 / kappa[covid_t]
  }
  # Grid sobre rho na região de razão sinal-ruído BAIXA (delta < 1):
  # rho < 0 gera ciclo PRÓ-cíclico (KMW); rho > 0 reproduz a patologia
  # anticíclica do BN clássico. Maximiza-se a amplitude do hiato no grid
  # (RI, Ap. 1g); na amostra brasileira o máximo ocorre na borda inferior,
  # com delta implícito ~0,28 — próximo do delta=0,24 de KMW para os EUA.
  grid_rho <- seq(-0.97, -0.05, by = 0.02)
  amps <- sapply(grid_rho, function(r) {
    s <- bn_solve(x, r, p, w)
    if (is.null(s)) NA else sd(s$cycle, na.rm = TRUE)
  })
  rho_star <- grid_rho[which.max(amps)]
  s <- bn_solve(x, rho_star, p, w)
  delta_imp <- s$sig2e / ((1 - rho_star)^2 * var(x))
  list(hiato = c(NA, s$cycle), rho = rho_star, delta = delta_imp)
}

bn_res <- bn_kmw(y, pib$data, p = 12)
pib$hiato_BN <- bn_res$hiato
cat(sprintf("    rho* = %.2f | delta implícito = %.3f | Hiato atual: %+.2f p.p.\n",
            bn_res$rho, bn_res$delta, tail(na.omit(pib$hiato_BN), 1)))
cat("    [Ref.: pró-cíclico (corr>0 c/ HP); BCB Tabela 1 dava -0,30 em 2024T1 no vintage de jun/24]\n")

# ════════════════════════════════════════════════════════════════════════════
# GRUPO II — MULTIVARIADO I: FUNÇÃO DE PRODUÇÃO (COMBINAÇÃO SIMPLES)
# ŷ_t = (1-α)·ĉ_t − α·û_t   |  α=0,6 (trabalho), 1-α=0,4 (capital)
# ĉ_t = hiato NUCI (HP), û_t = hiato desemprego (HP)
# ════════════════════════════════════════════════════════════════════════════
cat("\nMI.  Função de Produção — Combinação Simples...\n")
alpha <- 0.6

df_multi <- pib_raw %>% rename(pib = valor) %>%
  inner_join(nuci_raw  %>% rename(nuci   = valor), by = "data") %>%
  inner_join(desemp_raw %>% rename(desemp = valor), by = "data") %>%
  arrange(data)

cat(sprintf("    Período multivariado: %s a %s (%d obs.)\n",
            min(df_multi$data), max(df_multi$data), nrow(df_multi)))

t_ini_multi <- c(year(min(df_multi$data)), quarter(min(df_multi$data)))

nuci_ts   <- ts(df_multi$nuci,   start = t_ini_multi, frequency = 4)
desemp_ts <- ts(df_multi$desemp, start = t_ini_multi, frequency = 4)

hp_nuci   <- mFilter::hpfilter(nuci_ts,   freq = 1600, type = "lambda")
hp_desemp <- mFilter::hpfilter(desemp_ts, freq = 1600, type = "lambda")

df_multi <- df_multi %>%
  mutate(
    hiato_nuci    = as.numeric(hp_nuci$cycle),
    hiato_desemp  = as.numeric(hp_desemp$cycle),
    hiato_FP      = (1-alpha) * hiato_nuci - alpha * hiato_desemp,
    contrib_cap   = (1-alpha) * hiato_nuci,
    contrib_trab  = -alpha    * hiato_desemp
  )

cat(sprintf("    Hiato atual: %+.2f p.p. | DP: %.2f p.p.\n",
            tail(df_multi$hiato_FP, 1), sd(df_multi$hiato_FP)))

# ════════════════════════════════════════════════════════════════════════════
# GRUPO II — MULTIVARIADO V: COMPONENTES PRINCIPAIS (PCA)
# BCB RI Jun/2024, Apêndice 2e:
# Séries: ex-tendência PIB (hiato HP), NUCI (nível, dessaz.),
#         desemprego PNADC (sinal invertido), Novo Caged (hiato HP)
# Procedimento: padronizar → SVD → 1º componente principal
# ════════════════════════════════════════════════════════════════════════════
cat("\nMV.  PCA — Componentes Principais (BCB RI Jun/2024)...\n")

# Novo Caged: SGS 28763 — saldo mensal de empregos formais
# Agregar para trimestral (soma dos meses) e calcular hiato via HP
caged_trim <- caged_raw %>%
  rename(caged = valor) %>%
  mutate(
    ano = year(data), tri = quarter(data)
  ) %>%
  group_by(ano, tri) %>%
  summarise(caged = sum(caged, na.rm = TRUE), .groups = "drop") %>%
  mutate(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3 + 1))) %>%
  select(data, caged) %>%
  arrange(data)

# Período comum: interseção de todas as séries
# NUCI: desde 2001 | Desemp (PNADC): desde 2012 | Caged Novo: desde 2020
# Período PCA: 2012T1 em diante (limitado pelo início do PNADC)
df_pca_raw <- pib_raw %>% rename(pib = valor) %>%
  inner_join(nuci_raw   %>% rename(nuci   = valor), by = "data") %>%
  inner_join(desemp_raw %>% rename(desemp = valor), by = "data") %>%
  inner_join(caged_trim,                             by = "data") %>%
  arrange(data)

cat(sprintf("    Período PCA: %s a %s (%d obs.)\n",
            min(df_pca_raw$data), max(df_pca_raw$data), nrow(df_pca_raw)))

# Caged 28763 cobre desde 2010 — período limitado pelo PNADC (2012T1)
if (nrow(df_pca_raw) < 8) {
  cat("    AVISO: poucos dados — verificar séries disponíveis\n")
}

t_ini_pca <- c(year(min(df_pca_raw$data)), quarter(min(df_pca_raw$data)))

# 1. Hiato do PIB via HP (já calculado globalmente, mas recortar ao período PCA)
pib_pca_ts <- ts(log(df_pca_raw$pib), start = t_ini_pca, frequency = 4)
hp_pib_pca <- mFilter::hpfilter(pib_pca_ts, freq = 1600, type = "lambda")
hiato_pib_pca <- as.numeric(hp_pib_pca$cycle) * 100

# 2. NUCI: entra em nível (dessaz.) — série estacionária por natureza
nuci_pca <- df_pca_raw$nuci

# 3. Desemprego: sinal invertido (desemprego alto = atividade baixa)
desemp_pca <- -df_pca_raw$desemp

# 4. Caged: hiato via HP
caged_pca_ts <- ts(df_pca_raw$caged, start = t_ini_pca, frequency = 4)
hp_caged <- mFilter::hpfilter(caged_pca_ts, freq = 1600, type = "lambda")
hiato_caged <- as.numeric(hp_caged$cycle)

# 5. Montar matriz e padronizar (z-score)
X_pca <- cbind(hiato_pib_pca, nuci_pca, desemp_pca, hiato_caged)
colnames(X_pca) <- c("Hiato PIB (HP)", "NUCI", "-Desemprego", "Hiato Caged (HP)")

X_std <- scale(X_pca)   # padronização: mean=0, sd=1

# 6. PCA via SVD (equivalente a prcomp com scale=FALSE após padronização manual)
pca_fit <- prcomp(X_std, center = FALSE, scale. = FALSE)

# Variância explicada
var_exp <- pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100
cat(sprintf("    Variância explicada: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%, PC4=%.1f%%\n",
            var_exp[1], var_exp[2], var_exp[3], var_exp[4]))
cat(sprintf("    BCB reporta: PC1 = 71,9%% (referência)\n"))

# 7. Primeiro componente principal = hiato PCA
# Verificar sinal: PC1 deve correlacionar positivamente com hiato HP do PIB
pc1 <- pca_fit$x[, 1]
if (cor(pc1, hiato_pib_pca) < 0) pc1 <- -pc1

# Escalar para p.p. comparáveis com outros hiatos
# (PC1 está em unidades de desvio-padrão — multiplicar pelo sd do hiato PIB HP)
pc1_scaled <- pc1 * sd(hiato_pib_pca) / sd(pc1)

df_pca_raw$hiato_PCA <- pc1_scaled

cat(sprintf("    Loadings do PC1: PIB=%.3f, NUCI=%.3f, -Desemp=%.3f, Caged=%.3f\n",
            pca_fit$rotation[1,1] * sign(cor(pc1, hiato_pib_pca)),
            pca_fit$rotation[2,1] * sign(cor(pc1, hiato_pib_pca)),
            pca_fit$rotation[3,1] * sign(cor(pc1, hiato_pib_pca)),
            pca_fit$rotation[4,1] * sign(cor(pc1, hiato_pib_pca))))
cat(sprintf("    Hiato PCA atual: %+.2f p.p. | DP: %.2f p.p.\n",
            tail(df_pca_raw$hiato_PCA, 1), sd(df_pca_raw$hiato_PCA)))

# Robustez: NUCI como hiato HP (em vez de nível) — diagnóstico apenas
hiato_nuci_rob <- as.numeric(mFilter::hpfilter(
  ts(df_pca_raw$nuci, start=t_ini_pca, frequency=4), freq=1600)$cycle)
pca_rob <- prcomp(scale(cbind(hiato_pib_pca, hiato_nuci_rob,
                               desemp_pca, hiato_caged)),
                  center=FALSE, scale.=FALSE)
var_rob <- pca_rob$sdev[1]^2 / sum(pca_rob$sdev^2) * 100
cat(sprintf("    [Robustez] PC1 com NUCI em hiato HP: %.1f%% (vs %.1f%% com NUCI em nível)\n",
            var_rob, var_exp[1]))

# ════════════════════════════════════════════════════════════════════════════
# GRUPO II — MULTIVARIADO II: AREOSA (2008)
# Implementação fiel ao BCB RI Jun/2024 e WP 172:
# Resolve simultaneamente três filtros HP (λ comum) com a restrição da
# função de produção ĥ_Y = β2·ĥ_C − β1·ĥ_U imposta em todo t.
#
# PROPOSIÇÃO 1 (solução analítica fechada; verificada em verif_areosa.R):
# As CPOs do problema acoplado fatoram em (I+λD'D); usando a identidade
# λSA = I−S, onde S=(I+λD'D)^{-1} é o suavizador HP, a solução exata é
#     ĥ_Y^Areosa = (1/κ)·ĥ_Y^FP + ((β1²+β2²)/κ)·ĥ_Y^HP,  κ = 1+β1²+β2²
# i.e., MÉDIA PONDERADA EXATA do hiato da função de produção simples (MII.I)
# e do hiato HP do próprio PIB. Com β1=0,6, β2=0,4: 0,658·FP + 0,342·HP.
# Verificação numérica: sistema 2T×2T exato vs fórmula, max|dif|=4,6e-11;
# corr(Areosa,FP)=0,979; max|Areosa−FP|=1,0 p.p.
#
# Corolário: sob λ comum, Areosa não adiciona informação além de FP e HP;
# como combinação convexa, jamais expande a faixa mín–máx do thick modeling.
# A divergência MII.I vs MII.II no Gráfico 2 do BCB decorre, portanto, dos
# dados internos (PNADC retropolada, Alves & Fasolo WP 400), não do método.
# Referência: Areosa (2008), BCB WP 172; BCB RI Jun/2024, Apêndice 2b.
# ════════════════════════════════════════════════════════════════════════════
cat("\nMII.2 Areosa (2008) — HP triplo com restrição CD (solução exata)...\n")

df_areosa <- df_multi %>% arrange(data)
t_ar <- c(year(min(df_areosa$data)), quarter(min(df_areosa$data)))

beta1 <- 0.6; beta2 <- 0.4
kap_ar <- 1 + beta1^2 + beta2^2

# Hiato HP do PIB no MESMO período amostral do problema (2012T2+)
hp_pib_ar <- mFilter::hpfilter(ts(log(df_areosa$pib), start=t_ar, frequency=4),
                               freq=1600, type="lambda")
hiato_hp_pib_ar <- as.numeric(hp_pib_ar$cycle) * 100

# Solução exata (Proposição 1)
df_areosa$hiato_AR <- (1/kap_ar) * df_multi$hiato_FP +
                      ((beta1^2 + beta2^2)/kap_ar) * hiato_hp_pib_ar

cat(sprintf("    Pesos: %.3f·FP + %.3f·HP (κ=%.2f)\n",
            1/kap_ar, (beta1^2+beta2^2)/kap_ar, kap_ar))
cat(sprintf("    Corr(Areosa, FP): %.4f | Corr(Areosa, HP-PIB): %.4f\n",
            cor(df_areosa$hiato_AR, df_multi$hiato_FP),
            cor(df_areosa$hiato_AR, hiato_hp_pib_ar)))
cat(sprintf("    Hiato Areosa atual: %+.2f p.p. | DP: %.2f p.p.\n",
            tail(df_areosa$hiato_AR, 1), sd(df_areosa$hiato_AR)))
cat(sprintf("    max|Areosa − FP| na amostra: %.2f p.p.\n",
            max(abs(df_areosa$hiato_AR - df_multi$hiato_FP))))

# ════════════════════════════════════════════════════════════════════════════
# GRUPO II — MULTIVARIADO IV: JAROCINSKI & LENZA (2018) — versão linearizada
# Modelo de fator dinâmico: hiato g_t como fator comum de 4 séries de atividade
# Versão simplificada sem volatilidade estocástica (h_t = f_t = 0)
# ainda captura a ideia central: hiato consistente com a curva de Phillips
#
# Observáveis: PIB (y¹), NUCI (y²), -desemprego (y³), Caged (y⁴)
#   y_t^n = b0^n·g_t + b1^n·g_{t-1} + w_t^n + eps_t^n
# Curva de Phillips (sem vol. estoc.):
#   (π_t - z_t) = ag0·g_t + ag1·g_{t-1} + ap·(π_{t-1}-z_{t-1}) + eps_t^π
# Tendência inflação:
#   z_t = z_{t-1} + eps_t^z
# Hiato: g_t = φ₁·g_{t-1} + φ₂·g_{t-2} + eta_t^g
# Tendências: w_t^n = d^n + w_{t-1}^n + eta_t^n
#
# Dados: PIB(22109) + NUCI(28561) + Desemp(24369) + Caged(28763)
#        + núcleo IPCA (proxy: IPCA serviços SGS 10844)
#        + Focus 12m (proxy expectativas: SGS 13522)
# ════════════════════════════════════════════════════════════════════════════
cat("\nMII.4 Jarocinski & Lenza (2018) — fator dinâmico (versão linearizada)...\n")

# Baixar dados adicionais para JL
cat("    Baixando núcleo IPCA (serviços, SGS 10844)... ")
ipca_serv_raw <- tryCatch(get_sgs(10844, "01/01/2012"), error=function(e) NULL)

cat("    Baixando Focus 12m (SGS 13522)... ")
focus_raw <- tryCatch(get_sgs(13522, "01/01/2012"), error=function(e) NULL)

if (is.null(ipca_serv_raw) || is.null(focus_raw)) {
  cat("\n    AVISO: séries de inflação não disponíveis — usando versão sem curva de Phillips\n")
  usar_phillips <- FALSE
} else {
  cat("OK\n")
  usar_phillips <- TRUE
}

# Agregar IPCA serviços para trimestral (variação acumulada no trimestre)
if (usar_phillips) {
  ipca_trim <- ipca_serv_raw %>%
    rename(ipca = valor) %>%
    mutate(ano=year(data), tri=quarter(data)) %>%
    group_by(ano, tri) %>%
    summarise(ipca = sum(ipca, na.rm=TRUE), .groups="drop") %>%
    mutate(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3+1))) %>%
    select(data, ipca)

  focus_trim <- focus_raw %>%
    rename(focus = valor) %>%
    mutate(ano=year(data), tri=quarter(data)) %>%
    group_by(ano, tri) %>%
    summarise(focus = mean(focus, na.rm=TRUE), .groups="drop") %>%
    mutate(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3+1))) %>%
    select(data, focus)
}

# ── Montar dataset JL ─────────────────────────────────────────────────────────
# Base: df_pca_raw (PIB + NUCI + desemp + Caged, 2012T1 em diante)
df_jl <- df_pca_raw %>% select(data, pib, nuci, desemp, caged)

if (usar_phillips) {
  df_jl <- df_jl %>%
    left_join(ipca_trim, by="data") %>%
    left_join(focus_trim, by="data")
}

df_jl <- df_jl %>% arrange(data) %>% filter(!is.na(pib))
T_jl  <- nrow(df_jl)
t_jl  <- c(year(min(df_jl$data)), quarter(min(df_jl$data)))
cat(sprintf("    Período JL: %s a %s (%d obs.)\n",
            min(df_jl$data), max(df_jl$data), T_jl))

# ── Pré-processamento das séries de atividade ─────────────────────────────────
# Padronizar para mesma escala (z-score) antes de entrar no fator dinâmico
# Restrição: eps_t^1 = 0 → hiato coincide com ciclo do PIB (BCB, Apêndice 2d)
y1 <- scale(as.numeric(mFilter::hpfilter(
  ts(log(df_jl$pib), start=t_jl, freq=4), freq=1600)$cycle))[,1]   # hiato PIB (HP, padronizado)
y2 <- scale(df_jl$nuci)[,1]                                          # NUCI
y3 <- scale(-df_jl$desemp)[,1]                                       # -desemprego
y4 <- scale(as.numeric(mFilter::hpfilter(
  ts(df_jl$caged, start=t_jl, freq=4), freq=1600)$cycle))[,1]       # hiato Caged (HP)

Y_obs_jl <- cbind(y1, y2, y3, y4)

# ── Modelo de espaço de estados (versão linearizada, sem vol. estocástica) ────
# Estado: [g_t, g_{t-1}, w_t^1, w_t^2, w_t^3, w_t^4]  (6×1)
# Observação: y_t^n = b0^n·g_t + b1^n·g_{t-1} + w_t^n
# Restrição BCB: eps_t^1 = 0 → y1 = g_t (hiato = fator)
#   → b0^1 = 1, b1^1 = 0, w_t^1 = 0, sigma_eps^1 = 0

# Para simplificar e garantir identificação: fixar
# b0^n livres (loadings contemporâneos)
# b1^n = 0 (sem lag nos loadings para identificação)
# w_t^n = passeio aleatório com drift d^n
# Estimar via MLE: phi1, phi2, sigma_g, sigma_w^n, sigma_eps^n, b0^n

build_jl <- function(theta) {
  # Parâmetros
  # theta = [phi1, phi2, log_sg, log_sw2, log_sw3, log_sw4,
  #          log_se1, log_se2, log_se3, log_se4,
  #          b02, b03, b04]
  phi1  <- tanh(theta[1]) * 0.99
  phi2  <- tanh(theta[2]) * (1 - abs(phi1)) * 0.99   # garante estacionariedade
  sg    <- exp(theta[3])   # DP choque do hiato
  sw    <- exp(theta[4:6]) # DP choques tendências w2, w3, w4
  se    <- exp(theta[7:10])# DP erros de medida (se1 fixo pequeno pela restrição)
  b0    <- theta[11:13]    # loadings: b0^2, b0^3, b0^4

  # Dimensões
  # Estado: [g_t, g_{t-1}, w_t^2, w_t^3, w_t^4]  (5×1)
  # w_t^1 = 0 pela restrição BCB (y1 = g_t diretamente)

  # Transição: g_t = phi1·g_{t-1} + phi2·g_{t-2} + eta^g
  #            w_t^n = w_{t-1}^n + eta^n
  Tt_jl <- matrix(c(
    phi1, phi2, 0,  0,  0,
    1,    0,    0,  0,  0,
    0,    0,    1,  0,  0,
    0,    0,    0,  1,  0,
    0,    0,    0,  0,  1
  ), 5, 5, byrow=TRUE)

  # Choques de estado R (5×4): [g, w2, w3, w4]
  Rt_jl <- matrix(c(
    1, 0, 0, 0,
    0, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ), 5, 4, byrow=TRUE)

  Qt_jl <- diag(c(sg^2, sw^2))   # 4×4

  # Medida: [y1, y2, y3, y4] = Z·estado + eps
  # y1 = 1·g_t + 0·g_{t-1} + 0·w2 + 0·w3 + 0·w4
  # y2 = b02·g_t + 0·g_{t-1} + 1·w2 + 0·w3 + 0·w4
  # y3 = b03·g_t + 0·g_{t-1} + 0·w2 + 1·w3 + 0·w4
  # y4 = b04·g_t + 0·g_{t-1} + 0·w2 + 0·w3 + 1·w4
  Zt_jl <- matrix(c(
    1,     0, 0, 0, 0,
    b0[1], 0, 1, 0, 0,
    b0[2], 0, 0, 1, 0,
    b0[3], 0, 0, 0, 1
  ), 4, 5, byrow=TRUE)

  Ht_jl <- diag(se^2)   # 4×4, erros de medida independentes

  y_mv_jl <- ts(Y_obs_jl, start=t_jl, frequency=4)

  SSModel(y_mv_jl ~ -1 + SSMcustom(
    Z     = array(Zt_jl, dim=c(4,5,1)),
    T     = array(Tt_jl, dim=c(5,5,1)),
    R     = array(Rt_jl, dim=c(5,4,1)),
    Q     = array(Qt_jl, dim=c(4,4,1)),
    P1inf = diag(5), P1 = matrix(0,5,5)
  ), H = array(Ht_jl, dim=c(4,4,1)))
}

nll_jl <- function(theta) {
  tryCatch({
    mod <- build_jl(theta)
    ll  <- logLik(mod)
    if (!is.finite(ll)) return(1e10)
    -as.numeric(ll)
  }, error=function(e) 1e10)
}

# Inicialização
theta0_jl <- c(
  atanh(0.7),   # phi1
  atanh(0.1),   # phi2
  log(0.5),     # sg
  log(0.3), log(0.3), log(0.3),   # sw2, sw3, sw4
  log(0.3), log(0.3), log(0.3), log(0.3),  # se1..se4
  0.5, -0.5, 0.5  # b02, b03, b04
)

cat("    Otimizando MLE (Nelder-Mead)...")
opt_jl <- optim(theta0_jl, nll_jl, method="Nelder-Mead",
                control=list(maxit=10000, reltol=1e-12))
cat(sprintf(" val=%.2f | BFGS...", opt_jl$value))
opt_jl2 <- tryCatch(
  optim(opt_jl$par, nll_jl, method="BFGS",
        control=list(maxit=3000, reltol=1e-12)),
  error=function(e) opt_jl)
if (is.finite(opt_jl2$value) && opt_jl2$value < opt_jl$value) opt_jl <- opt_jl2
cat(sprintf(" val_final=%.4f\n", opt_jl$value))

# Extrair hiato (estado 1 = g_t)
mod_jl <- build_jl(opt_jl$par)
ks_jl  <- KFS(mod_jl, smoothing="state")
g_smooth <- as.numeric(ks_jl$alphahat[,1])

# Parâmetros estimados
phi1_est <- tanh(opt_jl$par[1]) * 0.99
phi2_est <- tanh(opt_jl$par[2]) * (1 - abs(phi1_est)) * 0.99
b0_est   <- opt_jl$par[11:13]
cat(sprintf("    φ₁=%.3f, φ₂=%.3f | loadings: b²=%.3f, b³=%.3f, b⁴=%.3f\n",
            phi1_est, phi2_est, b0_est[1], b0_est[2], b0_est[3]))

# O fator g_t está em desvios-padrão padronizados — reescalar para p.p. do PIB
# usando o sd do hiato HP do PIB como âncora
hiato_hp_ref <- as.numeric(mFilter::hpfilter(
  ts(log(df_jl$pib), start=t_jl, freq=4), freq=1600)$cycle) * 100
df_jl$hiato_JL <- g_smooth * sd(hiato_hp_ref) / sd(g_smooth)

# Verificar sinal: deve correlacionar positivamente com HP-PIB
if (cor(df_jl$hiato_JL, hiato_hp_ref) < 0) {
  df_jl$hiato_JL <- -df_jl$hiato_JL
  cat("    [sinal invertido — corrigido]\n")
}

cat(sprintf("    Hiato JL atual: %+.2f p.p. | DP: %.2f p.p.\n",
            tail(df_jl$hiato_JL, 1), sd(df_jl$hiato_JL)))
cat(sprintf("    Corr(JL, PCA): %.3f | Corr(JL, FP): %.3f\n",
            cor(df_jl$hiato_JL, df_pca_raw$hiato_PCA),
            cor(df_jl$hiato_JL, df_multi$hiato_FP)))

# ════════════════════════════════════════════════════════════════════════════
# GRUPO II — MULTIVARIADO III: CBO (Shackleton, 2018; BCB RI Jun/2024, Ap. 2c)
# Função de produção com potenciais via REGRESSÕES LINEARES POR PARTES:
# cada insumo é regredido no hiato de ocupação (Egap) + tendências de tempo
# segmentadas pelos picos do CODACE; o potencial zera os coeficientes do Egap.
#
# Dados (todos baixados por código — reprodutível):
#  - Desemprego retropolado: PME nova (FTP IBGE, tab. 177, mar/02-fev/16,
#    dessaz. via STL) empalmada com PNADC (SGS 24369) por ajuste de nível
#  - Estoque de capital: Ipea/DIMAC série DIMAC_CF_ELC_TOT12
#    (Souza Júnior & Cornelio, 2020), mensal 1980-2024
#  - NUCI FGV (28561), PIB (22109)
#
# Simplificações documentadas (vs. BCB):
#  (i)   Nairu via HP(λ=1600) sobre o desemprego retropolado
#        (BCB usa Nairu do modelo semiestrutural)
#  (ii)  PEA ≡ PEA*: hiato de trabalho vem somente do Egap; o componente
#        tendencial da PEA é absorvido pela PTF e suas tendências de tempo
#  (iii) K estendido após dez/2024 pela média da FBKF dos últimos 8 trimestres
#        (espelha o pro-rata do RI para 2024T1)
# ════════════════════════════════════════════════════════════════════════════
cat("\nMII.3 CBO — Função de produção com regressões por partes...\n")

# ── 1. Desemprego retropolado 2002T2+ ────────────────────────────────────────
cat("    Baixando PME nova (FTP IBGE)... ")
pme_zip <- file.path(tempdir(), "pme_pd.zip")
pme_ok <- tryCatch({
  download.file(paste0("https://ftp.ibge.gov.br/Trabalho_e_Rendimento/",
                       "Pesquisa_Mensal_de_Emprego/Tabelas/2016/pme_201602pd.zip"),
                pme_zip, mode = "wb", quiet = TRUE)
  unzip(pme_zip, files = "tab177022016.xls", exdir = tempdir())
  TRUE
}, error = function(e) FALSE)

if (pme_ok) {
  cat("OK\n")
  tab177 <- readxl::read_excel(file.path(tempdir(), "tab177022016.xls"),
                                col_names = FALSE, .name_repair = "minimal")
  col1 <- as.character(tab177[[1]])
  l_ini <- which(grepl("^Estimativas", col1))[1] + 1
  l_fim <- which(grepl("^Coeficientes", col1))[1] - 1
  meses_pt <- c("Janeiro","Fevereiro","Março","Abril","Maio","Junho",
                "Julho","Agosto","Setembro","Outubro","Novembro","Dezembro")
  reg <- list(); ano_cur <- NA
  for (i in l_ini:l_fim) {
    c0 <- trimws(col1[i])
    if (grepl("^[0-9]{4}$", c0)) { ano_cur <- as.integer(c0); next }
    m <- match(c0, meses_pt)
    if (!is.na(m) && !is.na(ano_cur)) {
      v <- suppressWarnings(as.numeric(gsub(",", ".", as.character(tab177[[2]][i]))))
      if (!is.na(v)) reg[[length(reg)+1]] <- tibble(
        data = as.Date(sprintf("%d-%02d-01", ano_cur, m)), pme = v)
    }
  }
  pme_m <- bind_rows(reg) %>% arrange(data)

  # Mensal → trimestral (apenas trimestres completos) → dessaz STL
  pme_t <- pme_m %>%
    mutate(ano = year(data), tri = quarter(data)) %>%
    group_by(ano, tri) %>%
    summarise(pme = mean(pme), nm = n(), .groups = "drop") %>%
    filter(nm == 3) %>%
    mutate(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3+1))) %>%
    arrange(data)
  stl_pme <- stl(ts(pme_t$pme, start = c(pme_t$ano[1], pme_t$tri[1]),
                    frequency = 4), s.window = "periodic", robust = TRUE)
  pme_t$pme_sa <- pme_t$pme - as.numeric(stl_pme$time.series[, "seasonal"])

  # Empalme com PNADC (ajuste de nível aditivo no período sobreposto)
  overlap <- pme_t %>% inner_join(desemp_raw %>% rename(pnadc = valor), by = "data")
  ajuste_nivel <- mean(overlap$pnadc - overlap$pme_sa)
  cat(sprintf("    Empalme PME→PNADC: overlap=%d tri | corr=%.3f | ajuste=%+.2f p.p.\n",
              nrow(overlap), cor(overlap$pnadc, overlap$pme_sa), ajuste_nivel))

  desemp_retro <- bind_rows(
    pme_t %>% filter(data < min(desemp_raw$data)) %>%
      transmute(data, desemp = pme_sa + ajuste_nivel),
    desemp_raw %>% rename(desemp = valor)
  ) %>% arrange(data)
} else {
  cat("FALHOU — CBO usará apenas PNADC (2012+)\n")
  desemp_retro <- desemp_raw %>% rename(desemp = valor)
}

# ── 2. Estoque de capital (Ipea/DIMAC) ───────────────────────────────────────
cat("    Baixando estoque de capital (Ipeadata DIMAC_CF_ELC_TOT12)... ")
cap_ok <- tryCatch({
  cap_raw <- jsonlite::fromJSON(paste0(
    "http://www.ipeadata.gov.br/api/odata4/",
    "ValoresSerie(SERCODIGO='DIMAC_CF_ELC_TOT12')"))$value
  TRUE
}, error = function(e) FALSE)

if (cap_ok) {
  cat(sprintf("OK (%d obs)\n", nrow(cap_raw)))
  cap_m <- tibble(data = as.Date(substr(cap_raw$VALDATA, 1, 10)),
                  K    = as.numeric(cap_raw$VALVALOR)) %>%
    filter(!is.na(K)) %>% arrange(data)
  # Estoque: valor do ÚLTIMO mês de cada trimestre
  cap_t <- cap_m %>%
    mutate(ano = year(data), tri = quarter(data)) %>%
    group_by(ano, tri) %>% slice_tail(n = 1) %>% ungroup() %>%
    transmute(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3+1)), K)
} else {
  cat("FALHOU\n")
}

if (pme_ok && cap_ok) {

  # ── 3. Base CBO (2002T2+) ──────────────────────────────────────────────────
  df_cbo <- pib_raw %>% rename(pib = valor) %>%
    inner_join(nuci_raw %>% rename(nuci = valor), by = "data") %>%
    inner_join(desemp_retro, by = "data") %>%
    left_join(cap_t, by = "data") %>%
    filter(data >= as.Date("2002-04-01")) %>%
    arrange(data)

  # Estender K após o fim da série (média da FBKF dos últimos 8 trimestres)
  if (any(is.na(df_cbo$K))) {
    K_obs <- df_cbo$K[!is.na(df_cbo$K)]
    fbkf_med <- mean(diff(tail(K_obs, 9)))
    idx_na <- which(is.na(df_cbo$K))
    for (i in idx_na) df_cbo$K[i] <- df_cbo$K[i-1] + fbkf_med
    cat(sprintf("    K estendido em %d trimestres (FBKF média = %.0f)\n",
                length(idx_na), fbkf_med))
  }
  T_cbo <- nrow(df_cbo)
  cat(sprintf("    Amostra CBO: %s a %s (%d obs.)\n",
              min(df_cbo$data), max(df_cbo$data), T_cbo))

  # ── 4. Egap: hiato da taxa de ocupação ─────────────────────────────────────
  # Nairu via HP (simplificação i); E = 1 − u/100
  # Nairu via HP super-suave (λ=25600): com λ=1600 a Nairu acompanha o
  # desemprego de perto e comprime o Egap; λ maior aproxima a rigidez da
  # Nairu do modelo semiestrutural do BCB
  hp_u_cbo <- mFilter::hpfilter(ts(df_cbo$desemp, frequency = 4),
                                freq = 25600, type = "lambda")
  df_cbo <- df_cbo %>% mutate(
    u_star = as.numeric(hp_u_cbo$trend),
    E      = 1 - desemp/100,
    E_star = 1 - u_star/100,
    Egap   = (E/E_star - 1) * 100,
    Egap_l = lag(Egap)
  )

  # ── 5. Tendências de tempo segmentadas (picos CODACE, nota 18 do RI) ───────
  # Rampa de 25/trimestre do pico anterior até o pico do ciclo; constante após
  picos <- as.Date(c("2002-10-01","2008-07-01","2014-01-01","2019-10-01"))
  make_T <- function(datas, pico_ant, pico) {
    idx <- as.integer(datas > pico_ant & datas <= pico)
    cumsum(idx) * 25 -> rampa
    rampa[datas > pico] <- max(rampa[datas <= pico])
    rampa
  }
  # Segmento ABERTO pós-2019T4: o CODACE não datou pico posterior, então o
  # fim da amostra atua como pico provisório (senão o potencial fica sem
  # tendência por 25 trimestres e o hiato cresce artificialmente — é também
  # uma fonte conhecida de revisão em tempo real do método CBO)
  df_cbo <- df_cbo %>% mutate(
    T2002 = make_T(data, as.Date("1900-01-01"), picos[1]),
    T2008 = make_T(data, picos[1], picos[2]),
    T2014 = make_T(data, picos[2], picos[3]),
    T2019 = make_T(data, picos[3], picos[4]),
    TPOS  = make_T(data, picos[4], max(data)),
    DCOV  = as.integer(data >= as.Date("2020-01-01"))
  )

  # ── 6. Regressões por partes (potencial = predição com Egap = 0) ───────────
  reg_piecewise <- function(dep, df) {
    fit <- lm(dep ~ Egap + Egap_l + T2002 + T2008 + T2014 + T2019 + TPOS + DCOV,
              data = df, na.action = na.exclude)
    df0 <- df %>% mutate(Egap = 0, Egap_l = 0)
    list(fit = fit, pot = predict(fit, newdata = df0))
  }

  # (a) FBKF e capital potencial (inventário perpétuo com FBKF*)
  df_cbo$FBKF <- c(NA, diff(df_cbo$K))
  r_fbkf <- reg_piecewise(log(df_cbo$FBKF), df_cbo)
  df_cbo$FBKF_star <- exp(r_fbkf$pot)
  K_star <- df_cbo$K
  for (i in 2:T_cbo) K_star[i] <- K_star[i-1] + df_cbo$FBKF_star[i]
  df_cbo$K_star <- K_star

  # (b) NUCI potencial
  r_nuci <- reg_piecewise(df_cbo$nuci, df_cbo)
  df_cbo$nuci_star <- r_nuci$pot

  # (c) Índices base 2002T4 e PTF
  i_base <- which(df_cbo$data == as.Date("2002-10-01"))
  df_cbo <- df_cbo %>% mutate(
    IOCUP      = 100 * E / E[i_base],                 # simplificação (ii): PEA ≡ PEA*
    IOCUP_star = 100 * E_star / E[i_base],
    IK         = 100 * K / K[i_base],
    IK_star    = 100 * K_star / K[i_base],
    IKN        = IK * nuci / 100,
    IKN_star   = IK_star * nuci_star / 100,
    lnA        = log(pib) - 0.6*log(IOCUP) - 0.4*log(IKN)
  )
  r_ptf <- reg_piecewise(df_cbo$lnA, df_cbo)
  df_cbo$lnA_star <- r_ptf$pot

  # (d) Produto potencial e hiato (+0,4·ln IKN*, consistente com Cobb-Douglas)
  df_cbo <- df_cbo %>% mutate(
    ln_pib_star = lnA_star + 0.6*log(IOCUP_star) + 0.4*log(IKN_star),
    hiato_CBO   = (log(pib) - ln_pib_star) * 100
  )

  cat(sprintf("    Coef. Egap: FBKF=%.3f | NUCI=%.3f | PTF=%.4f\n",
              coef(r_fbkf$fit)["Egap"], coef(r_nuci$fit)["Egap"],
              coef(r_ptf$fit)["Egap"]))
  cat(sprintf("    Hiato CBO atual: %+.2f p.p. | DP: %.2f p.p.\n",
              tail(df_cbo$hiato_CBO, 1), sd(df_cbo$hiato_CBO, na.rm = TRUE)))
  hp_ref_cbo <- df_cbo %>% inner_join(pib %>% select(data, hiato_HP), by = "data")
  cmp_fp <- df_cbo %>% inner_join(df_multi %>% select(data, hiato_FP), by = "data")
  cat(sprintf("    Corr(CBO, HP): %.3f | Corr(CBO, FP 2012+): %.3f\n",
              cor(hp_ref_cbo$hiato_CBO, hp_ref_cbo$hiato_HP, use = "complete.obs"),
              cor(cmp_fp$hiato_CBO, cmp_fp$hiato_FP, use = "complete.obs")))
} else {
  df_cbo <- NULL
  cat("    CBO não estimado (falha no download de dados)\n")
}

# ════════════════════════════════════════════════════════════════════════════
# ELEMENTOS COMUNS AOS GRÁFICOS
# ════════════════════════════════════════════════════════════════════════════
recessoes <- tibble(
  inicio = as.Date(c("1997-10-01","2001-04-01","2002-10-01",
                     "2008-10-01","2014-04-01","2019-10-01","2020-01-01")),
  fim    = as.Date(c("1999-01-01","2001-10-01","2003-07-01",
                     "2009-04-01","2016-10-01","2020-01-01","2020-10-01"))
)
df_qbr <- tibble(data  = quebras_datas,
                 label = c("2000T1","2008T4","2013T3","2020T2"))

# Tema padrão compartilhado
tema_padrao <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12.5),
    plot.subtitle    = element_text(size = 9, color = "#555555", margin = margin(b = 8)),
    plot.caption     = element_text(size = 7.5, color = "#777777", hjust = 0,
                                    margin = margin(t = 10)),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.4, "cm"),
    panel.grid.major = element_line(color = "grey91", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 20, 10, 10)
  )

escala_x_uni <- scale_x_date(
  breaks      = seq(as.Date("1996-01-01"), max(pib$data) + 365, by = "4 years"),
  date_labels = "%Y",
  expand      = expansion(mult = c(0.01, 0.02))
)
escala_y <- scale_y_continuous(
  limits = c(-5.5, 5.5), breaks = seq(-5, 5, by = 1),
  labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%")
)

# ════════════════════════════════════════════════════════════════════════════
# GRÁFICO 1 — 6 MÉTODOS UNIVARIADOS
# ════════════════════════════════════════════════════════════════════════════
metodo_uni <- c(
  hiato_TQ  = "I — Tend. Quadrática c/ Quebras",
  hiato_NP  = "II — Tend. Não-Paramétrica (LOESS)",
  hiato_HP  = "III — Hodrick-Prescott (λ=1600)",
  hiato_L1  = "IV — Tendência ℓ₁ (Kim et al., 2009)",
  hiato_HPM = "V — HP Modificada (Andrle/BCB)",
  hiato_BP  = "VI — Band-Pass (CF, 8-32T)",
  hiato_BN  = "VII — Beveridge-Nelson (Kamber et al., 2018)"
)
cores_uni <- c(
  "I — Tend. Quadrática c/ Quebras"              = "#1A5276",
  "II — Tend. Não-Paramétrica (LOESS)"           = "#C0392B",
  "III — Hodrick-Prescott (λ=1600)"              = "#1E8449",
  "IV — Tendência ℓ₁ (Kim et al., 2009)"        = "#7D3C98",
  "V — HP Modificada (Andrle/BCB)"               = "#D4AC0D",
  "VI — Band-Pass (CF, 8-32T)"                   = "#117A65",
  "VII — Beveridge-Nelson (Kamber et al., 2018)" = "#5D6D7E"
)

pib_long <- pib %>%
  select(data, all_of(names(metodo_uni))) %>%
  pivot_longer(-data, names_to = "metodo", values_to = "hiato") %>%
  mutate(metodo = factor(recode(metodo, !!!metodo_uni), levels = unname(metodo_uni)))

g_uni <- ggplot(pib_long, aes(x = data, y = hiato, color = metodo)) +
  geom_rect(data = recessoes, aes(xmin=inicio, xmax=fim, ymin=-Inf, ymax=Inf),
            fill = "grey85", alpha = 0.4, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(data = df_qbr, aes(xintercept = data),
             linetype = "dotted", color = "#E07B39", linewidth = 0.55,
             alpha = 0.8, inherit.aes = FALSE) +
  geom_label(data = df_qbr, aes(x = data, y = 5.5*0.86, label = label),
             hjust = -0.07, size = 2.4, color = "#E07B39", fill = "white",
             label.size = 0.15, fontface = "bold", inherit.aes = FALSE) +
  geom_line(linewidth = 0.75, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(values = cores_uni) +
  escala_x_uni + escala_y +
  labs(
    title    = "Hiato do Produto — Grupo I: Métodos Univariados (BCB RI Jun/2024)",
    subtitle = "I. TQ c/ Quebras · II. LOESS · III. HP · IV. ℓ₁ · V. HP Mod. · VI. Band-Pass · VII. Beveridge-Nelson",
    x = NULL, y = "Hiato (% do PIB potencial)", color = NULL,
    caption = paste0(
      "Fonte: BCB/SGS série 22109 (PIB dessaz., índice). Elaboração própria.\n",
      "I: Quebras Bai-Perron: 2000T1,2008T4,2013T3,2020T2 | II: LOESS bwidth=0,4 (Cleveland,1979) | ",
      "III: HP λ=1600 | IV: ℓ₁ λ=(1/2)⁵λ_max (Kim et al.,2009)\n",
      "V: HP Mod. ciclo AR(2) (Andrle,2013) | VI: CF Band-Pass 8-32T | ",
      "VII: BN-KMW AR(12), δ=1, WLS COVID (Morley et al.,2023) | Recessões: CODACE/FGV."
    )
  ) +
  tema_padrao +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

# ════════════════════════════════════════════════════════════════════════════
# GRÁFICO 2 — MULTIVARIADO I: FUNÇÃO DE PRODUÇÃO
# Mesmo padrão visual do gráfico univariado
# ════════════════════════════════════════════════════════════════════════════
# Recessões recortadas ao período multivariado
rec_multi <- recessoes %>% filter(fim >= min(df_multi$data))

escala_x_multi <- scale_x_date(
  breaks      = seq(as.Date("2012-01-01"), max(df_multi$data) + 365, by = "2 years"),
  date_labels = "%Y",
  expand      = expansion(mult = c(0.01, 0.02))
)

g_multi <- ggplot(df_multi, aes(x = data)) +
  geom_rect(data = rec_multi, aes(xmin=inicio, xmax=fim, ymin=-Inf, ymax=Inf),
            fill = "grey85", alpha = 0.4, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +

  # Contribuições empilhadas (barras)
  geom_col(aes(y = contrib_cap,  fill = "Capital — (1−α)·ĉ_t  [α=0,4]"),
           alpha = 0.65, width = 70) +
  geom_col(aes(y = contrib_trab, fill = "Trabalho — (−α)·û_t  [α=0,6]"),
           alpha = 0.65, width = 70) +

  # Hiato total FP (linha preta — mesmo estilo das linhas univariadas)
  geom_line(aes(y = hiato_FP, color = "MI — Função de Produção (comb. simples)"),
            linewidth = 0.9) +

  scale_fill_manual(values = c(
    "Capital — (1−α)·ĉ_t  [α=0,4]"  = "#2980B9",
    "Trabalho — (−α)·û_t  [α=0,6]"  = "#E74C3C"
  )) +
  scale_color_manual(values = c(
    "MI — Função de Produção (comb. simples)" = "#1C1C1C"
  )) +

  escala_x_multi + escala_y +
  labs(
    title    = "MII.1 — Função de Produção — Combinação Simples (BCB RI Jun/2024)",
    subtitle = "ŷ_t = (1−α)·ĉ_t − α·û_t  |  α=0,6 (trabalho), 1−α=0,4 (capital)  |  Hiatos via HP (λ=1600)",
    x = NULL, y = "Hiato (% do PIB potencial)",
    fill = "Contribuição:", color = NULL,
    caption = paste0(
      "Fonte: BCB/SGS. NUCI FGV (série 28561, dessaz.); Desemprego PNADC (série 24369, dessaz.).\n",
      "Hiatos de NUCI e desemprego estimados via filtro HP (λ=1600). ",
      "Pesos: α=0,6 e 1−α=0,4 (média 1999-2019, BCB RI Jun/2024, nota 11).\n",
      "Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao +
  guides(fill  = guide_legend(order = 1, nrow = 1),
         color = guide_legend(order = 2, nrow = 1))

# ── Gráfico PCA ─────────────────────────────────────────────────────────────
rec_pca <- recessoes %>% filter(fim >= min(df_pca_raw$data))

escala_x_pca <- scale_x_date(
  breaks      = seq(min(df_pca_raw$data), max(df_pca_raw$data) + 365, by = "1 year"),
  date_labels = "%Y",
  expand      = expansion(mult = c(0.01, 0.02))
)

g_pca <- ggplot(df_pca_raw, aes(x = data)) +
  geom_rect(data = rec_pca, aes(xmin=inicio, xmax=fim, ymin=-Inf, ymax=Inf),
            fill = "grey85", alpha = 0.4, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +

  # Área preenchida
  geom_ribbon(aes(ymin = pmin(hiato_PCA, 0), ymax = 0), fill = "#E74C3C", alpha = 0.2) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(hiato_PCA, 0)), fill = "#1E8449", alpha = 0.2) +

  geom_line(aes(y = hiato_PCA, color = "MII.V — Componentes Principais (PCA)"),
            linewidth = 0.9) +

  scale_color_manual(values = c("MII.V — Componentes Principais (PCA)" = "#8E44AD")) +
  escala_x_pca + escala_y +
  labs(
    title    = "MII.V — Componentes Principais — PCA (BCB RI Jun/2024)",
    subtitle = sprintf("1º CP explica %.1f%% da variância | Séries: Hiato PIB-HP, NUCI, −Desemprego, Hiato Caged-HP", var_exp[1]),
    x = NULL, y = "Hiato (% do PIB potencial)", color = NULL,
    caption = paste0(
      "Fonte: BCB/SGS. PIB (22109), NUCI FGV (28561), Desemprego PNADC (24369), Novo Caged (28763).\n",
      "Proc.: hiatos do PIB e Caged via HP (λ=1600); NUCI em nível; desemprego com sinal invertido. ",
      "Séries padronizadas (z-score) antes do PCA.\n",
      "Recessões (áreas cinzas): CODACE/FGV. Referência: BCB RI Jun/2024, Apêndice 2e."
    )
  ) +
  tema_padrao

# ── Gráfico Areosa ───────────────────────────────────────────────────────────
rec_areosa <- recessoes %>% filter(fim >= min(df_areosa$data))

g_areosa <- ggplot(df_areosa, aes(x = data)) +
  geom_rect(data = rec_areosa,
            aes(xmin=inicio, xmax=fim, ymin=-Inf, ymax=Inf),
            fill="grey85", alpha=0.4, inherit.aes=FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.4) +

  geom_ribbon(aes(ymin=pmin(hiato_AR,0), ymax=0), fill="#E74C3C", alpha=0.2) +
  geom_ribbon(aes(ymin=0, ymax=pmax(hiato_AR,0)), fill="#1E8449", alpha=0.2) +

  geom_line(aes(y=hiato_AR, color="MII.II — Areosa (2008)"),
            linewidth=0.9) +

  scale_color_manual(values=c("MII.II — Areosa (2008)" = "#E67E22")) +
  escala_x_multi + escala_y +
  labs(
    title    = "MII.II — Função de Produção — Areosa (2008) (BCB RI Jun/2024)",
    subtitle = "Solução exata da otimização conjunta (Proposição 1): ĥ_Y = 0,658·ĥ_Y^FP + 0,342·ĥ_Y^HP | λ=1600 comum",
    x=NULL, y="Hiato (% do PIB potencial)", color=NULL,
    caption=paste0(
      "Fonte: BCB/SGS. NUCI FGV (28561), Desemprego PNADC (24369), PIB (22109).\n",
      "Modelo: três filtros HP interligados via restrição da função de produção Cobb-Douglas (Areosa, 2008).\n",
      "Solução analítica exata do problema acoplado (CPOs fatoram em (I+λD'D); identidade λSA=I−S):\n",
      "hiato de Areosa = média ponderada do hiato FP simples (peso 1/κ=0,658) e do hiato HP do PIB (peso 0,342).\n",
      "Verificação numérica: sistema 2T×2T vs fórmula, erro máx 4,6e-11. Corr(Areosa,FP)=0,979.\n",
      "Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao

# ── Gráfico JL ───────────────────────────────────────────────────────────────
rec_jl <- recessoes %>% filter(fim >= min(df_jl$data))

escala_x_jl <- scale_x_date(
  breaks      = seq(min(df_jl$data), max(df_jl$data)+365, by="2 years"),
  date_labels = "%Y",
  expand      = expansion(mult=c(0.01,0.02))
)

g_jl <- ggplot(df_jl, aes(x=data)) +
  geom_rect(data=rec_jl, aes(xmin=inicio,xmax=fim,ymin=-Inf,ymax=Inf),
            fill="grey85", alpha=0.4, inherit.aes=FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.4) +
  geom_ribbon(aes(ymin=pmin(hiato_JL,0), ymax=0), fill="#E74C3C", alpha=0.2) +
  geom_ribbon(aes(ymin=0, ymax=pmax(hiato_JL,0)), fill="#1E8449", alpha=0.2) +
  geom_line(aes(y=hiato_JL, color="MII.IV — Jarocinski & Lenza (2018)"),
            linewidth=0.9) +
  scale_color_manual(values=c("MII.IV — Jarocinski & Lenza (2018)"="#2E86C1")) +
  escala_x_jl + escala_y +
  labs(
    title    = "MII.IV — Jarocinski & Lenza (2018) — Versão Linearizada (BCB RI Jun/2024)",
    subtitle = "Fator dinâmico comum: PIB, NUCI, -Desemprego, Caged | Hiato AR(2) | MLE via KFAS",
    x=NULL, y="Hiato (% do PIB potencial)", color=NULL,
    caption=paste0(
      "Fonte: BCB/SGS. PIB (22109), NUCI FGV (28561), Desemprego PNADC (24369), Novo Caged (28763).\n",
      "Versão linearizada sem volatilidade estocástica (h_t=0). ",
      "Hiato=fator AR(2) comum às séries de atividade, estimado por MLE via KFAS.\n",
      "Referência: Jarocinski & Lenza (2018), JMCB; BCB RI Jun/2024, Apêndice 2d.\n",
      "Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao

# ════════════════════════════════════════════════════════════════════════════
# GRÁFICO 6 — CONSOLIDADO (estilo Gráfico 3 do BCB RI Jun/2024)
# Faixa mín–máx + mediana + média de todos os métodos disponíveis em cada t
# Areosa excluído: combinação convexa de FP e HP (Proposição 1) — não pode
# expandir a faixa mín–máx por construção
# ════════════════════════════════════════════════════════════════════════════
cat("\nMontando gráfico consolidado (thick modeling)...\n")

hiatos_todos <- bind_rows(
  pib %>%
    select(data, hiato_TQ, hiato_NP, hiato_HP, hiato_L1, hiato_HPM, hiato_BP, hiato_BN) %>%
    pivot_longer(-data, names_to="metodo", values_to="hiato"),
  df_multi   %>% transmute(data, metodo="hiato_FP",  hiato=hiato_FP),
  df_pca_raw %>% transmute(data, metodo="hiato_PCA", hiato=hiato_PCA),
  df_jl      %>% transmute(data, metodo="hiato_JL",  hiato=hiato_JL),
  if (!is.null(df_cbo)) df_cbo %>%
    transmute(data, metodo="hiato_CBO", hiato=hiato_CBO) else NULL
) %>% filter(!is.na(hiato))

df_consol <- hiatos_todos %>%
  group_by(data) %>%
  summarise(h_min = min(hiato), h_max = max(hiato),
            h_med = median(hiato), h_mean = mean(hiato),
            n_met = n(), .groups="drop")

cat(sprintf("    Métodos por trimestre: %d (1996-2012) a %d (2012+)\n",
            min(df_consol$n_met), max(df_consol$n_met)))
cat(sprintf("    Amplitude da faixa em %s: %.2f p.p. (de %+.2f a %+.2f)\n",
            max(df_consol$data),
            tail(df_consol$h_max,1) - tail(df_consol$h_min,1),
            tail(df_consol$h_min,1), tail(df_consol$h_max,1)))

g_consol <- ggplot(df_consol, aes(x=data)) +
  geom_rect(data=recessoes, aes(xmin=inicio,xmax=fim,ymin=-Inf,ymax=Inf),
            fill="grey85", alpha=0.4, inherit.aes=FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.4) +
  geom_ribbon(aes(ymin=h_min, ymax=h_max,
                  fill="Intervalo mín–máx (10 métodos)"), alpha=0.35) +
  geom_line(aes(y=h_mean, color="Média"),   linewidth=0.6, linetype="dashed") +
  geom_line(aes(y=h_med,  color="Mediana"), linewidth=1.0) +
  scale_fill_manual(values=c("Intervalo mín–máx (10 métodos)"="#85C1E9")) +
  scale_color_manual(values=c("Mediana"="#1A5276", "Média"="#C0392B")) +
  escala_x_uni + escala_y +
  labs(
    title    = "Hiato do Produto — Síntese dos Métodos (estilo Gráfico 3, BCB RI Jun/2024)",
    subtitle = "Faixa mín–máx, mediana e média | 7 univariados (1996+) + CBO (2002+) + FP, PCA e JL (2012+) | Areosa: comb. convexa (excl.)",
    x=NULL, y="Hiato (% do PIB potencial)", fill=NULL, color=NULL,
    caption = paste0(
      "Fonte: BCB/SGS. Elaboração própria. A faixa agrega, em cada trimestre, todos os métodos com estimativa disponível:\n",
      "7 univariados desde 1996T1; CBO desde 2002T2; Função de Produção, PCA e Jarocinski-Lenza desde 2012T2. ",
      "A amplitude da faixa é a medida visual da incerteza de mensuração do hiato (thick modeling).\n",
      "Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao

# ════════════════════════════════════════════════════════════════════════════
# EXERCÍCIO DE TEMPO REAL (pseudo-real-time — conexão com WP 203, Cusinato et al.)
# A cada vintage τ (trimestral, desde 2006T1), estima com dados ATÉ τ e guarda
# o hiato do último ponto ("tempo real"). Compara com a amostra completa ("final").
# Sem revisões de dados (série atual do SGS): isola o VIÉS DE FINAL DE AMOSTRA.
# Indicadores no espírito do WP 203: RAM, RRQM, correlação, troca de sinal.
# ════════════════════════════════════════════════════════════════════════════
cat("\nExercício de tempo real (HP e TQ, vintages trimestrais desde 2006T1)...\n")

idx_vint <- 41:n          # 2006T1 em diante (≥40 obs por vintage)
rt_hp <- rep(NA_real_, n)
rt_tq <- rep(NA_real_, n)

for (tau in idx_vint) {
  # HP em tempo real
  hp_v <- mFilter::hpfilter(ts(y[1:tau], start=c(1996,1), frequency=4),
                            freq=1600, type="lambda")
  rt_hp[tau] <- as.numeric(hp_v$cycle)[tau] * 100

  # TQ em tempo real: apenas quebras com ≥6 obs no segmento final
  qb_v <- idx_q[idx_q <= tau - 6]
  kv   <- length(qb_v)
  ini_v <- c(1, qb_v); fim_v <- c(qb_v - 1, tau)
  Xv <- matrix(0, tau, 3*(kv+1))
  for (j in 1:(kv+1)) {
    tl <- ini_v[j]:fim_v[j]
    Xv[tl, 3*(j-1)+1] <- 1
    Xv[tl, 3*(j-1)+2] <- tl
    Xv[tl, 3*(j-1)+3] <- tl^2
  }
  rt_tq[tau] <- residuals(lm(y[1:tau] ~ Xv - 1))[tau] * 100
}

df_rt <- bind_rows(
  tibble(data=pib$data, metodo="III — Hodrick-Prescott (λ=1600)",
         tempo_real=rt_hp, final=pib$hiato_HP),
  tibble(data=pib$data, metodo="I — Tend. Quadrática c/ Quebras",
         tempo_real=rt_tq, final=pib$hiato_TQ)
) %>% filter(!is.na(tempo_real))

cat("\n=== REVISÕES (final − tempo real), indicadores estilo WP 203 ===\n")
cat(sprintf("%-38s %6s %6s %14s %12s\n",
            "Método","RAM","RRQM","corr(RT,final)","troca sinal"))
cat(strrep("-",80),"\n")
for (m in unique(df_rt$metodo)) {
  d   <- df_rt %>% filter(metodo == m)
  rev <- d$final - d$tempo_real
  cat(sprintf("%-38s %6.2f %6.2f %14.3f %11.1f%%\n",
              m, mean(abs(rev)), sqrt(mean(rev^2)),
              cor(d$tempo_real, d$final),
              100*mean(sign(d$tempo_real) != sign(d$final))))
}

g_rt <- ggplot(df_rt, aes(x=data)) +
  geom_rect(data=recessoes %>% filter(fim >= as.Date("2006-01-01")),
            aes(xmin=inicio,xmax=fim,ymin=-Inf,ymax=Inf),
            fill="grey85", alpha=0.4, inherit.aes=FALSE) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.4) +
  geom_line(aes(y=final,      color="Final (amostra completa)"),
            linewidth=0.85) +
  geom_line(aes(y=tempo_real, color="Tempo real (último ponto de cada vintage)"),
            linewidth=0.7, linetype="longdash") +
  facet_wrap(~metodo, ncol=1) +
  scale_color_manual(values=c(
    "Final (amostra completa)"                       = "#1A5276",
    "Tempo real (último ponto de cada vintage)"      = "#C0392B")) +
  scale_x_date(breaks=seq(as.Date("2006-01-01"), max(pib$data)+365, by="2 years"),
               date_labels="%Y", expand=expansion(mult=c(0.01,0.02))) +
  scale_y_continuous(limits=c(-5.5,5.5), breaks=seq(-4,4,by=2),
                     labels=function(x) paste0(ifelse(x>0,"+",""),x,"%")) +
  labs(
    title    = "Hiato em Tempo Real vs. Final — Viés de Final de Amostra (WP 203)",
    subtitle = "Estimativas recursivas trimestrais desde 2006T1 | Pseudo-tempo real: sem revisões de dados, isola o efeito de borda",
    x=NULL, y="Hiato (% do PIB potencial)", color=NULL,
    caption = paste0(
      "Fonte: BCB/SGS 22109. Elaboração própria, no espírito de Cusinato, Minella & Pôrto Jr. (WP 203 / Empirical Economics 2013).\n",
      "Em cada vintage τ, o método é reestimado com dados até τ e guarda-se o hiato de τ (linha tracejada). ",
      "A linha cheia é a estimativa com a amostra completa.\n",
      "TQ: quebras Bai-Perron entram apenas quando há ≥6 trimestres no segmento final. ",
      "Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao +
  theme(strip.text = element_text(face="bold", size=10))

# ── Gráfico CBO ──────────────────────────────────────────────────────────────
if (!is.null(df_cbo)) {
  rec_cbo <- recessoes %>% filter(fim >= min(df_cbo$data))
  escala_x_cbo <- scale_x_date(
    breaks = seq(as.Date("2002-01-01"), max(df_cbo$data)+365, by="2 years"),
    date_labels = "%Y", expand = expansion(mult=c(0.01,0.02)))

  g_cbo <- ggplot(df_cbo, aes(x=data)) +
    geom_rect(data=rec_cbo, aes(xmin=inicio,xmax=fim,ymin=-Inf,ymax=Inf),
              fill="grey85", alpha=0.4, inherit.aes=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.4) +
    geom_ribbon(aes(ymin=pmin(hiato_CBO,0), ymax=0), fill="#E74C3C", alpha=0.2) +
    geom_ribbon(aes(ymin=0, ymax=pmax(hiato_CBO,0)), fill="#1E8449", alpha=0.2) +
    geom_line(aes(y=hiato_CBO, color="MII.III — CBO (Shackleton, 2018)"),
              linewidth=0.9, na.rm=TRUE) +
    scale_color_manual(values=c("MII.III — CBO (Shackleton, 2018)"="#B03A2E")) +
    escala_x_cbo + escala_y +
    labs(
      title    = "MII.III — Função de Produção CBO (BCB RI Jun/2024)",
      subtitle = "Regressões por partes (picos CODACE) | Egap zerado nos potenciais | Único multivariado com amostra desde 2002",
      x=NULL, y="Hiato (% do PIB potencial)", color=NULL,
      caption=paste0(
        "Fontes: PME nova (IBGE, tab.177, dessaz. STL, empalmada com PNADC); ",
        "estoque de capital Ipea/DIMAC (Souza Júnior & Cornelio, 2020); NUCI FGV; PIB BCB/SGS 22109.\n",
        "Tendências de tempo segmentadas nos picos CODACE (2002T4, 2008T3, 2014T1, 2019T4), ",
        "segmento aberto pós-2019T4 (pico provisório = fim da amostra) + dummy COVID. Nairu via HP (simplificação ante o modelo semiestrutural do BCB).\n",
        "PEA \u2261 PEA* (hiato de trabalho via Egap). Recessões (áreas cinzas): CODACE/FGV. ",
        "Referência: Shackleton (2018); BCB RI Jun/2024, Apêndice 2c."
      )
    ) +
    tema_padrao
}

# ════════════════════════════════════════════════════════════════════════════
# CENÁRIO 2 — CONTEÚDO PREDITIVO PARA INFLAÇÃO (Phillips fora-da-amostra)
# Orphanides & van Norden (2005); Stock & Watson (1999)
# Pergunta: o hiato em TEMPO REAL ajuda a prever a inflação acumulada
# 4 trimestres à frente, além do que a própria inflação passada prevê?
#   Modelo:    pi4_{t+4} = a + b·pi4_t + g·hiato_t^RT + e   (janela expansiva)
#   Benchmark: pi4_{t+4} = a + b·pi4_t + e                  (AR puro)
# Métrica: RMSE relativo ao AR (<1 → hiato agrega); teste DM (Newey-West).
# Inclui a MEDIANA dos hiatos RT — o teste do thick modeling.
# ════════════════════════════════════════════════════════════════════════════
cat("\nCenário 2 — Previsão de inflação (Phillips fora-da-amostra)...\n")

# ── IPCA trimestral e acumulado em 4T ────────────────────────────────────────
ipca_raw <- get_sgs(433, "01/01/1995")
ipca_tri <- ipca_raw %>% rename(ipca = valor) %>%
  mutate(ano = year(data), tri = quarter(data), fator = 1 + ipca/100) %>%
  group_by(ano, tri) %>%
  summarise(fator = prod(fator), nm = n(), .groups = "drop") %>%
  filter(nm == 3) %>%
  mutate(data = as.Date(sprintf("%d-%02d-01", ano, (tri-1)*3+1))) %>%
  arrange(data) %>%
  mutate(pi4 = (fator * lag(fator) * lag(fator,2) * lag(fator,3) - 1) * 100)

# ── Hiatos em tempo real: 6 univariados (HP e TQ já computados) ─────────────
cat("    Vintages NP/L1/BP/BN (pode levar ~1-2 min) ")
rt_mat <- matrix(NA_real_, n, 6, dimnames = list(NULL, c("TQ","HP","NP","L1","BP","BN")))
rt_mat[, "HP"] <- rt_hp
rt_mat[, "TQ"] <- rt_tq
for (tau in idx_vint) {
  yv <- y[1:tau]
  # NP (LOESS)
  rt_mat[tau, "NP"] <- (yv[tau] - loess_bcb(yv)[tau]) * 100
  # L1 (lambda recalibrado por vintage)
  Dv <- diff(diag(tau), differences = 2)
  lam_max_v <- max(abs(solve(Dv %*% t(Dv), Dv %*% yv)))
  lam_v <- (0.5^5) * lam_max_v
  rm_v  <- lam_v / (max(abs(Dv %*% yv)) / 5)
  rt_mat[tau, "L1"] <- (yv[tau] -
    as.numeric(l1_filter(yv, lam_v, rho_mult = rm_v))[tau]) * 100
  # BP (CF assimétrico)
  rt_mat[tau, "BP"] <- cf_filter(yv)[tau] * 100
  # BN (KMW; AR(12) exige vintage razoável)
  bn_v <- tryCatch(bn_kmw(yv, pib$data[1:tau], p = 12), error = function(e) NULL)
  if (!is.null(bn_v)) rt_mat[tau, "BN"] <- tail(na.omit(bn_v$hiato), 1)
  if (tau %% 10 == 0) cat(".")
}
cat(" OK\n")

# Mediana dos hiatos RT (thick modeling)
rt_med <- apply(rt_mat, 1, function(r)
  if (all(is.na(r))) NA_real_ else median(r, na.rm = TRUE))

# ════════════════════════════════════════════════════════════════════════════
# CENÁRIO 1 COMPLETO — REVISÕES DE TODOS OS UNIVARIADOS + MEDIANA
# (o rt_mat acabou de ser computado; HPM fica fora: MLE por vintage é
#  computacionalmente proibitivo — documentado como limitação)
# ════════════════════════════════════════════════════════════════════════════
finais_uni <- pib %>% select(data, TQ = hiato_TQ, NP = hiato_NP, HP = hiato_HP,
                             L1 = hiato_L1, BP = hiato_BP, BN = hiato_BN)
finais_uni$MED <- apply(finais_uni %>% select(-data), 1, median, na.rm = TRUE)

ordem_rt <- c("TQ","NP","HP","L1","BP","BN","MED")
df_rt_full <- bind_rows(lapply(ordem_rt, function(m) {
  rt_v <- if (m == "MED") rt_med else rt_mat[, m]
  tibble(data = pib$data, metodo = m,
         tempo_real = rt_v, final = finais_uni[[m]])
})) %>% filter(!is.na(tempo_real), !is.na(final))

cat("\n=== CENÁRIO 1 — REVISÕES (final − tempo real): UNIVARIADOS + MEDIANA ===\n")
cat(sprintf("%-6s %6s %6s %14s %13s\n",
            "Hiato","RAM","RRQM","corr(RT,final)","troca sinal"))
cat(strrep("-", 52), "\n")
for (m in ordem_rt) {
  d <- df_rt_full %>% filter(metodo == m)
  rev <- d$final - d$tempo_real
  cat(sprintf("%-6s %6.2f %6.2f %14.3f %12.1f%%\n", m,
              mean(abs(rev)), sqrt(mean(rev^2)),
              cor(d$tempo_real, d$final),
              100 * mean(sign(d$tempo_real) != sign(d$final))))
}
cat("Nota: vintages trimestrais 2006T1+; pseudo-tempo real (isola viés de final de amostra).\n")
cat("      HPM omitido (MLE/KFAS por vintage computacionalmente proibitivo).\n")

# g7 atualizado: RT vs final para os 6 métodos + mediana (sobrescreve g_rt)
rotulos_rt <- c(TQ  = "I — Tend. Quadrática c/ Quebras",
                NP  = "II — LOESS",
                HP  = "III — Hodrick-Prescott",
                L1  = "IV — Tendência ℓ₁",
                BP  = "VI — Band-Pass (CF)",
                BN  = "VII — Beveridge-Nelson (KMW)",
                MED = "Mediana (thick modeling)")
df_rt_plot <- df_rt_full %>%
  mutate(metodo = factor(rotulos_rt[metodo], levels = unname(rotulos_rt)))

g_rt <- ggplot(df_rt_plot, aes(x = data)) +
  geom_rect(data = recessoes %>% filter(fim >= as.Date("2006-01-01")),
            aes(xmin = inicio, xmax = fim, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.4, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_line(aes(y = final, color = "Final (amostra completa)"), linewidth = 0.7) +
  geom_line(aes(y = tempo_real, color = "Tempo real (último ponto do vintage)"),
            linewidth = 0.55, linetype = "longdash") +
  facet_wrap(~metodo, ncol = 2) +
  scale_color_manual(values = c(
    "Final (amostra completa)"               = "#1A5276",
    "Tempo real (último ponto do vintage)"   = "#C0392B")) +
  scale_x_date(breaks = seq(as.Date("2006-01-01"), max(pib$data)+365, by = "4 years"),
               date_labels = "%Y", expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(-5.5, 5.5), breaks = seq(-4, 4, by = 2),
                     labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%")) +
  labs(
    title    = "Hiato em Tempo Real vs. Final — Todos os Univariados + Mediana (WP 203)",
    subtitle = "Estimativas recursivas trimestrais desde 2006T1 | Pseudo-tempo real: sem revisões de dados, isola o viés de final de amostra",
    x = NULL, y = "Hiato (% do PIB potencial)", color = NULL,
    caption = paste0(
      "Fonte: BCB/SGS 22109. Elaboração própria, no espírito de Cusinato, Minella & Pôrto Jr. (WP 203).\n",
      "Em cada vintage τ, cada método é reestimado com dados até τ; a linha tracejada liga os últimos pontos dos vintages. ",
      "TQ: quebras entram com ≥6 trimestres no segmento final; ℓ₁: λ recalibrado por vintage; ",
      "BN: ρ reotimizado por vintage.\n",
      "HPM omitido (custo computacional). Recessões (áreas cinzas): CODACE/FGV."
    )
  ) +
  tema_padrao +
  theme(strip.text = element_text(face = "bold", size = 9))

df_phil <- tibble(data = pib$data) %>%
  bind_cols(as_tibble(rt_mat)) %>%
  mutate(MED = rt_med) %>%
  left_join(ipca_tri %>% select(data, pi4), by = "data") %>%
  mutate(pi4_fwd = lead(pi4, 4))

metodos_phil <- c("TQ","NP","HP","L1","BP","BN","MED")

# ── Previsão recursiva (janela expansiva) ────────────────────────────────────
# Origem s: estima com pares cujo alvo já realizou (s'+4 <= s); prevê pi4_{s+4}
min_est  <- 16
origens  <- which(!is.na(df_phil$pi4_fwd) & !is.na(df_phil$pi4) &
                  !is.na(df_phil$MED) &
                  pib$data >= as.Date("2011-01-01"))

prev_err <- matrix(NA_real_, length(origens), length(metodos_phil) + 1,
                   dimnames = list(NULL, c("AR", metodos_phil)))
for (k in seq_along(origens)) {
  s <- origens[k]
  est <- df_phil[1:max(1, s-4), ] %>%
    filter(!is.na(pi4_fwd), !is.na(pi4), !is.na(MED))
  if (nrow(est) < min_est) next
  # Benchmark AR
  fAR <- lm(pi4_fwd ~ pi4, data = est)
  prev_err[k, "AR"] <- df_phil$pi4_fwd[s] -
    predict(fAR, newdata = df_phil[s, ])
  # Com hiato RT
  for (m in metodos_phil) {
    estm <- est %>% filter(!is.na(.data[[m]]))
    if (nrow(estm) < min_est || is.na(df_phil[[m]][s])) next
    fm <- lm(reformulate(c("pi4", m), "pi4_fwd"), data = estm)
    prev_err[k, m] <- df_phil$pi4_fwd[s] - predict(fm, newdata = df_phil[s, ])
  }
}

# ── Métricas: RMSE relativo + Diebold-Mariano (Newey-West, h-1 lags) ─────────
dm_pval <- function(e_b, e_m, h = 4) {
  ok <- !is.na(e_b) & !is.na(e_m)
  d  <- e_b[ok]^2 - e_m[ok]^2
  nT <- length(d); if (nT < 10) return(NA)
  db <- mean(d); dc <- d - db
  s  <- mean(dc^2)
  for (l in 1:(h-1))
    s <- s + 2 * (1 - l/h) * mean(dc[-(1:l)] * dc[1:(nT-l)])
  2 * (1 - pnorm(abs(db / sqrt(s/nT))))
}

rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
rmse_AR <- rmse(prev_err[, "AR"])
n_prev  <- sum(!is.na(prev_err[, "AR"]))

cat(sprintf("\n=== CENÁRIO 2 — PHILLIPS FORA-DA-AMOSTRA (h=4, %d previsões) ===\n", n_prev))
cat(sprintf("Benchmark AR: RMSE = %.3f p.p. (inflação acumulada 4T)\n\n", rmse_AR))
cat(sprintf("%-8s %8s %10s %12s\n", "Hiato", "RMSE", "RMSE/AR", "DM p-valor"))
cat(strrep("-", 44), "\n")
tab_phil <- tibble(metodo = metodos_phil) %>%
  mutate(RMSE  = sapply(metodo, function(m) rmse(prev_err[, m])),
         rel   = RMSE / rmse_AR,
         dm_p  = sapply(metodo, function(m)
                  dm_pval(prev_err[, "AR"], prev_err[, m])))
for (i in seq_len(nrow(tab_phil)))
  cat(sprintf("%-8s %8.3f %10.3f %12.3f\n",
              tab_phil$metodo[i], tab_phil$RMSE[i],
              tab_phil$rel[i], tab_phil$dm_p[i]))
cat("\nRMSE/AR < 1: hiato agrega conteúdo preditivo. DM: H0 = mesmo RMSE (poder baixo, ~",
    n_prev, "obs).\n")

# ── Robustez: subamostra PRÉ-COVID (origens até 2019T4) ─────────────────────
# O surto inflacionário 2021-22 foi dominado por choques de oferta globais
# (inflação sem hiato por construção) e pode dominar a amostra completa
idx_pre <- which(pib$data[origens] <= as.Date("2019-10-01"))
if (length(idx_pre) >= 20) {
  rmse_AR_pre <- rmse(prev_err[idx_pre, "AR"])
  cat(sprintf("\n--- Subamostra pré-COVID (origens 2011T1-2019T4, %d previsões) ---\n",
              sum(!is.na(prev_err[idx_pre, "AR"]))))
  cat(sprintf("%-8s %10s %12s\n", "Hiato", "RMSE/AR", "DM p-valor"))
  for (m in metodos_phil) {
    cat(sprintf("%-8s %10.3f %12.3f\n", m,
                rmse(prev_err[idx_pre, m]) / rmse_AR_pre,
                dm_pval(prev_err[idx_pre, "AR"], prev_err[idx_pre, m])))
  }
}

# ── Gráfico g9 ────────────────────────────────────────────────────────────────
df_g9 <- tab_phil %>%
  mutate(metodo = factor(metodo, levels = metodos_phil),
         destaque = ifelse(metodo == "MED", "Mediana (thick modeling)",
                           "Método individual"))
g_phil <- ggplot(df_g9, aes(x = metodo, y = rel, fill = destaque)) +
  geom_col(width = 0.62, alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#C0392B", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", rel)), vjust = -0.55, size = 3.1) +
  annotate("text", x = 0.7, y = 1.012, label = "Benchmark AR (=1)",
           hjust = 0, size = 3, color = "#C0392B") +
  scale_fill_manual(values = c("Método individual" = "#5DADE2",
                               "Mediana (thick modeling)" = "#1A5276")) +
  scale_y_continuous(limits = c(0, max(1.12, max(df_g9$rel)*1.1)),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Cenário 2 — Conteúdo Preditivo para Inflação (Phillips fora-da-amostra)",
    subtitle = sprintf("RMSE relativo ao AR | π4(t+4) ~ π4(t) + hiato RT | janela expansiva, %d previsões (2011+)", n_prev),
    x = NULL, y = "RMSE / RMSE do benchmark AR", fill = NULL,
    caption = paste0(
      "Fonte: IPCA (SGS 433); hiatos em TEMPO REAL (último ponto de cada vintage recursivo desde 2006T1). ",
      "Elaboração própria, no espírito de Orphanides & van Norden (2005) e Stock & Watson (1999).\n",
      "Barras abaixo da linha tracejada: o hiato em tempo real melhora a previsão da inflação acumulada ",
      "4 trimestres à frente. MED = mediana dos 6 hiatos univariados em tempo real (thick modeling)."
    )
  ) +
  tema_padrao + theme(legend.position = "bottom")

# ════════════════════════════════════════════════════════════════════════════
# CENÁRIO 4 — EVENTO COVID: revisão retroativa do hiato pré-pandemia
# Congela cada método em 2019T4 e compara com a estimativa final:
# quanto cada método REESCREVEU a história de 2017-2019 após observar 2020?
# Métodos two-sided (HP, LOESS) devem contaminar o passado; o BN (one-sided)
# e a TQ (segmento fechado pela quebra 2020T2) devem revisar pouco.
# ════════════════════════════════════════════════════════════════════════════
cat("\nCenário 4 — Evento COVID (vintage congelado em 2019T4)...\n")

tau19 <- which(pib$data == as.Date("2019-10-01"))
y19   <- y[1:tau19]

# Séries COMPLETAS de cada método no vintage 2019T4
hiato_2019 <- list()
# TQ: quebras disponíveis até 2019T4 (2000T1, 2008T4, 2013T3)
qb19 <- idx_q[idx_q <= tau19 - 6]; kv <- length(qb19)
ini19 <- c(1, qb19); fim19 <- c(qb19 - 1, tau19)
X19 <- matrix(0, tau19, 3*(kv+1))
for (j in 1:(kv+1)) {
  tl <- ini19[j]:fim19[j]
  X19[tl, 3*(j-1)+1] <- 1; X19[tl, 3*(j-1)+2] <- tl; X19[tl, 3*(j-1)+3] <- tl^2
}
hiato_2019$TQ <- residuals(lm(y19 ~ X19 - 1)) * 100
# NP, HP, L1, BP, BN
hiato_2019$NP <- (y19 - loess_bcb(y19)) * 100
hiato_2019$HP <- as.numeric(mFilter::hpfilter(
  ts(y19, start = c(1996,1), frequency = 4), freq = 1600)$cycle) * 100
D19 <- diff(diag(tau19), differences = 2)
lam_max19 <- max(abs(solve(D19 %*% t(D19), D19 %*% y19)))
lam19 <- (0.5^5) * lam_max19
hiato_2019$L1 <- (y19 - as.numeric(l1_filter(
  y19, lam19, rho_mult = lam19 / (max(abs(D19 %*% y19)) / 5)))) * 100
hiato_2019$BP <- cf_filter(y19) * 100
hiato_2019$BN <- bn_kmw(y19, pib$data[1:tau19], p = 12)$hiato
# Mediana
mat19 <- do.call(cbind, hiato_2019)
hiato_2019$MED <- apply(mat19, 1, median, na.rm = TRUE)

# Janela de avaliação: 2017T1–2019T4
jan_aval <- which(pib$data >= as.Date("2017-01-01") & pib$data <= as.Date("2019-10-01"))

cat("\n=== CENÁRIO 4 — REVISÃO RETROATIVA DO HIATO 2017-2019 (final − vintage 2019T4) ===\n")
cat(sprintf("%-6s %12s %8s %10s %14s\n",
            "Hiato","rev. média","RAM","|rev| máx","rev. em 2019T4"))
cat(strrep("-", 56), "\n")
df_covid <- list()
for (m in ordem_rt) {
  fin  <- finais_uni[[m]][jan_aval]
  cong <- hiato_2019[[m]][jan_aval]
  rev  <- fin - cong
  cat(sprintf("%-6s %+11.2f %8.2f %10.2f %+13.2f\n", m,
              mean(rev, na.rm=TRUE), mean(abs(rev), na.rm=TRUE),
              max(abs(rev), na.rm=TRUE), rev[length(rev)]))
  df_covid[[m]] <- tibble(data = pib$data[jan_aval], metodo = m,
                          congelado = cong, final = fin)
}
cat("rev > 0: o método passou a ver a economia pré-COVID MAIS aquecida após observar 2020.\n")

# Gráfico g10
rotulos_cov <- c(TQ="I — Tend. Quadrática", NP="II — LOESS", HP="III — HP",
                 L1="IV — Tendência ℓ₁", BP="VI — Band-Pass",
                 BN="VII — Beveridge-Nelson", MED="Mediana")
df_cov_plot <- bind_rows(df_covid) %>%
  mutate(metodo = factor(rotulos_cov[metodo], levels = unname(rotulos_cov)))

g_covid <- ggplot(df_cov_plot, aes(x = data)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_ribbon(aes(ymin = pmin(congelado, final), ymax = pmax(congelado, final)),
              fill = "#F5B7B1", alpha = 0.5) +
  geom_line(aes(y = congelado, color = "Vintage 2019T4 (sem COVID na amostra)"),
            linewidth = 0.75) +
  geom_line(aes(y = final, color = "Final (com COVID na amostra)"),
            linewidth = 0.75, linetype = "longdash") +
  facet_wrap(~metodo, ncol = 2) +
  scale_color_manual(values = c(
    "Vintage 2019T4 (sem COVID na amostra)" = "#1A5276",
    "Final (com COVID na amostra)"          = "#C0392B")) +
  scale_x_date(date_labels = "%Y", breaks = seq(as.Date("2017-01-01"),
               as.Date("2020-01-01"), by = "1 year")) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%")) +
  labs(
    title    = "Evento COVID — Quanto Cada Método Reescreveu a História Pré-Pandemia?",
    subtitle = "Hiato 2017-2019 estimado no vintage 2019T4 vs. com a amostra completa | área rosa = revisão retroativa",
    x = NULL, y = "Hiato (% do PIB potencial)", color = NULL,
    caption = paste0(
      "Fonte: BCB/SGS 22109. Elaboração própria. O choque de 2020 não deveria alterar a leitura do estado da economia ",
      "em 2017-2019;\na revisão retroativa mede a contaminação do passado pelo final da amostra em cada método. ",
      "TQ no vintage 2019T4 usa quebras até 2013T3.\n",
      "BN é one-sided por construção (esperança condicional com informação passada) — referência de robustez."
    )
  ) +
  tema_padrao +
  theme(strip.text = element_text(face = "bold", size = 9))

# ════════════════════════════════════════════════════════════════════════════
# TABELA-SÍNTESE — CENÁRIO × MÉTODO (a resposta à pergunta do TCC)
# C1 Política monetária em tempo real:  RRQM das revisões (menor = melhor)
# C2 Previsão de inflação:              RMSE/AR pré-COVID   (menor = melhor)
# C3 Análise histórica:                 % recessões CODACE com hiato em queda
#                                       na série final       (maior = melhor)
# C4 Robustez a rupturas (COVID):       RAM da revisão retroativa 2017-19
#                                       (menor = melhor)
# ════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat("              TABELA-SÍNTESE: CENÁRIO × MÉTODO\n")
cat("════════════════════════════════════════════════════════════════════\n")

# C1: RRQM das revisões
c1 <- sapply(ordem_rt, function(m) {
  d <- df_rt_full %>% filter(metodo == m)
  sqrt(mean((d$final - d$tempo_real)^2))
})

# C2: RMSE/AR pré-COVID (recomputado defensivamente)
idx_pre2 <- which(pib$data[origens] <= as.Date("2019-10-01"))
rmse_AR_p <- rmse(prev_err[idx_pre2, "AR"])
c2 <- sapply(ordem_rt, function(m) rmse(prev_err[idx_pre2, m]) / rmse_AR_p)

# C3: captura de recessões CODACE (série FINAL, amostra completa)
fase_rec <- sapply(pib$data, function(d)
  any(d >= recessoes$inicio & d <= recessoes$fim))
c3 <- sapply(ordem_rt, function(m) {
  h  <- finais_uni[[m]]
  dh <- c(NA, diff(h))
  ok <- fase_rec & !is.na(dh)
  100 * mean(dh[ok] < 0)
})

# C4: RAM da revisão retroativa COVID
c4 <- sapply(ordem_rt, function(m) {
  mean(abs(finais_uni[[m]][jan_aval] - hiato_2019[[m]][jan_aval]), na.rm = TRUE)
})

sintese <- tibble(metodo = ordem_rt,
                  C1_RRQM = c1, C2_RMSErel = c2, C3_CODACE = c3, C4_RAMcov = c4) %>%
  mutate(r1 = rank(C1_RRQM), r2 = rank(C2_RMSErel),
         r3 = rank(-C3_CODACE), r4 = rank(C4_RAMcov),
         rank_medio = (r1 + r2 + r3 + r4) / 4)

cat(sprintf("%-6s | %-14s | %-14s | %-14s | %-14s\n",
            "", "C1 Tempo real", "C2 Prev.infl.", "C3 Histórica", "C4 Rupturas"))
cat(sprintf("%-6s | %-14s | %-14s | %-14s | %-14s\n",
            "Hiato", "RRQM (rank)", "RMSE/AR (rank)", "%CODACE (rank)", "RAM (rank)"))
cat(strrep("-", 76), "\n")
for (i in seq_len(nrow(sintese))) {
  s <- sintese[i, ]
  cat(sprintf("%-6s | %6.2f    (%d°) | %7.3f   (%d°) | %6.1f%%   (%d°) | %6.2f    (%d°)\n",
              s$metodo, s$C1_RRQM, s$r1, s$C2_RMSErel, s$r2,
              s$C3_CODACE, s$r3, s$C4_RAMcov, s$r4))
}
cat(strrep("-", 76), "\n")
venc <- sintese$metodo[c(which.min(sintese$C1_RRQM), which.min(sintese$C2_RMSErel),
                          which.max(sintese$C3_CODACE), which.min(sintese$C4_RAMcov))]
cat(sprintf("VENCEDOR POR CENÁRIO:  C1: %s | C2: %s | C3: %s | C4: %s\n",
            venc[1], venc[2], venc[3], venc[4]))
cat(sprintf("Melhor rank médio: %s (%.2f)\n",
            sintese$metodo[which.min(sintese$rank_medio)], min(sintese$rank_medio)))

# ── Heatmap de rankings (g11) ────────────────────────────────────────────────
df_heat <- sintese %>%
  select(metodo, r1, r2, r3, r4) %>%
  pivot_longer(-metodo, names_to = "cenario", values_to = "rank") %>%
  mutate(cenario = recode(cenario,
           r1 = "C1\nTempo real\n(RRQM)",
           r2 = "C2\nPrevisão inflação\n(RMSE/AR pré-COVID)",
           r3 = "C3\nAnálise histórica\n(% recessões CODACE)",
           r4 = "C4\nRupturas/COVID\n(RAM retroativa)"),
         metodo = factor(metodo, levels = rev(ordem_rt)))

g_sintese <- ggplot(df_heat, aes(x = cenario, y = metodo, fill = rank)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = paste0(rank, "°")), size = 4.2, fontface = "bold",
            color = "grey15") +
  scale_fill_gradient2(low = "#1E8449", mid = "#F9E79F", high = "#C0392B",
                       midpoint = 4, name = "Ranking\n(1 = melhor)") +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Síntese — Qual Método para Qual Cenário?",
    subtitle = "Ranking de cada método de hiato em quatro cenários de uso | 1° = melhor desempenho no critério",
    x = NULL, y = NULL,
    caption = paste0(
      "Elaboração própria. C1: RRQM das revisões em pseudo-tempo real (vintages 2006T1+). ",
      "C2: RMSE relativo ao AR na previsão da inflação acumulada 4T à frente (pré-COVID).\n",
      "C3: %% de trimestres de recessão CODACE com hiato final em queda (1996+). ",
      "C4: RAM da revisão retroativa do hiato 2017-19 após o choque COVID.\n",
      "MED = mediana dos 6 univariados (thick modeling)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "#555555"),
        plot.caption = element_text(size = 7.5, color = "#777777", hjust = 0),
        axis.text = element_text(size = 9.5, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right")

# ════════════════════════════════════════════════════════════════════════════
# SALVAR E EXIBIR
# ════════════════════════════════════════════════════════════════════════════
dir_saida <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                      error = function(e) getwd())
if (is.null(dir_saida) || nchar(dir_saida) < 2) dir_saida <- getwd()

arq_uni    <- file.path(dir_saida, "g1_univariados.png")
arq_fp     <- file.path(dir_saida, "g2_multi_funcao_producao.png")
arq_pca    <- file.path(dir_saida, "g3_multi_pca.png")
arq_areosa <- file.path(dir_saida, "g4_multi_areosa.png")
arq_jl     <- file.path(dir_saida, "g5_multi_jl.png")
arq_consol <- file.path(dir_saida, "g6_consolidado.png")
arq_rt     <- file.path(dir_saida, "g7_tempo_real.png")

ggsave(arq_uni,    g_uni,    width=13, height=7,   dpi=180, bg="white")
ggsave(arq_fp,     g_multi,  width=12, height=6.5, dpi=180, bg="white")
ggsave(arq_pca,    g_pca,    width=12, height=6.5, dpi=180, bg="white")
ggsave(arq_areosa, g_areosa, width=12, height=6.5, dpi=180, bg="white")
ggsave(arq_jl,     g_jl,     width=12, height=6.5, dpi=180, bg="white")
ggsave(arq_consol, g_consol, width=13, height=7,   dpi=180, bg="white")
ggsave(arq_rt,     g_rt,     width=12, height=11,  dpi=180, bg="white")
if (!is.null(df_cbo)) {
  arq_cbo <- file.path(dir_saida, "g8_multi_cbo.png")
  ggsave(arq_cbo, g_cbo, width=12, height=6.5, dpi=180, bg="white")
  cat(sprintf("  %s\n", arq_cbo))
}
arq_phil <- file.path(dir_saida, "g9_phillips_oos.png")
ggsave(arq_phil, g_phil, width=11, height=6.5, dpi=180, bg="white")
cat(sprintf("  %s\n", arq_phil))
arq_cov <- file.path(dir_saida, "g10_evento_covid.png")
ggsave(arq_cov, g_covid, width=12, height=11, dpi=180, bg="white")
cat(sprintf("  %s\n", arq_cov))
arq_sint <- file.path(dir_saida, "g11_sintese.png")
ggsave(arq_sint, g_sintese, width=11.5, height=6.5, dpi=180, bg="white")
cat(sprintf("  %s\n", arq_sint))
cat(sprintf("\nGráficos salvos:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n",
            arq_uni, arq_fp, arq_pca, arq_areosa, arq_jl, arq_consol, arq_rt))

print(g_uni)
print(g_multi)
print(g_pca)
print(g_areosa)
print(g_jl)
print(g_consol)
print(g_rt)
if (!is.null(df_cbo)) print(g_cbo)
print(g_phil)
print(g_covid)
print(g_sintese)

# ════════════════════════════════════════════════════════════════════════════
# ESTATÍSTICAS COMPARATIVAS
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== GRUPO I — UNIVARIADOS ===\n")
cat(sprintf("%-45s %8s %8s %8s\n","Método","Atual","Média","DP"))
cat(strrep("-", 73), "\n")
for (nm in names(metodo_uni)) {
  x <- pib[[nm]]
  cat(sprintf("%-45s %+7.2f%% %+7.2f%% %7.2f%%\n",
              metodo_uni[nm], tail(na.omit(x),1), mean(x,na.rm=TRUE), sd(x,na.rm=TRUE)))
}

cat("\n=== GRUPO II — MULTIVARIADO I ===\n")
cat(sprintf("Período: %s a %s\n", min(df_multi$data), max(df_multi$data)))
cat(sprintf("Hiato FP atual:  %+.2f p.p.\n", tail(df_multi$hiato_FP, 1)))
cat(sprintf("Média:           %+.2f p.p.\n", mean(df_multi$hiato_FP)))
cat(sprintf("DP:              %.2f p.p.\n",  sd(df_multi$hiato_FP)))

cat("\n=== GRUPO II — MULTIVARIADO V: PCA ===\n")
cat(sprintf("Período: %s a %s (%d obs.)\n",
            min(df_pca_raw$data), max(df_pca_raw$data), nrow(df_pca_raw)))
cat(sprintf("Variância explicada PC1: %.1f%%\n", var_exp[1]))
cat(sprintf("Hiato PCA atual:  %+.2f p.p.\n", tail(df_pca_raw$hiato_PCA, 1)))
cat(sprintf("Média:            %+.2f p.p.\n", mean(df_pca_raw$hiato_PCA)))
cat(sprintf("DP:               %.2f p.p.\n",  sd(df_pca_raw$hiato_PCA)))

cat("\n=== GRUPO II — MULTIVARIADO II: AREOSA (2008) ===\n")
cat(sprintf("Período: %s a %s\n", min(df_areosa$data), max(df_areosa$data)))
cat(sprintf("Hiato Areosa atual: %+.2f p.p.\n", tail(df_areosa$hiato_AR, 1)))
cat(sprintf("Média:              %+.2f p.p.\n",  mean(df_areosa$hiato_AR)))
cat(sprintf("DP:                 %.2f p.p.\n",   sd(df_areosa$hiato_AR)))
cat(sprintf("Corr(Areosa, FP):   %.3f\n",
            cor(df_areosa$hiato_AR, df_multi$hiato_FP)))

cat("\n=== GRUPO II — MULTIVARIADO IV: JAROCINSKI & LENZA ===\n")
cat(sprintf("Período: %s a %s (%d obs.)\n",
            min(df_jl$data), max(df_jl$data), nrow(df_jl)))
cat(sprintf("Hiato JL atual:  %+.2f p.p.\n", tail(df_jl$hiato_JL,1)))
cat(sprintf("Média:           %+.2f p.p.\n",  mean(df_jl$hiato_JL)))
cat(sprintf("DP:              %.2f p.p.\n",   sd(df_jl$hiato_JL)))
cat(sprintf("Corr(JL, PCA):   %.3f\n", cor(df_jl$hiato_JL, df_pca_raw$hiato_PCA)))
cat(sprintf("Corr(JL, FP):    %.3f\n", cor(df_jl$hiato_JL, df_multi$hiato_FP)))

if (!is.null(df_cbo)) {
  cat("\n=== GRUPO II — MULTIVARIADO III: CBO ===\n")
  cat(sprintf("Período: %s a %s (%d obs.)\n",
              min(df_cbo$data), max(df_cbo$data), nrow(df_cbo)))
  cat(sprintf("Hiato CBO atual: %+.2f p.p.\n", tail(df_cbo$hiato_CBO,1)))
  cat(sprintf("Média:           %+.2f p.p.\n", mean(df_cbo$hiato_CBO, na.rm=TRUE)))
  cat(sprintf("DP:              %.2f p.p.\n",  sd(df_cbo$hiato_CBO, na.rm=TRUE)))
}

# ════════════════════════════════════════════════════════════════════════════
# CORRELAÇÕES CRUZADAS — período comum (2012T2+), estilo Tabela 2 do RI
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== CORRELAÇÕES CRUZADAS — UNIVARIADOS × MULTIVARIADOS (período comum) ===\n")
df_corr_comum <- pib %>%
  select(data, TQ=hiato_TQ, NP=hiato_NP, HP=hiato_HP, L1=hiato_L1,
         HPM=hiato_HPM, BP=hiato_BP, BN=hiato_BN) %>%
  inner_join(df_multi   %>% select(data, FP  = hiato_FP),  by="data") %>%
  inner_join(df_pca_raw %>% select(data, PCA = hiato_PCA), by="data") %>%
  inner_join(df_jl      %>% select(data, JL  = hiato_JL),  by="data")
if (!is.null(df_cbo))
  df_corr_comum <- df_corr_comum %>%
    inner_join(df_cbo %>% select(data, CBO = hiato_CBO), by="data")

cat(sprintf("Período: %s a %s (%d obs.)\n\n",
            min(df_corr_comum$data), max(df_corr_comum$data), nrow(df_corr_comum)))
print(round(cor(df_corr_comum %>% select(-data)), 3))
cat("\nNota: Areosa omitido — comb. convexa 0,658·FP+0,342·HP (Proposição 1; corr 0,979 c/ FP),\n",
    "     não altera a faixa mín–máx. Referência de comparação: Tabela 2, BCB RI Jun/2024.\n")
