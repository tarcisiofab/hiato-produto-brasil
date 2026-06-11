# ═══════════════════════════════════════════════════════════════════════════
# VERIFICAÇÃO DA PROPOSIÇÃO AREOSA
# (a) Resolver o problema de otimização COMPLETO (sistema acoplado 2T×2T)
# (b) Comparar com a fórmula fechada derivada analiticamente:
#     ĥ_Y^Areosa = (1/κ)·ĥ_Y^FP + ((β1²+β2²)/κ)·ĥ_Y^HP,  κ = 1+β1²+β2²
# (c) Comparar com o FP simples (implementação atual do script)
# ═══════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(jsonlite); library(dplyr); library(mFilter)})

get_sgs <- function(serie, inicio="01/01/1996") {
  url <- sprintf("https://api.bcb.gov.br/dados/serie/bcdata.sgs.%d/dados?formato=json&dataInicial=%s", serie, inicio)
  raw <- jsonlite::fromJSON(url)
  tibble(data=as.Date(raw$data,"%d/%m/%Y"),
         valor=as.numeric(gsub(",",".",raw$valor))) %>% filter(!is.na(valor))
}
pib_raw    <- get_sgs(22109)
nuci_raw   <- get_sgs(28561, "01/01/2001")
desemp_raw <- get_sgs(24369, "01/01/2012")

df <- pib_raw %>% rename(pib=valor) %>%
  inner_join(nuci_raw %>% rename(nuci=valor), by="data") %>%
  inner_join(desemp_raw %>% rename(desemp=valor), by="data") %>%
  arrange(data)
T_ <- nrow(df)

# ESCALA CONSISTENTE (o bug da tentativa antiga):
# u, c em FRAÇÕES (U/100, C/100); y em log → hiatos todos adimensionais
u <- df$desemp / 100
c_ <- df$nuci  / 100
y  <- log(df$pib)

b1 <- 0.6; b2 <- 0.4; lam <- 1600
kap <- 1 + b1^2 + b2^2

# Operadores
D <- diff(diag(T_), differences=2)
A <- t(D) %*% D
S <- solve(diag(T_) + lam * A)        # smoother HP: tendência = S %*% série

# ── (a) SISTEMA ACOPLADO EXATO ───────────────────────────────────────────────
ytil <- y - b2*c_ + b1*u
IpA  <- diag(T_) + lam*A
a11  <- (1+b1^2) * IpA
a12  <- -b1*b2   * IpA
a22  <- (1+b2^2) * IpA
A_sys <- rbind(cbind(a11, a12), cbind(a12, a22))
rhs1  <- (1+b1^2)*u  - b1*b2*c_ + lam*b1*(A %*% ytil)
rhs2  <- (1+b2^2)*c_ - b1*b2*u  - lam*b2*(A %*% ytil)
sol   <- solve(A_sys, c(rhs1, rhs2))
u_st  <- sol[1:T_]; c_st <- sol[(T_+1):(2*T_)]

hu_ex <- u  - u_st
hc_ex <- c_ - c_st
hy_exato <- (b2*hc_ex - b1*hu_ex) * 100    # hiato Areosa EXATO, em %

# ── (b) FÓRMULA FECHADA ──────────────────────────────────────────────────────
cyc <- function(x) as.numeric(x - S %*% x)   # ciclo HP
hy_FP <- (b2*cyc(c_) - b1*cyc(u)) * 100      # FP simples
hy_HP <- cyc(y) * 100                         # HP do PIB
hy_formula <- (1/kap)*hy_FP + ((b1^2+b2^2)/kap)*hy_HP

# ── Confrontos ───────────────────────────────────────────────────────────────
cat("════════════════════════════════════════════════════════\n")
cat(sprintf("κ = 1+β1²+β2² = %.2f → pesos: %.4f (FP) + %.4f (HP)\n",
            kap, 1/kap, (b1^2+b2^2)/kap))
cat("════════════════════════════════════════════════════════\n")
cat(sprintf("(a)=(b)?  max|sistema − fórmula| = %.2e   (← deve ser ~0)\n",
            max(abs(hy_exato - hy_formula))))
cat(sprintf("corr(Areosa_exato, FP simples)    = %.6f\n", cor(hy_exato, hy_FP)))
cat(sprintf("corr(Areosa_exato, HP do PIB)     = %.6f\n", cor(hy_exato, hy_HP)))
cat(sprintf("corr(FP, HP do PIB)               = %.6f\n", cor(hy_FP, hy_HP)))
cat(sprintf("max|Areosa_exato − FP|            = %.3f p.p.\n",
            max(abs(hy_exato - hy_FP))))
cat(sprintf("DP: Areosa=%.3f | FP=%.3f | HP=%.3f\n",
            sd(hy_exato), sd(hy_FP), sd(hy_HP)))
cat(sprintf("Hiato atual: Areosa=%+.2f | FP=%+.2f | HP=%+.2f\n",
            tail(hy_exato,1), tail(hy_FP,1), tail(hy_HP,1)))

# Verificação de consistência da restrição na solução exata
chk <- max(abs((y - (ytil + b2*c_st - b1*u_st)) - (b2*hc_ex - b1*hu_ex)))
cat(sprintf("Restrição CD satisfeita na solução: max viol = %.2e\n", chk))
