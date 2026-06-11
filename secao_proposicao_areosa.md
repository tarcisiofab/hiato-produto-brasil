# Uma solução fechada para o método de Areosa (2008)

> *Seção redigida para o TCC — numeração de seções, equações e referências a ajustar ao template final. A verificação numérica reportada na Tabela 1 é reproduzível pelo script `verif_areosa.R`.*

## Motivação

Entre os métodos multivariados empregados pelo Banco Central do Brasil na estimação do hiato do produto (BCB, 2024, Apêndice 2b), o método de Areosa (2008) — doravante MII.II — distingue-se por estimar *conjuntamente* os componentes potenciais do desemprego, da utilização da capacidade instalada (Nuci) e do produto, mediante três filtros de Hodrick-Prescott interligados pela restrição de uma função de produção Cobb-Douglas. À primeira vista, a otimização conjunta sugere que o método extrai informação adicional relativamente à combinação simples dos hiatos dos fatores (método MII.I), na qual os potenciais de cada série são estimados por filtros HP independentes.

Esta seção demonstra que essa impressão é, em sentido preciso, falsa: sob parâmetro de suavização comum, o problema de Areosa admite **solução analítica fechada**, e o hiato do produto resultante é uma **média ponderada exata, com pesos fixos**, de dois estimadores que já integram o conjunto de métodos do BCB — o hiato da função de produção simples (MII.I) e o hiato HP do próprio PIB (método III). O resultado tem três consequências práticas, formalizadas nos Corolários 1 a 3: o método não adiciona informação ao conjunto; não pode, por construção, ampliar a faixa de estimativas do *thick modeling*; e a divergência entre MII.I e MII.II reportada pelo BCB deve ser atribuída integralmente aos insumos de dados, e não à estrutura de estimação.

## Notação e o problema

Sejam $u, c, y \in \mathbb{R}^{T}$ os vetores das séries observadas de desemprego, Nuci e log-PIB, respectivamente, em unidades tais que a restrição da função de produção seja dimensionalmente consistente (na implementação: $u = U/100$ e $c = C/100$ em frações; $y$ em logaritmo natural). Sejam $u^{*}, c^{*}, y^{*}$ os respectivos componentes potenciais, e defina os hiatos $h_{u} \equiv u - u^{*}$, $h_{c} \equiv c - c^{*}$ e $h_{y} \equiv y - y^{*}$.

Seja $D \in \mathbb{R}^{(T-2)\times T}$ a matriz de segundas diferenças, com linhas $(\,\dots, 1, -2, 1, \dots)$, e defina $A \equiv D'D$. O filtro de Hodrick-Prescott com parâmetro $\lambda$ aplicado a uma série $x$ entrega a tendência $S x$ e o ciclo $(I - S)x$, em que

$$S \equiv (I + \lambda A)^{-1}$$

é o operador linear (suavizador) do filtro. Duas propriedades elementares de $S$ serão usadas: a **linearidade** e a identidade

$$\lambda S A = I - S, \tag{1}$$

obtida diretamente de $S(I + \lambda A) = I$.

O problema de Areosa (2008), na formulação do BCB (2024, Apêndice 2b), é o programa quadrático

$$\min_{u^{*},\,c^{*},\,y^{*}}\;\; \lVert u - u^{*}\rVert^{2} + \lVert c - c^{*}\rVert^{2} + \lVert y - y^{*}\rVert^{2} + \lambda\left(\lVert D u^{*}\rVert^{2} + \lVert D c^{*}\rVert^{2} + \lVert D y^{*}\rVert^{2}\right) \tag{2}$$

sujeito à restrição da função de produção, imposta em todos os períodos:

$$h_{y} = \beta_{2}\, h_{c} - \beta_{1}\, h_{u}, \qquad \beta_{1}, \beta_{2} > 0. \tag{3}$$

As hipóteses mantidas são:

- **(H1)** o parâmetro de suavização $\lambda$ é comum às três séries;
- **(H2)** os três termos de fidelidade em (2) têm pesos unitários nas unidades adotadas;
- **(H3)** a restrição (3) vale exatamente em todo $t$, com pesos $\beta_{1}, \beta_{2}$ fixos.

Defina, por fim, os dois estimadores de referência:

$$\hat{h}^{FP} \equiv \beta_{2}(I-S)c - \beta_{1}(I-S)u \qquad \text{(função de produção simples, MII.I)}$$
$$\hat{h}^{HP} \equiv (I-S)y \qquad \text{(hiato HP do PIB, método III)}$$

## Proposição 1

**Proposição 1.** *Sob (H1)–(H3), o problema (2)–(3) tem solução única, e o hiato do produto resultante é*

$$\hat{h}^{Areosa} \;=\; \frac{1}{\kappa}\,\hat{h}^{FP} \;+\; \frac{\beta_{1}^{2}+\beta_{2}^{2}}{\kappa}\,\hat{h}^{HP}, \qquad \kappa \equiv 1 + \beta_{1}^{2} + \beta_{2}^{2}. \tag{4}$$

*Em particular, com os pesos do BCB ($\beta_{1} = 0{,}6$; $\beta_{2} = 0{,}4$), tem-se $\kappa = 1{,}52$ e*

$$\hat{h}^{Areosa} = 0{,}6579\,\hat{h}^{FP} + 0{,}3421\,\hat{h}^{HP}.$$

*Adicionalmente, os potenciais ótimos dos fatores são*

$$u^{*} = S u + \frac{\beta_{1}}{\kappa}(I-S)\tilde{y}, \qquad c^{*} = S c - \frac{\beta_{2}}{\kappa}(I-S)\tilde{y}, \tag{5}$$

*em que $\tilde{y} \equiv y - \beta_{2} c + \beta_{1} u$ — isto é, o filtro HP univariado de cada fator mais uma correção proporcional ao ciclo HP da série composta $\tilde{y}$.*

### Demonstração

**Passo 1 (eliminação de $y^{*}$).** A restrição (3) determina $y^{*}$ como função afim de $(u^{*}, c^{*})$:

$$y^{*} = y - \beta_{2}(c - c^{*}) + \beta_{1}(u - u^{*}) = \tilde{y} + \beta_{2} c^{*} - \beta_{1} u^{*}.$$

Substituindo em (2), o objetivo torna-se função apenas de $(u^{*}, c^{*})$:

$$J(u^{*}, c^{*}) = \lVert h_{u}\rVert^{2} + \lVert h_{c}\rVert^{2} + \lVert \beta_{2} h_{c} - \beta_{1} h_{u}\rVert^{2} + \lambda\lVert D u^{*}\rVert^{2} + \lambda\lVert D c^{*}\rVert^{2} + \lambda\big\lVert D\big(\tilde{y} + \beta_{2} c^{*} - \beta_{1} u^{*}\big)\big\rVert^{2}.$$

Note que o termo $\lVert y - y^{*}\rVert^{2}$ **não desaparece** com a substituição: ele equivale a $\lVert \beta_{2} h_{c} - \beta_{1} h_{u}\rVert^{2}$, o quadrado do próprio hiato do produto — fonte do acoplamento entre as equações de $u^{*}$ e $c^{*}$.

**Passo 2 (condições de primeira ordem).** Usando $\partial h_{u}/\partial u^{*} = -I$, $\partial h_{y}/\partial u^{*} = +\beta_{1} I$ e $\partial y^{*}/\partial u^{*} = -\beta_{1} I$ (e análogos para $c^{*}$), as CPOs são

$$-h_{u} + \beta_{1} h_{y} + \lambda A u^{*} - \lambda \beta_{1} A y^{*} = 0, \tag{6a}$$
$$-h_{c} - \beta_{2} h_{y} + \lambda A c^{*} + \lambda \beta_{2} A y^{*} = 0. \tag{6b}$$

**Passo 3 (o sistema fatora).** Substituindo $h_{u}$, $h_{c}$, $h_{y}$ e $y^{*}$ por suas expressões em $(u^{*}, c^{*})$ e agrupando termos, todos os blocos de coeficientes fatoram no operador $(I + \lambda A)$:

$$\underbrace{\begin{pmatrix} (1+\beta_{1}^{2}) & -\beta_{1}\beta_{2} \\ -\beta_{1}\beta_{2} & (1+\beta_{2}^{2}) \end{pmatrix}}_{\textstyle B} \otimes\, (I + \lambda A) \begin{pmatrix} u^{*} \\ c^{*} \end{pmatrix} = \begin{pmatrix} (1+\beta_{1}^{2})\,u - \beta_{1}\beta_{2}\,c + \lambda\beta_{1} A \tilde{y} \\ (1+\beta_{2}^{2})\,c - \beta_{1}\beta_{2}\,u - \lambda\beta_{2} A \tilde{y} \end{pmatrix}. \tag{7}$$

Como $\det B = (1+\beta_{1}^{2})(1+\beta_{2}^{2}) - \beta_{1}^{2}\beta_{2}^{2} = 1 + \beta_{1}^{2} + \beta_{2}^{2} = \kappa > 0$ e $(I+\lambda A) \succ 0$, a Hessiana $2\,[B \otimes (I+\lambda A)]$ é definida positiva: $J$ é estritamente convexo, a solução existe, é única, e as CPOs são necessárias e suficientes.

**Passo 4 (inversão e identidade do HP).** Multiplicando (7) por $B^{-1} \otimes S$, com

$$B^{-1} = \frac{1}{\kappa}\begin{pmatrix} 1+\beta_{2}^{2} & \beta_{1}\beta_{2} \\ \beta_{1}\beta_{2} & 1+\beta_{1}^{2} \end{pmatrix},$$

e aplicando a identidade (1) — que converte cada termo $\lambda S A \tilde{y}$ em $(I - S)\tilde{y}$, o **ciclo HP** de $\tilde{y}$ — os coeficientes colapsam: na equação de $u^{*}$, o coeficiente de $S u$ é $[(1+\beta_{2}^{2})(1+\beta_{1}^{2}) - \beta_{1}^{2}\beta_{2}^{2}]/\kappa = 1$; o de $S c$ é $[-(1+\beta_{2}^{2})\beta_{1}\beta_{2} + \beta_{1}\beta_{2}(1+\beta_{2}^{2})]/\kappa = 0$; e o de $(I-S)\tilde{y}$ é $[(1+\beta_{2}^{2})\beta_{1} - \beta_{1}\beta_{2}^{2}]/\kappa = \beta_{1}/\kappa$. Cálculo análogo vale para $c^{*}$. Obtém-se (5).

**Passo 5 (o hiato do produto).** De (5),

$$h_{u} = (I-S)u - \frac{\beta_{1}}{\kappa}(I-S)\tilde{y}, \qquad h_{c} = (I-S)c + \frac{\beta_{2}}{\kappa}(I-S)\tilde{y},$$

e, pela restrição (3),

$$\hat{h}^{Areosa} = \beta_{2} h_{c} - \beta_{1} h_{u} = \hat{h}^{FP} + \frac{\beta_{1}^{2}+\beta_{2}^{2}}{\kappa}\,(I-S)\tilde{y}.$$

Pela linearidade do operador $(I-S)$ e pela definição de $\tilde{y}$,

$$(I-S)\tilde{y} = (I-S)y - \beta_{2}(I-S)c + \beta_{1}(I-S)u = \hat{h}^{HP} - \hat{h}^{FP},$$

donde

$$\hat{h}^{Areosa} = \left(1 - \frac{\beta_{1}^{2}+\beta_{2}^{2}}{\kappa}\right)\hat{h}^{FP} + \frac{\beta_{1}^{2}+\beta_{2}^{2}}{\kappa}\,\hat{h}^{HP} = \frac{1}{\kappa}\,\hat{h}^{FP} + \frac{\beta_{1}^{2}+\beta_{2}^{2}}{\kappa}\,\hat{h}^{HP},$$

usando $\kappa - (\beta_{1}^{2}+\beta_{2}^{2}) = 1$. $\blacksquare$

## Corolários

**Corolário 1 (conteúdo informacional).** *Sob (H1)–(H3), o estimador de Areosa é combinação linear, com pesos fixos e conhecidos ex-ante, de dois estimadores já pertencentes ao conjunto de métodos do BCB. Sua adição a um conjunto que contenha MII.I e o método III não acrescenta informação.*

**Corolário 2 (thick modeling).** *Como os pesos em (4) são positivos e somam um, $\hat{h}^{Areosa}_{t}$ é combinação convexa de $\hat{h}^{FP}_{t}$ e $\hat{h}^{HP}_{t}$ em cada $t$; logo $\min\{\hat{h}^{FP}_{t}, \hat{h}^{HP}_{t}\} \le \hat{h}^{Areosa}_{t} \le \max\{\hat{h}^{FP}_{t}, \hat{h}^{HP}_{t}\}$. A inclusão do método de Areosa jamais amplia a faixa mínimo–máximo de um conjunto de estimativas que já contenha MII.I e o filtro HP do PIB.*

**Corolário 3 (origem da divergência MII.I vs. MII.II).** *Qualquer divergência entre os métodos MII.I e MII.II além da implicada por (4) — como a reportada no Gráfico 2 do boxe do BCB (2024) — não pode decorrer da estrutura de otimização conjunta sob $\lambda$ comum. Decorre, necessariamente, de diferenças nos insumos (em particular, a série de desemprego da PNAD Contínua retropolada internamente; Alves e Fasolo, BCB Working Paper 400) ou do relaxamento de (H1)–(H2) ($\lambda$ ou pesos de fidelidade específicos por série).*

## Observações

**(i) Papel das hipóteses.** A fatoração no Passo 3 depende crucialmente de (H1): com $\lambda$ distintos por série, os blocos envolvem operadores $(I + \lambda_{x} A)$ diferentes, a estrutura de Kronecker se perde e a forma fechada (4) deixa de valer — embora o sistema permaneça linear e solúvel numericamente. (H2) entra na mesma posição: pesos de fidelidade $w_{u}, w_{c}, w_{y}$ distintos alteram a matriz $B$ e, com ela, os pesos da combinação. A escala das séries só importa por intermédio de (H2): a derivação nunca usa as unidades, de modo que (4) vale para qualquer convenção, desde que a mesma seja adotada no cômputo de $\hat{h}^{FP}$.

**(ii) Casos-limite e leitura econômica.** O peso do hiato HP do PIB, $(\beta_{1}^{2}+\beta_{2}^{2})/\kappa$, é crescente em $\beta_{1}$ e $\beta_{2}$. Quando $\beta_{1}, \beta_{2} \to 0$, o método colapsa na função de produção simples; quando $\beta_{1}, \beta_{2} \to \infty$, colapsa no filtro HP do próprio PIB. A intuição é direta: os $\beta$ governam o quanto o bloco do produto — fidelidade $\lVert y - y^{*}\rVert^{2}$ e suavidade $\lambda\lVert D y^{*}\rVert^{2}$, que isoladamente definem um filtro HP em $y$ — pesa na função-perda relativamente aos blocos dos fatores. O método de Areosa, portanto, **interpola** entre um estimador baseado nos fatores e um estimador baseado no produto, com a posição da interpolação fixada pelos coeficientes da função de produção.

**(iii) Sobre a implementação.** A equação (4) reduz o custo computacional do método a dois filtros HP adicionais, dispensando a solução do sistema $2T \times 2T$ — útil, por exemplo, em exercícios recursivos de tempo real.

## Verificação numérica

A Proposição 1 foi verificada confrontando a fórmula (4) com a solução numérica direta do sistema (7), nos dados utilizados neste trabalho (PIB: SGS 22109; Nuci: SGS 28561; desemprego PNADC: SGS 24369; amostra trimestral 2012T2–2026T1; $\lambda = 1600$; $\beta_{1}=0{,}6$, $\beta_{2}=0{,}4$). Os resultados constam da Tabela 1.

**Tabela 1 — Verificação numérica da Proposição 1**

| Quantidade | Valor |
|---|---|
| $\max_{t}\,\lvert \text{sistema (7)} - \text{fórmula (4)} \rvert$ | $4{,}55 \times 10^{-11}$ p.p. |
| Violação máxima da restrição (3) na solução | $1{,}32 \times 10^{-15}$ |
| $\text{corr}(\hat{h}^{Areosa}, \hat{h}^{FP})$ | $0{,}979$ |
| $\text{corr}(\hat{h}^{Areosa}, \hat{h}^{HP})$ | $0{,}963$ |
| $\max_{t}\,\lvert \hat{h}^{Areosa}_{t} - \hat{h}^{FP}_{t}\rvert$ | $1{,}01$ p.p. |
| Desvio-padrão: Areosa / FP / HP | $1{,}57$ / $1{,}39$ / $2{,}03$ p.p. |

*Fonte: elaboração própria (script `verif_areosa.R`). A coincidência entre o sistema completo e a fórmula fechada à precisão de máquina confirma a derivação; as demais linhas quantificam a relação do estimador de Areosa com seus dois componentes.*

A primeira linha estabelece que a fórmula fechada reproduz a solução exata do problema completo à precisão de máquina. As linhas seguintes ilustram o conteúdo econômico da proposição: o estimador de Areosa correlaciona-se em 0,979 com a função de produção simples e difere dela em até 1 ponto percentual — exatamente a contribuição do componente HP do PIB com peso 0,342 —, posicionando-se, em variabilidade, entre os dois estimadores que o compõem.

## Implicação para o conjunto de métodos do BCB

O boxe do BCB (2024) apresenta MII.I e MII.II como métodos distintos e exibe trajetórias visivelmente diferentes para os dois (Gráfico 2 do boxe). A Proposição 1 mostra que, sob $\lambda$ comum e com os mesmos insumos, os dois métodos não podem divergir além da combinação (4) — cuja distância ao MII.I é limitada pelo termo $0{,}342\,(\hat{h}^{HP} - \hat{h}^{FP})$. Conclui-se que a divergência adicional observada nos resultados oficiais é informativa sobre os **dados**, não sobre os **métodos**: ela reflete, sobretudo, o uso interno da série de desemprego retropolada (amostra desde 2003) em lugar da PNAD Contínua pública (desde 2012). Para o argumento central deste trabalho, o episódio ilustra uma lição metodológica mais geral: parte da dispersão entre estimativas de hiato usualmente atribuída a diferenças de método pode, na realidade, originar-se de diferenças de insumo — uma distinção que apenas a replicação independente, como a aqui realizada, permite estabelecer.
