# Dívida técnica — o que precisa ser feito direito

Escrito em 2026-08-14, ao fim da sessão que executou as Fases 0-3 do
`ROADMAP.md`.

**Critério deste documento**: não é uma lista de melhorias desejáveis. É a
lista do que, se não for feito, **força remendo depois** — cada item está aqui
porque contorná-lo cria trabalho maior do que resolvê-lo.

Itens marcados **[criada nesta sessão]** são dívida que eu introduzi, não
herdada. Estão separados de propósito: são as que menos merecem indulgência.

---

## 1. Bloqueiam correção, não só elegância

### 1.1 ~~Fluxo de saída inconsistente~~ — RESOLVIDO em 2026-08-14

`projectToDivergenceFree()` forma o fluxo de saída a partir de `velocityStar`
(antes da correção de pressão); `netBoundaryFlux()` e `maxFaceDivergence()`
leem de `velocity_` (depois). A projeção equilibra um fluxo e o diagnóstico
mede outro, diferindo por `−dt·∇p·A` nas células de saída.

**Evidência**: desbalanço de 13,2% da vazão, idêntico até a quinta casa em
t = 3, 8 e 20 — constante, não transiente.

**Por que não dá para contornar**: qualquer escoamento externo (cilindro,
aerofólio, duto) tem o mesmo erro sistemático. Afrouxar a tolerância do teste
é exatamente o paliativo a evitar — esconde um erro que se propaga para todo
caso com entrada e saída.

**Resolvido**: a projeção já estava correta — `a_b` na diagonal e `a_b·p_saída`
no lado direito *é* a correção de pressão no contorno. O erro estava só no
diagnóstico, que recomputava o fluxo pelo gradiente de mínimos quadrados em
vez do gradiente compacto que a projeção usou. Passou a registrar o fluxo que
a projeção de fato impôs. **Desbalanço: 1,3e-01 → 6,3e-14.**

Vale registrar o padrão, porque é o segundo caso idêntico no projeto: na Fase 1
a "divergência ~0,2" da cavidade estruturada também era propriedade do
diagnóstico, não do escoamento. **Medir o operador que de fato se resolve.**

### 1.2 ~~Poisson da pressão sem correção não-ortogonal no solver de NS~~ — RESOLVIDO em 2026-08-16

`UnstructuredDiffusionSolver` tinha correção não-ortogonal por correção
diferida; `UnstructuredCavitySolver3D::applyPoissonOperator()` **não tinha**.
Eram o mesmo operador matemático, com tratamentos diferentes.

**Por que importava**: a Fase 2.2 mediu que, sem essa correção, o erro
**estagna** (ordem 0,10 em vez de convergir). O solver de NS estava com a
versão que sabidamente não converge, na equação mais importante que resolve.

**Resolvido**: não existe mais um segundo Laplaciano. `applyPoissonOperator()`
foi apagado; a projeção chama `solveDeferredCorrection()` da base — a mesma
função que a difusão usa, com a mesma correção — passando a célula fixada
quando não há saída. Foi a mudança de uma linha que o item 2.1 prometia.

**Medido — divergência máxima por faces na cavidade, contra o número de
correctores por passo** (topologia idêntica em todos os casos):

| correctores | divergência | |
|---|---|---|
| 0 (antes) | 1,150e-01 | sem correção nenhuma |
| 1 | **diverge** (2,9e+88) | a correção nunca alcança a pressão que corrige |
| 2 | 2,284e-02 | mínimo estável |
| 3 | 1,322e-02 | |
| **4** | **7,129e-03** | escolhido: +0,8 s numa suíte de 25 s |
| 8 | 8,560e-04 | |
| 16 | 1,286e-04 | |

Duas coisas se leem nessa tabela. O resíduo cai pela metade a cada corrector —
geométrico, que é o que diz que a correção diferida está convergindo e não
brigando com alguma coisa —, então **esse diagnóstico reporta quantos
correctores foram pagos, não uma propriedade do esquema**. E um corrector não
é apenas impreciso, é instável, que é por que a contagem não é simplesmente
minimizada. Escolhido 4: 16× melhor que a versão sem correção, por 3% de
tempo de suíte.

**Consequência que precisou ser perseguida**: com a correção ligada, o balanço
de massa do canal saltou de 6e-14 para 2,5e-03. Não era regressão da física —
a projeção passou a impor um fluxo com um termo a mais, e `boundaryFlux_`
ainda registrava só a parte compacta `a_b(p_saída − p_P)`. Passou a registrar
`(∇p)_f·A` inteiro, que é a mesma quantidade pela identidade
`A = a_b·d + A_nonorth`. Voltou a 2,1e-13. **Terceira vez que o mesmo erro
aparece neste projeto** (Fase 1, item 1.1, agora aqui): medir o operador que
de fato se resolve.

### 1.3 ~~Gradiente zero silencioso em célula de estêncil deficiente~~ — RESOLVIDO em 2026-08-16

`UnstructuredCavitySolver3D::scalarGradients()` devolvia gradiente **zero**
quando a matriz de mínimos quadrados era singular. `UnstructuredDiffusionSolver`
caía para Green-Gauss no mesmo caso.

**Por que não dava para contornar**: gradiente zero não é "sem resposta", é
uma resposta errada e plausível. Numa malha de geometria real, células de
estêncil deficiente aparecem justamente perto de contornos complicados — onde
o gradiente importa mais.

**Resolvido**: `computeCellGradients()` na base faz mínimos quadrados onde o
estêncil suporta e Green-Gauss onde não suporta. **A escolha não é mais
oferecida a quem chama** — não existe API que devolva zero.

**Medido, e a medição é o achado**: `deficientStencilCount()` (novo
diagnóstico de qualidade de malha, ao lado de `maxNonOrthogonality()`) mostra
que isto **não era hipotético**:

| caso | células | estêncil deficiente |
|---|---|---|
| placa (difusão), n=4/6/8 | 415 / 1358 / 3184 | **0** |
| cavidade não-estruturada | 415 | **29 (7,0%)** |
| canal com saída | 177 | **21 (11,9%)** |

A difusão nunca caiu no recuo porque inclui os valores de Dirichlet do
contorno no ajuste, o que completa o posto. O solver de NS não inclui (parede
é gradiente nulo, não carrega informação), então suas células de canto ficam
com poucos vizinhos — e **7 a 12% das células vinham recebendo gradiente de
pressão zero em todo teste que já rodou**.

Efeito isolado, desligando só o recuo: u médio no topo da cavidade
0,06782 → 0,06929 (+2,2%); a divergência praticamente não se mexe
(7,158e-03 → 7,129e-03). O item mexe no campo, não no diagnóstico da
projeção — como esperado, já que a projeção zera o operador que resolve
independentemente de o gradiente das células deficientes estar certo.

### 1.4 ~~Face de saída arrasta viscosamente no operador de Helmholtz~~ — RESOLVIDO em 2026-08-16

`applyHelmholtzOperator()` montava a diagonal como
`V/dt + ν·laplacianDiagonal_[célula]` e depois pulava as faces de saída no
laço de contorno, com o comentário "zero-gradient: no viscous flux through an
outlet". Só que `laplacianDiagonal_` **já continha** o `a_b` da saída — somado
ali de propósito, porque é o que fixa o nível de pressão no Poisson.
Resultado: a célula de saída ganhava `ν·a_b` na diagonal **sem termo
correspondente no lado direito**, o que impõe `u = 0` na face de saída em vez
de gradiente nulo. O contrário do que o comentário afirmava.

**Resolvido**: a base mantém duas diagonais. `interiorDiagonal_` só as faces
interiores, `laplacianDiagonal_` interiores mais as de contorno que entram no
operador. O Poisson usa a segunda, o Helmholtz usa a primeira mais as paredes.
Uma saída é Dirichlet para pressão e gradiente nulo para velocidade; compartilhar
uma diagonal entre as duas era assumir que os dois operadores concordam sobre
quais faces são conexões, e eles não concordam.

**Medido, e o resultado contraria o que eu previa**: u médio nas células
junto à saída do canal, 1,00367 com o termo espúrio, 0,98698 sem ele — 1,7%
de diferença, com o balanço de massa idêntico em precisão de máquina nos dois
casos (−8,0e-15 vs −2,1e-13). Vale registrar que o deslocamento é **para
baixo**: o termo não é um freio simples cuja remoção acelera o escoamento,
porque a projeção responde ao preditor alterado redistribuindo pressão. A
direção não é o argumento para a correção — não haver fluxo viscoso através de
uma face de gradiente nulo é —, mas é o motivo pelo qual isto precisou ser
medido em vez de previsto.

**O que este item ensina, além dele mesmo**: um diagnóstico global que fecha
não prova que o campo está certo. O balanço de massa fechou em 1e-13 durante
todo o tempo em que a saída estava com condição de parede no preditor, porque
a projeção impõe o balanço independentemente do que a equação de momento fez
lá. O teste do canal passou a reportar o perfil de saída por isso.

---

## 2. Duplicação que vai divergir

### 2.1 ~~Dois solvers não-estruturados com a mesma geometria copiada~~ — RESOLVIDO em 2026-08-16

`UnstructuredDiffusionSolver` e `UnstructuredCavitySolver3D` tinham cópias
independentes de: coeficientes de face (`|A|²/(A·d)`), decomposição
over-relaxed, estêncil de mínimos quadrados, inversa 3×3 simétrica com guarda
de posto, e o laço de CG.

**Por que não dava para contornar**: é exatamente a situação que motivou
extrair `StaggeredCavityBase3D` das seis cavidades 3D antes dos Módulos 9-14 —
e a razão registrada lá vale igual aqui: **portar uma correção seis vezes é a
versão cara do problema**. O item 1.2 acima já era a primeira divergência
entre as duas cópias, com menos de um dia de vida.

**Resolvido**: `UnstructuredFvmBase` (`engine/solver/.../UnstructuredFvmBase.hpp`)
passa a ser dona da geometria de face, do estêncil de gradientes, da inversa
simétrica com guarda de posto, do CG matrix-free e do
`maxNonOrthogonality()`. As duas classes herdam dela. Código (sem comentários
nem linhas em branco): 836 → 720 linhas, com 388 linhas saindo dos dois
solvers e 272 passando a existir **uma** vez em vez de duas.

**Onde as duas cópias divergiam de verdade, a diferença virou argumento
nomeado em vez de divergência silenciosa**: `buildFaceGeometry()` recebe o
predicado que decide se uma face de contorno entra no operador (Dirichlet na
difusão, saída no NS — a mesma regra dita para condições diferentes), e
`buildGradientStencils()` recebe o que informa o valor prescrito de uma face
de contorno (a difusão informa; o NS não). Quem lê o cabeçalho vê as duas
diferenças; antes elas estavam escondidas em dois corpos de código que
ninguém compara.

**A extração é preservadora de comportamento, e isso foi projetado, não
torcido para acontecer**: cada laço compartilhado percorre as faces na mesma
ordem das duas versões originais e acumula nas mesmas somas na mesma
sequência. Foi por isso que `buildFaceGeometry()` recebe um predicado em vez
de deixar cada solver somar suas contribuições de contorno depois — montar
"interiores primeiro, contornos depois" seria igualmente correto e teria
perturbado toda soma no último bit, o que destrói a única checagem barata de
que um refactor não mudou nada.

**Medido**: a placa vs. série de Fourier saiu idêntica em todos os dígitos
impressos (0.98949 / 0.69966 / 0.53084, ordens 0,85 e 0,96); a topologia da
cavidade idêntica (topo +0.06850, fundo −0.01169). O único número que se
mexeu foi o balanço do canal, −6,347e-14 → −6,425e-14: ruído de arredondamento,
da única mudança aritmética real da extração — o estêncil do NS reconstruía
`d` como `unitD·|d|` e agora usa o `c_N − c_P` exato, como a difusão sempre
fez. Suíte inteira: 11/11.

**O que isso destrava**: 1.2 e 1.3 agora são correção em um lugar só, que era
o motivo de este item vir primeiro na ordem.

### 2.2 ~~Duas representações de malha que não se falam~~ — RESOLVIDO em 2026-08-16

`core::Mesh` (vértices + células) existe desde o Módulo 1 e é o que
`ScalarField`/`VectorField` usam. `TetrahedralMesh` foi criada na Fase 2.1 e
não é construída sobre ela.

**Por que importa**: os campos do engine (`ScalarField`, `VectorField`) não
podem ser usados com os solvers não-estruturados, então esses solvers
carregam `std::vector<double>` cru. Isso vai forçar conversões manuais em todo
ponto de integração — pós-processamento, persistência, visualização.

**O que fazer**: decidir qual é a representação canônica e fazer a outra
construir sobre ela. Não manter as duas.


**Resolvido: `TetrahedralMesh` passou a guardar um `core::Mesh` em vez de uma
segunda cópia dos vértices.** `core::Mesh` é a canônica porque é o que a
camada de campos já fala; o que a `TetrahedralMesh` acrescenta é
**conectividade** — quais células dividem uma face, e a geometria dessa face —
que é a parte de que o método de volumes finitos precisa e que uma lista de
vértices e células não expressa. Uma representação com uma camada em cima, em
vez de duas paralelas.

**As duas concordam exatamente, não aproximadamente**, e isso vale dizer
porque não é verdade em geral: `core::Mesh::cellCentroid()` é a média dos
vértices da célula, e para um tetraedro isso *é* o centroide. Para qualquer
outra forma de célula seria aproximação e a identidade precisaria ser
reexaminada. O teste compara os dois com `== 0.0`, não com tolerância.

**O dano que o item descrevia era "conversões manuais em todo ponto de
integração"**. `coreMesh()` devolve a malha sobre a qual `ScalarField` e
`VectorField` são definidos, célula por célula e na mesma ordem, então não há
conversão: o teste constrói os dois campos sobre ela, preenche a partir de um
solver que rodou 200 passos e confirma que carregam os números do solver, não
uma reamostragem deles.

### 2.3 ~~Pressão da saída fica de fora do ajuste de mínimos quadrados~~ — RESOLVIDO em 2026-08-16

`buildGradientStencils(useBoundaryValues)` era chamado com `false` no solver
de NS, então **nenhuma** face de contorno entrava no estêncil de gradiente.
Para parede isso é certo e o motivo está no código: gradiente nulo afirma que
a *derivada normal* some, e o ajuste leria isso como "o valor na face é igual
ao da célula", enviesando também as componentes tangenciais. Mas a saída tem
pressão prescrita — é Dirichlet, exatamente como as faces quentes e frias da
placa, que a difusão sempre incluiu.

**Resolvido**: uma linha (`true` em vez de `false`). O callback registrado por
`setBoundaryFaceValue()` já responde "esta face tem valor?", e só a saída diz
que sim — parede continua fora, pelo motivo certo.

**Medido, e o efeito é muito maior do que "completar alguns estêncils".**
Canal, medido a correctores suficientes para isolar o efeito do atraso da
correção:

| malha | estêncil deficiente | u médio na saída |
|---|---|---|
| C++ n=3 (suíte) | 21 → **14** | 0,98698 → **0,99916** |
| Python n=3 jitter 0,25 | 22 → **15** | 0,91759 → **0,96881** |
| Python n=4 jitter 0,25 | 30 → **19** | 1,02152 → **0,99941** |

As três malhas movem o `u` da saída **em direção a 1,0**, que é a velocidade
de bulk que a conservação de massa impõe num canal — e nas duas mais finas
chega a 0,999. Nenhuma delas foi ajustada para isso; é o gradiente de pressão
na saída deixando de ser inventado.

**E dissolveu o item 5.4.** O canal precisava de 64 correctores não-ortogonais
para fechar o balanço de massa em 1e-13; com a saída no ajuste ele fecha em
1e-13 com **dois**, em toda malha testada até não-ortogonalidade 2,24:

| correctores | antes (naoOrtog 1,651) | depois |
|---|---|---|
| 2 | — | +3,1e-13 |
| 4 | +2,7e-04 | +3,6e-13 |
| 64 | −7,1e-14 | +9,9e-14 |

A correção diferida não estava convergindo devagar por causa da malha: estava
sendo construída a partir de um gradiente ruim exatamente onde importava. **A
lição é a de sempre neste projeto, num terceiro disfarce** — o número que
parecia medir a malha estava medindo um erro nosso. Item 1.3 (gradiente zero
silencioso), 2.3 (saída fora do ajuste) e o sintoma registrado como 5.4 eram
a mesma causa vista de três ângulos.

**O que continua valendo do 5.4**: a cavidade *fechada* — que não tem face de
Dirichlet nenhuma e fixa uma célula de referência no lugar — ainda mostra a
queda geométrica com o número de correctores (2,28e-02 com 2, 7,13e-03 com 4,
1,29e-04 com 16). É um problema genuinamente mais difícil para essa iteração,
e é a razão do default ser 4 e não 2.

---

## 3. Acurácia declarada mas não atingida

### 3.1 ~~Convecção upwind de primeira ordem~~ — RESOLVIDO em 2026-08-16

Escolhida deliberadamente por estabilidade (diferença central é
incondicionalmente instável acima de Re de célula 2, e a Fase 1 mediu a
cavidade estruturada em 16,7). Mas é primeira ordem, e **difusão numérica de
upwind cresce com o refinamento da malha em direções oblíquas ao escoamento**
— justamente o caso em tetraedros.

**Primeiro a régua, depois o esquema** — na ordem, porque construir o esquema
antes seria não saber o que ele comprou. `UnstructuredScalarTransportSolver`
resolve `div(U φ) = Γ lap(φ) + S` com **U prescrito**, que é o que torna a
medição possível: no solver de Navier-Stokes a velocidade é a incógnita, então
uma solução manufaturada lá mede a equação de momento, a projeção e o esquema
de convecção ao mesmo tempo.

Caso: mesma `φ = sin(πx) cos(πy) exp(z)` do item 3.2, carregada por
`U = (1,1,1)/√3` — uniforme, logo de divergência discreta exatamente nula, e
**oblíqua a toda face de tetraedro**, que é onde a difusão do upwind é pior.
Γ = 0,01 põe o Péclet de célula entre 12,5 e 4,4.

| n | upwind 1ª ordem | ordem | linear-upwind limitado | ordem |
|---|---|---|---|---|
| 4 | 2,378150e-01 | — | 1,313639e-01 | — |
| 6 | 1,567389e-01 | **1,028** | 5,683272e-02 | **2,066** |
| 8 | 1,115700e-01 | **1,182** | 3,361285e-02 | **1,826** |
| 10 | 8,621935e-02 | 1,155 | 2,360342e-02 | 1,584 |

**A premissa do item estava certa: upwind é ordem 1.** E o esquema limitado é
segunda ordem na região suave, com **3,7× menos erro na malha mais fina**, sem
custo de malha nenhum.

**A queda de 2,07 para 1,58 é esperada e é o preço da limitação**, não um
defeito escondido: um limitador TVD corta em extremos suaves, e esta solução
tem extremos interiores em x e y — aquelas células caem para primeira ordem e
tomam uma fatia crescente da norma conforme a malha refina. Registrado porque
seria fácil (e errado) apresentar "segunda ordem" sem a ressalva.

**O esquema**: reconstrução linear a partir da célula a montante, escrita como
`φ_f = φ_C + ψ (φ_central − φ_C)`, de modo que ψ = 0 é upwind exato e ψ = 1 é
central exato. `φ_central` é a interpolação ponderada por distância que o
resto da base já usa, não média simples — em malha torta as duas diferem em
primeira ordem. Razão do limitador `r = 2(∇φ_C · d)/(φ_D − φ_C) − 1` e
limitador de van Leer, escolhido liso de propósito: um limitador com quina
(minmod, superbee) torna o resíduo da iteração estacionária não-diferenciável
e pode travá-la antes de convergir. Vive na base, então o solver de NS e o de
transporte usam o mesmo código — disciplina do item 2.1.

**Aplicado ao momento do Navier-Stokes**, por componente (o limitador se
constrói do gradiente da grandeza convectada, e momento tem três). Efeito
medido: u médio na saída do canal 0,99916 → **1,00020** contra a velocidade de
bulk 1,0 que a conservação de massa impõe; topo da cavidade +0,06929 →
+0,06973. Pequenos e na direção física esperada, como se espera de uma malha
grosseira onde a topologia já estava certa.

**Efeito colateral necessário**: o march explícito ao estado estacionário era
inviável — 150 mil passos e 246 s para n = 8, porque o passo é limitado pela
convecção enquanto a relaxação acontece na escala difusiva. Passo **local por
célula** resolveu: 150.235 → 780 passos, com o estado estacionário convergido
idêntico. O operador continua simétrico (V/dt entra só na diagonal), então o
CG da base segue valendo — que é o motivo de essa aceleração estar disponível
aqui. O caso de verificação inteiro custa 1,2 s na suíte.

**O que este item não afirma**: que o solver de NS é de segunda ordem em
convecção. Mede o esquema de face isoladamente, com U prescrito. A ordem do NS
completo depende também do acoplamento pressão-velocidade, cuja
incompatibilidade entre operadores é a suspeita aberta do item 4.3.

### 3.2 ~~FVM não-estruturado em ordem ~1, não 2~~ — RESOLVIDO em 2026-08-16: é ordem 2

Medido na Fase 2.2: 0,91 → 0,98 na norma global; 1,06 → 1,63 excluindo os
cantos singulares. Segunda ordem existia na região suave, mas não tinha sido
demonstrada num caso sem singularidade — era inferência, e estava marcada
como tal.

**Resolvido, e a inferência estava certa.** Solução manufaturada suave em todo
o cubo fechado:

    phi(x,y,z) = sin(pi x) cos(pi y) exp(z)
    lap(phi) = (1 - 2 pi^2) phi   =>   S = (2 pi^2 - 1) phi

com `phi` imposto nas seis faces. Tudo o que sobra entre o campo calculado e
`phi` é erro de discretização — nada recordado de tabela, nada ajustado.

| n | células | rms | ordem | naoOrtog |
|---|---|---|---|---|
| 4 | 415 | 2,338802e-02 | — | 1,48 |
| 6 | 1358 | 1,033868e-02 | **2,013** | 1,45 |
| 8 | 3184 | 5,775569e-03 | **2,024** | 1,63 |
| 10 | 6204 | 3,675188e-03 | **2,026** | 1,75 |
| 12 | 10667 | 2,577107e-03 | 1,947 | 1,64 |

**Segunda ordem, em malha cuja não-ortogonalidade fica entre 1,45 e 1,75 e não
melhora com refino.** O ~1 da placa era o teto da singularidade de canto, não
o esquema. As duas primeiras linhas estão na suíte; n = 10 e 12 rodam fora
dela, porque tetraedralizar n = 10 custa 36 s contra 0,16 s de solve.

**Três coisas verificadas antes de acreditar no número**, cada uma capaz de
limitar a ordem medida exatamente no valor sob teste:

- A fonte é integrada pela regra do ponto médio, `S(centroide)·V` — exata para
  fonte linear e O(h²) caso contrário. Fica na mesma ordem sob teste, não
  abaixo dela, então não pode se disfarçar de esquema de primeira ordem.
- O erro é medido contra `phi` no centroide, enquanto a incógnita de volumes
  finitos é a média da célula. Diferem por O(h²), o que significa que **esta
  norma não poderia demonstrar ordem acima de 2** — mas distingue 1 de 2, que
  é a pergunta.
- A correção diferida foi levada à convergência, não ao teto de varreduras: em
  n = 8 o erro é 5,77434e-03 com 20 varreduras, 5,775568e-03 com 60 e
  5,7755687e-03 convergido. A iteração externa parou de importar quatro
  dígitos antes da discretização.

**O que isto exigiu de API**, e ambas são funcionalidade legítima que faltava:
`setSourceTerm()` (Laplace vira Poisson — sem fonte, as únicas soluções
exatas são as harmônicas, família estreita que deixa a integração de volume
sem teste) e uma sobrecarga de `setDirichletBoundary()` com valor variável
pela face, porque a condição de contorno de uma solução manufaturada *é* a
solução exata amostrada na superfície. Ambas com binding Python.

**Efeito colateral que valeu o passo**: as malhas de teste passaram a ser
construídas uma vez e compartilhadas (`cubeLatticeMesh`). A tetraedralização
domina — 0,24 s / 1,99 s / 10,26 s / 36,2 s para n = 4/6/8/10, contra
0,00 s / 0,02 s / 0,06 s / 0,16 s de solve — e quatro testes queriam as
mesmas resoluções. A suíte inteira **caiu de 29,3 s para 28,7 s mesmo com o
caso novo dentro**.

**O que este item não afirma**: que o *solver de Navier-Stokes* é de segunda
ordem. Isto mede o Laplaciano com a correção não-ortogonal e os gradientes de
mínimos quadrados — a parte compartilhada. A convecção continua upwind de
primeira ordem (item 3.1), e agora existe uma régua confiável para medir o que
um esquema de alta ordem melhoraria.

### 3.3 ~~Interpolação de face sem correção de skewness~~ — MEDIDO em 2026-08-16; a conclusão nula não generaliza, e a atribuição não fecha

Implementada e medida como **nula** na malha de rede com jitter. Mas essa
malha é benigna: gerada de uma rede regular perturbada. Malha de geometria
real tem skewness muito pior.

**O que fazer**: repetir a medição numa malha de geometria real antes de
concluir que não importa. O resultado nulo é válido só para o caso testado.

---


**Medido, com a régua que o item 3.2 construiu** — a solução manufaturada mede
ordem sem precisar implementar nada. Malhas graduadas geometricamente, que é o
que produz skewness alta, contra a rede uniforme com jitter em que o resultado
nulo foi obtido:

| malha | skew máx | naoOrtog | ordem | convergiu? |
|---|---|---|---|---|
| uniforme + jitter (a original) | 1,15 | 1,63 | **2,02** | sim, 99 varreduras |
| graduada 6× | 3,17 | 4,90 | **1,77** | sim, **282** varreduras |
| graduada 20× | 10,05 | 15,08 | — | **não: diverge** |

**A conclusão nula não generaliza**: em skewness 3,17 a ordem cai de 2,02 para
1,77, num solve genuinamente convergido (o teto de 200 varreduras não bastava;
com 2000 ele fecha em 282 com o mesmo rms). O item estava certo em dizer que o
nulo valia só para o caso testado.

**Mas a atribuição não fecha, e isso precisa ser dito.** As malhas graduadas
têm skewness *e* não-ortogonalidade altas juntas (1,63 → 4,90 → 15,08), então
esta medição não separa as duas causas. Concluir "é a skewness" seria
exatamente o erro que os itens 1.1, 1.2 e 4.3 já custaram. Separar exige uma
malha de skewness alta e não-ortogonalidade baixa, que não é trivial de
construir.

---

**Atualização 2026-08-25: metade da atribuição fecha — não-ortogonalidade
sozinha NÃO custa ordem.**

A malha que faltava é obtida por **cisalhamento afim**, não por gradação.
Um cisalhamento leva toda célula numa paralelepípeda congruente: a linha
centroide-a-centroide continua cruzando cada face no centro (skewness
fica baixa) enquanto a normal da face se inclina em relação a essa linha
(não-ortogonalidade sobe). É exatamente a separação que faltava.

Medido com a mesma régua do item 3.2 (solução manufaturada
`sin(πx)sin(πy)sin(πz)`, fonte exata, n = 4 → 6), em
`python/research/skewness_vs_nonorthogonality.py`:

| cisalhamento | skew (n=4 / n=6) | naoOrtog (n=4 / n=6) | ordem |
|---|---|---|---|
| 0,0 | 0,44 / 0,41 | 1,54 / 1,59 | **2,72** |
| 0,6 | 0,48 / 0,51 | 2,15 / 3,29 | **2,06** |
| 1,2 | 1,02 / 0,99 | 3,95 / **9,43** | **2,15** |

**Com a skewness segurada em torno de 1,0, a não-ortogonalidade sobe de
1,54 a 9,43 e a ordem não se move de ~2,1.** A correção diferida dá conta
dela — o que é consistente com o próprio achado deste item de que a
correção não-ortogonal converge até não-ortogonalidade 563 com o
amortecimento certo. As varreduras crescem (74 → 179 → 546), ou seja, ela
trabalha mais, mas entrega a mesma ordem.

**A janela de isolação tem teto, e vale registrar.** Em cisalhamento 1,8 e
2,4 a tetraedralização deixa de preservar a estrutura afim e a skewness
sobe junto (2,35 e 2,20), com a ordem caindo para 0,81 e 1,02 — de volta
ao regime em que as duas andam juntas e nada se separa. A tentativa
complementar (gradação axial, que deveria subir skewness segurando
não-ortogonalidade) **não isolou**: a não-ortogonalidade subiu a 16,33
junto. Está no script, registrada como o que é.

**O que isto conclui e o que não conclui.** Conclui, com medição: a
não-ortogonalidade **não é** a causa da perda de ordem, na faixa onde as
duas podem ser separadas. Não conclui que a causa seja a skewness — isso
seria eliminação por sobra, que é o erro que este item se propôs a evitar.
Uma malha de skewness alta com não-ortogonalidade baixa continua sendo o
experimento que falta, e a gradação axial não a produz.

**O achado que veio de graça, e é maior que o item.** A primeira leitura desta
medição deu "ordem 0,62" para a malha graduada 20× — e era **lixo**: aquele
solve estava divergindo, e o número saiu de onde o teto de varreduras por acaso
o interrompeu. Com teto 2000 a mudança por varredura chega a 7,14e+07 e o
campo devolvido a 1,2e+08. **A correção diferida é um ponto fixo que não é
sempre contração**, e `solveDeferredCorrection()` devolvia o resultado
divergente como se fosse solução.

Passou a recusar: se a mudança por varredura cresce mais de três ordens de
grandeza acima do melhor ponto já alcançado, levanta. Comparado contra o
*melhor* e não contra o primeiro, porque a correção legitimamente piora por
uma ou duas varreduras antes de assentar. É a mesma família do 4.3 — mas aqui
sem nenhuma célula deficiente, então não é o recuo Green-Gauss: é a correção
não-ortogonal em si, a não-ortogonalidade 15.

**O que fazer**: para atribuir, uma malha que isole skewness de
não-ortogonalidade. Para usar malha de geometria real, o limite prático agora
é a divergência da correção acima de não-ortogonalidade ~15, não a
interpolação de face.


**Atualização 2026-08-16: o bloqueio prático caiu.** A correção diferida
ganhou sub-relaxação adaptativa — o fator multiplica a correção, portanto
multiplica também o fator de contração do ponto fixo, e é reduzido pela metade
sempre que a mudança por varredura cresce dez vezes acima do seu melhor valor.
Adaptativo e não derivado da malha porque **três métricas estáticas já
falharam** em separar uma malha que roda de uma que não roda
(não-ortogonalidade, a razão |A_nonorth|/|A_orth| e a contagem de estêncils
deficientes).

| malha graduada | não-ortogonalidade | antes | depois |
|---|---|---|---|
| 6× | 4,9 | 282 varreduras | 282 varreduras |
| 20× | 15,1 | **divergia** | 344 varreduras, rms 8,61e-02 |
| 50× | **563,1** | divergia | 175 varreduras, rms 1,13e-01 |

O limite prático que este item registrava — "a divergência da correção acima
de não-ortogonalidade ~15" — subiu para além de 563. A atribuição
skewness-versus-não-ortogonalidade continua aberta.

---

## 4. Estabilidade sobre margem nula

### 4.1 ~~Solvers estruturados em CFL exatamente 1,0000~~ — RESOLVIDO em 2026-08-16

`LidDrivenCavitySolver2D::stableTimeStep()` devolvia o limite convectivo sem
fator de segurança, e a cavidade a Re=400 roda com Re de célula 16,7 — muito
acima do limite 2 da diferença central.

**Evidência original**: durante a Fase 1, uma perturbação modesta no esquema
levou esse caso a NaN em ~500 passos. Ele estava à beira o tempo todo.

**O que se encontrou ao olhar de perto foi pior que o registrado.** O fator de
segurança existia — morava em **49 cópias no chamador**. Testes, o app do
viewer e os primers dos solvers turbulentos escreviam cada um seu
`0.3 * solver.stableTimeStep()`. Uma função chamada `stableTimeStep()` que
devolve um passo instável é uma armadilha, e ela **já tinha sido acionada**:
`engine/analysis` e `engine/persistence` chamam
`solver.step(solver.stableTimeStep())` sem fator nenhum, tendo lido o nome e
acreditado. Essas duas suítes rodavam em CFL 1,0000.

**E uma das sete cópias da fórmula estava simplesmente errada.**
`TransientDiffusionSolver::stableTimeStep()` devolvia `1/Σ(1/h²)`, sem o fator
2 de `1/(2ν Σ1/h²)` — ou seja, **o dobro do limite de estabilidade**, com um
comentário mandando o chamador dividir por dois. Medido em vez de argumentado:
perfil senoidal em 41 células, passo exatamente no valor devolvido, **NaN no
passo 654**.

É a mesma "duplicação que vai divergir" do item 2.1, do lado estruturado, e
com a divergência já materializada.

**Resolvido**: `explicitStableTimeStep()` — uma fórmula, o fator dentro dela,
os sete solvers chamando-a, e as 49 margens privadas dos chamadores removidas.

**Preservador de comportamento onde já havia margem, e isso é verificável**:
0,3 é o valor em que 49 chamadores independentes tinham convergido, então
movê-lo para dentro deixa esses casos idênticos. A saída da suíte de solvers
saiu **bit a bit igual** antes e depois. O que mudou foi exatamente o que
devia: os chamadores que confiavam no nome passaram a receber margem, e o
solver de difusão transiente passou a devolver um passo que é de fato estável.

**Margem medida, não declarada**: o mesmo perfil senoidal agora sobrevive a
**3,33×** o valor devolvido e diverge em **3,40×**. 3,33 = 1/0,3 exatamente —
o valor devolvido está no fator de segurança que a constante declara, e o
ponto de falha está onde a teoria diz.

**O que este item não resolve, e o texto original já dava a escolha**: "fator
de segurança **e/ou** upwind na convecção 2D". O fator entrou; a diferença
central a Re de célula 16,7 continua. Essa metade é propriedade do *esquema*,
não do passo, e passo menor não cura — oscilação de diferença central acima de
Re 2 não é instabilidade temporal. O item 3.1 agora produziu um esquema
limitado medido em segunda ordem para o lado não-estruturado; portá-lo para os
solvers estruturados é o trabalho que fecha essa metade, e fica registrado
como **4.4**.

### 4.2 ~~Convecção ainda explícita no solver não-estruturado~~ — RESOLVIDO em 2026-08-16

A difusão virou implícita e destravou o passo. A convecção não — e passou a
ser ela o limite, o que o item 5.3 mediu: dt caindo de 7,3e-04 a 7,7e-05 entre
n=3 e n=6, e o estudo de convergência custando 309 s.

**Resolvido pela metade que importa: a parte de saída da convecção foi para o
implícito.** Não é o implícito completo, e a escolha tem uma razão precisa. Um
operador upwind totalmente implícito é **não-simétrico** — os dois
coeficientes de vizinho diferem conforme a direção do escoamento — e exigiria
BiCGSTAB ou GMRES. Mover só o fluxo de *saída* torna a adição puramente
diagonal, então o operador continua simétrico positivo-definido e nada no
solve linear muda.

E é incondicionalmente limitado, que é o ponto. Escrevendo a parte de primeira
ordem da atualização,

    u* = [ (V/dt) u^n + Σ_entrada F u_N^n ] / [ V/dt + Σ_saída F ]

e usando Σ_entrada F = Σ_saída F para um fluxo de face de divergência nula,
`u*` é **combinação convexa** de `u^n` e dos vizinhos, para qualquer dt.
Nenhum passo consegue fazê-la ultrapassar. A parte do limitador fica
explícita, na mesma estrutura de correção diferida do termo não-ortogonal.

**Medido na cavidade n=4, marchando a t=8 em múltiplos do passo devolvido:**

| fator | passos | u topo | divergência | estável |
|---|---|---|---|---|
| 1× | 16356 | +0,069704 | 1,3e-11 | sim |
| 10× | 1635 | +0,072180 | 1,8e-10 | sim |
| 25× | 654 | +0,073303 | 2,9e-10 | sim |
| 100× | 163 | +0,073137 | 1,4e-09 | sim |

**Estável a cem vezes o passo devolvido**, com a resposta deslocando 5% — que
é *acurácia*, não estabilidade, e está bem dentro dos 44% que a própria malha
custa nessa resolução (item 5.3). O limite convectivo deixou de existir.

`stableTimeStep()` continua devolvendo o valor conservador de antes, e o
cabeçalho passou a dizer o que ele significa agora: **menos** que o limite, não
mais — o nome segue honesto, no sentido exato que o item 4.1 cobra. O que
mudou é que o número virou sugestão de acurácia em vez de barreira.

**Efeito colateral que não era o alvo**: o caso do item 4.3 (canal em malha
jitter ±0,45) ia a NaN no passo 20; com a convecção semi-implícita ele vai ao
passo **292**. Não sumiu — em jitter ±0,65 o raio espectral em torno do
repouso é 44,05 e o caso levanta no passo 25 —, mas o guarda que o 4.3
instalou passou a disparar tão mais tarde que o teste Python que o exercitava
**parou de falhar por conta própria**, e precisou de uma malha pior. Vale
registrar como aconteceu: o teste quebrou porque o solver melhorou, não porque
o guarda quebrou.


### 4.3 ~~NaN em malha muito distorcida~~ — RESOLVIDO em 2026-08-25, para domínio fechado e com saída, na sétima tentativa

Reproduz em segundos pelos bindings: rede cúbica n=3 com jitter de ±0,45/n
(contra os ±0,25/n dos testes), canal com entrada e saída. O campo vai a
**NaN no passo 20 de 1333**.

**O que foi medido, e o que cada medição eliminou.** A instrução deste item
era instrumentar antes de propor qualquer correção; foi o que se fez, e a
resposta fácil estava errada em todas as tentativas:

| hipótese | teste | resultado |
|---|---|---|
| NaN súbito (divisão por zero) | traçar \|u\|max por passo | **não** — cresce exponencialmente desde o passo 1 |
| passo de tempo (CFL) | dt × 0,25 e × 0,1 | **não** — morre no *mesmo passo* (20→21), não no mesmo tempo |
| face descartada (`A·d ≤ 0`) | contar em Python | **não** — zero em toda malha; cos(A,d) mínimo é 0,436 em jitter 0,45 contra 0,464 em jitter 0,25 |
| estêncil deficiente | `deficientStencilCount()` | **não** — 15 nas duas malhas |
| a saída | rodar cavidade fechada na mesma malha | **não** — morre também (passo 31) |
| convecção / Reynolds | entrada 0,1; ν = 1,0 (Re=1) | **não** — morre nos dois (passo 29) |
| correctores de pressão | 2, 4 e 16 | **não** — morre nos três |

**O que sobrou, e é um resultado positivo, não uma eliminação.** Rodando a
cavidade fechada com tampa a 10⁻⁶ e o *mesmo* dt do caso a 1,0 — amplitude
seis ordens de grandeza menor, onde a convecção (quadrática em u) é
irrelevante:

| malha | tampa 1,0 | tampa 10⁻⁶ |
|---|---|---|
| jitter 0,25 | razão 1,0293 por passo | razão 1,0294 |
| jitter 0,45 | NaN no passo 31 | **razão 2,0524 por passo** |

A razão de crescimento é a mesma em 10⁶ de amplitude. **A instabilidade é
linear, com fator ≈ 2,05 por passo, propriedade da malha e do operador — não
da física, não do passo, não da convecção.** Em jitter 0,45 ela sobrevive 40
passos com tampa 10⁻⁶ apenas porque começou seis ordens abaixo; é a mesma
instabilidade.

**A suspeita restante, registrada como suspeita e não como conclusão**: o
único elemento do passo que é linear, independente de dt e independente de
amplitude é a incompatibilidade entre os dois operadores da projeção — a
pressão é resolvida contra o operador *de face* (compacto mais correção
não-ortogonal), e a velocidade é corrigida com o gradiente *de célula* por
mínimos quadrados. A difusão é implícita (incondicionalmente estável) e a
convecção foi excluída pela medição de amplitude. Confirmar isso exige medir
o raio espectral do operador de um passo, o que precisa de uma API para impor
um campo de velocidade inicial — não existe hoje.

**Contido, não resolvido.** `step()` passa a recusar em vez de propagar: se o
campo perde finitude, levanta com uma mensagem que diz o que foi medido,
inclusive que reduzir o passo não adianta — o primeiro instinto de quem vê
uma explosão, e aqui só desperdiça tempo. Antes disso, todo número devolvido
depois era NaN, **inclusive os diagnósticos que um chamador consultaria para
perceber o problema**.

**Por que não há portão de qualidade de malha na construção**: porque não há
critério a-priori honesto para usar, e isso foi medido, não suposto. Uma malha
com não-ortogonalidade 2,24 (jitter 0,35) roda bem; esta, com 2,07, diverge.
Nenhuma quantidade que esta classe calcula hoje separa as duas. Inventar um
limiar seria fingir um critério.

**O que fazer**: expor um modo de carregar estado (velocidade e pressão), que
serve tanto para checkpoint quanto para medir o raio espectral do operador de
um passo diretamente — potência iterada em Python, poucos minutos. Só depois
propor correção; o candidato natural é tornar consistentes os dois operadores
da projeção, mas isso é reescrever o acoplamento pressão-velocidade e não se
faz sobre uma suspeita.

---


---

**Atualização: a causa foi encontrada, e ela é o item 1.3.**

O que faltava era impor um campo inicial arbitrário para medir o operador de
um passo. `loadState()` passou a existir — a mesma superfície que o
`StaggeredCavityBase3D` já tinha para checkpoint — e com ela o raio espectral
sai por potência iterada em torno do repouso, que é ponto fixo exato quando as
paredes estão imóveis.

| malha | raio espectral |
|---|---|
| jitter 0,25 (roda) | 0,9840 |
| jitter 0,35 (roda) | 0,9799 |
| jitter 0,45 (morre) | **2,0611** |

O 2,0611 bate com o 2,0524 medido antes por simulação não-linear a amplitude
10⁻⁶ — dois caminhos independentes, mesma resposta.

**A suspeita registrada estava errada.** Não era a incompatibilidade entre o
operador de face e o gradiente de célula. Três medições apontaram para outro
lugar:

- **Viscosidade não amortece**: a ν=100, com difusão implícita, ρ ainda é
  1,36. Não é momento nem difusão.
- **Mais correctores pioram**: 2 → 1,58, 4 → 2,06, 16 → **6,47**. A correção
  diferida não converge devagar ali; ela amplifica.
- **O modo é local**: 8 de 178 células concentram 96,3% da energia, e a
  dominante é a que tem **2 vizinhos interiores** e cai no recuo Green-Gauss.

E a atribuição, desligando só o recuo (o estado anterior ao item 1.3):

| malha | sem recuo | com recuo |
|---|---|---|
| jitter 0,25 | 0,9839 | 0,9840 |
| jitter 0,35 | 0,9794 | 0,9799 |
| jitter 0,45 | **0,9776** | **2,0611** |

**A correção do item 1.3 criou esta instabilidade.** O argumento do 1.3
continua de pé — gradiente zero era resposta errada, e removê-lo mexeu 2,2% no
campo da cavidade — mas Green-Gauss numa célula com um ou dois vizinhos
interiores não é estimativa, é chute: sua magnitude escala com |A|/V e nada a
limita. A energia do modo instável em células de recuo é 3–7% nas malhas que
rodam e **88,8%** na que morre.

**Corrigido pela metade, e a metade que falta está medida.** O recuo passou a
ser limitado pela maior inclinação que a célula de fato enxerga — nenhuma
diferença de vizinho dividida pela distância. Um gradiente maior que toda
diferença de que foi construído não é extrapolação, é invenção. Efeito:

- ρ em jitter 0,45: **2,0611 → 0,9785**, praticamente o valor sem recuo, então
  o limite devolveu a estabilidade sem voltar ao zero;
- cavidade fechada em jitter 0,45: de **NaN no passo 31** para 1333 passos
  completos;
- suíte: 12/12, com a cavidade mexendo na quinta casa (u topo 0,06973 →
  0,06972) e o canal na quarta (u saída 1,00020 → 1,00071) — o limite quase
  não morde nas malhas da suíte, que é o que se queria.

**O que ainda cai**: o canal com entrada e saída em jitter 0,45, agora no
passo 25 em vez de 20. O raio espectral que medi é em torno do **repouso**, e
o canal linearizado em torno do estado desenvolvido é outro operador. Medir
ali é o próximo passo, e `loadState()` já o permite: carregar um estado
desenvolvido, perturbar e iterar.


**Atualização 2026-08-16: melhorou muito, não fechou, e o que resta está
localizado.**

Duas mudanças desta rodada mexeram no item. A convecção semi-implícita (4.2)
levou o canal em jitter ±0,45 do passo 20 ao passo 292. E a correção diferida
ganhou sub-relaxação adaptativa, **persistente entre passos** — um projeção
que permite quatro varreduras por passo nunca aprende o amortecimento se o
fator reinicia a cada passo, então ele é carregado, e `pressureRelaxation()`
o expõe porque um valor abaixo de 1 é uma afirmação sobre a malha.

Isso resolveu o lado da *difusão* de forma decisiva (ver 3.3: converge agora
até não-ortogonalidade 563, contra 15 antes). **Não resolveu o solver de
Navier-Stokes em malha muito distorcida**, e a medição diz por quê: em jitter
±0,55 e acima o que dispara é o guarda de finitude do `step()`, não o da
correção — ou seja, a amplificação restante **não está na correção de
pressão**, que é a única coisa que a sub-relaxação toca.

| jitter | raio espectral (4 correctores) | canal até t=2 |
|---|---|---|
| 0,25 | 0,981 | ok |
| 0,45 | 0,974 | levanta no passo 292 |
| 0,55 | 5,19 | levanta no passo 13 |
| 0,65 | 44,05 | levanta no passo 25 |

**O que fazer**: o raio espectral já é mensurável em torno do repouso, e o
próximo corte é decompor o passo — medir o operador com a projeção desligada,
o que exige poder pular a projeção. Isso separaria o preditor de momento da
correção de velocidade por gradiente de célula, que é o único elemento do
passo ainda não isolado.


**A causa foi isolada, e é a suspeita original — que eu registrei ao abrir o
item, depois declarei errada, e que volta confirmada.**

O que faltava era rodar o passo *sem uma de suas peças*. `stepWith()` faz
isso: um passo com convecção, viscosidade ou projeção desligadas. Não é modo
físico — um passo sem projeção não conserva massa — é instrumento de medida, e
nada no solver o chama.

Raio espectral do operador de um passo, em torno do repouso:

| malha | completo | **sem projeção** | sem convecção | sem viscosidade | só projeção |
|---|---|---|---|---|---|
| jitter 0,25 | 0,9791 | 0,9839 | 0,9791 | 1,0007 | 1,0007 |
| jitter 0,45 | 0,9729 | 0,9803 | 0,9729 | 1,0000 | 1,0000 |
| jitter 0,55 | **5,1924** | **0,9880** | 5,1924 | 5,2363 | 5,2363 |
| jitter 0,65 | **2,4259** | **0,9918** | 2,4259 | 2,4356 | 2,4356 |

Três leituras, e juntas fecham:

- **Sem a projeção, ρ < 1 em toda malha.** O preditor de momento — convecção,
  difusão, contornos — é estável mesmo onde o passo completo explode.
- **Sem convecção, ρ é idêntico ao completo.** Em torno do repouso a convecção
  não contribui: já se supunha, agora está medido.
- **Só a projeção reproduz o número inteiro.** Ela carrega toda a amplificação.

E dentro da projeção, o solve de Poisson **já está sub-relaxado** e isso não
mudou o resultado. O que sobra é a **correção de velocidade**:
`velocity = velocityStar − ∇p·dt` usa o gradiente de *célula* por mínimos
quadrados, enquanto a pressão foi resolvida contra o operador *de face*. São
dois operadores diferentes, e é o único elemento que a decomposição não
consegue desligar separadamente — porque ele *é* a projeção.

Vale registrar o percurso: essa suspeita foi anotada ao abrir o item;
declarada errada quando a medição do recuo Green-Gauss explicou o jitter 0,45;
e volta confirmada para as malhas piores. **As duas coisas eram verdade ao
mesmo tempo** — o recuo ilimitado dominava numa faixa de malha, a
incompatibilidade da projeção domina na seguinte. Foi por isso que eliminá-la
cedo pareceu justificado e não era.

**O que fecha**: corrigir a velocidade com um gradiente consistente com o
operador que a projeção resolveu — reconstruir a correção de célula a partir
dos fluxos de face que foram zerados, em vez do gradiente de mínimos
quadrados. É reescrever o acoplamento pressão-velocidade, agora sobre causa
medida em vez de suspeita.


**Atualização: a correção que a isolação sugeria foi implementada e medida —
e não funciona.** Registro porque é achado, não silêncio: a hipótese era
sólida e a medição a derrubou, o que é exatamente o padrão que os itens 3.2 e
3.1 já mostraram (inferência plausível, decidida por medição em vez de
argumento).

**A hipótese.** A correção de velocidade usa `computeCellGradients(pressure_)`
— mínimos quadrados sobre diferenças célula-a-célula ao longo de `d` — para
corrigir `velocity_[cell] = velocityStar[cell] - ∇p[cell]·dt`. Mas o Poisson
foi resolvido contra o operador **de face**. Teoria de métodos de projeção diz
que a composição divergência→Poisson→gradiente só é um projetor não-expansivo
(ρ ≤ 1 garantido) quando divergência e gradiente são **adjuntos discretos**; a
divergência daqui é construída de fluxos de face e o gradiente de mínimos
quadrados não é o adjunto dela.

**A correção testada**: reconstruir a correção por célula a partir dos
próprios fluxos de face corrigidos, ao estilo `fvc::reconstruct` do OpenFOAM —
ajustar por mínimos quadrados o vetor `v` que melhor reproduz
`v·n_f ≈ deltaFlux_f/|A_f|` em cada face da célula, com `deltaFlux_f` sendo a
correção exata que a projeção aplicou àquele fluxo (0 numa parede, o termo de
gradiente consistente em face interior ou de saída). A matriz
`M = Σ_f n_f⊗n_f` depende só da geometria e foi extraída como estêncil
compartilhado; a extração do inversor 3×3 comum (`invertSymmetric3x3`) que a
viabilizou **foi mantida** — verificada bit-idêntica em toda a suíte e é
limpeza de código legítima independente do resultado abaixo.

**Medido — raio espectral, antes e depois:**

| malha | antes | depois |
|---|---|---|
| jitter 0,25 | 0,9791 | 0,9826 |
| jitter 0,45 | 0,9729 | 0,9751 |
| jitter 0,55 | **5,1924** | 4,4704 |
| jitter 0,65 | **2,4259** | **32,9884** |
| jitter 0,75 | — | 12,3937 |

**Não funciona.** Melhora em jitter 0,55, mas **piora em treze vezes** em
jitter 0,65. Verificado que não é bug óbvio nem estêncil mal condicionado: zero
células inválidas na malha 0,65, e as células que dominam o modo instável têm
número de condição do estêncil de fluxo entre 2 e 29 — tão bem-condicionado
quanto o estêncil de gradiente original. **A reconstrução por fluxo de face
não é universalmente melhor, é apenas um operador diferente**, com sua própria
estrutura de acoplamento entre células que pode amplificar em topologias que a
versão anterior não amplificava.

**Revertido**: a mudança na projeção foi descartada por inteiro. O código de
`UnstructuredCavitySolver3D` volta ao estado anterior; só a extração do
inversor 3×3 compartilhado permanece, verificada bit-idêntica.

**O que isso ensina para a próxima tentativa**: a hipótese da adjunção
continua matematicamente correta como explicação de *por que* ρ pode exceder
1, mas "usar fluxos de face em vez de gradiente de célula" não é a única forma
de restaurar a adjunção, e a forma testada introduz seu próprio modo de
amplificação. Duas direções ficam abertas e não testadas: (1) um peso
diferente no ajuste de mínimos quadrados do fluxo (o `1/|A_f|²` — variantes
existem na literatura), ou (2) sub-relaxar especificamente a correção de
velocidade, do mesmo jeito que a correção não-ortogonal do Poisson já é
sub-relaxada — mais barato de testar e mais próximo do que já funcionou
alhures neste item.


**Segunda tentativa, e ela é mais informativa que a primeira porque falhou de
um jeito diferente.**

A primeira tentativa (reconstrução por fluxo de face) trocava o operador de
correção por um adjunto ao operador de face — falhou ao piorar malhas já
ruins. A segunda tentativa manteve o operador original (gradiente de mínimos
quadrados) e tentou **iterá-lo com amortecimento**, exatamente o padrão que já
funciona para a correção não-ortogonal do Poisson: um laço de correctores
extras, cada um resolvendo uma correção de pressão contra o resíduo de
divergência da velocidade *já corrigida*, aceito só se o resíduo encolhe,
desfeito e re-tentado com metade da relaxação caso contrário — relaxação
persistente entre passos.

**Não funciona, e falha pior**: numa malha que **já convergia bem** (n=4 da
suíte, 415 células, antes chegando a divergência ~7e-03 em 16 mil passos), o
laço de correctores extras esgota o teto de 8 tentativas em **todo passo**,
sem nunca atingir a tolerância, e a divergência cresce monotonicamente —
4,6e-02 no passo 0, 1,1e-01 no passo 29, sem sinal de estabilizar. É pior que
não fazer nada.

**Por que isto é mais informativo que o primeiro fracasso**: o padrão
"amortecer até parar de crescer" funciona para a correção não-ortogonal
porque aquele mapa, perto do ponto fixo, tem um fator de contração
essencialmente escalar — sub-relaxar por um fator α reduz o fator de
crescimento pelo mesmo α, e "parou de crescer" é teste suficiente. **A
composição interpolação-para-face → Poisson → gradiente-de-célula não parece
ter essa estrutura**: iterá-la sem melhorar monotonicamente sugere autovalores
complexos ou um subespaço onde reduzir α não move na direção certa — halving
ingênuo baseado em "o resíduo cresceu, tente metade" não é a estratégia
correta para esse tipo de operador. Confirma, por outro ângulo, que o problema
não está na magnitude da correção (que sub-relaxação resolveria) mas na
**direção**: o operador de correção aponta para um lugar que o operador de
divergência não reconhece como progresso.

**Revertido por inteiro.**

**O que isso deixa para a próxima tentativa**: nem trocar o operador nem
amortecer o existente resolveu sozinho. O caminho que resta e não foi tentado
é medir os **autovalores** (não só o raio espectral) do operador
projeção-sozinha nas malhas ruins — se de fato há um par complexo dominante,
a estratégia certa não é escalar a correção, é girar a direção (um passo tipo
GMRES/Anderson na correção, não Richardson simples). Isso é uma ferramenta de
diagnóstico nova (não existe hoje: só a magnitude ρ é medida, via potência
iterada), e vale mais que uma terceira tentativa às cegas.


**Correção do que este documento afirmou antes**: "sugere autovalores
complexos" estava errado, e a medição que o corrige está registrada abaixo.
Isso importa porque orientou mal a segunda tentativa.

**A medição que corrige**: rastreando o ângulo entre iterados sucessivos da
potência iterada (não só a magnitude), `cos(ângulo) → 1,00000` em toda malha
testada — a direção do estado converge para um vetor fixo, não gira. Autovalor
dominante **real**, não complexo. A composição interpolação-para-face →
Poisson → gradiente-de-célula não tem a estrutura que a hipótese anterior
propunha.

**Terceira tentativa, mais cuidadosa que a segunda**: em vez de amortecer a
partir de α=1 (que a medição isolada mostrou divergir 2,5–5,8× num único
corretor extra, mesmo em malha boa), medi α diretamente antes de tentar de
novo:

| α | resíduo depois / antes (malha boa, vários passos) |
|---|---|
| 1,0 | 2,5×–5,8× |
| 0,4 | 1,2×–2,4× |
| 0,2 | 0,8×–1,3× |
| 0,15 | 0,8×–1,1× |

Contração real, ainda que fraca, em torno de α ≈ 0,15–0,2 — medida um
corretor de cada vez, a partir do estado que o primeiro corretor já produziu.
Implementei o laço com α inicial 0,15 (não 1,0) e até 16 correctores (a
contração fraca, ~0,85×, precisa de mais iterações que a correção
não-ortogonal do Poisson, que corta pela metade a cada uma).

**Também não funciona, e de um jeito que ensina algo específico**: mesmo na
malha *boa* (jitter 0,25, que sempre convergiu bem), em marcha real no tempo
(sem o protocolo sintético de potência iterada) a divergência **cresce**:
8,5e-02 no passo 0 → 2,6e-01 no passo 14, sempre esgotando os 16 correctores,
α nunca saindo de 0,15 — ou seja, o laço nunca detecta que está piorando o
suficiente para amortecer mais, mas está.

**Por que a medição isolada (que era boa) não previu isto**: o teste de um
corretor mede a contração a partir de um resíduo *genérico*, cuja direção é
dominada pelas componentes que o primeiro corretor já resolveu bem. **Cadear
muitos correctores filtra essas componentes e isola progressivamente a que
sobrou** — exatamente o mesmo mecanismo, em espírito, de por que a potência
iterada revela o autovalor dominante em vez da média: iteração repetida
converge para o pior caso, não para o caso típico. Se essa componente
remanescente tiver, sob α=0,15, um fator de amplificação *acima* de 1 — o que
a medição de um único corretor, dominada por outras direções, não conseguiria
enxergar — o laço a persegue e amplifica precisamente por insistir nela passo
após passo.

**O que isto estabelece, e é a conclusão real deste item**: **relaxação
escalar não pode consertar isto**, porque o operador tem comportamento
diferente por direção no espaço de resíduos — nenhum α único serve para
todas. Um remédio de fato precisa de um método que adapte a *direção* da
correção, não só sua magnitude: um Krylov genuíno (GMRES ou BiCGSTAB, que este
projeto já tem implementado para `ImplicitConvectionDiffusionSolver1D`)
aplicado ao sistema acoplado pressão-velocidade não-simétrico, em vez de CG
mais um escalar por fora. É trabalho de escopo maior — um solver novo, não um
parâmetro — e não foi tentado aqui.

**As três tentativas, juntas, delimitam o problema com precisão**: não é a
malha em si (item 3.3 mostrou a correção não-ortogonal convergindo até
não-ortogonalidade 563 com o amortecimento certo — aquele problema *é*
escalar). É especificamente o acoplamento pressão-velocidade desta projeção,
que não é.


**Quarta tentativa: o Krylov genuíno que a terceira tentativa pedia — e ele
resolve o problema errado.** Este é o achado mais importante do item até
agora, porque não é "a correção não funcionou": é "a correção estava mirando
a métrica errada o tempo todo", o que muda o diagnóstico de toda a
investigação anterior.

**O que foi construído.** `UnstructuredFvmBase` ganhou um `gmres()` genérico
— GMRES(m) reiniciado, matriz-livre, extraído fielmente do
`ImplicitConvectionDiffusionSolver1D::solveGmresOnce()` já testado (Gram-Schmidt
modificado, rotações de Givens incrementais). `UnstructuredCavitySolver3D`
ganhou um segundo corretor, depois do corretor CG original (mantido
bit-idêntico): mede o resíduo de divergência remanescente após a primeira
correção e, se não-desprezível, resolve por GMRES um `correctionOperator`
construído para ser a parte genuinamente linear de `projectionDivergenceRhs`
— a divergência que uma correção adicional de pressão `x` introduziria via
`-dt·∇x`. Dois bugs de implementação foram encontrados e corrigidos *durante*
esta mesma tentativa, o que já é informativo por si:

1. **Erro de sinal**: a primeira versão resolvia `correctionOperator(δp) =
   +resíduo` em vez de `= -resíduo` — a correção *somava* à divergência em vez
   de cancelá-la, e como isso realimenta a cada passo de tempo, o campo
   explodia para ~1e60 em poucos passos mesmo em malha boa. Corrigido.
2. **Operador que mentia sobre ser linear**: `computeCellGradients()` aplica
   um corte não-linear (limita o recuo Green-Gauss à maior inclinação que a
   célula enxerga — exatamente o limite do item 1.3/4.3 anterior) nas células
   de estêncil deficiente. Usar essa função dentro do operador do GMRES viola
   a premissa fundamental de um método de Krylov: com o sinal corrigido, a
   malha boa (415 células, 29 estênceis deficientes) que sempre convergiu bem
   **regrediu** — divergência 1,86e-01 contra o teto de 5e-2 que ela sempre
   cumpriu. Corrigido expondo uma variante `computeCellGradients(field,
   clampFallback=false)` usada só dentro do operador de correção, com o
   incremento de velocidade aplicado depois do mesmo jeito (não recomputando
   o gradiente clampeado do total).

**Depois dos dois bugs corrigidos, medi de novo — e o número não mudou em
nada perceptível** (divergência 1,860e-01, idêntica às três primeiras casas
decimais, antes e depois do segundo bug). Isso é o sinal de que havia um
terceiro problema, mais fundamental que os outros dois, e a instrumentação
direta (script Python chamando `last_coupling_iterations()`,
`last_coupling_residual_before/after()`, `max_face_divergence()` a cada
passo) o expôs:

| passo | iterações GMRES | resíduo antes (alvo do GMRES) | resíduo depois | `max_face_divergence()` (métrica real) |
|---|---|---|---|---|
| 0 | 200 (teto) | 0,0729 | 2,29e-09 | 0,4358 |
| 5 | 200 (teto) | 0,0784 | 2,17e-09 | 0,2923 |
| 14 | 200 (teto) | 0,0849 | 2,19e-10 | 0,2239 |

**O GMRES resolve exatamente o que foi pedido — o resíduo que ele mira cai de
~0,08 para ~1e-9 a cada passo, um fator de 10⁷ — e a divergência real
(`max_face_divergence()`) mal se move.** A razão é estrutural, não numérica:
`projectionDivergenceRhs()` (o alvo do GMRES, herdado sem alteração da
construção original do lado direito de Poisson) chama
`faceMassFluxes(velocity, pressure_, 0.0)` — **dt = 0 fixo no termo
Rhie-Chow**, decisão original e correta para montar o RHS a partir do
preditor sem pressão. Mas `max_face_divergence()` — a métrica que de fato
importa, e a mesma que a convecção usa para o fluxo de massa por face — chama
`faceMassFluxes(velocity_, pressure_, lastDt_)`, com o **dt real** no termo
Rhie-Chow. São duas definições de fluxo diferentes por construção, não por
acidente: o corretor CG original **já zera o resíduo dt=0 até 1e-10** (é
literalmente o critério de parada da correção diferida) — por isso o GMRES,
mirando esse mesmo resíduo, não tinha quase nada a fazer, e resolveu-o ainda
melhor sem que isso tivesse qualquer efeito sobre a divergência que a
instabilidade do item 4.3 de fato mede.

**Revertido por inteiro** (as cinco alterações — `gmres()`,
`projectionDivergenceRhs()`, o `correctionOperator`, os três getters de
diagnóstico e seus bindings Python) via cópia do conteúdo commitado por cima
do working tree, com `git diff` confirmando zero diferença de conteúdo antes
do commit deste registro. Suíte de volta a 12/12.

**O que isto muda para a próxima tentativa.** A hipótese da não-adjunção
(gradiente de célula vs. divergência de face) segue de pé, mas incompleta: há
uma segunda incompatibilidade, independente dela, entre **duas definições de
fluxo** dentro do próprio código — uma com o termo Rhie-Chow em dt=0 (usada
para montar e resolver o Poisson) e outra com dt real (usada para medir a
divergência e para a convecção). Um GMRES apontado corretamente precisa
resolver contra a *segunda* definição, não a primeira — o que por sua vez
exige decidir se o RHS original do Poisson (que usa dt=0 de propósito, por um
motivo documentado e válido: o preditor não carrega correção de pressão
nenhuma) também deveria mudar, ou se apenas o resíduo pós-correção deve ser
reavaliado com dt real antes de alimentar o Krylov. Qualquer um dos dois
caminhos é uma mudança na formulação do Poisson desta projeção, não mais um
ajuste no corretor extra — escopo maior do que as quatro tentativas feitas
até aqui, e não tentado.

**As quatro tentativas, juntas, delimitam o problema com mais precisão que
as três**: não é a malha em si (item 3.3). Não é falta de um método de
Krylov (quarta tentativa: o Krylov funciona perfeitamente para o problema que
recebeu). É o acoplamento pressão-velocidade desta projeção **e** uma
inconsistência de definição de fluxo dentro dele, ambas ainda sem correção.


**Quinta tentativa: mirar a definição de fluxo certa — resolve para domínio
fechado, expõe uma nova causa em domínio com saída.**

A quarta tentativa mostrou exatamente o que faltava: o GMRES precisa mirar
`faceMassFluxes(·, ·, dt real)` — a mesma definição que `maxFaceDivergence()`
usa — não `faceMassFluxes(·, ·, 0)`, que o primeiro corretor já zera por
construção. `projectionDivergenceRhs()` virou `divergenceOfFlux(velocity,
pressure, fluxDt, outerDt)`, com `fluxDt` e `outerDt` independentes:
`fluxDt=0` reproduz o RHS original do primeiro corretor (inalterado,
bit-idêntico); `fluxDt=dt` reproduz exatamente `lastDt_`, o que
`max_face_divergence()` mede. O `correctionOperator` agora varia **os dois**
argumentos de `faceMassFluxes` (velocidade *e* pressão) — a quarta tentativa
variava só a velocidade e deixava a pressão em `pressure_`, descartando
silenciosamente a própria contribuição do termo dependente de `fluxDt` à
linearização.

**Medido — raio espectral em torno do repouso, malha fechada:**

| malha | antes (quarta tentativa / linha de base) | depois (quinta tentativa) |
|---|---|---|
| jitter 0,25 | 0,981 | **0,127** |
| jitter 0,45 | 0,974 | **0,748** |
| jitter 0,55 | 5,19 | **0,476** |
| jitter 0,65 | 44,05 | **0,580** |
| jitter 0,75 | 12,39 | **0,270** |

Toda malha testada, incluindo as duas que mais divergiam, cai bem abaixo de
1. Cavidade n=4 (415 células): divergência por faces caiu de dentro da
margem que já passava para **6,19e-05**. `cavidade fechada` nos bindings
Python: **4,35e-13** — precisão de máquina, quando antes desta rodada era
3,49e-01 na mesma malha.

**Dois bugs de implementação apareceram e foram corrigidos no caminho**: o
mesmo erro de sinal da quarta tentativa (`correctionOperator(δp) = -resíduo`,
não `= +resíduo`) e o mesmo gradiente clampeado sem linearidade
(`computeCellGradients(x, false)` para o operador, incremento aplicado do
mesmo jeito na velocidade). Ambos já eram conhecidos da tentativa anterior;
registrar que reapareceram é o motivo de tê-los corrigido rápido desta vez.

**O que não fecha: domínio com saída.** No canal, o GMRES estagna cedo — 63
de até 182 iterações possíveis num único ciclo de reinício, **não** por
esgotar o orçamento (aumentar o reinício de 30 para 200 e o teto de 200 para
2000 iterações não mudou o resultado em nada perceptível). É quebra
(`breakdown`) ou pivô numericamente nulo, não estagnação por falta de
direções de Krylov. Rastreado a uma causa específica:
`computeCellGradients(x, false)` — usado tanto para construir o operador
quanto para aplicar a correção — busca o valor de contorno de uma face de
saída no estêncil de mínimos quadrados a partir de `boundaryValueCache_`,
que foi fixado em `outletPressure_` na construção do solver. Para o campo
físico de pressão isso é exatamente certo; para `x`, um campo de *correção*
que deveria responder a si mesmo (idealmente contorno zero, não
`outletPressure_`), o estêncil usa um valor constante e alheio ao próprio
`x`, o que deixa o operador mal-condicionado exatamente onde o posto de um
domínio com saída normalmente vem. Medido, não só suposto: o desbalanço de
massa do canal, que fechava a 6e-14 antes desta rodada, ficou em torno de
0,01–0,1% e não melhorou com mais iterações — consistente com o operador
convergindo para o resíduo mínimo que essa quase-singularidade permite, não
com falta de orçamento.

**Decisão tomada por medição, não por suposição**: a correção de acoplamento
agora só roda quando `pinnedCell != kNoPinnedCell` (domínio fechado, sem
saída — exatamente onde a quinta tentativa foi provada). Num domínio com
saída o passo volta a ser bit-idêntico ao comportamento anterior a todo este
item — verificado: canal fecha a 6,475e-14 nos testes C++ e a ~4e-13 nos
testes Python em 2, 4 e 16 correctores (era 6e-14/1e-13 antes de qualquer
tentativa deste item), e a malha muito distorcida (jitter 0,65, com saída)
continua recusando em vez de propagar NaN no passo 25, exatamente como antes.
Suíte: 12/12.

**O que fica para a próxima tentativa**: uma variante de
`computeCellGradients` cujo estêncil, numa face de saída, use o valor de
contorno consistente com o campo que está sendo diferenciado -- zero para um
campo de correção, `outletPressure_` só para a pressão física -- em vez de
sempre `boundaryValueCache_`. Com isso o `correctionOperator` recupera posto
pleno numa saída do mesmo jeito que o Laplaciano original recupera via
`boundaryEntersOperator`, e a quinta tentativa deveria se generalizar sem
precisar de uma sexta reformulação.

**Custo**: a suíte C++ foi de ~26s para ~55-60s -- o GMRES roda a cada passo
em toda malha fechada, mesmo as que já convergiam bem, porque medir o
resíduo real primeiro (e pular quando desprezível) ainda não foi ligado a
essa decisão; hoje só a malha ter ou não saída decide se ele roda. Aceitável
para a suíte atual; vale revisitar se um caso maior tornar isso o gargalo.

**Sexta tentativa: um defeito real do operador, corrigido — e ainda assim
o caso com saída não fecha.**

A quinta tentativa deixou registrado que o GMRES estagnava em domínio com
saída por uma direção quase-nula, e apontou a suspeita:
`computeCellGradients` alimenta células adjacentes à saída com o valor
*prescrito* `outletPressure_` como valor do vizinho de contorno — mesmo
quando o campo sendo diferenciado é uma **correção** de pressão, cujo valor
num contorno de Dirichlet é necessariamente zero (o campo base já satisfaz
a pressão prescrita exatamente; qualquer correção adicional ali violaria
essa condição).

**A suspeita estava certa, e o defeito é real**: usar o valor prescrito
para um campo de correção torna o gradiente **afim, não linear** — o
`outletPressure_` entra como termo constante. É exatamente a mesma classe
de erro que o clamp de recuo Green-Gauss (quarta tentativa), pela mesma
razão: um método de Krylov recebe um operador que mente sobre ser linear e
converge com confiança para a resposta errada.

Corrigido com `homogeneousBoundaryValues`, um parâmetro que substitui 0
pelo valor prescrito em toda face de contorno que carregue um, propagado
por `computeCellGradients`, `greenGaussGradients`, `faceMassFluxes` e
`divergenceOfFlux`. **A correção fica, independentemente do que ela
conserta**, porque é correta por si só — e é comprovadamente inócua em
domínio fechado, onde nenhuma face de contorno carrega pressão prescrita
(só saídas carregam, ver `setBoundaryFaceValue` no construtor). Verificado:
cavidade 6,484e-05 e canal 6,475e-14, idênticos aos valores anteriores à
mudança.

**Mediu-se a melhora, e ela não basta.** No canal, o resíduo por passo
depois da correção caiu de **0,047 para 0,033** (~30%). Mas o GMRES ainda
para em **64 iterações contra um teto de 728** — quebra numérica, não
orçamento esgotado, então sobra pelo menos mais uma direção quase-nula.

**E aplicar a correção nesse estado é regressão, não vitória parcial**: o
desbalanço de massa do canal vai de **6,5e-14** (correção pulada) para
**1,2e-4** (correção aplicada mas estagnada). Por isso o portão para
domínio fechado **permanece**, agora por um motivo mais estreito e melhor
entendido do que "o caso com saída não funciona".

**E o diagnóstico seguinte foi feito na mesma sessão, corrigindo o que o
parágrafo acima presumia.** A suposição era que sobrava uma direção
quase-nula e que o GMRES quebrava numericamente. Instrumentei os dois
caminhos de quebra (Givens e "happy breakdown") e **nenhum dispara**: o
GMRES *converge*, à sua própria tolerância relativa de 1e-8, em 64
iterações. Não há singularidade.

Uma segunda hipótese — que o piso de resíduo vinha do clamp em células de
estêncil deficiente, onde o operador usa o gradiente sem clamp e o resíduo
recomputado usa o clampeado — também foi **descartada por medição**: dois
domínios fechados com **zero** estêncils deficientes (911 e 3493 células)
continuam reduzindo o resíduo por um fator de 100–370 em vez de ir a
precisão de máquina.

**A causa real, encontrada relendo o código com essa pista**: numa face de
contorno as duas quantidades divergem por construção.

| quem | fluxo numa saída |
|---|---|
| `divergenceOfFlux` (o que o GMRES zera) | `velocity[cell] · areaVector` |
| `maxFaceDivergence` (o que reporta conservação) | `boundaryFlux_[i]`, a fórmula com mistura de gradiente |

É **a mesma classe de erro que a quarta tentativa encontrou** entre
`fluxDt=0` e `fluxDt=dt`, um nível para fora: *medir o operador que de
fato se resolve* — agora pela terceira vez neste item, em faces de
contorno em vez de faces internas. Reconciliar as duas definições de fluxo
de saída é o próximo passo concreto, e é uma pergunta mais estreita que
qualquer formulação anterior desta lacuna.

**Sétima tentativa: RESOLVIDO. A causa era a face de saída não receber a
correção de Rhie-Chow que toda face interna já recebia.**

A sexta tentativa terminou com o diagnóstico correto e a correção errada.
Duas medições fecharam a questão:

1. **Instrumentei os dois caminhos de quebra do GMRES e nenhum dispara.**
   O método *converge*, à sua tolerância relativa de 1e-8, em 64 iterações.
   Isso derrubou a afirmação anterior deste documento de que havia uma
   direção quase-nula sobrevivente: não há singularidade nenhuma.
2. **A hipótese seguinte — de que o piso vinha do clamp em células de
   estêncil deficiente — também caiu por medição**: dois domínios fechados
   com **zero** estêncils deficientes mantinham o mesmo piso.

Se o GMRES converge e ainda assim o resíduo medido não cai, então o
operador resolvido e a quantidade medida **não são o mesmo operador**. E
numa face de saída eles provavelmente divergiam:

| quem | fluxo numa saída |
|---|---|
| `divergenceOfFlux` (o que o GMRES zerava) | `velocity[cell] · A` |
| `maxFaceDivergence` (o que reporta massa) | `boundaryFlux_[i]`, com mistura de gradiente |

**A correção, e por que ela é pequena.** Escrevendo o fluxo *que a projeção
de fato impõe* em termos da velocidade já corrigida:

```
imposto = velocityStar·A − dt (grad p)_f·A
        = velocity·A + dt (grad p)_P·A − dt (grad p)_f·A
        = velocity·A + dt · unitD ((grad p)_P·unitD − dPdn) · A
```

que é **exatamente o termo de correção que o laço de faces internas já
aplicava**, com a pressão prescrita do contorno no lugar do vizinho. Ou
seja: não faltava um operador novo, faltava aplicar o operador existente
numa face onde ele havia sido esquecido. E como o termo é proporcional a
`fluxDt`, o caminho do lado direito do Poisson (`fluxDt = 0`) fica
inalterado **por construção**, não por um guarda que alguém precise
lembrar — verificado: canal continua em 6,475e-14 nesse caminho.

**Medido no canal com saída (177 células):**

| | antes | depois |
|---|---|---|
| divergência por faces | 3,0e-02 | **1,1e-11** |
| fluxo líquido de contorno | 1,2e-04 | **1,2e-13** |
| resíduo de acoplamento por passo | 0,033 | 2,7e-11 |

E nos testes Python o canal fecha em **1,24e-13 com 2, 4 e 16
correctores** — o número deixou de depender da contagem de correctores,
que é o sinal de que a correção não está mais compensando um erro.

**Varredura de jitter em malha com saída, raio espectral em torno do
repouso pelo protocolo do próprio item:**

| jitter | não-ortog. | ρ antes | ρ depois | 400 passos |
|---|---|---|---|---|
| 0,45 | 2,07 | 0,974 | 0,712 | ok |
| 0,65 | 3,53 | **44,05** | **0,286** | ok |
| 0,75 | 3,01 | 12,39 | 0,208 | ok |
| 0,85 | 10,04 | — | 0,495 | ok |
| 0,95 | 7,40 | — | 0,580 | **morre no passo 41** |

**A instabilidade linear que abriu este item acabou.** Era o que a medição
original descreveu — fator ≈2,05 por passo em torno do repouso, propriedade
da malha e do operador — e hoje ρ < 1 em toda malha testada.

**Um erro de protocolo corrigido no caminho, que vale registrar**: a
primeira versão desta varredura mediu ρ com a *entrada ligada a 1,0* e
obteve valores absurdos (5.000 a 42.000) enquanto a marcha real rodava
400 passos com divergência 1e-06. Com entrada não-nula o repouso **não é
ponto fixo**, então a razão de normas mede o termo forçante, não o
operador. É a mesma armadilha de "medir a coisa errada" que este item
armou três vezes; desta vez a contradição entre as duas medidas a expôs
em minutos.

**O que continua aberto, e é outra coisa** — agora caracterizado.

Jitter 0,95 ainda morre no passo 41 numa marcha *dirigida*, apesar de
ρ = 0,58 em repouso. Não é não-ortogonalidade (0,85 roda 400 passos com
não-ortogonalidade 10,04 enquanto 0,95 morre com 7,40). Decomposto em
`python/research/jitter095_failure.py`:

| teste | resultado |
|---|---|
| passo de tempo × 0,25 e × 0,05 | **não salva** (64 e 50 passos contra 41) |
| viscosidade × 3, × 10, × 30 | **não salva** (59 e 53 passos) |
| sem convecção (`stepWith`) | sobrevive 237 passos, ~6× mais |
| sem projeção (`stepWith`) | sobrevive os 400 |
| entrada 1e-6 e 1e-3 | **sobrevive os 400** |
| entrada 0,01 / 0,05 | **sobrevive** / morre em 50 |

**É não-linear, e isso o separa definitivamente deste item.** A
instabilidade que abriu o 4.3 foi *provada linear* pela própria medição de
amplitude: crescia com o mesmo fator por passo com a tampa a 1,0 e a
10⁻⁶. Esta faz o oposto — tem um limiar nítido entre entrada 0,01 e 0,05.
Um crescimento que só aparece acima de uma amplitude crítica é a assinatura
de uma instabilidade **subcrítica não-linear**, não de um operador linear
com raio espectral acima de 1.

Que não seja passo de tempo nem viscosidade descarta as duas explicações
mais baratas: não é CFL e não é balanço convecção-difusão (Re de célula).
Convecção participa (removê-la estende a sobrevida 6×) mas não é
suficiente sozinha, já que sem ela ainda morre no passo 237.

**Não está explicado, e é deliberadamente registrado como item separado**
e não como resquício deste — a lição do próprio 4.3, três vezes, é que
juntar dois fenômenos sob um rótulo atrasa os dois. O guarda de finitude
continua sendo o piso que transforma isso em recusa em vez de campo de
NaN, e o teste dele subiu de jitter 0,65 para 0,95 pela terceira vez,
sempre porque o solver melhorou.

### 4.4 ~~Diferença central na convecção estruturada~~ — RESOLVIDO em 2026-08-16 no solver 2D; medido em 2026-08-26 nos seis 3D, sem justificar trocar o default

Os solvers estruturados usam diferença central na convecção, e a cavidade a
Re=400 roda com Re de célula 12,5 (medido, não estimado). Acima de 2,
diferença central produz oscilação — propriedade do *esquema*, não do passo, e
o fator de segurança do item 4.1 não a toca.

**O defeito deixou de ser citação e virou medição.** Toda a literatura afirma
isso, e este projeto vinha repetindo a afirmação em comentários desde a Fase 1
sem nunca tê-la medido — exatamente o estado em que estavam a ordem do
Laplaciano (item 3.2, a inferência estava certa) e a difusão do upwind (item
3.1, também certa). Só medir distingue.

Caso com resposta em forma fechada: convecção-difusão 1D estacionária sem
fonte, cuja solução exata é derivada da própria EDO no cabeçalho do solver. É
**monótona** entre os dois valores de Dirichlet, então qualquer excesso fora
dessa faixa é o esquema falhando — princípio do máximo, que é teorema e não
tolerância. É o que torna este caso decisivo onde a cavidade não é: a cavidade
não tem solução exata, então um campo suspeito lá sempre dá margem a discussão.

| Pe de célula | upwind | central | limitado |
|---|---|---|---|
| 0,1 | 3,59e-03 | 2,07e-04 | **1,71e-04** |
| 0,5 | 2,21e-02 | 6,83e-03 | **3,10e-03** |
| 2,0 | 2,91e-02 | 5,87e-02 | **7,17e-03** |
| 8,3 | 2,83e-02 | 6,36e-01 · **excesso 3,17** | **1,54e-02** · excesso 0 |

Três coisas se leem daí. Central **passa a ser pior que upwind** exatamente ao
cruzar Pe 2, que é onde a teoria diz. Em Pe 8,3 ele põe a solução a **3,17×
a faixa dos contornos fora dela** — não é imprecisão, é resposta sem sentido
físico. E o esquema limitado é o melhor dos três em **todo** Péclet, com
excesso zero em todos: 1,8× melhor que upwind e 41× melhor que central em
Pe 8,3.

**Resolvido em `ImplicitConvectionDiffusionSolver1D`**, que ganhou
`setConvectionScheme()` com os três esquemas. A matriz montada continua sendo
a de upwind qualquer que seja a escolha, e a diferença vai para o lado direito
por **correção diferida** em varreduras externas — a mesma mecânica que a
`UnstructuredFvmBase` usa para o termo não-ortogonal, pelas mesmas duas
razões: a matriz preserva a propriedade de M-matriz que faz o Krylov se
comportar, e um esquema *não-linear* (qualquer limitador) não cabe num
operador linear. Converge em 1 a 3 varreduras. O default segue
`FirstOrderUpwind`, então o comportamento existente da classe é idêntico.

**O limitador passou a ser um só** (`ConvectionLimiter.hpp`), compartilhado
entre o lado estruturado e o não-estruturado. Vale registrar por que as duas
razões do limitador são a mesma: a clássica estruturada é
`r = (φ_C − φ_CC)/(φ_D − φ_C)`, que precisa de uma célula a-montante-da-montante
que uma malha tetraédrica não tem; a não-estruturada usa o gradiente da célula
a montante. Substituindo o gradiente por diferença central numa malha
uniforme, `∇φ_C·d = (φ_D − φ_CC)/2`, a segunda vira a primeira **exatamente**.
Não são duas convenções.

**O que falta, e por que parei aqui**: portar para `LidDrivenCavitySolver2D` e
`StaggeredCavityBase3D` (que seis fechamentos de turbulência herdam). Isso é
reescrever o preditor de forma não-conservativa `u·∇u` para forma conservativa
de fluxo de face, mudando o resultado de sete solvers validados — e **a suíte
atual não conseguiria mostrar o ganho**: seus casos rodam a Re de célula 0,31
a 3,12, onde central ainda não tem excesso nenhum. Trocar sete resultados
validados para corrigir um defeito que os testes existentes não exercitam é a
troca errada de fazer às cegas.

**O que isso exige antes**: um caso de cavidade a Re=400 na suíte — que é
onde o Re de célula 12,5 aparece — com um diagnóstico que o defeito mova. O
índice de checkerboard da pressão foi medido e cresce com o Re de célula
(3,4e-03 a Re=10, 9,8e-03 a Re=100, 1,96e-02 a Re=400, caindo para 8,5e-03 ao
refinar a malha no mesmo Re), então é candidato — mas checkerboard em malha
colocada é sintoma de acoplamento pressão-velocidade, não necessariamente do
esquema de convecção, e atribuí-lo sem separar as duas causas seria repetir o
erro que os itens 1.1 e 1.2 já custaram caro.

---


**Medido em 2026-08-16 o que faltava para decidir o port — e ele não se
justifica ainda.** A pergunta era se a diferença central causa dano
*mensurável* na cavidade estruturada no regime em que ela é teoricamente
ilimitada. Cavidade 2D a Re=400, marchada a t=12, comparando o perfil de u na
linha vertical central entre resoluções:

| n | Re de célula | E cinética | desvio do perfil vs. n anterior |
|---|---|---|---|
| 16 | 25,00 | 6,978e-03 | — |
| 32 | 12,50 | 1,337e-02 | 0,0445 |
| 64 | 6,25 | 1,619e-02 | 0,0180 |

**Ela converge.** O desvio cai pela metade a cada refino e o campo não
apresenta o colapso que o caso 1D exibe acima de Pe 2 — lá, com a resposta
exata conhecida, a diferença central põe a solução a 3,17× a faixa dos
contornos fora dela; aqui, a Pe 12,5 e 25, não há excesso visível. A diferença
provável é que o campo da cavidade é suave em relação à malha, enquanto o caso
1D impõe uma camada limite exponencial íngreme.

**Por que isso decide contra o port por ora**: trocar o preditor de sete
solvers validados — `LidDrivenCavitySolver2D` mais os seis que herdam de
`StaggeredCavityBase3D` — exige uma medição que mostre o ganho, e a que existe
mostra que o dano no regime real ainda não foi quantificado. Atribuir o erro
de 21% em energia cinética entre n=32 e n=64 ao *esquema* seria confundi-lo
com falta de resolução, que é exatamente o erro que os itens 1.1, 1.2 e 4.3
custaram caro.

**O que fecharia**: implementar o esquema conservativo limitado no
`LidDrivenCavitySolver2D` (os fluxos de face já existem em
`rhieChowFaceU/rhieChowFaceV`, então a forma conservativa está ao alcance),
mantendo o default atual, e comparar os dois contra a solução em n=128 no
mesmo Re. Se o limitado ficar mais perto da referência em n=32, o port se
justifica e o default vira dado; se não, o item vira nota de rodapé.

---


**Fechado para o `LidDrivenCavitySolver2D`, e o default só mudou depois de
medido.** A pergunta que faltava era se o esquema limitado de fato ganha *na
cavidade*, e não só no caso 1D com resposta exata. Referência: a própria
cavidade a Re=400 em n=128 (Re de célula 3,1). Desvio rms de u na linha
vertical central, marchando a t=12:

| | n=16 (ReCel 25,0) | n=32 (ReCel 12,5) |
|---|---|---|
| central | 0,066457 | 0,023450 |
| upwind 1ª ordem | 0,072535 | 0,040897 |
| **limitado** | **0,054422** | **0,021399** |

O limitado é o melhor nas duas resoluções: 18% melhor que a central em n=16,
9% em n=32, e 25% a 48% melhor que o upwind. **Não é só mais seguro, é mais
exato** — e isso, somado à limitação que a central não oferece, é o que
justificou mexer no default de um solver validado. A forma conservativa usa
os mesmos fluxos de face que a projeção equilibra, porque usar um fluxo aqui e
outro lá é exatamente o erro que os itens 1.1 e 1.2 custaram.

`ConvectionScheme::Central` continua disponível: todo número que este solver
produziu antes foi medido com ela, e reproduzi-los tem de continuar possível.

**Um efeito que só o teste de persistência pegaria**: o preditor conservativo
usa `lastDt_`, e um solver recém-carregado não tem passo anterior — retomar de
checkpoint deixou de reproduzir bit a bit uma corrida ininterrupta. O recuo é
usar o dt corrente quando não há anterior, o que restaura a igualdade exata.
Vale registrar que foi um teste de *outra camada* que pegou.

**O que fica**: os seis fechamentos 3D construídos sobre
`StaggeredCavityBase3D` continuam com diferença central. Eles compartilham um
preditor diferente, e portá-lo é a sua própria medição — a régua e o esquema
agora existem, então é trabalho de repetir o que está acima, não de inventar.

---

**Atualização 2026-08-26: portado para os seis 3D, com o mesmo cuidado do
2D — mas a medição não repetiu o resultado do 2D, e o default não mudou.**

`StaggeredCavityBase3D::ConvectionScheme` (`Central`/`FirstOrderUpwind`/
`LimitedLinearUpwind`) foi adicionado ao preditor conservativo compartilhado
pelos seis solvers de cavidade 3D, reaproveitando `ConvectionLimiter.hpp`
(o mesmo limitador de van Leer do 2D e do lado não-estruturado). `Central`
continua default: **cada um dos nove termos de fluxo do preditor (self e
cruzado, para u/v/w) foi reescrito para que o caso `Central` devolva o
valor central original, sem nenhuma aritmética nova** — não apenas igual,
literalmente o mesmo `double`, o que faz o default reproduzir todo número
que este preditor já produziu bit a bit. Verificado, não assumido: as 13
suítes passam sem alteração de nenhum resultado numérico, incluindo os
testes dos seis solvers que herdam este preditor.

Uma diferença de geometria em relação ao 2D exigiu uma decisão que o 2D não
precisou tomar: a direção própria de cada componente (i para u, j para v,
k para w) não tem célula-fantasma além das duas faces físicas de contorno
— ao contrário das direções tangenciais, que já espelham corretamente para
qualquer distância. O ponto "mais distante" do limitador nessa direção
própria é então fixado (`std::clamp`) na face de contorno mais próxima em
vez de lido fora do domínio, o que é seguro (nunca acessa memória fora do
array) e apenas reduz o limitador ao upwind na única célula adjacente a
cada parede — degradação padrão e localizada, não um defeito.

**Medição, mesma régua do 2D**: cavidade a Re=400 (viscosidade =
lidVelocity·L/Re), marchando a t=12, desvio rms de u na linha vertical
central. O 2D usou n=128 como referência (Re de célula 3,1); um equivalente
3D não cabe em tempo razoável (custo escala com n⁴, não n²), então a
referência aqui é n=32 (Re de célula 12,5) — mais grosseira que a do 2D em
termos absolutos. Script preservado em
`python/research/convection_scheme_3d_cavity.py`.

| esquema | n=8 (ReCel 50,0) | n=16 (ReCel 25,0) |
|---|---|---|
| central | 0,081547 | **0,041718** |
| upwind 1ª ordem | 0,079147 | 0,047443 |
| limitado (van Leer) | **0,079598** | 0,044765 |

**Ao contrário do 2D, não há vencedor claro — e a razão é mensurável, não
uma desculpa.** Central é o melhor em n=16; upwind e limitado disputam o
melhor em n=8, todos dentro de 3% um do outro. O que decide isso não é o
esquema: `peak|u|` entre n=8, n=16 e n=32 sobe 0,222 → 0,433 → 0,650 —
ainda mudando ~50% a cada duplicação de resolução, muito longe de
convergido. Nessas resoluções o erro de malha grosseira domina completamente
sobre a escolha de esquema, então comparar contra uma referência (n=32) que
ela mesma não convergiu não separa "esquema melhor" de "ruído de malha
grosseira". A diferença central-versus-limitado entre malhas consecutivas
também não cai monotonicamente (rms central-vs-limitado: 0,00935 → 0,01072
→ 0,00385 em n=8/16/32) — mais um sinal de que a comparação ainda não está
no regime assintótico onde ela seria decisiva.

**Decisão**: `Central` permanece default nos seis solvers 3D. Diferente do
2D, onde a medição decisiva (n=128, cell Re 3,1) justificou trocar o
default, aqui a medição afordável é honesta mas **inconclusiva** — trocar o
default sem essa justificativa repetiria exatamente o erro que os itens
1.1, 1.2 e 4.3 já custaram (mudar comportamento validado sem medição que
mostre ganho). `FirstOrderUpwind` e `LimitedLinearUpwind` ficam disponíveis
como opção explícita — quem precisa do limite formal (boundedness acima de
Re de célula 2) já pode escolher — mas nenhum dos seis fechamentos muda de
comportamento por padrão. Reabrir esta decisão exigiria uma referência 3D
genuinamente fina, o que é uma questão de orçamento de tempo computacional,
não de código.

---

## 5. Verificação incompleta

### 5.1 ~~O CI nunca rodou~~ — RESOLVIDO em 2026-08-27

O workflow existia e a metade local do portão estava provada (regressão
deliberada reprovou). Mas nada tinha sido enviado ao GitHub, então o CI
nunca tinha executado.

**Por que importava**: um workflow que nunca rodou é uma hipótese, não uma
rede de proteção. Erros de YAML, de caminho e de ambiente só aparecem na
primeira execução — e foi exatamente isso que aconteceu.

**O push já tinha sido feito numa sessão anterior — e o CI falhou nas 6
execuções seguintes, sempre no mesmo passo, sem que ninguém notasse.**
GitHub agora exige login para ler o log de uma execução mesmo em
repositório público, então a falha ficou invisível para qualquer coisa
sem uma sessão de navegador autenticada — inclusive para este agente.
Tentativas de reproduzir localmente os comandos exatos do workflow (Python
limpo do sistema, com e sem CUDA disponível) não reproduziram nada: o
build local sempre passava.

**A saída foi tornar o próprio erro público.** O passo `Configure` passou a
capturar a saída do `cmake` e, ao falhar, ecoar as últimas linhas como uma
única anotação `::error::` — que a API REST pública de check-runs devolve
sem autenticação, ao contrário do log em si. Isso revelou a causa real na
primeira tentativa:

```
CMake Error at CMakeLists.txt:39 (add_subdirectory):
  add_subdirectory given source "engine/testing" which is not an existing
  directory.
```

**A causa raiz não era de ambiente — era `engine/testing/` nunca ter sido
versionado, desde o primeiro commit deste repositório.** A regra `Testing/`
em `.gitignore` (pensada para o diretório de log que o próprio CTest gera)
casava também com `engine/testing/`: um padrão sem barra inicial casa em
qualquer profundidade, e em filesystem case-insensitive (Windows)
`Testing/` e `testing/` são o mesmo padrão. `engine/testing/` — a biblioteca
`aether_testing` com a macro `AETHER_CHECK` que as 13 suítes C++ usam —
ficou de fora do controle de versão inteiro esse tempo, e isso nunca
apareceu localmente porque o diretório já existia em disco antes mesmo do
primeiro commit. Só um clone de verdade — exatamente o que o runner do CI
faz a cada execução — expõe o problema.

**Corrigido** em duas partes: a regra virou `/Testing/` (ancorada na raiz,
não casa mais `engine/testing/`), e `engine/testing/` entrou no
repositório pela primeira vez. Verificado com um clone de verdade (não só
localmente confiado): `git clone` para um diretório novo, `cmake -S . -B
build` com os mesmos argumentos do workflow, configure limpo. **O CI
passou pela primeira vez na história deste projeto** na execução seguinte
(run #8, ~10 min, build completo + suíte de testes).

O passo de diagnóstico (as anotações `::error::`/`::notice::` no
`Configure`) foi deixado no lugar — descrito no próprio `ci.yml` como
"remover assim que Configure passar pela primeira vez", mas vale manter
mais uma rodada: é barato, e a próxima falha real (se houver) se beneficia
da mesma visibilidade sem depender de login.

### 5.2 ~~Solvers não-estruturados sem bindings Python~~ — RESOLVIDO em 2026-08-16

Todas as outras camadas do engine tinham bindings — é o invariante
arquitetural "hybrid C++/Python" desde o Módulo 1. `UnstructuredDiffusionSolver`,
`UnstructuredCavitySolver3D` e `TetrahedralMesh` não tinham.

**Por que não dava para contornar**: sem bindings, todo experimento com esses
solvers exige escrever e compilar C++, o que é exatamente o atrito que a
camada Python existe para eliminar. A investigação dos 13,2% teria sido
minutos em Python em vez de vários ciclos de build.

**Resolvido**: `TetrahedralMesh` e `TetrahedralFace` em
`aether_mesh_bindings.cpp`, os dois solvers em `aether_solver_bindings.cpp`,
exportados pelo pacote `aether`. Seletores de contorno, velocidade de parede
e predicado de saída atravessam como *callables* Python.

**E os bindings passaram a ser exercitados, não só compilados.**
`python/tests/test_unstructured_bindings.py` está registrado no ctest (12
suítes agora, era 11). Isso importa porque compilar um binding pega erro de
assinatura e nada mais: ele pode compilar e ainda devolver ponteiro pendurado,
perder um argumento default ou trocar a ordem do construtor. O teste checa o
que a suíte C++ estruturalmente não pode — a superfície do binding — mais duas
afirmações físicas que saem de graça (princípio do máximo do Laplaciano,
balanço de massa do canal). Inclui, de propósito, o caso `del mesh; gc.collect()`
seguido de `step()`: sem `keep_alive` isso é use-after-free, e largar a malha
depois de construir o solver é a coisa natural para um script fazer.

**O primeiro achado veio na primeira execução, que é o argumento inteiro do
item.** O teste falhou no balanço de massa do canal: 2,7e-04 da vazão, onde a
suíte C++ fecha em 1e-13. Três scripts Python e nenhum rebuild depois:

- **não é face descartada** — contei em Python as faces com `A·d ≤ 0`, que a
  montagem pula em silêncio: zero, nas duas malhas;
- **não é transiente** — t = 2, 4 e 8 todos em 1e-04;
- **é o número de correctores não-ortogonais**, e depende da malha:

| correctores | naoOrtog 0,707 | naoOrtog 1,651 |
|---|---|---|
| 4 | 2,6e-13 | 2,7e-04 |
| 16 | 2,7e-12 | 1,3e-07 |
| 64 | 2,7e-12 | 7,1e-14 |

Isso virou API: `pressureCorrectors` passou a ser argumento de construtor (com
o default documentado como "o que a suíte calibra", não como suficiência), e
`lastPressureChange()` é o diagnóstico que diz se a projeção convergiu ou só
acabou os correctores — mesmo papel que `lastOuterChange()` na difusão. Os
dois casos C++ reportam 2,9e-11 e 0,0e+00, ou seja, **estão convergidos** com
4; foi a malha que eu gerei em Python, mais torta, que não estava. Sem o
diagnóstico, os três números seriam indistinguíveis.

O teste passou a medir a *relação* em vez de um limiar: mais correctores têm
de fechar melhor, e o suficiente tem de chegar à precisão de máquina. Um
limiar fixo no default teria codificado a malha que aquele teste por acaso
constrói.

**Também encontrado ao escrever os bindings**: o callback de valor de
contorno era chamado de dentro do laço de Green-Gauss e da correção
não-ortogonal — milhares de vezes por passo, contradizendo a afirmação da
própria base de que nenhum laço interno paga despacho. Em C++ isso era só
ineficiência; com uma lambda Python seria o GIL milhares de vezes por passo,
o que tornaria os bindings inúteis para os experimentos que eles existem para
permitir. Passou a ser avaliado uma vez por face e cacheado. Números da suíte
idênticos; o canal de 5016 passos roda em 1,4 s a partir de Python (0,28
ms/passo).

### 5.3 ~~Caso da cavidade não-estruturada em malha grosseira por custo~~ — RESOLVIDO em 2026-08-16

n=3 (177 células) e t=4 foram escolhidos para caber no tempo de suíte, não
pela física. A afirmação de topologia é válida nessa resolução; qualquer
afirmação quantitativa não seria.

**O que fazer**: com a convecção implícita (4.2) e os bindings (5.2), rodar a
convergência de malha de verdade fora da suíte de testes.

---


**Resolvido: a convergência foi rodada, fora da suíte, e o item deixou de ser
afirmação.** Mesma cavidade, mesma família de malha, marchada a t = 8:

| n | células | u topo | u fundo | E cinética | Δ topo | Δ E |
|---|---|---|---|---|---|---|
| 3 (suíte) | 177 | +0,050504 | −0,011337 | 1,830e-03 | — | — |
| 4 | 415 | +0,069719 | −0,011931 | 2,716e-03 | +38,1% | +48,4% |
| 5 | 790 | +0,081655 | −0,011524 | 3,208e-03 | +17,1% | +18,1% |
| 6 | 1358 | +0,088595 | −0,011102 | 3,586e-03 | +8,5% | +11,8% |
| 7 | 2146 | +0,090057 | −0,011816 | 3,837e-03 | **+1,7%** | +7,0% |

**O item dizia "qualquer afirmação quantitativa não seria válida"; agora se
sabe o tamanho disso**: a malha da suíte está 44% abaixo do valor mais fino, e
o n = 4 do teste de topologia, 23% abaixo. O *sinal* está assentado em toda
resolução — que é exatamente a afirmação que o teste faz, e a razão de ser um
teste de topologia. A energia cinética ainda se move 7% entre as duas malhas
mais finas, então nem n = 7 está convergido.

A divergência por faces fica em 1e-12 em toda resolução: a projeção não é o
que limita.

**E a medição decidiu o item 4.2.** Custou 309 s em n = 7 contra 8 s em n = 4,
e o custo é o *passo*: dt cai de 7,3e-04 para 7,7e-05 entre n = 3 e n = 6 — e
então **sobe** para 8,3e-05 em n = 7, porque o limite é posto pelo sliver que
a tetraedralização por acaso produziu, não pela resolução. O item 4.2 dizia
"não urgente enquanto o limite convectivo for folgado". Ele não é mais.

### 5.4 ~~Projeção com contagem fixa de correctores~~ — RESOLVIDO em 2026-08-16

Registrado quando o canal precisava de 64 correctores não-ortogonais para
fechar o balanço de massa, e lido como "o número necessário é propriedade da
malha". **Não era**: com a pressão da saída no ajuste de mínimos quadrados
(item 2.3) o mesmo canal fecha em 1e-13 com dois correctores, em toda malha
testada até não-ortogonalidade 2,24. A correção estava sendo construída a
partir de um gradiente ruim, não lutando com a geometria.

**O que sobra**: a cavidade fechada, sem nenhuma face de Dirichlet, ainda cai
geometricamente com o número de correctores (2,28e-02 / 7,13e-03 / 1,29e-04
para 2 / 4 / 16). Uma contagem fixa continua sendo um botão em vez de um
critério, e `solveDeferredCorrection()` já para por tolerância — falta ela ser
relativa à escala da pressão em vez de 1e-10 absoluto, que quase nunca dispara
cedo.

**Por que não é urgente**: o diagnóstico existe (`lastPressureChange()`) e o
botão existe (argumento de construtor), então hoje é escolha informada e não
erro silencioso. E a evidência que fazia isso parecer urgente já foi explicada
por outra coisa — o que é, por si, motivo para não agir sobre a próxima
suspeita antes de medi-la.


**Resolvido, e a objeção que restava caiu junto com o 4.3.** O critério de
parada passou a ser **relativo à escala do campo**, com a forma absoluta
mantida só como piso para o campo identicamente nulo — um estado de repouso é
solução legítima e não tem escala. `tolerance` passou a significar a mesma
coisa em qualquer unidade, que é o que quem passa 1e-10 espera.

Isso só é seguro porque o item 4.3 mudou o quadro. Antes, **mais correctores
amplificavam**: 2 → 1,58, 4 → 2,06, 16 → 6,47 de raio espectral na malha
distorcida, então um critério que iterasse até convergir seria pior que o teto
fixo. Com o recuo Green-Gauss limitado, o raio espectral ficou **independente
da contagem**:

| malha | 2 corr. | 4 corr. | 16 corr. |
|---|---|---|---|
| jitter 0,25 | 0,9819 | 0,9819 | 0,9819 |
| jitter 0,35 | 0,9773 | 0,9773 | 0,9773 |
| jitter 0,45 | 0,9764 | 0,9760 | 0,9760 |

A amplificação por corrector era inteiramente o recuo ilimitado. Iterar voltou
a ser seguro.

**Medido**: a solução manufaturada cai de 90/84/105 varreduras para 87/80/99
com o erro **idêntico em todos os dígitos** (2,338802e-02 / 1,033868e-02 /
5,775569e-03). O ganho é modesto ali porque a escala do campo é ~2,7; ele
cresce com a magnitude do campo, que é exatamente o caso da pressão.

### 5.5 ~~GMRES duplicado entre o solver 1D e a base não-estruturada~~ — RESOLVIDO em 2026-08-25

Registrado quando a segunda cópia foi escrita, e não escondido: o item 4.3
precisou de GMRES para o acoplamento pressão-velocidade e o extraiu de
`ImplicitConvectionDiffusionSolver1D::solveGmresOnce()`, deixando duas
implementações completas do mesmo método na árvore.

**A correção não foi fazer uma chamar a outra**, porque nenhuma das duas
era a dona certa: um Krylov matrix-free não sabe nada sobre malha
tetraédrica nem sobre convecção-difusão 1D — é álgebra linear. Deixá-lo
numa base de FVM não-estruturado obrigaria qualquer chamador futuro a
herdar de algo sem relação só para alcançá-lo. Mesma lição do item 2.1, um
nível acima: **o lugar de uma coisa compartilhada se decide pelo que ela
é**, não por qual chamador precisou dela primeiro.

Vive agora em `engine/solver/include/aether/solver/KrylovSolvers.hpp`.

Duas políticas que a rotina compartilhada deliberadamente não assume
viraram parâmetros opcionais em vez de ficarem embutidas:

- `onResidual` — chamado a cada iteração com o resíduo relativo. É a
  observabilidade que o solver 1D precisa para seu `residualHistory()`
  (usado pelos testes de monotonicidade do GMRES reiniciado), sem que a
  rotina tome posição sobre o que fazer com o número.
- `brokeDown` — reportado em vez de dobrado no valor de retorno.
  "Convergiu", "esgotou as iterações" e "quebrou" são três desfechos
  diferentes, e chamadores diferentes se importam com desfechos
  diferentes. O solver 1D converte para o `-1` que sua API sempre
  devolveu; o não-estruturado ignora e checa o resíduo, como faria para CG.

**Saldo: 305 linhas removidas contra 28 adicionadas**, com a suíte em 13/13
— incluindo os testes que comparam a contagem de iterações do GMRES 1D com
valores literais, que passariam a falhar ao menor desvio numérico.

---

## 6. Continua sendo pesquisa, não dívida

Estes **não** são remendo — são trabalho genuinamente difícil, registrados
para não serem confundidos com descuido:

- ~~**Recuperação combinatória de facetas** (caso do quadrilátero
  coplanar)~~ — **FEITO em 2026-08-27**, ver abaixo: o bloqueio do cilindro
  fecha para o caso convexo.
- **Tetraedralização restrita preservando superfície importada, caso
  genuinamente não-convexo**: o que resta depois do item acima é
  especificamente uma obstrução tipo Schönhardt — contorno com estrutura
  reflexa severa o bastante para nenhuma troca de diagonal local nem
  subdivisão de Steiner limitada recuperar. Nunca construído nem encontrado
  neste código; citado como teto teórico, não como bug observado.
- **Validação contra literatura**: evitada por princípio (recordar tabelas de
  memória é risco real). Quando credibilidade externa importar, obter os dados
  da fonte real vira tarefa legítima.

**Atualização 2026-08-27: recuperação combinatória por troca de diagonal
implementada — o cilindro deixa de estar bloqueado, para geometria
convexa.**

`DelaunayTetrahedralization3D::recoverFacets()` ganhou uma primeira fase
combinatória (`tryFlipCoplanarQuadDiagonal()`), tentada antes da subdivisão
por ponto de Steiner: quando duas facetas pedidas e ausentes compartilham
uma aresta e descrevem, juntas, um quadrilátero coplano cuja *outra*
diagonal já está presente como duas faces de tetraedro com um ápice comum,
troca a diagonal no lugar — dois tetraedros removidos, dois outros
cobrindo exatamente o mesmo volume adicionados, nenhum ponto novo. A
convexidade do quadrilátero (precondição para a troca não sobrepor nem
deixar buraco) é verificada explicitamente via `orientation3D`, não
assumida a partir do sucesso de `makeTetrahedron()`.

**Medido, não suposto, em dois níveis:**

- **O caso adversarial (cubo de 8 pontos, sem nenhum ponto interior
  próximo)** ainda não fecha por completo — e agora sabe-se exatamente por
  quê: dos 6 quadriláculos com diagonal errada, a troca recupera 2; os
  outros 4 têm a aresta bloqueadora compartilhada por **três** tetraedros
  (um leque aberto passando por dois vértices do lado oposto do cubo), não
  dois — retriangular um leque de tamanho arbitrário é o caso geral
  genuinamente difícil que a própria classe já documentava. Volume e
  propriedade de Delaunay conferidos antes/depois: preservados
  exatamente.
- **O caso que de fato importa — geometria convexa real — fecha por
  completo.** Um cilindro aproximado por prisma poligonal (tampas
  triangularizadas em leque, faces laterais retangulares cada uma dividida
  à mão), a exata forma que produz esta ambiguidade, testado em três
  tamanhos (10, 24 e 60 lados, altura e raio deliberadamente não-redondos
  para descartar simetria acidental): **100% das facetas pedidas
  recuperadas em todos os três**, só pela troca combinatória, sem nenhum
  ponto de Steiner, com o volume exato do prisma preservado. Um cubo de 8
  pontos, sem estrutura local próxima, é próximo do pior caso para esta
  técnica; geometria realista tem estrutura local mais densa e raramente
  cai nele.

**O que isso mudou na frase "único bloqueio restante do cilindro"**: deixou
de ser verdade no nível da tetraedralização — medido, com um cilindro de
verdade, não um substituto escolhido para evitar o problema (como o
icosaedro que fechou o portão da Fase 3).

**Atualização 2026-08-28: o pipeline completo fecha também, não só a
tetraedralização isolada.** `mesh_flow_around_object` +
`freestream_boundary` + `UnstructuredCavitySolver3D` rodados de ponta a
ponta sobre um cilindro (prisma de 12 lados, 48 triângulos) **importado
via round-trip real por arquivo STL binário** (`save_stl_binary` +
`load_stl`, não mantido em memória) — o gap nomeado explicitamente na
atualização anterior. Medido, registrado em
`python/tests/test_pipeline.py::test_cylinder_from_stl_recovers_every_facet_and_conserves_mass`
(ctest, roda a cada build): **0 de 48 facetas perdidas** (recuperação
inteira via `tryFlipCoplanarQuadDiagonal`, nenhum ponto de Steiner),
volume esculpido exato (erro relativo 1,6e-15), escoamento passante com
desbalanço de massa de 0,0014% da vazão de entrada, sem NaN, velocidade
máxima limitada. Essa é a prova de integração ponta a ponta que faltava —
não mais "o próximo passo natural", já fechado.

---

## Ordem sugerida

Por dependência, não por tamanho:

1. ~~**1.1** (fluxo de saída)~~ — FEITO
2. ~~**2.1** (extrair base compartilhada)~~ — FEITO
3. ~~**1.2, 1.3 e 1.4**~~ — FEITO, cada um com sua medição isolada
4. ~~**5.2** (bindings)~~ — FEITO, e já pagou: encontrou 5.4 na primeira execução
5. ~~**2.3** (pressão da saída no estêncil)~~ — FEITO, e dissolveu o 5.4
6. ~~**4.3** (NaN em malha distorcida)~~ — RESOLVIDO (2026-08-25) para
   domínio fechado **e** com saída: raio espectral em jitter 0,65 foi de
   44,05 para 0,286, canal fecha em 1,2e-13. Sétima tentativa; a causa
   final foi a face de saída não receber a correção de Rhie-Chow que toda
   face interna já recebia
7. ~~**3.2** (solução manufaturada)~~ — FEITO: é ordem 2, a inferência estava certa
8. ~~**3.1** (esquema de alta ordem)~~ — FEITO: upwind medido em ordem 1, limitado em 2
9. ~~**4.1** (margem nos estruturados)~~ — FEITO; abriu o 4.4
9b. ~~**5.4**, ~~**5.3**~~, ~~**3.3**~~ (medido), ~~**4.2**~~ — FEITOS nesta rodada
10. ~~**4.4**~~ — FEITO: resolvido no 2D e no solver 1D; portado (opção, não
    default) para os seis 3D e medido em Re=400 — inconclusivo, decisão de
    manter `Central` como default está registrada e justificada
11. ~~**5.1**~~ — FEITO: o push já tinha sido feito, mas o CI falhava 100%
    das vezes (`engine/testing/` nunca versionado, por um `.gitignore`
    ambíguo) sem que ninguém pudesse ver por quê, já que o GitHub exige
    login até para ler o log em repositório público. Corrigido e
    **verificado rodando** — primeira vez que o CI passa neste projeto

O item **6** (tetraedralização restrita) ficou por último não por
prioridade, mas porque era o único que não se resolvia com disciplina — se
resolvia com pesquisa. Sua metade tratável (recuperação combinatória,
caso do quadrilátero coplano) fechou em 2026-08-27; o que resta —
contorno genuinamente não-convexo, tipo Schönhardt — continua sendo
pesquisa de verdade.
