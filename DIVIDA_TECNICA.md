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

### 2.2 Duas representações de malha que não se falam

`core::Mesh` (vértices + células) existe desde o Módulo 1 e é o que
`ScalarField`/`VectorField` usam. `TetrahedralMesh` foi criada na Fase 2.1 e
não é construída sobre ela.

**Por que importa**: os campos do engine (`ScalarField`, `VectorField`) não
podem ser usados com os solvers não-estruturados, então esses solvers
carregam `std::vector<double>` cru. Isso vai forçar conversões manuais em todo
ponto de integração — pós-processamento, persistência, visualização.

**O que fazer**: decidir qual é a representação canônica e fazer a outra
construir sobre ela. Não manter as duas.

### 2.3 Pressão da saída fica de fora do ajuste de mínimos quadrados [notada ao resolver 1.3]

`buildGradientStencils(useBoundaryValues)` é chamado com `false` no solver de
NS, então **nenhuma** face de contorno entra no estêncil de gradiente. Para
parede isso é certo e o motivo está no código: gradiente nulo afirma que a
*derivada normal* some, e o ajuste leria isso como "o valor na face é igual ao
da célula", enviesando também as componentes tangenciais. Mas a saída tem
pressão prescrita — é Dirichlet, exatamente como as faces quentes e frias da
placa, que a difusão inclui.

**Por que importa**: é a causa direta dos 21 de 177 estêncils deficientes do
canal (item 1.3). Incluir a saída completaria o posto de boa parte deles, e
essas células são justamente as que decidem quanto o escoamento leva embora.
Enquanto ficarem de fora, o recuo Green-Gauss — de primeira ordem — carrega
12% do domínio.

**Por que não foi feito junto**: é mudança de comportamento com medição
própria a fazer, e empilhá-la sobre três outras no mesmo passo tira a
atribuição de cada uma. É uma linha (`true` em vez de `false`) mais a medição
do perfil de saída antes e depois.

---

## 3. Acurácia declarada mas não atingida

### 3.1 Convecção upwind de primeira ordem [criada nesta sessão]

Escolhida deliberadamente por estabilidade (diferença central é
incondicionalmente instável acima de Re de célula 2, e a Fase 1 mediu a
cavidade estruturada em 16,7). Mas é primeira ordem, e **difusão numérica de
upwind cresce com o refinamento da malha em direções oblíquas ao escoamento**
— justamente o caso em tetraedros.

**Por que não dá para contornar**: qualquer resultado quantitativo de
escoamento (arrasto, sustentação, perda de carga) fica dominado pela difusão
do esquema, não pela física. Serve para topologia, não para número.

**O que fazer**: esquema limitado de alta ordem (TVD / linear-upwind com
limitador), que é o padrão exatamente por esse motivo.

### 3.2 FVM não-estruturado em ordem ~1, não 2

Medido na Fase 2.2: 0,91 → 0,98 na norma global; 1,06 → 1,63 excluindo os
cantos singulares. Segunda ordem existe na região suave, mas não foi
demonstrada num caso sem singularidade.

**O que fazer**: um caso de verificação com solução suave em todo o domínio
(solução manufaturada) para medir a ordem sem o teto da singularidade. Sem
isso, "é de segunda ordem" continua sendo inferência, não medição.

### 3.3 Interpolação de face sem correção de skewness

Implementada e medida como **nula** na malha de rede com jitter. Mas essa
malha é benigna: gerada de uma rede regular perturbada. Malha de geometria
real tem skewness muito pior.

**O que fazer**: repetir a medição numa malha de geometria real antes de
concluir que não importa. O resultado nulo é válido só para o caso testado.

---

## 4. Estabilidade sobre margem nula

### 4.1 Solvers estruturados em CFL exatamente 1,0000

`LidDrivenCavitySolver2D::stableTimeStep()` devolve o limite convectivo sem
fator de segurança, e a cavidade a Re=400 roda com Re de célula 16,7 — muito
acima do limite 2 da diferença central.

**Evidência**: durante a Fase 1, uma perturbação modesta no esquema levou esse
caso a NaN em ~500 passos. Ele estava à beira o tempo todo.

**Por que não dá para contornar**: qualquer mudança futura no esquema tem
chance de derrubá-lo, e o modo de falha é NaN — não um aviso.

**O que fazer**: fator de segurança em `stableTimeStep()` (como o solver
não-estruturado já tem) e/ou upwind na convecção 2D.

### 4.2 Convecção ainda explícita no solver não-estruturado [criada nesta sessão]

A difusão virou implícita e destravou o passo. A convecção não — e agora é ela
o limite.

**O que fazer**: se malhas finas com escoamento rápido passarem a importar,
tratamento implícito ou semi-implícito da convecção. Não urgente enquanto o
limite convectivo for folgado.

---

## 5. Verificação incompleta

### 5.1 O CI nunca rodou

O workflow existe e a metade local do portão está provada (regressão
deliberada reprovou). Mas nada foi enviado ao GitHub, então o CI nunca
executou.

**Por que importa**: um workflow que nunca rodou é uma hipótese, não uma rede
de proteção. Erros de YAML, de caminho e de ambiente só aparecem na primeira
execução.

**O que fazer**: um push. O resto já está pronto.

### 5.2 Solvers não-estruturados sem bindings Python [criada nesta sessão]

Todas as outras camadas do engine têm bindings — é o invariante arquitetural
"hybrid C++/Python" desde o Módulo 1. `UnstructuredDiffusionSolver`,
`UnstructuredCavitySolver3D` e `TetrahedralMesh` não têm.

**Por que não dá para contornar**: sem bindings, todo experimento com esses
solvers exige escrever e compilar C++, o que é exatamente o atrito que a
camada Python existe para eliminar. A investigação dos 13,2% teria sido
minutos em Python em vez de vários ciclos de build.

**O que fazer**: bindings, como toda outra camada tem.

### 5.3 Caso da cavidade não-estruturada em malha grosseira por custo

n=3 (177 células) e t=4 foram escolhidos para caber no tempo de suíte, não
pela física. A afirmação de topologia é válida nessa resolução; qualquer
afirmação quantitativa não seria.

**O que fazer**: com a convecção implícita (4.2) e os bindings (5.2), rodar a
convergência de malha de verdade fora da suíte de testes.

---

## 6. Continua sendo pesquisa, não dívida

Estes **não** são remendo — são trabalho genuinamente difícil, registrados
para não serem confundidos com descuido:

- **Tetraedralização restrita preservando superfície importada**: o Módulo 3
  cobre a fatia tratável; recuperação de facetas em geometria arbitrária
  esbarra na obstrução de Schönhardt. É o único bloqueio restante do cilindro.
- **Recuperação combinatória de facetas** (caso do quadrilátero coplanar).
- **Validação contra literatura**: evitada por princípio (recordar tabelas de
  memória é risco real). Quando credibilidade externa importar, obter os dados
  da fonte real vira tarefa legítima.

---

## Ordem sugerida

Por dependência, não por tamanho:

1. ~~**1.1** (fluxo de saída)~~ — FEITO
2. ~~**2.1** (extrair base compartilhada)~~ — FEITO
3. ~~**1.2, 1.3 e 1.4**~~ — FEITO, cada um com sua medição isolada
4. **5.2** (bindings) — destrava todo o resto da investigação
5. **2.3** (pressão da saída no estêncil) — uma linha mais a medição; com os
   bindings prontos vira um experimento em minutos em vez de um ciclo de build
6. **3.2** (solução manufaturada) — mede o que hoje é inferido
7. **3.1** (esquema de alta ordem) — só depois que 3.2 der uma régua confiável
8. **4.1** (margem nos estruturados) — independente, pode entrar quando quiser
9. **5.1** (push) — um comando

O item **6** (tetraedralização restrita) fica por último não por prioridade,
mas porque é o único que não se resolve com disciplina — se resolve com
pesquisa.
