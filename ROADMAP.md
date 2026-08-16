# Roadmap — Aether CFD Engine

Escrito em 2026-08-14, depois de fechar os 14 módulos do plano original.

**Princípio deste roadmap**: cada passo só começa depois do anterior estar
concluído *e testado*. Por isso cada fase abaixo declara explicitamente seu
**portão de conclusão** — a verificação objetiva que precisa passar antes da
fase seguinte começar. Sem portão cumprido, não se avança.

---

## Onde o projeto está

Os 14 módulos do plano original têm primeira etapa implementada e validada:
math core, geometria (STL/OBJ), malha, física incompressível, solvers
lineares (CG/PCG/Multigrid/GMRES/BiCGSTAB), turbulência (mixing-length,
k-ε, k-ω SST, LES, DES), pós-processamento, visualização GL 3.3, UI em modo
imediato, GPU (CUDA), persistência, IA (otimizador/analítico/surrogate/
conversacional), API REST e sistema de plugins.

11 suítes C++ passando (~7s), pacote Python com 58 símbolos.

## A lacuna estrutural central

**A malha não-estruturada não alimenta nenhum solver.** Verificado
diretamente no código: `DelaunayTetrahedralization3D`, a CDT, os pontos de
Steiner e os predicados robustos — trabalho sofisticado e validado — são
consumidos apenas pelo `MarchingCubes3D` (pós-processamento) e pelos
bindings Python. Nenhum solver os referencia. Os cinco solvers consomem
exclusivamente `StructuredGrid3D`, uma grade cartesiana uniforme.

Consequência prática: **o engine hoje só simula caixas.** Importa um STL,
tetraedraliza, e depois resolve num cubo. Duas metades bem construídas que
nunca se encontraram.

Essa é a diferença entre um engine didático muito bom e um engine de CFD
utilizável, e por isso é o eixo das Fases 2-3 abaixo.

---

## Fase 0 — Rede de proteção — PARCIAL, remoto suspenso (2026-08-14)

**Objetivo original**: CI que roda as 11 suítes a cada push.

**O que foi feito**: `.github/workflows/ci.yml` existe e está verificado até
onde dá sem GitHub. O build limpo com `-DAETHER_ENABLE_CUDA=OFF` (o ambiente
que um runner teria, e um caminho que nunca havia sido exercitado desde o
Módulo 10) roda com 0 erros, 10 suítes registradas — `aether_gpu_tests`
corretamente ausente — e 10 módulos de binding compilados.

**Metade do portão cumprida**: uma regressão numérica deliberada
(normalização do `checkerboardIndex` de `2.0*fieldRms` para `1.0*fieldRms`)
fez `aether_analysis_tests` falhar e o `ctest` sair com código 8 — ou seja,
o push *seria* reprovado. Revertida, volta a 100% e código 0. Isso prova que
a suíte detecta regressão sutil; é a parte que de fato importa.

**Metade suspensa, por decisão explícita**: o projeto não usa GitHub por ora
(13 commits nunca enviados; o repositório remoto responde 404). Sem push, o
CI não roda, então o portão "ser visto reprovando no CI" fica **em aberto,
não cumprido** — registrado assim de propósito, em vez de marcado como
concluído.

**Como retomar**: basta um `git push`; o workflow já está no repositório e
dispara sozinho. Enquanto isso, a proteção equivalente é local e manual —
rodar `ctest --test-dir build -C Release --no-tests=error` antes de commitar.
Um hook de pre-commit automatizaria isso sem GitHub nenhum, se a fricção de
lembrar virar problema.

---

## Fase 1 — PREMISSA REFUTADA (2026-08-14). Rhie-Chow não é necessário aqui

**O que a fase supunha**: que o `checkerboardIndex` de 0,020357 na cavidade
2D media o desacoplo par-ímpar de grade colocada, e que Rhie-Chow o
eliminaria. **Essa premissa estava errada, e foi escrita aqui sem ter sido
verificada.**

**Como foi refutada — estudo de refinamento de malha** (Re=100, tempo
simulado fixo em 12,0 para comparar o mesmo estado físico):

| n | h | residualRms | índice | ordem p |
|---|---|---|---|---|
| 16 | 0,0625 | 1,713e-03 | 0,028979 | — |
| 24 | 0,0417 | 8,930e-04 | 0,014144 | 1,61 |
| 32 | 0,0313 | 5,427e-04 | 0,008393 | 1,73 |
| 48 | 0,0208 | 2,571e-04 | 0,003914 | 1,84 |
| 64 | 0,0156 | 1,490e-04 | 0,002257 | 1,90 |

**Um checkerboard genuíno é um artefato de amplitude fixa na escala da
grade: não converge.** O que se mede aqui converge em O(h²) — ordem 1,90 e
subindo em direção a 2. Logo é *estrutura física resolvida* perto dos cantos
da tampa, ou seja, erro de discretização se comportando exatamente como
deveria. Não há patologia a curar.

**O que foi tentado e rejeitado**: uma projeção incremental com o termo
Rhie-Chow no lado direito do Poisson. Ela *parecia* funcionar — o índice caía
61%, de 0,020357 para 0,007961. Olhando numerador e denominador separados,
porém:

| | residualRms | fieldRms | índice | max\|p\| |
|---|---|---|---|---|
| original | 1,40e-03 | 3,45e-02 | 0,020357 | 0,148 |
| "melhorada" | 1,22e-01 | 7,68e+00 | 0,007961 | 28,6 |

O conteúdo absoluto de checkerboard ficou **87× pior**; o denominador cresceu
223×, e só por isso a razão caiu. A análise de amplificação explica: com o
preditor sem termo de pressão, o modo de pressão suave amplifica por
1/sin²(k·h/2) a cada passo — a pressão inflou de 0,15 para 28,6 e o caso
Re=400 divergiu para NaN em ~500 passos, tendo sido estável antes. Revertido.

**Duas lições que valem mais que a fase**:
1. **`checkerboardIndex` sozinho não serve como portão**: por ser normalizado,
   baixa quando a pressão infla. Um portão que o use precisa fixar também o
   numerador absoluto, ou a métrica é gameável — inclusive de boa-fé.
2. **Convergência de malha é o que distingue artefato de erro físico.** Foi o
   único teste que respondeu a pergunta, e devia ter vindo antes da
   implementação, não depois.

**O que ficou de real, e é valioso**: `maxFaceDivergence()`. A equação de
Poisson usa o Laplaciano **compacto**, mas a divergência e o gradiente usavam
o stencil **largo** — compor dois largos não reproduz o compacto. Medindo a
divergência pelos fluxos de face Rhie-Chow (a álgebra dá exatamente
`wide_div(u*) − dt·compact_lap(p)`, que é o que o Poisson zera):

| caso | div stencil largo | div por faces |
|---|---|---|
| 20×16 Re=100 | 2,652e-01 | 1,012e-12 |
| 20×20 Re=100 | 2,659e-01 | 1,156e-12 |
| 32×32 Re=100 | 1,806e-01 | 1,994e-13 |
| 24×24 Re=400 | 2,827e-01 | 7,253e-13 |

**O solver sempre conservou massa quase exatamente.** A "divergência ~0,2"
documentada nesta classe desde que foi escrita era propriedade do
diagnóstico, não do escoamento. Coberto por
`testLidDrivenCavityFaceDivergenceIsAtSolverTolerance`.

**Fragilidade pré-existente registrada de passagem**: a cavidade a Re=400 roda
com CFL = 1,0000 exato (sem margem) e Re de célula = 16,7 — muito acima do
limite 2 da diferença central para convecção. Não foi introduzido por esta
fase, mas é por isso que ela tolerou tão mal uma perturbação. Vale um fator
de segurança em `stableTimeStep()` ou upwind na convecção, se casos de Re mais
alto passarem a importar.

---

## Fase 2 — A ponte malha→solver (difusão primeiro)

**Objetivo**: um solver de difusão FVM que roda sobre malha de tetraedros
gerada pelo próprio `DelaunayTetrahedralization3D`.

**Por que aqui**: é o menor passo que atravessa a lacuna central. Difusão
é a física mais simples do projeto e **já tem validação analítica
independente** — a série de Fourier 2D que o `SteadyDiffusionSolver`
estruturado valida hoje. Resolver o mesmo problema, com a mesma solução
exata, em malha não-estruturada isola o que mudou (a discretização
geométrica) do que não mudou (a física). Mesma disciplina de sempre: prove
a máquina nova num problema cuja resposta já se conhece.

**O que precisa ser construído**: conectividade por faces sobre
`core::Mesh` (hoje é só vértices + células, sem vizinhança de face),
cálculo de área/normal/distância entre centroides por face, e tratamento
de não-ortogonalidade — que em tetraedros é a regra, não a exceção, e é
onde mora a dificuldade real desta fase.

**Portão de conclusão**:
1. O mesmo problema de placa (3 lados frios, 1 quente) resolvido em malha
   de tetraedros bate com a série de Fourier dentro de uma tolerância
   **medida antes de ser escolhida** (como foi feito no caso estruturado).
2. Refinar a malha reduz o erro na ordem esperada — convergência de malha
   demonstrada, não assumida.
3. As 11 suítes continuam passando.

### Estado em 2026-08-14: 2.1 e 2.2 CUMPRIDOS

**2.1 — `TetrahedralMesh` (conectividade de faces): CUMPRIDO.** Validado por
identidades exatas, não tolerâncias sobre física: fechamento de célula
(teorema da divergência discreto) em 2,944e-17; contagem 4×136 = 544 =
2×266 + 12 exata em inteiros; volume total 1,000000000000000; antissimetria
exata entre os dois lados de cada face interna.

Um bug real foi pego pelo teste escrito exatamente para isso: a tabela de
winding das quatro faces estava inteiramente invertida, com todas as normais
apontando para dentro. **Nem o fechamento nem a antissimetria pegam isso** —
um conjunto globalmente invertido ainda soma zero na célula e ainda nega
corretamente entre os lados. Só a checagem owner→neighbour pega.

**2.2 — `UnstructuredDiffusionSolver`: PARCIAL.** Item 1 do portão cumprido;
item 2 cumprido em "o erro cai", não em "na ordem esperada".

A primeira versão, só com a parte ortogonal da decomposição over-relaxed,
**estagnou** — e o experimento foi montado justamente para descobrir isso
antes de assumir:

| n | células | rmsErro (só ortogonal) | ordem |
|---|---|---|---|
| 4 | 415 | 2,879 | — |
| 6 | 1358 | 2,433 | 0,42 |
| 8 | 3184 | 2,365 | **0,10** |

Platô, não convergência: a não-ortogonalidade medida (~1,5) é propriedade da
*forma* da célula, então refinar não a reduz e o termo omitido vira piso de
erro. Com a correção não-ortogonal por correção diferida (gradientes
Green-Gauss, termo explícito no lado direito, matriz segue SPD):

| n | rmsErro (com correção) | ordem |
|---|---|---|
| 4 | 1,604 | — |
| 6 | 1,056 | 1,03 |
| 8 | **0,895** | 0,58 |

Erro ~2,6× menor e convergindo, mas com a ordem **caindo** (1,03 → 0,58) —
sinal de outro piso adiante. O suspeito principal era o gradiente
Green-Gauss, que só tem primeira ordem em malha distorcida. Trocado por
mínimos quadrados ponderados por 1/|d|² (matriz normal 3×3 pré-computada por
célula, com Green-Gauss de reserva para as raras células de estêncil
deficiente):

| n | rmsErro (mín. quadrados) | ordem |
|---|---|---|
| 4 | 0,999 | — |
| 6 | 0,691 | 0,91 |
| 8 | **0,522** | **0,98** |

Erro mais 1,7× menor e — o que mais importa — a ordem agora é **estável em
~0,95** em vez de degradar. Isso é o que diz que o erro restante é
discretização honesta, não termo negligenciado.

**Resumo das três versões medidas**, cada uma respondendo a pergunta aberta
pela anterior em vez de ser assumida:

| versão | n=4 | n=6 | n=8 | ordem |
|---|---|---|---|---|
| só ortogonal | 2,879 | 2,433 | 2,365 | 0,42 → 0,10 |
| + correção não-ortogonal | 1,604 | 1,056 | 0,895 | 1,03 → 0,58 |
| + gradiente mín. quadrados | 0,999 | 0,691 | **0,522** | 0,91 → **0,98** |

**A quarta versão foi um resultado nulo honesto.** Interpolação ponderada por
distância mais substituição da componente normal do gradiente médio pela
diferença compacta: 0,989 / 0,700 / 0,531 — dentro do ruído da terceira. A
interpolação de face **não** era o termo que limitava a ordem, ao contrário da
hipótese que motivou a mudança. Mantida por ser a formulação mais correta,
não por ter melhorado número nenhum.

**O que limitava era o próprio problema.** A placa é descontínua nos dois
cantos superiores, onde a borda quente encontra uma fria, e a solução exata
não tem gradiente limitado ali — **nenhum esquema converge na sua ordem
formal numa norma que inclua essas células**. Excluindo-as:

| n | rmsErro (todas) | ordem | sem os cantos | ordem |
|---|---|---|---|---|
| 4 | 0,989 | — | 0,901 | — |
| 6 | 0,700 | 0,85 | 0,585 | 1,06 |
| 8 | 0,531 | 0,96 | **0,367** | **1,63** |

Ordem 1,63 e subindo rumo a 2: **o esquema é de segunda ordem onde a solução
exata é suave**, que é exatamente a afirmação que um FVM deve poder fazer. A
versão estruturada deste mesmo teste evita a armadilha amostrando só pontos
interiores — eu caí nela primeiro e só saí medindo.

O laço externo de correção diferida foi descartado como causa antes disso:
sua mudança final mede 1e-5 a 1e-8, desprezível frente ao erro de
discretização (ele para no teto de varreduras só porque a tolerância pedida é
bem mais estrita que o necessário).

**Portão cumprido**: bate com a série de Fourier dentro de tolerância medida,
e a convergência de malha está demonstrada na região suave, com a taxa menor
da norma global explicada pela singularidade e não pelo esquema.

---

## Fase 3 — Navier-Stokes em malha não-estruturada

**Objetivo**: escoamento incompressível sobre geometria arbitrária
importada de STL/OBJ.

**Por que aqui**: só faz sentido depois da Fase 2 provar que o FVM
não-estruturado está correto na física simples, e depois da Fase 1 porque
malha não-estruturada colocada precisa de Rhie-Chow **mais** que a
cartesiana, não menos (a não-ortogonalidade agrava o desacoplo).

**Portão de conclusão**:
1. Cavidade com tampa deslizante resolvida em malha de tetraedros reproduz
   a mesma topologia de vórtice que a versão estruturada já valida
   (reversão de fluxo no fundo, divergência limitada).
2. Um caso com geometria genuinamente não-cartesiana (ex.: escoamento em
   torno de um cilindro, a partir de STL) roda estável e conserva massa.
3. As 11 suítes continuam passando.

**Marco**: cumprido este portão, o engine deixa de simular apenas caixas.

---

### Estado em 2026-08-14: item 1 do portão CUMPRIDO

`UnstructuredCavitySolver3D` — Navier-Stokes incompressível colocado sobre
`TetrahedralMesh`. A equação de Poisson da pressão **é** o operador validado
na Fase 2.2, reusado sem alteração; o que é novo é transporte de momento e a
projeção que o acopla à pressão.

Escolhas e por quê: armazenamento colocado (grade deslocada não tem análogo
natural em tetraedros — é por isso que códigos não-estruturados reais são
colocados e usam fluxos de face tipo Rhie-Chow); convecção **upwind**, porque
diferença central é incondicionalmente instável acima de Re de célula 2 e a
Fase 1 já mediu isso mordendo; passo explícito com fator de segurança 0,4,
justamente porque a Fase 1 achou a cavidade estruturada rodando em CFL
1,0000 sem margem.

**Topologia de vórtice reproduzida** (177 células, Re=10, t=4):
u médio no topo **+0,0439** (arrasto viscoso direto da tampa), no fundo
**−0,0099** (escoamento de retorno, que conservação de massa numa caixa
fechada torna obrigatório). Divergência por faces 3,3e-04, limitada.

**Um custo real medido, não estimado**: o passo explícito é limitado pela
*menor* célula, e uma tetraedralização Delaunay de rede com jitter sempre
produz slivers — então refinar encolhe o dt bem mais rápido do que adiciona
células. A n=4 com t=8 este único teste levou **6m57s**, contra ~21s para
todo o resto da suíte. Malha e tempo simulado foram fixados no que a
afirmação precisa (topologia é resolvível em malha grossa), não no que
pareceria impressionante. Suíte total hoje: 59s.

**O que falta para fechar a Fase 3**: item 2 do portão — um caso genuinamente
não-cartesiano (escoamento em torno de cilindro a partir de STL). E, para
tornar prático rodar malhas finas por tempos longos, **difusão implícita**:
o limite difusivo explícito escala com o quadrado da menor célula, o que é a
restrição real aqui, não a capacidade do esquema.

## Fase 4 — GPU para valer

**Objetivo**: (4.1) solve de CG inteiramente residente na GPU; (4.2)
preditor de momento do `StaggeredCavityBase3D` portado para CUDA.

**Por que aqui**: hoje o `PoissonOperatorCuda` faz cópia host→device,
lança o kernel e copia de volta a cada chamada — validado bit a bit, mas
paga PCIe por iteração, então não acelera nada de verdade. Ganho real exige
o laço de CG inteiro (produtos internos inclusive) residente na GPU. E isso
só passa a importar quando as malhas ficam grandes — ou seja, depois da
Fase 3. Portar o preditor de momento vem depois porque se beneficia da
infraestrutura de dados residentes do 4.1, e beneficia os seis fechamentos
3D de uma vez (foi exatamente por isso que a base compartilhada foi
extraída antes dos módulos 9-14).

**Portão de conclusão**:
1. Speedup **medido** contra a versão CPU no mesmo caso, não estimado.
2. Resultado numericamente equivalente ao da CPU dentro da tolerância do
   método iterativo (aqui não cabe exigir bit a bit: a ordem de somatório
   do produto interno paralelo difere legitimamente da sequencial — e essa
   diferença precisa ser medida e explicada, não tolerada por padrão).
3. As 11 suítes continuam passando.

---

## Fase 5 — Geometria CAD real

**Objetivo**: importar STEP/IGES.

**Por que aqui**: só tem valor depois que malha não-estruturada alimenta
solver (Fase 3) — antes disso, importar CAD sofisticado para simular uma
caixa não entrega nada.

**Decisão pendente, sua**: STEP/IGES exige um kernel CAD de verdade
(OpenCASCADE é o caminho realista). Seria a **segunda** dependência externa
pesada do projeto, depois do CUDA — e diferente do CUDA, que era
insubstituível, aqui é uma escolha entre depender de uma biblioteca grande
ou não suportar CAD. Vale decidir quando chegar a hora, não agora.

**Portão de conclusão**: um arquivo STEP real importa, tetraedraliza e
resolve, com o volume da geometria batendo com o do CAD de origem.

---

## Fase 6 — Física nova

Em ordem de dependência crescente:

- **6.1 Transferência de calor conjugada** (sólido + fluido acoplados) — o
  passo mais curto, porque difusão e Navier-Stokes já existem; falta o
  acoplamento na interface.
- **6.2 Escoamento compressível** — exige tratamento de choque e um
  esquema de fluxo próprio; mudança de regime, não extensão.
- **6.3 Multifásico** — o maior dos três (rastreamento de interface,
  tensão superficial).

**Portão de conclusão de cada um**: uma solução analítica ou uma
propriedade de conservação exata que o caso satisfaça — o mesmo padrão
usado em todo solver deste projeto, nunca uma tabela de literatura
recordada de memória.

---

## Fase 7 — Produto

Itens independentes entre si, ordenáveis conforme a necessidade:

- **UI**: janelas móveis/dockáveis, árvore de cena (Módulo 9).
- **Plugins**: pontos de extensão para solvers e fechamentos de turbulência
  (hoje só diagnósticos de campo); descoberta por diretório; sandboxing —
  este último só importa se algum dia rodar plugin de terceiro não confiável.
- **REST**: autenticação, mais rotas, documento OpenAPI (Módulo 13).
- **IA**: surrogate de campo inteiro (hoje é escalar→escalar); operadores
  neurais / PINNs, se a Fase 3 gerar dados suficientes para treinar.
- **AMR** (refinamento adaptativo de malha) — nomeado no Módulo 3 original
  e nunca construído; depende inteiramente da Fase 2/3, porque refinar
  adaptativamente uma grade cartesiana uniforme não faz sentido.

---

## Explicitamente adiado, não esquecido

- **DNS**: computacionalmente impraticável nas escalas de malha atuais.
  Reavaliar só depois da Fase 4 (GPU) e com malhas muito maiores.
- **Backend Vulkan** (Módulo 8): esforço grande e isolado; o pipeline GL 3.3
  atende tudo que o projeto renderiza hoje. Só se o gargalo virar
  renderização, o que não é o caso.
- **Recuperação combinatória de facetas na CDT** (caso do quadrilátero
  coplanar) e o caso Schönhardt-difícil: a CDT atual cobre a fatia tratável.
  Vira prioridade se a Fase 3 esbarrar nesses casos com geometria real —
  o que é plausível, e aí sobe de prioridade por evidência, não por
  antecipação.
- **Validação contra benchmarks de literatura** (Ghia et al. e afins): o
  projeto evitou isso por princípio, porque recordar tabelas de memória é
  risco real de imprecisão. Continua correto. Quando credibilidade externa
  importar, o caminho honesto é obter os dados da fonte real, não de
  memória — e nesse momento vira uma tarefa legítima, não um atalho.

---

## Nota sobre sequenciamento

As Fases 0 e 1 são pequenas e de baixo risco: dá para concluí-las
rapidamente e ganhar confiança no processo de portões.

A Fase 2 é onde o projeto muda de patamar, e é substancialmente maior que
tudo que veio antes — FVM não-estruturado com não-ortogonalidade é um
problema genuinamente difícil, não uma extensão do que já existe. Vale
tratar com o mesmo cuidado que o k-ε recebeu (que precisou de três bugs
reais até chegar num resultado fisicamente sensato) e reservar tempo de
depuração de verdade.
