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

**Difusão implícita: FEITA (2026-08-14).** O termo viscoso passou a ser
resolvido como `(V/dt + νL)u* = ...` — o mesmo Laplaciano com a diagonal
deslocada, ainda SPD, mesmo CG. Com isso o `stableTimeStep()` perde o limite
difusivo inteiro, que era o que escalava com o **quadrado** da menor célula.
Validação cruzada na mesma malha (n=4): explícito dá u topo +0,06686 / fundo
−0,01156; implícito dá +0,06850 / −0,01169 — dois tratamentos temporais
independentes chegando ao mesmo escoamento dentro de ~2%. O caso caiu de
~6m37s para segundos.

**O caso do cilindro tem DOIS bloqueios, não um** (levantado em 2026-08-14
verificando o código, não assumido):

1. **O solver não tem condição de saída.** O construtor recebe apenas
   `wallVelocity` e `buildFaces()` trata toda face de contorno como parede; o
   operador de Poisson não recebe contribuição de contorno e o fluxo de massa
   por faces de contorno é assumido zero. **Entrada já funciona** — é uma
   parede com velocidade prescrita não-nula. **Saída não existe**: precisa de
   pressão Dirichlet contribuindo para a diagonal do Poisson e para o lado
   direito, velocidade zero-gradiente, e — o ponto crítico — massa poder
   deixar o domínio, o que hoje é estruturalmente impossível. Contido,
   independente de malha, ~80 linhas.
2. **Tetraedralização restrita preservando a superfície do cilindro** — o
   Módulo 3 registra isso resolvido só na fatia tratável desde que foi escrito.

**Fazer o (1) primeiro**: é o menor, não depende do (2), e sem ele nenhuma
malha resolve o problema.

**O que falta para fechar a Fase 3**: item 2 do portão — um caso genuinamente
não-cartesiano (escoamento em torno de cilindro a partir de STL). E, para
tornar prático rodar malhas finas por tempos longos, **difusão implícita**:
o limite difusivo explícito escala com o quadrado da menor célula, o que é a
restrição real aqui, não a capacidade do esquema.

### Estado em 2026-08-20: pipeline geometria→malha→solver conectado pela primeira vez; dois bloqueios novos, precisos, ainda não corrigidos

**O bloqueio (1) registrado acima (sem condição de saída) já não existe** —
`UnstructuredCavitySolver3D` ganhou `is_outlet`/`outletPressure` em trabalho
posterior a esta entrada (ver DIVIDA_TECNICA.md 1.1, 2.3). Esta nota do
roadmap ficou desatualizada nesse ponto; corrigido aqui. O que **ainda**
falta do lado do solver: domínio com saída em malha muito distorcida segue
com um raio espectral instável (DIVIDA_TECNICA.md 4.3, correção restrita a
domínio fechado por ora).

**O bloqueio (2) (tetraedralização restrita preservando a superfície
importada) foi exercitado de ponta a ponta pela primeira vez**, via novo
`python/aether/pipeline.py::mesh_flow_around_object()`: importa/constrói uma
`TriangleMesh`, tetraedraliza pontos do objeto + cantos de uma caixa
delimitadora, recupera as facetas do objeto (`recover_facets`), esculpe o
interior do objeto (`remove_region`) e devolve um `TetrahedralMesh` pronto
para `UnstructuredCavitySolver3D` — o encadeamento que **não existia
nenhum código ligando** antes desta entrada (confirmado por grep: `load_stl`/
`load_obj` não eram usados em lugar nenhum fora dos próprios bindings).

Testado num icosaedro real (12 vértices, 20 triângulos, convexo, estanque) —
não a rede jitterada canônica que toda a suíte usa até aqui:
- 0 facetas não recuperadas.
- Volume esculpido bate com caixa-menos-objeto **essencialmente exato**
  (1061,111279737648 vs. 1061,1112797376475).
- Classificação de face por posição (`classify_boundary_face`) separa
  corretamente as 20 faces do objeto das 12 faces da caixa (2 por lado, como
  esperado de uma caixa cortada por um tetraedro por face).

**Dois bloqueios novos, medidos, não presentes na malha jitterada canônica
que toda a validação anterior usou**:

1. **Malha sem refino interior explode o solver** — não é o item 4.3
   (não-ortogonalidade máxima medida em 0,80, zero estêncil deficiente,
   ambos bem dentro do que já roda estável em outros lugares). É razão de
   volume entre células: **80x** entre a menor e a maior, porque só os
   vértices do objeto e os 8 cantos da caixa entram como pontos —
   praticamente toda a caixa vazia vira poucos tetraedros enormes. Instalado
   e medido diretamente: divergência por faces foi a **3,9e6** em 50 passos.
   Uma malha real de CFD nunca tetraedraliza só o contorno; sempre distribui
   pontos de fundo no volume, e é exatamente essa etapa que falta aqui.
2. **Refino por pontos de Steiner não é seguro perto de paredes
   recuperadas.** Tentativa: inserir o centroide de cada tetraedro grande
   como ponto de Steiner (`insert_steiner_point`, já validado para outros
   usos) até reduzir a razão de volume. Resultado medido: **5 das 20 facetas
   do objeto, antes recuperadas, desapareceram** depois do refino
   (`missing_facets()` confirmou). `insertSteinerPoint()` não tem noção
   nenhuma de "parede protegida" — sua retriangulação local pode atravessar
   uma faceta que `remove_region()` já havia respeitado. Isto é uma lacuna
   real e específica na ferramenta de refino, não um bug do pipeline novo.

**Próximo passo concreto para fechar o item 2 do portão**: refino
consciente de parede, de uma das duas formas — (a) distribuir os pontos de
fundo **antes** de `tetrahedralize()` (grade de fundo graduada por
distância ao objeto, mais pontos perto da superfície) em vez de refinar
depois via Steiner, evitando o problema por completo; ou (b) verificar
`missing_facets()` após cada inserção de Steiner e re-recuperar
imediatamente as que se perderam. (a) é mais simples de implementar e mais
parecido com o que geradores de malha reais fazem; (b) reaproveita mais
código já validado. Nenhum dos dois foi tentado ainda.

Suíte completa (12 suítes, incluindo o novo binding
`TriangleMesh::triangle()` que este trabalho precisou adicionar): 12/12.

### Estado em 2026-08-20 (continuação, mesmo dia): opção (a) implementada — melhora sete ordens de grandeza, não fecha ainda; dois achados novos, mais estreitos

`mesh_flow_around_object()` ganhou uma rede de fundo jitterada (mesma
convenção "camada de contorno exatamente no plano, interior com jitter" já
usada em `build_jittered_lattice()` — reaproveitada, não inventada) gerada
**antes** de `tetrahedralize()`, com espaçamento derivado automaticamente do
tamanho médio de triângulo da própria superfície importada, e uma zona de
exclusão esférica ao redor do objeto para não amontoar pontos de fundo perto
da geometria fina. Isto evita por construção o problema do refino por
Steiner (achado 2 da entrada acima): a rede de fundo nunca toca uma faceta
já recuperada, porque é adicionada antes de qualquer faceta existir.

**Medido no mesmo icosaedro de teste**: 50 → **2197 células**, razão de
volume ainda em 97x (quase igual a antes — ver por quê abaixo), mas
**velocidade máxima caiu de ~3,3e12 para ~44** e divergência por faces de
3,9e6 para 56,6 — sete ordens de grandeza de melhora, e a malha em si segue
exata (volume esculpido bate, 0 facetas perdidas). **Ainda instável**, não
fechado.

**Dois achados novos, mais estreitos que os da entrada anterior**:

1. **A razão de volume de 97x não está espalhada pela malha — está
   concentrada numa cauda pequena de células na interface objeto↔fundo.**
   Percentis medidos: 1%→0,059, 5%→0,17, mediana→0,46, 95%→0,88, 99%→1,11,
   mas o mínimo absoluto é 0,014 — bem abaixo do percentil 1%. A maior parte
   da malha está bem dimensionada; um punhado de células na costura entre a
   densidade esparsa do icosaedro (só 12 vértices) e a rede de fundo mais
   fina é que puxa a razão. **Não-ortogonalidade máxima subiu para 2,76**
   (contra 0,80 na malha sem refino, e acima de tudo que a quinta tentativa
   do item 4.3 validou — jitter até 0,75) — o sintoma agora parece
   genuinamente aparentado ao item 4.3, só que disparado pela transição de
   densidade superfície↔fundo, não por distorção de rede regular.
2. **`recover_facets()` não escala para uma superfície mais fina.** Testado
   com um icosaedro subdividido uma vez (42 vértices, 80 triângulos, ainda
   um objeto simples e convexo): a chamada não retornou em 60s, contra bem
   menos de 1s com 12 vértices. Consistente com a complexidade documentada
   na própria classe (`missing_facets`: O(facetas × contagem_de_tetraedros))
   — cara o bastante para travar antes mesmo de medir se a superfície mais
   fina resolveria o achado 1.

**O que isto muda para o próximo passo**: aumentar a resolução do objeto
(que resolveria a transição de densidade do achado 1) esbarra imediatamente
no achado 2. As duas lacunas precisam ser atacadas juntas, não em sequência
— ou a densidade de fundo perto da zona de exclusão precisa se aproximar
suavemente da densidade da própria superfície (grade graduada de verdade,
não a exclusão binária atual) para que o objeto não precise ficar mais fino
só para casar com o fundo. Nenhuma das duas tentada ainda.

Suíte C++ inalterada nesta rodada (só `pipeline.py`, Python puro) — não
re-executada por não haver mudança de binding.

### Estado em 2026-08-20 (terceira rodada): **item 2 do portão da Fase 3 CUMPRIDO**; e uma correção a esta própria seção

**Primeiro, a correção — o achado 2 da entrada anterior estava errado.**
Aquela entrada afirmou que `recover_facets()` não escalava, com base em um
caso que não retornava em 60s. Perfilando as fases separadamente, em vez de
atribuir o tempo à última chamada feita: `missing_facets` = **0,00s**,
`recover_facets` = **0,00s**, `tetrahedralize` = **131s**. A recuperação de
facetas nunca foi o gargalo. O muro é `tetrahedralize()`, e ele é
**O(N²)** — medido, não estimado: dobrar N multiplica o tempo por 4,29 →
4,21 → 4,03 (convergindo a 4,0) e t/N² fica constante em 17,6 → 20,0.
Mesma família de erro que o item 4.3 já registrou duas vezes (atribuir
causa antes de isolar), agora numa terceira forma: atribuir custo à etapa
mais próxima do sintoma.

**Filtro de ponto flutuante nos predicados exatos: 29x mais rápido, sem
mudar nenhuma resposta.** A causa do O(N²) é `findCavity()` testar a
circumsfera de *todo* tetraedro por ponto inserido, e cada teste ia direto
para aritmética BigInt exata sem filtro. `RobustPredicates` ganhou o
filtro padrão da literatura (Shewchuk; é o que CGAL/TetGen/Qhull fazem, e
o que o próprio comentário do header já citava): avalia o determinante em
`double` junto com um limite rigoroso do próprio erro de arredondamento, e
só cai no caminho exato quando |det| não excede esse limite. Constantes
deliberadamente várias vezes mais conservadoras que o mínimo teórico —
superestimar só manda casos a mais para o caminho exato (mais lento,
nunca errado), subestimar devolveria sinal errado.

Medido: N=1600 caiu de **51,3s para 1,75s**. Continua O(N²) (o filtro
mexe na constante, não na complexidade), mas a constante caiu 29x — o caso
de 4874 pontos que nunca terminava agora roda. Suíte C++ inteira também
ficou mais rápida (58s → 44s).

A afirmação "não muda nenhuma resposta" é verificada, não assumida:
`orientation3DExact`/`inSphere3DExact` expõem o caminho não-filtrado e
`mesh_tests` compara os dois em **120 mil** avaliações — aleatórias mais
famílias deliberadamente degeneradas (quatro pontos exatamente coplanares,
cinco exatamente co-esféricos, slivers) — exigindo sinais idênticos. 35 mil
delas deram sinal exatamente zero, então a metade adversarial do espaço de
entrada foi de fato exercitada, não só a fácil.

**Gradação da rede de fundo: testada, medida, e mantida desligada por
padrão.** A hipótese era que a costura de densidade objeto↔fundo explicava
a instabilidade. A/B direto na mesma geometria e espaçamento: graduada
(3 níveis) dá 886 células de 168 pontos, razão de volume 965 e
não-ortogonalidade 2,52; uniforme dá 3493 células de 722 pontos, razão
**88** e não-ortogonalidade **2,06**. Graduada é mais barata e
mensuravelmente **pior** nas duas métricas de qualidade, porque cada
limiar de decimação cria seu próprio salto abrupto de densidade — trocou
uma costura por três. `max_level` default passou a 0 (sem gradação); o
código fica, documentado, porque engrossar o campo distante é a ideia
certa — só precisa de variação suave de tamanho (refino de Delaunay
contra uma função de tamanho), não decimação de rede.

**E o bloqueio de verdade não era nada disso: era o arnês de teste.** O
teste de repouso e o raio espectral — os dois diagnósticos que este
projeto aplica a todo solver — deram, nas duas malhas: repouso
**exatamente 0,000e+00** e raio espectral **0,245 / 0,279**, ambos bem
abaixo de 1. As malhas eram estáveis o tempo todo. O que estava errado era
o meu script: eu prescrevia velocidade em **x** na face `x_min`, cuja
normal é x — ou seja, soprando massa para dentro de uma caixa **selada**,
que não é um caso difícil e sim um caso insolúvel. O solver reportou isso
desde a primeira execução via `net_boundary_flux()` = −26 (deveria ser 0),
e eu não estava lendo o diagnóstico que já existia.

Com a tampa tangencial (função nova `driving_wall_velocity()`, que
**recusa** uma direção não-tangencial em vez de corrigi-la em silêncio),
no icosaedro, 3493 células, 400 passos:

| quantidade | medido |
|---|---|
| erro relativo de volume da malha | **1,5e-15** |
| fluxo líquido de contorno (caixa selada) | **0,000e+00** |
| divergência por faces | **6,1e-09** |
| velocidade máxima (tampa a 1,0) | 0,463 |
| u médio no topo | **+0,02992** (segue a tampa) |
| u médio no fundo | **−0,00312** (retorno, conservação de massa) |

Isto é **o item 2 do portão da Fase 3**: um caso genuinamente
não-cartesiano, com geometria importada e um objeto esculpido dentro do
domínio, rodando estável, conservando massa exatamente e reproduzindo a
topologia de vórtice primário que todo solver de cavidade deste projeto
valida. Suíte: 12/12.

**O que fica em aberto**, agora com prioridade medida em vez de suposta:
`tetrahedralize()` continua O(N²) mesmo 29x mais rápido, o que limita
malhas reais (4874 pontos ≈ 131s antes do filtro, ~4,5s depois; 10⁵ pontos
segue impraticável). A correção conhecida é localização de ponto +
propagação por adjacência a partir do tetraedro que contém o ponto, em vez
de varrer todos — troca O(N) por O(tamanho da cavidade) por inserção. Essa
é a próxima alavanca real, e agora é a maior.

### Estado em 2026-08-29: **a alavanca foi puxada** — localização de ponto + adjacência, O(N²) → quase-linear

Com as Fases 2 e 3 fechadas (item acima) e a linha da beta encerrada, esta
era a alavanca técnica mais citada neste documento como "a maior
restante" — puxada agora, escopo estreito de propósito: só o laço
principal de `tetrahedralize()` (que sempre insere pontos garantidamente
interiores ao super-tetraedro) foi reescrito. `insertSteinerPoint()`
(usado por `recoverFacets()`/refino) continua com o `findCavity()` de
varredura completa, intocado — não é o caminho medido como gargalo (o
profiling da entrada anterior já tinha achado `tetrahedralize`=131s contra
`recover_facets`=0,00s no mesmo caso).

**O mecanismo**, implementado em `DelaunayTetrahedralization3D.cpp`, sem
tocar o `.hpp` público nem o struct `Tetrahedron` nem os bindings Python
(confirmado antes por exploração: nenhum consumidor externo constrói
`Tetrahedron` por agregado, só lê `.vertices`):

- Uma lista de trabalho local à função (`WorkTet`: tetraedro + 4 vizinhos
  + vivo/morto) substitui a varredura de `tetrahedra_` inteira. Tetraedros
  mortos são marcados (`alive=false`), nunca fisicamente removidos —
  elimina todo o problema de reindexação por elemento que "moveu".
- Um mapa de faces expostas (`std::map` por `sortedFace()`, não
  `unordered_map` — evita escrever um hash só pra isso, e seu tamanho
  acompanha o perímetro da região exposta, não N) liga/expõe cada face em
  O(log k) amortizado, no mesmo padrão "primeiro toque expõe, segundo
  toque resolve" que `TetrahedralMesh::fromCells()` já usa em lote —
  aplicado aqui de forma incremental.
- **Localização de ponto**: uma caminhada por adjacência a partir do
  último tetraedro criado, usando o mesmo teste de substituição de
  `pointInTetrahedron()`/`isFaceVisible()` já existente. Nunca revisita um
  tetraedro na mesma chamada (marcado por geração, não por um vetor limpo
  a cada ponto — isso por si só reintroduziria custo O(N) por ponto) — o
  que torna a caminhada autolimitada a `work.size()` passos e fecha de
  vez a falha clássica de uma caminhada "primeira face violada" ciclar
  para sempre.
- **Expansão da cavidade**: BFS por adjacência a partir do tetraedro
  localizado (garantidamente "ruim" — o interior de um tetraedro está
  sempre dentro da própria circunsfera, já que seus 4 vértices estão
  sobre ela e uma bola é convexa), em vez de testar `inCircumsphere` em
  todos. O(tamanho da cavidade), não O(N).
- Rede de segurança: se a caminhada não convergir (não esperado para um
  ponto interior ao super-tetraedro, mas corretude nunca deve depender
  disso), cai para uma varredura exaustiva só daquele ponto — perde
  velocidade, nunca corretude.

**`tetrahedralizeReference()` — o algoritmo antigo, preservado como
oráculo**, não deletado: mesmo padrão que `orientation3DExact()`/
`inSphere3DExact()` já usam ao lado das versões filtradas em
`RobustPredicates.hpp`. Isso deu um teste de validação cruzada direto
(`testFastTetrahedralizeMatchesReferenceImplementation`, novo em
`mesh_tests.cpp`) que roda os dois algoritmos sobre os mesmos pontos, na
mesma ordem, e exige o **mesmo conjunto exato** de tetraedros — não só a
mesma contagem ou o mesmo volume — em todo caso geometricamente difícil
que esta classe já tinha um teste dedicado (o sliver de 7 pontos que
expôs o bug do super-tetraedro pequeno demais, o cubo esparso de 8
pontos, o octaedro oco duplo, o prisma denso de 24 lados) mais duas redes
jitteradas novas (4×4×4 e 6×6×6). Os dois algoritmos concordam
exatamente em todos os 8 casos.

**Medido, não estimado** (mesma construção de rede jitterada de
`buildCubeLatticeTetrahedralization` em `solver_tests.cpp`, cronometrando
`tetrahedralize()` isolado contra `tetrahedralizeReference()` isolado, no
mesmo processo):

| n | pontos | tetraedros | rápido | oráculo O(N²) | razão |
|---|---|---|---|---|---|
| 4 | 125 | 415 | 0,004s | 0,021s | 5,3x |
| 6 | 343 | 1358 | 0,015s | 0,255s | 17,0x |
| 8 | 729 | 3184 | 0,044s | 1,616s | 36,7x |
| 10 | 1331 | 6204 | 0,089s | 7,191s | 80,8x |
| 12 | 2197 | 10667 | 0,173s | — | — |
| 14 | 3375 | 16876 | 0,311s | — | — |

(n=12/14 não cronometrados contra o oráculo — extrapolando o crescimento
quadrático já medido, seriam da ordem de minutos, por nenhum ganho de
informação novo.) A razão de tempo do caminho rápido entre linhas
consecutivas (3,75x / 2,93x / 2,02x / 1,94x / 1,80x) acompanha de perto a
razão de **pontos** (2,74x / 2,13x / 1,83x / 1,65x / 1,54x), não o
quadrado dela — confirma escala quase-linear, não apenas "mais rápida
numa constante". A vantagem cresce com N por construção (um lado é O(N²),
o outro não), então o ganho real em malhas de 10⁴-10⁵ pontos — o alvo que
motivou esta fase — é maior ainda que os 80x já medidos em N=1331.

**Honestidade sobre o que isso muda hoje, não só o que muda em teoria**:
a suíte `ctest` completa continua em ~156s, praticamente inalterada
(`aether_solver_tests` foi de 122,74s pra 117,58s — dentro do ruído de
variação de máquina). Verificado por quê: `solver_tests.cpp` nunca chama
`cubeLatticeMesh()` com n maior que 8 (grep confirma: só 3, 4, 6 e 8
aparecem nos call sites reais), e nessa faixa a economia absoluta é de
segundos, não minutos — o ganho medido acima é real e vale para os
tamanhos que interessam (o problema declarado no início desta seção,
"10⁵ pontos é impraticável"), mas essa suíte específica ainda não exercita
malhas grandes o bastante para o número mudar visivelmente no agregado.
O benefício é para a frente: é exatamente o que faltava pra Fase 4 (GPU)
valer a pena em malhas de tamanho real, e pra física nova (Fase 6) sobre
geometria CAD importada não esbarrar no mesmo teto.

Suíte: **13/13**.

### Estado em 2026-08-25 (rodada autônoma): sete itens fechados além das cinco fases

Trabalho feito sem supervisão, com a instrução de pular o que exigisse
decisão. Resumo do que entrou e por quê — cada um tem seu próprio commit
com a medição que o sustenta.

**Capacidades novas**

- **Persistência de malha tetraédrica.** Era a lacuna que tornava um
  checkpoint não-estruturado inutilizável sozinho: os campos podiam ser
  salvos, a malha contra a qual eles são indexados não. Guardado como
  conectividade e não como nuvem de pontos, porque o resultado de Delaunay
  com empates co-esféricos não é único — re-tetraedralizar no load poderia
  produzir células diferentes e descasar os campos. Exigiu
  `TetrahedralMesh::fromCells()`, que era lacuna e não escolha: até então o
  único caminho para uma malha era o gerador Delaunay.
- **Transporte de temperatura** (Fase 6.1), passivo e Boussinesq, como dois
  ajustes separáveis. Passivo vem primeiro porque é exatamente checável
  (princípio do máximo); empuxo não tem checagem igualmente gratuita.
- **Escoamento passante no pipeline** (`freestream_boundary`) — entrada e
  saída, desbloqueado pela Fase C. É o caso de engenharia mais comum, e a
  malha não mudou: só as condições de contorno.

**Investigações que fecharam perguntas em aberto**

- **Item 3.3, metade da atribuição.** Cisalhamento afim isola o que a
  gradação não isolava: com skewness segurada em ~1,0, a não-ortogonalidade
  sobe de 1,54 a 9,43 e a ordem não sai de ~2,1. **Não-ortogonalidade não é
  a causa** da perda de ordem. Não conclui que seja skewness — isso seria
  eliminação por sobra.
- **Jitter 0,95 caracterizado.** É **não-linear** (limiar nítido entre
  entrada 0,01 e 0,05), logo categoricamente distinto da instabilidade
  linear que o item 4.3 fechou. Não é passo de tempo (÷20 não salva) nem
  viscosidade (×30 não salva).

**Dívida e usabilidade**

- **Item 5.5 fechado**: GMRES desduplicado para `KrylovSolvers.hpp`, 305
  linhas a menos. A correção não foi uma classe chamar a outra — nenhuma
  das duas era a dona certa de um Krylov matrix-free.
- **`docs/getting-started.md` e `examples/flow_around_object.py`**: a
  primeira documentação voltada a usar em vez de a continuar
  desenvolvendo, incluindo uma seção do que o motor **não** faz.

**O que ficou deliberadamente de fora**

- **STEP/IGES** — depende da decisão sobre OpenCASCADE, que é sua.
- **Tetraedralização O(N²)** — é a maior alavanca técnica restante
  (localização de ponto + propagação por adjacência), mas é refatoração
  profunda na classe com o histórico de bugs mais sutil do projeto, e um
  resultado parcial precisaria de avaliação. Documentada, não tentada.
- **Wheel/`pip install`**, **AMR**, **GPU residente** — escopo próprio.

### Estado em 2026-08-25: Fase C concluída — o item 4.3 fecha na sétima tentativa

**A correção é pequena e o diagnóstico é que era caro.** Faltava aplicar,
numa face de saída, exatamente o termo de correção de Rhie-Chow que o laço
de faces internas já aplicava — a projeção *impõe* um fluxo de saída
construído com o gradiente de face misturado, mas o resíduo que o GMRES
mirava usava a extrapolação simples `velocity·A`. O GMRES zerava uma
quantidade enquanto `maxFaceDivergence()` media outra.

Duas medições chegaram lá, ambas derrubando uma hipótese minha:

1. Instrumentei os dois caminhos de quebra do GMRES: **nenhum dispara**. O
   método converge. Isso derrubou a afirmação (que eu próprio havia escrito
   no código) de que sobrava uma direção quase-nula.
2. A hipótese seguinte — clamp em células de estêncil deficiente — caiu
   igualmente: dois domínios com **zero** estêncils deficientes mantinham o
   mesmo piso de resíduo.

Restando "o GMRES converge mas o resíduo medido não cai", só sobra uma
explicação possível: o operador resolvido e a quantidade medida não são o
mesmo operador.

**Medido, no canal com saída**: divergência por faces de 3,0e-02 para
**1,1e-11**; fluxo líquido de contorno de 1,2e-04 para **1,2e-13**. Nos
testes Python o canal fecha em 1,24e-13 com 2, 4 e 16 correctores — deixou
de depender da contagem de correctores, que é o sinal de que a correção
parou de compensar um erro.

**Raio espectral em torno do repouso, malha com saída**: em jitter 0,65 foi
de **44,05 para 0,286**; toda malha testada até jitter 0,95 fica abaixo de
1. A instabilidade linear que abriu o item 4.3 acabou.

O termo é proporcional a `fluxDt`, então o caminho do lado direito do
Poisson (`fluxDt = 0`) fica inalterado **por construção** — verificado:
6,475e-14, o mesmo de antes.

**Um erro de protocolo meu, corrigido no caminho**: a primeira varredura
mediu o raio espectral com a entrada ligada a 1,0 e deu valores absurdos
(5.000 a 42.000) enquanto a marcha real rodava 400 passos com divergência
1e-06. Com entrada não-nula o repouso não é ponto fixo, então a razão de
normas mede o forçamento, não o operador. Foi a contradição entre as duas
medidas que expôs isso.

**O que fica aberto, e é outro fenômeno**: jitter 0,95 ainda morre no passo
41 numa marcha dirigida, apesar de ρ = 0,58 em repouso — e não é
não-ortogonalidade (0,85 roda 400 passos com 10,04 enquanto 0,95 morre com
7,40). Registrado como distinto, não como resquício. O teste do guarda de
finitude subiu de jitter 0,65 para 0,95, terceira vez que sobe porque o
solver melhorou.

Suíte: **13/13**.

### Estado em 2026-08-25: Fase D — turbulência no solver não-estruturado (comprimento de mistura)

Primeiro fechamento de turbulência deste projeto a rodar sobre algo que
não é uma grade estruturada.

**Comprimento de mistura primeiro, pelo mesmo motivo de sempre.** Prandtl é
algébrico — não carrega equações de transporte próprias — então exercita a
única peça genuinamente nova (viscosidade efetiva variando no espaço dentro
do operador de momento) sem trazer junto os campos k/epsilon acoplados e
não-lineares, cuja numérica precisou de três correções de bug reais quando
foi construída para o canal 1D. Fechamentos de duas equações nesta malha
são o próximo passo natural, não uma lacuna esquecida.

**Feito em duas etapas, e a primeira era a que podia quebrar tudo.** Antes
de acoplar modelo nenhum, o termo viscoso passou a aceitar viscosidade *por
face* — e o caminho laminar foi mantido **bit-idêntico**, não apenas
próximo: `nu*(c1+c2+c3)` e `nu*c1+nu*c2+nu*c3` diferem em ponto flutuante,
então o caso uniforme continua sendo calculado pela expressão antiga, atrás
de um `if`. Verificado: cavidade 6,484e-05, canal 6,475e-14, erros da
solução manufaturada idênticos aos anteriores.

**O operador continua SPD**, o que importa porque é o que permite manter o
mesmo Gradiente Conjugado: os dois lados de uma face usam a *mesma*
viscosidade de face, então o par fora da diagonal permanece igual.

**Viscosidade molecular pura na face de parede**, não o nu_t da célula
vizinha — o comprimento de mistura é exatamente zero numa parede por
definição, e este projeto já pagou por errar isso uma vez:
`MixingLengthChannelFlowSolver1D` espelhava o nu_t da célula adjacente
através da parede e sua velocidade de atrito saiu 38% fora do balanço de
momento exato.

**Validação, na disciplina de sempre — propriedades do modelo, nunca tabela
de literatura.** O caso em repouso é uma afirmação incomumente forte aqui
justamente por o fechamento ser algébrico: sem equação de transporte e sem
memória, campo de velocidade exatamente zero dá taxa de deformação
exatamente zero e portanto nu_t **exatamente zero** — não pequeno, zero.
Um fechamento de duas equações não poderia ser checado assim, porque k e
epsilon difundem em direção aos seus valores de parede mesmo em repouso.

**E uma checagem que nenhuma das outras implica**: que o fechamento de fato
*alcança* o operador de momento. Um nu_t calculado todo passo e nunca
consultado passaria por nu_t máximo positivo, divergência limitada e
repouso perfeito, mudando nada. A única forma de saber é rodar o mesmo caso
laminar e ver os campos divergirem — medido, `|u_turb − u_lam|` máximo
**4,7e-03** contra nu_t máximo 6,8e-03 sobre nu 0,1, ou seja, os ~7% que a
viscosidade adicional prevê.

Rodando o pipeline inteiro (geometria importada → malha → solver
turbulento), 911 células: nu_t de **1,8e-06** junto às paredes a
**2,8e-02** no campo distante, divergência 1,3e-10, fluxo líquido
exatamente 0.

Suíte: **13/13**.

### Estado em 2026-08-25 (mesma sessão): Fase E — atrito de instalação removido; wheel segue em aberto

Dois atritos, e o maior deles não era o build.

**O pacote acha os próprios módulos compilados.** Até aqui, todo ponto de
entrada — script, teste, sessão interativa — precisava apontar `PYTHONPATH`
para o diretório de saída do build à mão, e esse diretório *muda com o
gerador*: multi-config (Visual Studio) acrescenta um nível `Release/` ou
`Debug/` que single-config não tem. `aether._ensure_extensions_importable()`
procura nos três lugares que este repositório de fato usa. Um módulo já
importável é deixado em paz e a função retorna na hora: `PYTHONPATH`
explícito, instalação em venv ou um wheel futuro significam que o chamador
já disse onde as extensões estão, e sobrepor isso seria pior que o atrito
removido.

**`build.py` — um comando do clone ao engine funcionando**: configura,
compila em Release e roda as 13 suítes. Existe em vez de virar uma linha do
README por causa de uma armadilha específica: `--config Release` vai no
comando de *build*, não no de configuração, porque gerador multi-config
escolhe a configuração na hora de compilar e ignora `CMAKE_BUILD_TYPE`
silenciosamente. Errar isso produz um erro que não nomeia a própria causa.
Roda `ctest` com `--no-tests=error`, não o padrão: uma suíte que parasse de
registrar seus testes reportaria sucesso, que é o único modo de falha que
uma execução de testes não pode ter.

**O que a Fase E deliberadamente NÃO fecha**: não há wheel nem
`pip install aether`. Isso exige as extensões C++ compiladas por plataforma
e por versão de Python — scikit-build-core na cadeia de build, ou CI
produzindo wheels por alvo. É trabalho real e separado; um script de build
não o cobre, e dizer que cobre seria o tipo de exagero que o resto desta
documentação evita. Registrado como lacuna, não como concluído.

### Estado em 2026-08-25: **Fase B concluída** — o resultado não-estruturado saiu do terminal

O item era "visualização do resultado da malha não-estruturada", e a
primeira decisão foi *não* construir mais um modo de viewer.

**Por que exportar em vez de renderizar, e a assimetria que decide.** Todo
solver estruturado deste projeto tem seu modo no `unified_viewer`, e
estavam certos: um campo estruturado é um array denso, então um mapa de
calor ou uma malha de setas o mostra. Um resultado tetraédrico
não-estruturado precisa de corte, limiar, semeadura de linhas de corrente
e iso-superfície sobre malha irregular — que é exatamente o conjunto de
operações que um pós-processador geral existe para oferecer, e que ParaView
e VisIt já fazem melhor do que este repositório faria dentro de um app
Win32/GL. Havia também um obstáculo de arquitetura real: o pipeline é
Python (por decisão registrada — orquestração fora do núcleo numérico) e o
viewer é C++, então um modo de viewer duplicaria a geração de malha, que é
precisamente o custo que o item 2.1 já mediu.

**`engine/postprocessing/VtkWriter`** (Módulo 7, onde saída de campo
pertence) escreve `TetrahedralMesh` mais campos por célula em VTK legado
ASCII. Formato legado e não XML `.vtu` de propósito: layout de texto
estável e completamente especificado, sem escritor XML, sem base64, sem
biblioteca de compressão e sem negociação de versão de esquema — e todo
programa que lê VTK lê esse. Campos são **por célula**, que é onde as
incógnitas de volumes finitos deste engine de fato vivem; escrevê-los por
vértice exigiria inventar uma interpolação e depois apresentá-la como se
fosse a solução.

Validado sem depender de biblioteca externa: o teste **relê o próprio
arquivo** e confere cada contagem do cabeçalho, cada coordenada, cada
índice de conectividade e cada valor de campo contra a malha que o
escritor recebeu — igualdade **exata**, não tolerância, já que 17 dígitos
significativos round-trip um double bit a bit. Também checa que entrada
inconsistente é recusada (campo de tamanho errado, nome com espaço — este
último produziria um arquivo que outra ferramenta *misparseia* em vez de
rejeitar, já que o formato separa tokens por espaço).

`export_result_vtk()` fecha o caminho no pipeline: velocidade, pressão e
`speed` derivado num arquivo que abre direto no ParaView. O teste do
pipeline confere as **911 pressões** contra o solver uma a uma, e exige que
o campo não seja trivialmente nulo — a primeira versão dessa asserção
olhava só a célula 0, que é a célula *fixada* e vale exatamente 0, logo
comparava 0 com 0 e não podia falhar.

Suíte: **13/13 em 55s**.

### Estado em 2026-08-20 (quarta rodada): checagens de pré-voo e relato de convergência — **Fase A do roadmap de usabilidade concluída**

Os dois itens que faltavam da Fase A existiam no plano antes do erro da
rodada anterior; depois dele ficou óbvio *qual* checagem faltava, então
ambos foram construídos a partir do que a falha real ensinou, não do que
parecia razoável no papel.

**`check_closed_domain_conservation()` — a checagem cuja ausência custou
mais caro.** Integra a velocidade prescrita sobre toda face de contorno e
diz se o fluxo líquido é zero, como um domínio fechado exige. Custa uma
passada de geometria pura, sem marchar um único passo. É exatamente o que
teria pego, em segundos, o erro que consumiu boa parte da rodada anterior.
Escrever o teste rendeu um diagnóstico mais afiado que o líquido: numa
caixa selada cada parede deve ser impermeável **face a face**, e uma tampa
tangencial dá `max_face_flux` exatamente 0 — então o relatório expõe esse
número também, e é ele que aponta *qual* parede está vazando quando os
totais por acaso se cancelam.

**`measure_mesh_stability()` — mede em vez de inventar um limiar.**
DIVIDA_TECNICA.md 4.3 afirma que não há critério a-priori honesto de
qualidade de malha (não-ortogonalidade 2,24 roda, 2,07 diverge). Então esta
sonda não inventa limiar nenhum sobre métrica-proxy: roda os dois
diagnósticos que de fato decidem — repouso exato e raio espectral por
potência iterada em torno dele — em dezenas de passos, contra os milhares
de uma marcha real. As métricas de qualidade (não-ortogonalidade, estêncil
deficiente, razão de volume) entram no relatório como contexto para quem
lê, e **deliberadamente não** entram no veredito `is_stable`, porque o
item 4.3 já mediu que elas não separam malha que roda de malha que não.

**`run_to_steady_state()` / `RunReport` — o relato que faltava.** Marcha
até o campo parar de mudar (tolerância *relativa*, pela mesma razão que o
item 5.4 registra: um limiar absoluto é frouxo demais ou apertado demais
conforme a escala do campo, que o chamador não deveria precisar saber de
antemão), detecta divergência suave (campo crescendo, não só explodindo) e
devolve divergência por faces, fluxo de contorno e velocidade máxima num
objeto — com `summary()` para quando um humano está lendo.

**Registrado no ctest como `aether_pipeline_tests`** (13 suítes agora), e
não como script solto: é o único caminho do engine que parte de superfície
*importada* em vez da rede sintética, então uma regressão ali não
apareceria em nenhum outro lugar. O teste inclui a verificação decisiva de
que o checador **pega** a tampa não-tangencial — um checador que sempre
aprova passaria por um teste que só olha o caso bom.

Custo mantido honesto: a primeira versão do teste levava 143s (duas
marchas de 400 passos duplicadas, malha fina à toa). Fundidas numa marcha
só e com malha mais grossa — nenhuma das asserções é sobre resolução, são
identidades exatas e checagens de sinal — caiu para **5,7s**. Suíte
completa: **13/13 em 55s**.

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

### Estado em 2026-08-29: **4.1 fechado** — CG residente, classe nova e autocontida, sem tocar `engine/solver`

Confirmado antes de começar: GPU real nesta máquina (RTX 5080, 16GB,
compute capability 12.0) e CUDA Toolkit 13.3 — trabalho compilável e
mensurável de verdade, não teórico.

**Escopo deliberadamente estreito, mesmo princípio da correção de
`tetrahedralize()`**: uma classe nova em `engine/gpu`
(`ConjugateGradientSolverCuda`), testada e medida de forma independente —
**sem** alterar `engine/solver`/`StaggeredCavityBase3D` (a base dos 6
fechamentos de turbulência 3D). Ligar essa classe ao
`projectToDivergenceFree()` real fica pra um passo seguinte, explícito —
mesmo raciocínio que deixou `insertSteinerPoint()` intocado na correção
anterior.

**Mecanismo**: `Impl` (pImpl) possui 5 buffers de device de tamanho
`nx*ny*nz` (`p`, `rhs`, `residual`, `direction`, `ad`) alocados uma vez no
construtor e reaproveitados por toda chamada de `solve()` — nenhum
`cudaMalloc`/`cudaFree` por iteração nem por chamada. O laço de CG inteiro
(inicialização + todas as iterações) roda no device; só dois escalares de
8 bytes cruzam o barramento por iteração (os dois produtos internos
reduzidos, pra computar `alpha`/`beta`/checar convergência do lado do
host) — os campos O(n) nunca saem do device entre iterações, que é a
propriedade de residência que de fato importa. A fórmula do stencil
(7 pontos, célula 0 fixada, espelho de Neumann nas bordas) foi extraída
pra `detail/PoissonStencilCuda.cuh`, compartilhada com o `PoissonOperatorCuda`
já existente — refatoração comportamento-preservando, confirmada
re-rodando o teste bit-a-bit dele sem alteração.

**Redução em duas passadas, shuffle de warp, sem atomics — decisão
deliberada**. Uma redução de passe único por `atomicAdd` soma na ordem em
que os blocos terminam — não-determinística execução a execução na mesma
GPU, o que tornaria "medir e explicar a diferença" (a exigência literal
do item 2 deste portão) um alvo móvel. A versão em duas passadas tem
ordem de soma fixa pela indexação de bloco, então a resposta da própria
GPU é reprodutível execução a execução — a única fonte de divergência
GPU-vs-CPU que sobra é a diferença legítima entre "várias somas parciais
paralelas combinadas" e "uma soma sequencial".

**Validação, três checagens independentes, nenhuma bit-a-bit** (o próprio
item 2 deste portão diz que essa régua não cabe aqui):

- **A — reverificação independente na CPU**: recalcula
  `residual = rhs - A*p` do zero a partir da pressão final devolvida pela
  GPU (reaproveitando `cpuReferencePoisson()` já existente em
  `gpu_tests.cpp`). Medido: 9,953e-11 — mesma ordem de grandeza da
  tolerância de convergência (1e-10) que a própria GPU usa internamente.
- **B — invariante da célula 0, caixa-preta**: `pressure[0]` tem que
  igualar `initialGuess[0]` bit a bit, sempre — pega imediatamente se
  algum dos três kernels de atualização esquecer sua guarda `idx==0`.
- **C — trajetória GPU vs. referência de CG na CPU** (função nova em
  `gpu_tests.cpp`, mesmo laço estrutural de `projectToDivergenceFree`):
  medido, os dois lados convergem em **exatamente as mesmas 117
  iterações**, com `maxDiff` entre os campos finais de **2,398e-14** — a
  diferença de ordem de soma acabou sendo desprezível nesta grade, bem
  abaixo do limiar de 1e-10 fixado (com folga de ~4 ordens de grandeza).

**Medido, não estimado — onde a residência realmente compensa**:

| n | células | iterações | GPU | CPU (referência) | razão |
|---|---|---|---|---|---|
| 16 | 4.096 | 166 | 0,0255s | 0,0024s | GPU **10,6x mais lenta** |
| 24 | 13.824 | 260 | 0,0385s | 0,0133s | GPU **2,9x mais lenta** |
| 32 | 32.768 | 352/420 | 0,0532s | 0,0490s | ~empate |
| 48 | 110.592 | 620/648 | 0,0992s | 0,2464s | GPU **2,5x mais rápida** |
| 64 | 262.144 | 871/899 | 0,1434s | 1,0196s | GPU **7,1x mais rápida** |
| 96 | 884.736 | 1331/1426 | 0,2776s | 5,7957s | GPU **20,9x mais rápida** |

**Honestidade sobre a forma da curva, não só o número final**: abaixo de
~33 mil células a GPU é mais lenta, não mais rápida — overhead fixo de
lançamento de kernel (6 lançamentos por iteração, centenas de iterações)
domina quando o trabalho por lançamento é pequeno demais. O cruzamento
real fica perto de n=32-40. Só teria sido "aceleração de verdade" no
sentido enganoso se essa tabela tivesse começado em n=48 — está completa
de propósito, do jeito que dá errado até o jeito que dá certo, porque é
exatamente essa curva que diz a partir de que tamanho de malha vale a
pena acionar o caminho GPU (informação que a futura integração com
`StaggeredCavityBase3D` vai precisar: escolher a rota por tamanho de
grade, não sempre GPU). A vantagem cresce rápido e continua crescendo:
o tempo de GPU de n=64→96 (3,375x mais células, 1,53x mais iterações)
subiu só 1,94x, sinal de estar limitada por banda/computação, não por
overhead de lançamento, nesse regime.

**O que fica explicitamente fora deste item**: 4.2 (preditor de momento
em CUDA) e a integração real com `projectToDivergenceFree()` — ambos
dependem desta infraestrutura residente, mas não foram tentados aqui,
por decisão de escopo, não por dificuldade encontrada.

Suíte: **13/13** (continua 13, não 11 — o número já cresceu desde que
este portão foi escrito; nenhuma suíte nova foi adicionada por este
item, a classe nova mora no `aether_gpu_tests` já existente).

### Estado em 2026-08-31: **4.2 fechado** — preditor de momento em CUDA, bit a bit idêntico à CPU

Mesmo escopo estreito da 4.1: classe nova e autocontida
(`MomentumPredictorCuda`), sem tocar `StaggeredCavityBase3D` nem os 6
fechamentos que a usam.

**Simplificação central, verificada linha a linha contra
`StaggeredCavityBase3D.cpp:173-364` antes de portar, não suposta**:
`computeMomentumPredictor` suporta 3 esquemas de convecção via
`schemeTransportValue(...)`, mas essa função, para
`ConvectionScheme::Central` — o único que qualquer um dos 6 fechamentos
de fato constrói hoje — retorna o argumento `centralValue` sem tocar em
mais nada. Toda a busca de ponto distante (`far0`/`far1`) que só existe
pra alimentar o limitador vira computação morta pra esse caso e foi
eliminada por completo, sem nenhuma mudança de comportamento — reduz o
kernel de "limitador de 3 vias + busca de ponto duplo" pra "média central
simples". v e w confirmados como permutações cíclicas exatas de u nessa
forma reduzida, sem surpresa.

**Turbulência (nu_t) totalmente suportada**, ao contrário do esquema de
convecção: os 5 fechamentos turbulentos precisam disso pra dar resultado
certo. Buffer de `nut` no device sempre real (nunca nulo) e zero-
preenchido no caso laminar — remove um branch de três kernels em vez de
um só. **Armadilha real e específica do buffer reaproveitado entre
chamadas**, pega por teste dedicado
(`testMomentumPredictorCudaReZerosNutBetweenCalls`): se uma chamada
turbulenta suja o buffer e a próxima é laminar na mesma instância, sem
`cudaMemset` explícito o campo turbulento anterior vazaria silenciosamente
pro resultado seguinte — nenhum crash, só número errado.

**Preservação de face de contorno por construção, não por caso
especial**: `predict()` copia `u→uStar`/`v→vStar`/`w→wStar` inteiros via
`cudaMemcpy` device-to-device antes de lançar qualquer kernel — os
kernels só escrevem faces interiores, então uma face de contorno nunca é
sequer tocada por eles. Testado à parte
(`testMomentumPredictorCudaPreservesBoundaryFaces`), isolado da fórmula
do interior.

**Validação, e por que bit a bit é a régua certa aqui** (diferente da
4.1): este kernel não tem redução nenhuma — cada célula de saída é uma
expressão fixa e independente de ordem sobre leituras de vizinhos, com
`-fmad=false` já valendo pro alvo inteiro. Medido, os três casos batem
**exatamente**, `maxDiff=0,000e+00`, zero incompatibilidades:

- Caso laminar, grade 6×5×4 (não múltipla do bloco de 8×8×8).
- Caso turbulento, mesma grade, campo `nut` não-nulo — confirma que
  `gammaE`/`gammaW`/`gammaTransverse` portam corretamente, não só o
  caminho laminar.
- Caso laminar rodado **depois** de um turbulento, na mesma instância —
  confirma que a re-zeragem do buffer reaproveitado funciona de verdade,
  não só em teoria.

**Mitigação do maior risco real, copiar-colar entre os três kernels
quase-simétricos**: `vMomentumKernel`/`wMomentumKernel` foram escritos
por transcrição direta do código-fonte já lido e verificado (não por
rotação cega de nomes de variável), e o resultado bit-a-bit contra a
referência de CPU independente é exatamente o que teria denunciado
qualquer deslize — nenhum encontrado.

**O que fica explicitamente fora deste item**: a integração real com
`projectToDivergenceFree()`/`computeMomentumPredictor()` do solver de
verdade, e a medição de speedup que só faz sentido uma vez que isso rode
dentro de um laço de passos de tempo de verdade — ambos dependem da
infraestrutura que 4.1+4.2 já entregam, mas não foram tentados aqui, por
decisão de escopo, não por dificuldade encontrada. Com 4.1 e 4.2 fechados,
esse é o próximo passo natural de toda a Fase 4.

Suíte: **13/13**.

### Estado em 2026-08-31: **Fase 4 fechada** — integração real em `StaggeredCavityBase3D`, portão cumprido

Com 4.1 (CG residente) e 4.2 (preditor de momento) prontos e validados
isoladamente, faltava exatamente o que o próprio portão desta fase exige:
"speedup medido contra a versão CPU no mesmo caso, não estimado" e
"resultado numericamente equivalente ao da CPU" — nenhum dos dois pode
ser medido de verdade fora de um solver real rodando muitos passos de
tempo. Este item fecha isso.

**Escopo, mesmo princípio das duas fases anteriores**: `useGpu=false`
opt-in em `StaggeredCavityBase3D`, comportamento padrão preservado bit a
bit para todo código existente — confirmado rodando a suíte completa
antes e depois sem nenhuma mudança de resultado. Só uma das seis
subclasses ganhou o parâmetro na própria API pública nesta fase,
`StaggeredLidDrivenCavitySolver3D` (a laminar, já era a única cujo
construtor expõe `convection` explicitamente, pelo mesmo motivo: "é a que
mede coisa"). As outras cinco continuam de fora, deliberadamente.

**Segurança de ABI, o ponto mais delicado desta integração**:
`StaggeredCavityBase3D` é base de 6 subclasses em várias unidades de
compilação, então um membro cuja presença dependesse de uma macro de
pré-processador seria um risco real de layout inconsistente entre elas.
Resolvido com pImpl: exatamente um `std::unique_ptr<GpuState> gpuState_`
sempre presente, com `GpuState` incompleto no header — o layout declarado
é sempre um ponteiro, idêntico em qualquer unidade de tradução, com ou
sem `AETHER_HAVE_GPU` definida em qualquer lugar. A definição real de
`GpuState` mora inteiramente em `StaggeredCavityBase3D.cpp`, a única
unidade que precisa saber qual dos dois ramos é real. Efeito colateral
documentado, não escondido: a classe virou não-copiável (já não era
movível, por já ter destrutor declarado) — confirmado por grep que nada
no projeto hoje copia ou move-atribui um `StaggeredCavityBase3D`/
subclasse.

**CMake**: `engine/gpu` precisou passar a ser adicionado antes de
`engine/solver` no `CMakeLists.txt` raiz (só a chamada
`add_subdirectory(engine/gpu)` em si — o bloco de detecção de CUDA já
rodava cedo), pra que `if(TARGET aether_gpu)` em
`engine/solver/CMakeLists.txt` enxergue o alvo. Mesmo condicional que
`bindings/python/CMakeLists.txt` já usa pra `aether_gpu_py` — não é
padrão novo pro projeto, só a primeira vez que uma cadeia de *duas*
bibliotecas estáticas (`aether_solver_tests` → `aether_solver` →
`aether_gpu`) precisa propagar o link de runtime CUDA. Verificado antes
de escrever qualquer ramo de GPU no `.cpp`: a cadeia linkou de primeira.

**Os dois pontos de ramificação**: `computeMomentumPredictor()` ramifica
por completo no topo (GPU calcula os três campos e retorna, ou cai pro
corpo de CPU inalterado); `projectToDivergenceFree()` ramifica só a etapa
do meio — montar `rhs` e corrigir `u_`/`v_`/`w_` a partir de `p_` ficam
incondicionais e idênticos nos dois caminhos, só o solve de `p_` em si
troca entre `gpuState_->cg->solve(...)` e o laço de CG da CPU já
existente.

**Validação — três testes de correção mais um de comparação direta**:
repouso com `useGpu=true` (zero continua exatamente zero, já que não há
nada a somar em lugar nenhum), conservação de massa/topologia de vórtice
com `useGpu=true` (mesmos limiares físicos já medidos pro caminho de
CPU), e a checagem que de fato fecha o portão — duas instâncias com
condições iniciais idênticas, uma de cada caminho, mesmos passos,
`maxDiff` medido e impresso antes de qualquer `AETHER_CHECK`.

**Achado honesto sobre a divergência entre os dois caminhos, medido em
duas contagens de passo pra caracterizar a tendência, não só um número
solto**: em n=16/Re=10, após 20 passos, `maxDiff` u=7,120e-04
v=2,555e-04 w=4,062e-04 p=9,979e-03; após 100 passos,
u=3,453e-03 v=8,491e-04 w=4,503e-03 p=2,643e-02 — um aumento de ~4,85x
para 5x mais passos, ou seja, **crescimento aproximadamente linear no
número de passos, não exponencial**. Isso é a assinatura de cada passo
injetando um viés sistemático pequeno e constante (a ordem de soma da
redução paralela do CG residente) que se acumula aditivamente, não de
uma instabilidade numérica amplificando uma única perturbação — e os
testes de repouso/conservação de massa já confirmam que as duas
trajetórias são fisicamente válidas de forma independente uma da outra.
Limiares do teste travados com ~2x de folga sobre o valor medido em 20
passos.

**Medido, não estimado — o número que fecha o portão**, mesmos tamanhos
da tabela da Fase 4.1, agora com preditor de momento + CG rodando juntos
por passo, dentro do solver de verdade:

| n | células | passos | CPU/passo | GPU/passo | razão |
|---|---|---|---|---|---|
| 16 | 4.096 | 10 | 0,0029s | 0,0275s | GPU **9,1x mais lenta** |
| 24 | 13.824 | 10 | 0,0160s | 0,0415s | GPU **2,6x mais lenta** |
| 32 | 32.768 | 8 | 0,0498s | 0,0579s | GPU **1,16x mais lenta** |
| 48 | 110.592 | 6 | 0,2615s | 0,0940s | GPU **2,78x mais rápida** |
| 64 | 262.144 | 5 | 1,1169s | 0,1334s | GPU **8,37x mais rápida** |
| 96 | 884.736 | 4 | 5,5696s | 0,2642s | GPU **21,08x mais rápida** |

**A confirmação mais importante desta tabela não é a velocidade em si —
é que ela bate quase exatamente com a tabela isolada da Fase 4.1** (que
media só o CG, sem preditor de momento nenhum): 20,9x medido lá contra
21,08x medido aqui em n=96, mesma ordem de grandeza em cada tamanho
intermediário. Isso confirma que o preditor de momento não desloca o
ponto de cruzamento (perto de n=32-40, ~33-65 mil células) — o CG
continua dominando o custo por passo, já que precisa de centenas a
milhares de iterações contra um único lançamento de kernel por
componente do preditor. O portão da Fase 4 está cumprido: speedup medido
(não estimado) no caso real, resultado numericamente equivalente dentro
de uma tolerância medida e explicada.

**O que fica explicitamente fora**: as outras cinco subclasses
turbulentas ainda não têm `useGpu` na própria API — herdam `gpuActive()`
(sempre `false`) mas não têm como ligar o caminho GPU ainda. Estender é
mecânico (mesmo parâmetro à direita, mesmo repasse posicional), fica pra
quando a necessidade aparecer.

Suíte: **13/13**.

### Estado em 2026-09-01: `useGpu` estendido às cinco subclasses turbulentas

Conveniência mecânica, não um novo marco de fase -- o portão da Fase 4 já
estava cumprido acima. `MixingLengthLidDrivenCavitySolver3D`,
`KEpsilonLidDrivenCavitySolver3D`, `KOmegaSSTLidDrivenCavitySolver3D`,
`SmagorinskyLesLidDrivenCavitySolver3D` e `DesSstLidDrivenCavitySolver3D`
ganharam o mesmo parâmetro `useGpu=false` à direita, repassado
posicionalmente ao construtor da base -- exatamente como previsto acima,
porque o preditor de momento da base já repassava para o caminho GPU
qualquer campo que `setEddyViscosityField()` tivesse registrado
(`StaggeredCavityBase3D.cpp`); nenhuma das cinco precisou de código
específico de GPU, só expor o parâmetro.

O que genuinamente precisava de verificação era a fiação, não o
mecanismo (já validado na classe laminar): que o nu_t de cada fechamento
realmente chega ao caminho GPU e produz um resultado fisicamente
consistente. Cinco testes novos, mesmo padrão de
`testStaggeredLidDrivenCavity3DGpuMatchesCpuWithinMeasuredTolerance`
(instâncias CPU/GPU idênticas, mesmos passos, maxDiff medido e impresso
antes do `AETHER_CHECK`), n=16, Re=100:

| Fechamento         | maxDiff u   | v          | w          | p          |
|--------------------|-------------|------------|------------|------------|
| Mixing length       | 1.011e-02   | 2.063e-03  | 6.819e-03  | 1.347e-02  |
| Smagorinsky LES      | 9.466e-03   | 1.795e-03  | 6.370e-03  | 1.246e-02  |
| k-epsilon            | 4.042e-02   | 1.768e-02  | 6.207e-02  | 4.052e-02  |
| k-omega SST          | 4.077e-02   | 1.848e-02  | 6.720e-02  | 4.145e-02  |
| DES-SST              | 4.105e-02   | 1.855e-02  | 6.853e-02  | 4.181e-02  |

**Achado, consistente e explicado**: os dois fechamentos algébricos
(mixing length, Smagorinsky -- nu_t/nu_sgs sem equação de transporte
própria) ficam na mesma faixa da classe laminar. Os três com equações de
transporte (k-epsilon, k-omega SST, DES-SST, este último construído sobre
o SST) divergem cerca de 4x mais: k/omega/epsilon realimentam nu_t na
quantidade de movimento a cada passo, então uma divergência CPU/GPU
pequena na velocidade se acumula por esse laço de realimentação em vez de
ficar uma diferença isolada -- um padrão consistente entre os três, não
um outlier de um deles.

Suíte: **13/13**.

---

## Fase 5 — Geometria CAD real

**Objetivo**: importar STEP/IGES.

**Por que aqui**: só tem valor depois que malha não-estruturada alimenta
solver (Fase 3) — antes disso, importar CAD sofisticado para simular uma
caixa não entrega nada.

**Concluído em 2026-08-25, parcialmente e por decisão explícita de escopo**:
importação de STEP (ISO 10303-21) para o subconjunto **BREP facetado** —
`FACETED_BREP`/`MANIFOLD_SOLID_BREP` cujas faces são todas `POLY_LOOP`
planos, sem furos. Sem OpenCASCADE, sem dependência nova: um tokenizer
Part 21 de ~450 linhas (`engine/geometry/src/StepIO.cpp`) mais o
triangulador planar que o motor já tinha (`PolygonTriangulation2D`,
reaproveitado sem alterar). Cada entidade EXPRESS usada
(`CARTESIAN_POINT`, `POLY_LOOP`, `FACE_BOUND`/`FACE_OUTER_BOUND`, `FACE`,
`CLOSED_SHELL`, `(FACETED_)MANIFOLD_SOLID_BREP`) foi conferida contra o
schema publicado da ISO 10303-42 antes de ser implementada, não assumida de
memória — a mesma disciplina que este projeto já aplica a física.

**O que continua em aberto, e por quê isso não é o portão original
completo**: essa fatia cobre só a parte de um arquivo STEP que já chega
triangulada. Uma face curva (cilíndrica, B-spline, `ADVANCED_FACE`) ou a
representação tesselada do AP242 (`TRIANGULATED_SURFACE_SET`, instanciada
via sintaxe de entidade complexa do Part 21) não são interpretadas — o
carregador **relata isso explicitamente** em
`StepLoadResult::unsupported_features` em vez de produzir uma malha errada
ou vazia em silêncio (ver o header de `StepIO.hpp` para o motivo exato de
`TRIANGULATED_SURFACE_SET` ter ficado de fora: a sintaxe da instância real
não foi verificada com confiança suficiente para implementar). Um arquivo
com geometria curva de verdade ainda exige um kernel CAD — a decisão sobre
OpenCASCADE (**a segunda dependência externa pesada do projeto**, depois do
CUDA) continua seus, e sem prazo.
**IGES não foi tocado** — nem a fatia facetada: seu formato de entidades é
diferente o bastante do Part 21 do STEP que não haveria reaproveitamento
direto do parser escrito aqui.
Placement de assembly (`AXIS2_PLACEMENT_3D` via item mapeado/NAUO) também
não é aplicado — um arquivo multi-sólido carrega a união bruta da geometria
de cada sólido, sem transformar nenhum deles.

**Portão de conclusão original, ainda parcialmente aberto**: um arquivo
STEP real importa, tetraedraliza e resolve, com o volume da geometria
batendo com o do CAD de origem. Verificado para o caso facetado (teste
`testStepLoadsFacetedTetrahedron` em `geometry_tests.cpp`, mais o binding
Python em `test_unstructured_bindings.py`); o caso curvo/IGES continua
dependendo da decisão acima.

### Estado em 2026-09-04: **decisão tomada — OpenCASCADE trazido**; geometria curva fecha o portão original

A decisão que era "sua, e sem prazo" foi tomada: trazer o OpenCASCADE
(OCCT 8.0.1, via vcpkg em modo manifesto -- `vcpkg.json` na raiz do
repositório, o mesmo modelo reprodutível que já documenta o CUDA Toolkit
como pré-requisito externo em vez de vendorizado). Build medido: 18
minutos via `vcpkg install` num ambiente sem cache binário.

**Detectado, não exigido -- mesmo padrão do CUDA.** `loadStep()` e
`StepLoadResult` não mudaram de assinatura; `AETHER_HAVE_OPENCASCADE`
(privado, nunca num header público -- a mesma regra que `AETHER_HAVE_GPU`
já segue) seleciona qual das duas implementações responde. Sem o
toolchain do vcpkg configurado, o build continua idêntico a antes desta
entrada: o tokenizer próprio de ~450 linhas continua sendo a
implementação, sem nenhuma dependência nova. **`engine/geometry/src/StepIO.cpp`
não teve uma linha do tokenizer removida** -- a implementação antiga
inteira permanece, ramo `#else` do mesmo arquivo.

**Com OpenCASCADE**, `loadStep()` delega a `STEPControl_Reader` +
`BRepMesh_IncrementalMesh` (tesselação com deflexão relativa 1e-3 /
angular 0,3 rad) para QUALQUER face, curva ou plana, e um `loadIges()`
novo (sempre declarado, mesma razão de estabilidade de ABI -- lança
`std::runtime_error` num build sem OpenCASCADE, já que não existe fatia
facetada de IGES para cair de volta) faz o mesmo via
`IGESControl_Reader`. Não existe mais uma categoria "sintaxe válida,
geometria não suportada" nesse caminho, então um arquivo que o kernel
genuinamente não resolve lança em vez de popular
`unsupported_features`.

**Validação -- um cilindro real, volume exato**: em vez de tentar
escrever à mão um Part 21 com face curva (impraticável, ao contrário do
tetraedro facetado), o teste novo (`testStepLoadsCurvedCylinderThroughOpenCascade`,
`geometry_tests.cpp`) constrói a geometria de referência *com o próprio
OpenCASCADE* (`BRepPrimAPI_MakeCylinder`), escreve um STEP real via
`STEPControl_Writer`, lê de volta pela API pública `loadStep()`, e compara
o volume tesselado contra `pi*r^2*h`: **erro relativo medido 8,305e-04**
(raio 2, altura 5, 352 triângulos) -- um portão fechado com round-trip de
verdade por geometria curva, não um caso plano disfarçado.

**Achado honesto, encontrado medindo em vez de supondo**: os dois testes
existentes do tokenizer (`testStepLoadsFacetedTetrahedron` e
`testStepReportsNonFacetedBoundaryAsUnsupported`, mais o equivalente em
`test_unstructured_bindings.py`) hand-escrevem Part 21 mínimo -- sem o
`PRODUCT`/`SHAPE_DEFINITION_REPRESENTATION`/`GEOMETRIC_REPRESENTATION_CONTEXT`
que uma exportação CAD real sempre carrega -- porque o tokenizer procura
`FACETED_BREP` pelo nome do tipo diretamente e não precisa de nada disso.
O `STEPControl_Reader` do OpenCASCADE **recusa transferir qualquer forma**
desse mesmo arquivo mínimo ("file transferred no shape"), porque sua
mecânica de transferência procura por uma entidade-raiz que este arquivo
não tem. Não é defeito de nenhum dos dois lados -- é exatamente a
diferença entre um parser deliberadamente permissivo para um subconjunto
restrito e um kernel real que insiste em um arquivo completo e
padrão-conforme, que é o que o torna confiável em entrada do mundo real.
Os dois testes do tokenizer passaram a compilar só na configuração sem
OpenCASCADE (`#ifndef AETHER_HAVE_OPENCASCADE`); o teste Python usa uma
nova função consultável em tempo de execução,
`aether.step_io_has_opencascade()`, para adequar as verificações à
implementação realmente ativa (uma macro de compilação não é visível do
lado Python) -- mantendo bit a bit as mesmas verificações de antes na
configuração sem OpenCASCADE, e cobrindo só o que genuinamente vale para
as duas (lista Python, sobrevivência do `mesh` após o `StepLoadResult`
ser coletado) quando ativo.

**Correção incidental, necessária, documentada**: o toolchain do vcpkg
envolve toda chamada `find_package` com sua própria macro, que não aceita
um caminho com barras invertidas cruas -- quebrava a detecção existente
do pybind11 (que lê `python -m pybind11 --cmakedir`, sempre em barras
invertidas no Windows) assim que o toolchain entra em jogo. Corrigido com
`file(TO_CMAKE_PATH ...)` antes do `find_package`, sem efeito em builds
sem o toolchain do vcpkg.

**Deliberadamente fora do escopo**: placement de assembly
(`AXIS2_PLACEMENT_3D`/NAUO -- `STEPCAFControl_Reader`/`XCAFDoc` do
próprio OpenCASCADE seria a extensão natural), reinicialização/qualidade
de malha configurável (deflexão fixa por ora). Nenhum dos dois é exigido
pelo portão original, que só pede importar+tetraedralizar+bater volume
num sólido.

Suíte: **13/13**, nas duas configurações (com e sem o toolchain do
vcpkg). Com isso, o portão original da Fase 5 está fechado para o caso
curvo -- só falta o que nunca foi prometido por ele (assembly).

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

### Estado em 2026-08-31: **6.1 fechada pelo escopo formal** — correção a uma rotulagem de 2026-08-25

Antes de abrir a 6.2, uma releitura desta seção achou uma divergência: o
item registrado em 2026-08-25 como fechando "Fase 6.1" (seção de sete
itens daquela rodada autônoma) é transporte de temperatura passivo e
Boussinesq em `UnstructuredCavitySolver3D` — inteiramente dentro do
fluido, sem nenhuma região sólida. A 6.1 como definida acima é outra
coisa: sólido e fluido **acoplados na interface**. Busca em todo o
repositório por qualquer acoplamento sólido-fluido deu zero resultados
fora deste texto (as únicas ocorrências de "conjugat" fora do ROADMAP são
"Conjugate Gradient", o solver linear, sem relação). A 6.1 formal nunca
tinha sido construída; o item de 2026-08-25 fechou uma capacidade real,
mas mais estreita do que o nome sugeria.

**O que fecha agora.** `UnstructuredDiffusionSolver::setConductivity()` --
condutividade térmica por célula, selecionada por posição na mesma
convenção de `wallVelocity`/`isOutlet`/`setDirichletBoundary`. Uma região
"sólida" e uma "fluida" são só dois valores diferentes de k; nenhuma
classe nova, porque a maquinaria (condição de contorno, termo fonte,
correção não-ortogonal deferida) já existia inteira em
`UnstructuredDiffusionSolver` -- faltava só o coeficiente variar por
célula, exatamente a lacuna que o texto de 6.1 já apontava ("difusão e
Navier-Stokes já existem; falta o acoplamento na interface").

O acoplamento em si é a condutividade de face como média harmônica
ponderada pela distância -- `k_face = 1/((1-w)/k_P + w/k_N)`, com `w` a
mesma `ownerWeight` que a decomposição de face já carrega -- o
coeficiente de dois pontos que a continuidade de fluxo através de uma
interface exige, e que se reduz exatamente ao operador antigo quando
`k_P == k_N`. `UnstructuredFvmBase` ganhou três métodos novos e paralelos
(`applyWeightedLaplacian`, `weightedNonOrthogonalCorrection`,
`solveWeightedDeferredCorrection`) em vez de mudar os três existentes:
risco zero para os outros dois consumidores da base
(`UnstructuredScalarTransportSolver`, a solve de pressão de
`UnstructuredCavitySolver3D`), que continuam chamando exatamente o que já
chamavam.

**Validação, a mesma regra de sempre**: a parede composta -- duas lâminas
de condutividade k1=1, k2=4 em série, regime permanente, sem fonte. Como
nada gera ou absorve calor, o fluxo é a mesma constante nas duas lâminas,
`q = (Tesq-Tdir)/(L1/k1+L2/k2)`, e cada lâmina é linear com inclinação
`-q/k` -- temperatura contínua, inclinação descontínua na interface pela
razão k2/k1. Derivada de um balanço de fluxo, não de uma tabela.

Medido em `testConjugateHeatTransferMatchesCompositeWallSolution`, n=8
(3184 células): erro RMS ponderado por volume 5.09e-01, erro máximo
4.18, contra um salto de temperatura de 80 -- 0.64% de erro RMS relativo.
Um segundo teste (`...UniformConductivityMatchesUnweightedSolve`) confirma
que condutividade uniforme 1.0 reproduz `solveConjugateGradient()` sem
`setConductivity()` dentro de 1.7e-12 -- o caminho novo não perturbou o
antigo.

**Achado honesto, não escondido**: o erro máximo (4.18) é bem mais largo
que o RMS (0.51) porque a interface não é conforme à malha -- os
centroides das células cruzam x=0.5 na malha jitterizada em vez de uma
camada de faces sentar exatamente ali, e o gradiente de mínimos quadrados
usado na correção não-ortogonal mistura valores dos dois lados bem na
célula onde a solução tem uma quina genuína (inclinação descontínua). Essa
mistura amplia o erro exatamente nas faces onde esta malha já é menos
ortogonal (~1.5-1.75, medido em outros itens deste projeto) -- um custo
real e explicado do tratamento mais simples possível da interface, não um
defeito no acoplamento em si: o erro cresce proporcionalmente ao tamanho
do salto de condutividade (medido também com k2=2: erro cai pela metade
nos mesmos pontos), a assinatura de um artefato de reconstrução de
gradiente, não de um bug. Uma malha conforme à interface, ou um gradiente
que exclua o vizinho do outro material do próprio ajuste de mínimos
quadrados, eliminaria isso -- não feito aqui porque nenhum dos dois é "o
passo mais curto" que 6.1 pedia.

**Deliberadamente fora do escopo desta etapa**: nenhuma velocidade em
lugar nenhum -- isola o que é novo (o acoplamento na interface) do que já
existe (advecção de temperatura, `EnergyModel` de
`UnstructuredCavitySolver3D`). Acoplar isto a um escoamento real ao redor
de um obstáculo sólido é uma extensão natural e separada, registrada aqui
como próximo passo explícito, não como lacuna escondida -- exatamente
para não repetir o que motivou esta correção.

Binding Python: `UnstructuredDiffusionSolver.set_conductivity(...)`, ao
lado de `set_source_term`. Suíte: **13/13**.

### Estado em 2026-08-31: **6.2 fechada** — Euler 1D, fluxo de Rusanov, choque de Sod

Módulo novo e autocontido: `CompressibleEulerSolver1D` (Euler 1D em forma
conservativa, gás ideal, fluxo numérico de Rusanov, Euler explícito no
tempo) -- a mudança de regime que o próprio texto desta fase já avisava,
não uma extensão de nenhum solver existente. Deliberadamente a fatia mais
simples que ainda captura choque de verdade sem oscilar: Rusanov não
precisa de decomposição de autovalores (Roe) nem separação de
velocidades de onda (HLLC), só o maior sinal local, o mesmo raciocínio
que já levou este projeto a upwind de primeira ordem antes de
linear-upwind limitado. MUSCL/Roe/HLLC de segunda ordem ficam registrados
aqui como extensão natural, não como parte deste portão.

**As duas formas de portão que esta fase permite, ambas fechadas**:

1. **Conservação exata, incondicional.** Com parede refletora nos dois
   extremos (fantasma: mesma densidade/pressão, velocidade negada -- força
   u=0 exatamente na face), o fluxo de massa e de energia é zero ali para
   qualquer condição inicial, sem exigir simetria nenhuma -- só o momento
   não é testado como invariante, porque uma parede empurra de volta
   sempre que a pressão nas duas pontas difere. Medido: deriva de massa
   1.11e-15, deriva de energia 5.77e-15 depois de t=1,0 (a condição de
   Sod, mas com paredes fechadas, dando tempo para as ondas baterem nas
   duas paredes várias vezes) -- ruído de ponto flutuante, não um limiar
   solto.

2. **Solução exata: o tubo de choque de Sod.** Entradas clássicas (não um
   resultado memorizado) -- esquerda rho=1,u=0,p=1; direita
   rho=0.125,u=0,p=0.1, descontinuidade em x=0.5, gamma=1.4, saída livre
   nas duas pontas, avançado até t=0.2. A solução exata é derivada, não
   copiada: um solver de Riemann exato (relação isentrópica do lado em
   rarefação, Rankine-Hugoniot do lado em choque, Newton-Raphson na
   pressão da região-estrela) implementado no próprio teste, o mesmo lugar
   onde a solução da parede composta foi derivada para 6.1.

   Medido, erro RMS ponderado por célula (mesma escolha da parede
   composta e da solução manufaturada, pela mesma razão: o esquema borra
   o choque e o contato por várias células, então um máximo pontual bem
   na descontinuidade não diz nada sobre o resto do domínio):

   | n   | rms densidade | rms pressão | ordem |
   |-----|----------------|-------------|-------|
   | 100 | 3.556e-02      | 3.732e-02   | --    |
   | 200 | 2.662e-02      | 2.576e-02   | 0.42  |

   **Achado honesto**: a ordem medida (0.42) é menor que o ~1 que um
   esquema de primeira ordem sugeriria pela parte suave sozinha. Não é
   defeito: a norma RMS sobre o domínio inteiro mistura um erro de região
   suave que cai em primeira ordem com um erro de choque/contato cuja
   *magnitude* fica quase constante enquanto só sua *largura* encolhe como
   O(h) -- propriedade conhecida de esquemas Godunov sobre dado
   genuinamente descontínuo, não algo que esta implementação especificamente
   devesse superar.

**Limite de CFL, medido, não suposto**: estável até cfl=1.05 no tubo de
Sod (n=100), diverge (NaN) em cfl=1.10 -- bate com o limite clássico
CFL<=1 de um esquema Godunov explícito. Default fixado em 0.4, ~2,5x de
margem abaixo do limite medido, a mesma folga conservadora que o 0.3 de
`ExplicitTimeStep.hpp` já mantém para os solvers incompressíveis.

**Deliberadamente fora do escopo desta etapa**: sem viscosidade (já existe
em todo o resto do projeto; o que faltava era choque e compressibilidade,
que é exatamente o que Euler isola), sem 2D/3D, sem fluxo de segunda
ordem, sem modo novo no `unified_viewer` -- nenhum exigido pelo portão
desta fase. Extensões naturais, registradas, não esquecidas.

Suíte: **13/13**.

### Estado em 2026-09-01: **6.3 fechada — e com ela, a Fase 6 inteira** — rastreamento de interface por level set

Módulo novo: `LevelSetAdvectionSolver2D` -- só a primeira das duas coisas
que "6.3" nomeia (rastreamento de interface), sob um campo de velocidade
**prescrito analiticamente, não resolvido**: nenhum Navier-Stokes de duas
fases, nenhuma tensão superficial. Multifásico completo é comparável em
escopo ao motor inteiro; esta é a fatia mínima que ainda demonstra algo
que uma advecção passiva comum (já presente em várias classes deste
projeto) nunca precisou provar -- uma forma reconhecível voltando à sua
posição exata depois de mover-se de verdade. Rotular o que já existia
como isto teria repetido o erro que a correção da 6.1 encontrou; por isso
o critério de validação abaixo exige mais do que "o campo advecta".

**Método: level set, não VOF.** Distância assinada `phi`, advectada em
forma não-conservativa `dphi/dt + u.grad(phi) = 0` -- sem reconstrução
geométrica de interface (PLIC), ao preço de precisar de reinicialização
periódica em escoamentos longos (não implementada aqui: só afeta a norma
do gradiente longe da interface, não a posição dela sob advecção pura,
que é tudo que o teste de uma rotação completa precisa). Esquema de
diferença finita upwind de primeira ordem (o clássico de Hamilton-Jacobi)
e Euler explícito no tempo -- o par mais simples, mesmo raciocínio de
6.2.

**Validação, as duas formas do mesmo experimento**: rotação de corpo
sólido de um círculo (raio 0.15 em (0.65,0.5), domínio [0,1]², rotação em
torno de (0.5,0.5), omega=2*pi dando período T=1). Sem divergência por
construção, então:

1. **Conservação exata** (teorema do transporte de Reynolds): a área
   onde phi>0 deveria ficar em pi*0.15² = 0,070686 durante todo o
   percurso.
2. **Solução exata**: rotação rígida por um período completo é a
   identidade -- a solução exata em t=T é phi0 de novo, sem precisar
   derivar nada à parte (diferente do tubo de choque de Sod, que
   precisou de um solver de Riemann).

**Achado honesto, medido, maior do que uma primeira suposição sugeriria**:
em n=128, a área cai de 0,070686 para 0,050171 ao longo de uma rotação
completa -- **29% perdidos para difusão numérica**, não um limiar solto.
Não é bug: confirmado dobrando a resolução (n=256), a perda cai pela
metade (14%, para 0,060715) -- a assinatura O(h) de um esquema
genuinamente de primeira ordem, não um esquema estagnado ou divergindo.
Upwind de primeira ordem é conhecido por ser particularmente difusivo
justamente em rotação de corpo sólido: a direção do escoamento relativo
aos eixos da malha varre todos os ângulos ao longo de um período, o que
maximiza a difusão numérica de sentido cruzado (a fraqueza clássica desta
família de esquema fora do alinhamento com a malha) -- não um caso de
canto. Exatamente por isso que ENO/WENO ou uma reconstrução geométrica
(PLIC) são a extensão natural registrada aqui, não parte deste portão.

**Limite de CFL, medido e mais restritivo do que a intuição de uma única
direção sugeriria**: estável até cfl=0,70, diverge (|phi| chegando a
1e24) em cfl=0,75. A atualização em 2D num único passo não fracionado
precisa do limite combinado `dt*(|u|/dx+|v|/dy)<=1`, que para |u|~|v| e
dx=dy permite só cerca da metade do que a fórmula de uma única direção
sugeriria em cfl=1 -- consistente com a quebra medida perto de 0,70-0,75,
não de 1,0. Default fixado em 0,4, ~1,75x de margem abaixo do limite
medido.

**Deliberadamente fora do escopo desta etapa**: reinicialização, ENO/WENO
ou PLIC, acoplamento a um Navier-Stokes de duas fases resolvido, tensão
superficial. Todas registradas como extensões naturais e separadas, não
como lacunas escondidas.

Suíte: **13/13**. Com isso, a Fase 6 do ROADMAP -- as três físicas novas,
em ordem de dependência crescente -- está inteiramente fechada.

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

- ~~**Wheel/`pip install`**~~ — FEITO em 2026-08-27: `pyproject.toml` +
  scikit-build-core orquestrando o `CMakeLists.txt` existente (nenhuma
  duplicação de lógica de build). `pip wheel .` produz uma wheel real com
  os 11 módulos compilados dentro do pacote `aether`; verificado instalando
  numa venv nova, fora do repositório, e rodando um solver de verdade a
  partir dela. Duas opções de CMake novas (`AETHER_BUILD_APPS`, junto da
  `AETHER_BUILD_TESTS` já existente) deixam a wheel pular o visualizador
  desktop e a suíte C++ — ambas `ON` por padrão, então nenhum build local
  ou de CI existente muda. O que ainda falta é distribuição, não
  empacotamento: publicar no PyPI e gerar wheels por plataforma/versão de
  Python via CI é trabalho separado.
- **DNS**: computacionalmente impraticável nas escalas de malha atuais.
  Reavaliar só depois da Fase 4 (GPU) e com malhas muito maiores.
- **Backend Vulkan** (Módulo 8): esforço grande e isolado; o pipeline GL 3.3
  atende tudo que o projeto renderiza hoje. Só se o gargalo virar
  renderização, o que não é o caso.
- ~~**Recuperação combinatória de facetas na CDT** (caso do quadrilátero
  coplanar)~~ — FEITO em 2026-08-27 (`tryFlipCoplanarQuadDiagonal()`,
  DIVIDA_TECNICA.md §6): medido recuperando 100% das facetas de um
  cilindro aproximado por prisma poligonal (10/24/60 lados), sem nenhum
  ponto de Steiner. O caso Schönhardt-difícil (contorno genuinamente
  não-convexo) continua em aberto — nunca construído neste código, citado
  como teto teórico. ~~Ainda falta validar o pipeline completo~~ — FEITO
  em 2026-08-28: `mesh_flow_around_object` + `freestream_boundary` +
  `UnstructuredCavitySolver3D` sobre um cilindro real importado via
  round-trip STL binário (não mantido em memória), registrado como teste
  permanente em `test_pipeline.py` (roda no ctest a cada build) — 0 de 48
  facetas perdidas, volume exato, desbalanço de massa de 0,0014% da
  entrada.
- ~~**Validação contra benchmarks de literatura** (Ghia et al. e afins)~~ —
  FEITA em 2026-08-28 (`python/research/ghia_1982_validation.py`), com os
  dados obtidos de fontes reais (duas transcrições independentes do artigo,
  a coluna Re=100 cruzada entre elas), não de memória — exatamente o
  cuidado que este item sempre disse que exigiria. **Resultado é honesto,
  não um "bateu": `LidDrivenCavitySolver2D` não reproduz a tabela de Ghia
  diretamente, e refinar a malha (n=64→128, 4x as células) não fecha a
  diferença** — o que descarta erro de discretização como explicação, já
  que esse encolheria com refino, e aponta pra uma diferença na *condição
  de contorno posta*, não um bug: o lid deste solver usa
  `lidVelocity·sin²(πx/L)`, suavizado nos cantos pra remover a
  singularidade de pressão/vorticidade de um lid literalmente descontínuo
  — decisão deliberada e documentada, não descoberta agora — enquanto a
  cavidade clássica de Ghia usa lid uniforme. A média desse taper sobre
  toda a largura é exatamente 0,5 (fato de forma fechada, verificado por
  integração numérica no próprio script) — metade do momento de um lid
  uniforme, uma redução bem maior que a das cavidades "regularizadas" da
  literatura (que só suavizam uma faixa estreita perto de cada canto,
  mantendo o lid em velocidade plena na maior parte da largura), o que
  explica por que o desvio medido aqui (~10-20% do próprio range do campo)
  é maior que o "pequeno, concentrado nos cantos" que essa literatura
  registra.

  **Atualização, mesmo dia: a ablação foi feita, e confirma a hipótese.**
  `LidDrivenCavitySolver2D` ganhou o parâmetro `taper_lid` (default `True`
  — reproduz bit a bit todo número já publicado deste solver, nenhum teste
  existente muda), especificamente pra rodar com lid uniforme (o problema
  que Ghia de fato propôs) em vez do lid regularizado. Resultado em Re=100,
  n=64: os erros caem de u rms/max 0,0618/0,1187 pra **0,0041/0,0087**, e
  v rms/max de 0,0318/0,0531 pra **0,0095/0,0175** — uma redução de 3x a
  15x, abaixo de 1% do próprio range do campo. A diferença era mesmo o
  perfil do lid, não um defeito de discretização ou física, e agora isso
  está medido, não só inferido. O default do solver continua sendo o lid
  regularizado (correto pra tudo que não é reproduzir o problema clássico
  de Ghia especificamente); a claim de credibilidade externa que este item
  buscava agora está feita de verdade, não aproximada.

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
