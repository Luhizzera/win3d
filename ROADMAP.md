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

## Fase 1 — Rhie-Chow nos solvers colocados 2D

**Objetivo**: interpolação de Rhie-Chow em `TaylorGreenVortexSolver2D`,
`LidDrivenCavitySolver2D` e nas três variantes 2D de turbulência.

**Por que aqui**: é a lacuna numérica documentada há mais tempo no projeto
(desde o primeiro solver de Navier-Stokes) e a única que **já tem métrica
pronta** — o `checkerboardIndex` do Módulo 12.2 mede exatamente o desacoplo
par-ímpar que Rhie-Chow existe para eliminar. Hoje o campo de pressão da
cavidade 2D mede **0,020357**. É um alvo pequeno, bem definido e com
resultado verificável antes/depois. Os solvers 3D staggered não precisam
disso (a grade deslocada evita o problema estruturalmente).

**Portão de conclusão**:
1. `checkerboardIndex` do campo de pressão cai mensuravelmente frente ao
   0,020357 de referência, medido no mesmo caso (20×20, Re=100, 300 passos).
2. Taylor-Green continua batendo com o decaimento exato dentro da tolerância
   atual — a correção não pode custar acurácia na solução analítica.
3. As 11 suítes continuam passando.

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
