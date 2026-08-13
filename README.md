# Aether CFD Engine (win3d)

Motor de simulação de fluidos em camadas, com núcleo numérico em C++ e
camada de orquestração/API/IA em Python.

## Layout do repositório

```
engine/core/          Módulo 1 - core matemático em C++ (Vector3, Tensor3x3,
                       Mesh, ScalarField, VectorField)
engine/geometry/       Módulo 2 (escopo STL+OBJ) - TriangleMesh, weld de
                       vértices, normais, área/volume, detecção de buracos,
                       import/export STL e OBJ (faces com mais de 3 vértices
                       são trianguladas em leque). STEP/IGES ficam de fora
                       por enquanto (exigem um kernel CAD completo, ex.
                       OpenCASCADE)
engine/mesh/           Módulo 3 - StructuredGrid3D (malha Cartesiana
                       uniforme), DelaunayTriangulation2D (Bowyer-Watson,
                       nuvem de pontos) e PolygonTriangulation2D (ear
                       clipping + flips de Delaunay restritos por contorno
                       = Triangulação de Delaunay Restrita real: preserva
                       cada aresta do polígono E otimiza qualidade via
                       critério do circumcírculo). DelaunayTetrahedralization3D
                       estende o mesmo núcleo de Bowyer-Watson para 3D
                       (circumsfera em vez de circumcírculo, cavidade
                       delimitada por faces triangulares em vez de arestas).
                       Inserção incremental de pontos de Steiner
                       *interiores* (`insertSteinerPoint()`, pra refinamento
                       de malha sem reconstruir do zero) já está pronta -
                       tetraedralização restrita por contorno (buracos/
                       múltiplos contornos, preservando uma malha de
                       superfície importada, e pontos de Steiner fora do
                       casco convexo atual) continua como próximo passo
engine/solver/         Módulo 5 - DiffusionProblem (base compartilhada:
                       malha, condições de contorno Dirichlet nas 6 faces,
                       peso 1/h^2 por eixo) com duas especializações:
                       SteadyDiffusionSolver (regime permanente - Poisson/
                       Laplace, Gauss-Seidel, Gradiente Conjugado matrix-free
                       (~100x menos iterações que GS) e CG pré-condicionado
                       por Jacobi) e MultigridPoissonSolver2D (V-cycle
                       geométrico, esquema de correção, convergência
                       independente da resolução da malha - medido em ~11
                       V-cycles contra ~244 iterações de CG e ~7425 de GS no
                       mesmo problema). TransientDiffusionSolver (regime
                       transiente - Euler
                       explícito, condicionalmente estável). Validado contra
                       soluções analíticas fechadas (1D linear, 2D via série
                       de Fourier, decaimento senoidal transiente, Poiseuille
                       plano). Também inclui TaylorGreenVortexSolver2D -
                       Módulo 4: primeiro solver de Navier-Stokes
                       incompressível de verdade (convecção + pressão-
                       velocidade acoplados via projeção de Chorin, grade
                       colocalizada periódica), validado contra o vórtice de
                       Taylor-Green (solução exata não-linear completa, não
                       um caso especial simplificado). LidDrivenCavitySolver2D
                       estende isso pra paredes sólidas - a cavidade com tampa
                       deslizante clássica, com tratamento de célula-fantasma
                       (Dirichlet pra velocidade, Neumann pra pressão) e perfil
                       de tampa regularizado (afunila a zero nos cantos, evita
                       a singularidade de canto do problema com tampa
                       descontínua). MixingLengthChannelFlowSolver1D - Módulo 6
                       (primeira etapa): fechamento de turbulência por
                       comprimento de mistura de Prandtl (algébrico, sem
                       equações de transporte extras) no escoamento turbulento
                       totalmente desenvolvido em canal, validado contra a
                       inclinação da lei da parede logarítmica (1/kappa).
                       KEpsilonChannelFlowSolver1D - Módulo 6 (segunda etapa):
                       k-epsilon padrão de duas equações, funções de parede de
                       equilíbrio, warm-start a partir do comprimento de
                       mistura, também validado contra a inclinação log-law.
                       KOmegaSSTChannelFlowSolver1D - Módulo 6 (terceira etapa):
                       k-omega SST de Menter (blend k-omega/k-epsilon via F1/F2,
                       limitador de Bradshaw). StaggeredNavierStokesSolver3D -
                       Navier-Stokes 3D em malha staggered (MAC), resolve o
                       checkerboard das versões colocalizadas (divergência
                       ~1e-13, praticamente precisão de máquina), validado
                       contra um escoamento de Beltrami genuinamente 3D,
                       periódico. StaggeredLidDrivenCavitySolver3D estende
                       isso pra paredes sólidas nas 6 faces + tampa
                       deslizante no topo (a cavidade 3D) - mesma
                       precisão de máquina na divergência (~2.5e-14) da
                       versão periódica.
                       MixingLengthLidDrivenCavitySolver2D - Módulo 6 (quarta
                       etapa): primeiro fechamento de turbulência acoplado a
                       um escoamento 2D real com convecção/recirculação (não
                       mais o canal 1D totalmente desenvolvido), comprimento
                       de mistura avaliado localmente a partir da taxa de
                       deformação 2D completa e da distância à parede mais
                       próxima. KEpsilonLidDrivenCavitySolver2D - Módulo 6
                       (quinta etapa): k-epsilon de duas equações acoplado à
                       mesma cavidade real, k/epsilon transportados
                       explicitamente a cada passo (advecção+difusão+
                       produção+destruição), tratamento de parede via
                       relação assintótica k~y² em vez de funções de parede
                       de equilíbrio. KOmegaSSTLidDrivenCavitySolver2D -
                       Módulo 6 (sexta etapa): mesma extensão pro SST,
                       condição de contorno de omega perto da parede
                       derivada da mesma relação assintótica via a
                       identidade epsilon=beta*·k·omega.
                       MixingLengthLidDrivenCavitySolver3D - primeiro
                       fechamento de turbulência acoplado a um escoamento
                       genuinamente 3D (a cavidade 3D em malha staggered,
                       não mais a cavidade 2D colocalizada), taxa de
                       deformação completa (6 componentes) e viscosidade
                       efetiva amostrada localmente na malha staggered.
                       KEpsilonLidDrivenCavitySolver3D - k-epsilon de duas
                       equações estendido pra mesma cavidade 3D real, k/epsilon
                       transportados em pontos cell-centered (mesma malha da
                       pressão/nu_t) via estêncil colocalizado de 6 vizinhos.
                       KOmegaSSTLidDrivenCavitySolver3D - SST estendido pra
                       mesma cavidade 3D real, limitador de Bradshaw usando a
                       magnitude completa da vorticidade 3D (vetor, não só o
                       componente z do caso 2D).
                       SmagorinskyLesLidDrivenCavitySolver3D - LES (não mais
                       RANS): resolve as escalas grandes e modela só o
                       submalha, com comprimento característico vindo da
                       MALHA (largura de filtro Delta), não da geometria
engine/postprocessing/ Módulo 7 - Streamline2D (integração real de linhas de
                       corrente via Runge-Kutta 4ª ordem, interpolação
                       bilinear do campo de velocidade cell-centered),
                       marchingSquares2D (extração de iso-contornos de
                       campos escalares 2D) e marchingCubes3D (extração de
                       iso-superfícies 3D via decomposição em tetraedros -
                       ver seção própria abaixo)
engine/testing/        Header-only: macro AETHER_CHECK para os testes (ver nota abaixo)
apps/common/           Módulo 8 - toolkit compartilhado de bootstrap OpenGL
                       3.3 core profile + WGL (contexto, shaders, VBO/VAO),
                       extraído do viewer original para não duplicar ~150
                       linhas de boilerplate (e o bug de janela que já
                       apareceu lá)
apps/unified_viewer/   Módulo 8 - os quatro visualizadores originais num
                       único executável com seleção de modo por argumento de
                       linha de comando (`mesh`/`heatmap`/`cavity`/`turbulence`),
                       mais o modo `cavity3d` (campo de vetores 3D da cavidade
                       staggered, o primeiro viewer de qualquer solver 3D) e
                       `isosurface` (iso-superfície real de nu_t via
                       marchingCubes3D, k-omega SST 3D)
bindings/python/       Bindings pybind11 do core C++ para Python
python/aether/         Pacote Python de orquestração (importa o core compilado)
python/research/       Experimentos isolados, fora da árvore da engine
  lbm_prototypes/       Protótipos LBM (validação de física, não é o solver oficial)
```

**Nota sobre os testes**: os executáveis de teste usam `AETHER_CHECK` (em
`engine/testing/`), não `assert()` do `<cassert>` - `assert()` vira no-op
quando `NDEBUG` está definido, que é o padrão do CMake em builds Release, o
que faria os testes "passarem" sempre, sem checar nada de verdade.

`engine/geometry` cobre importação STL (ASCII e binário) e agora também OBJ
(`ObjIO`) - linhas `v`/`f`, cada token de face aceita as formas `v`, `v/vt`,
`v/vt/vn` ou `v//vn` (só o índice de vértice é usado; textura/normal são
lidos e descartados), faces com mais de 3 vértices são trianguladas em
leque a partir do primeiro vértice (correto para faces planas convexas,
o caso comum de quads gerados por exportadores reais; não é um
triangulador geral de polígonos não-convexos - `PolygonTriangulation2D`
existe pra isso mas opera sobre um contorno 2D, não uma face 3D
arbitrariamente orientada, então incorporá-lo aqui exigiria uma etapa de
projeção que esta primeira versão não tenta). Diferente do `loadStl`, o
`loadObj` **não** chama `weldVertices()`: o formato OBJ já armazena uma
lista de vértices compartilhada referenciada por índice (ao contrário das
cópias independentes por triângulo que o STL grava em disco), então não há
nada pra soldar. Índices negativos/relativos do OBJ não são suportados.
STEP/IGES ficam de fora por enquanto - exigem um kernel CAD completo (ex.
OpenCASCADE), o que é uma decisão de dependência separada.

`DelaunayTriangulation2D` triangula uma nuvem de pontos no plano z=0 (ainda
não restringe a arestas de contorno de uma geometria real - isso é
trabalho futuro). Validado de duas formas independentes: a propriedade que
de fato define uma triangulação de Delaunay (nenhum ponto dentro do
circumcírculo de qualquer triângulo, checado por força bruta) e, num caso
convexo à parte (hexágono regular + centro), a contagem exata de
triângulos que a fórmula de Euler prevê (`2n-2-h`), que deve valer
independente de qual triangulação válida for escolhida entre pontos
co-circulares empatados.

`PolygonTriangulation2D` complementa isso: em vez de uma nuvem de pontos,
recebe um polígono simples (possivelmente côncavo), triangula por ear
clipping (preservando exatamente cada aresta de contorno) e em seguida
melhora o resultado via flips de aresta de Delaunay restritos (algoritmo
de Lawson, nunca invertendo uma aresta de contorno) - o resultado final é
uma Triangulação de Delaunay Restrita (CDT) de verdade. Ainda sem suporte
a buracos/múltiplos contornos, nem inserção de pontos de Steiner (então a
qualidade fica limitada aos vértices de entrada). Validado com um hexágono
em L (côncavo, pra exercitar de verdade o tratamento de vértice reflexo):
contagem exata de triângulos (`n-2`) e soma das áreas batendo exatamente
com a área do polígono (fórmula do shoelace), mais a checagem de que o
resultado é localmente Delaunay em toda aresta interna. Um segundo teste
usa um quadrilátero convexo específico onde o ear clipping puro escolhe
uma diagonal comprovadamente não-Delaunay - confirma que o passo de flip
realmente troca a diagonal (não é um no-op que só passa por acaso).

`DelaunayTetrahedralization3D` generaliza o mesmo algoritmo de Bowyer-Watson
para um ponto genuíno em 3D (não mais restrito ao plano z=0): o teste de
"ponto dentro do circumcírculo" vira "ponto dentro da circumsfera" (via o
determinante 4x4 do parabolóide levantado, mesma estrutura de sinais
alternados do caso 2D, generalizada uma dimensão acima), o triângulo-super
vira um tetraedro-super (4 vértices alternados de um cubo, grande o
suficiente para envolver todos os pontos de entrada) e a cavidade deixada
ao remover tetraedros "ruins" é delimitada por faces triangulares (em vez
de arestas) que pertencem a exatamente um tetraedro ruim. Validado de duas
formas: a propriedade de Delaunay por força bruta num conjunto de pontos
com jitter (grade 3x3x3), e um invariante exato de volume - um bipirâmide
de base quadrada com volume analítico conhecido (2.0), onde a soma dos
volumes (com sinal) de todos os tetraedros gerados precisa bater
exatamente com o volume do sólido, o análogo em 3D do teste de área do
shoelace em 2D. Cobre apenas tetraedralização não-restrita de nuvem de
pontos; suporte a buracos/contornos múltiplos (preservando faces de uma
malha de superfície importada) fica para depois, deliberadamente.

## Inserção de pontos de Steiner (refinamento incremental)

`insertSteinerPoint()` fecha a primeira metade do item que ficou em aberto
na task acima: insere um ponto adicional numa tetraedralização já
calculada, atualizando a malha incrementalmente pelo mesmo passo de
substituição de cavidade que o laço principal de `tetrahedralize()` já usa
por ponto - sem reconstruir a partir de um tetraedro-super do zero. É assim
que um gerador de malha de verdade adiciona pontos de Steiner pra
refinamento de qualidade (um circumcentro, o ponto médio de uma aresta, ou
qualquer outro ponto escolhido pra melhorar a forma dos tetraedros).

**Só suporta pontos *interiores*** (dentro do casco convexo já
tetraedralizado) - o caso comum de refinamento, já que esses pontos vêm de
dados já dentro da malha. Pontos exteriores precisariam de extensão do
casco/malha restrita por contorno, exatamente o trabalho que o comentário
da classe já sinaliza como adiado deliberadamente. Em vez de assumir
cegamente que o ponto é interior, a inserção verifica que a cavidade dos
tetraedros "ruins" que o novo ponto invalida fecha numa fronteira
estanque de verdade (toda aresta de fronteira compartilhada por exatamente
duas faces) - se o ponto está fora do casco, essa checagem falha e o ponto
é rejeitado sem modificar a malha, em vez de silenciosamente produzir uma
tetraedralização inválida.

**Validado com dois invariantes exatos, o mesmo estilo "medir, não supor"
já usado no resto da classe**: numa nuvem de pontos com jitter (grade 3x3x3,
mesma configuração do teste de Delaunay já existente), inserir o centroide
da grade (claramente interior) mantém a propriedade de Delaunay e conserva
o volume total somado dos tetraedros **exatamente** (diferença medida de
`~7e-15`, ruído de ponto flutuante) - já que inserir um ponto interior só
subdivide o volume existente, nunca adiciona nem remove nenhum. Um ponto
claramente fora da nuvem (bem distante da caixa delimitadora) é rejeitado
corretamente, com a malha inalterada (contagem de pontos/tetraedros
idêntica antes e depois).

Nenhum bug encontrado. Todos os 5 suites de teste C++ passam (`ctest -C
Debug`, 5/5) na primeira tentativa. Vinculado ao Python como
`aether.DelaunayTetrahedralization3D.insert_steiner_point`.

**Ainda em aberto**: a segunda metade do item original - tetraedralização
restrita por contorno (buracos, múltiplos contornos, preservando faces de
uma malha de superfície importada), que precisa recuperar um contorno
arbitrário e possivelmente não-convexo (pode exigir seus próprios pontos de
Steiner só pra ser tetraedralizável - a obstrução clássica do poliedro de
Schönhardt) - genuinamente mais difícil, deliberadamente adiada.

## Pontos de Steiner fora do casco convexo (extensão de casco)

`insertSteinerPoint()` ganhou suporte a pontos *exteriores* ao casco
convexo, fechando a lacuna deixada na etapa anterior - agora um único
método cobre tanto refinamento interior quanto extensão de casco.

**A derivação exigiu verificação numérica cuidadosa antes de qualquer
código**, porque a primeira tentativa "óbvia" estava errada. Usando o
truque padrão de "tetraedro fantasma" (toda face de casco tem um vizinho
implícito representando a região infinita além dela; um fantasma é "ruim"
quando o ponto é visível através da face, o limite natural de "dentro da
circumsfera" quando o vértice distante do fantasma vai pro infinito) - a
regra de visibilidade em si (substituir o vértice oposto pelo novo ponto
na mesma ordem de vértices do tetraedro; visível se o volume com sinal
resultante fica negativo) foi validada numericamente contra um exemplo
manual antes de implementar.

**Dois bugs reais foram encontrados e corrigidos durante o
desenvolvimento**, ambos via cross-check contra o `scipy.spatial.Delaunay`
(qhull) - a mesma disciplina de "medir contra uma implementação
independente" em vez de confiar só em invariantes internos:

1. A primeira versão marcava o tetraedro real *inteiro* como "ruim"
   sempre que possuía uma face visível, descartando-o e reconectando todas
   as suas faces não-visíveis ao novo ponto. Isso produzia uma
   tetraedralização com volume total *correto* (passava um teste de
   conservação de volume!) mas que não era a triangulação de Delaunay -
   confirmado comparando contra o `scipy` num exemplo construído à mão,
   que mostrou que o tetraedro original devia ficar intocado (seu próprio
   circumsfera não continha o novo ponto) e apenas *um* novo tetraedro
   devia ser adicionado ao lado dele. Corrigido tratando o "tetraedro
   real ruim" e o "fantasma ruim" como independentes: um fantasma ruim
   com dono real bom só adiciona um tetraedro-tampa, sem descartar o dono.
2. Depois de corrigir o primeiro bug, uma checagem de segurança
   ("fronteira precisa ser estanque") sobrevivente da primeira versão
   passou a rejeitar inserções exteriores *válidas* incorretamente -
   porque a fronteira de um caso exterior genuinamente não forma uma
   casca fechada por si só (o "buraco" da face visível descartada é
   exatamente o que os novos tetraedros preenchem). Substituída por uma
   checagem correta: cada novo tetraedro deve ter volume não-degenerado.

**Validação em larga escala, não só o exemplo manual**: 50-60 execuções
aleatórias comparando contra `scipy.spatial.Delaunay` (nuvens de 4-12
pontos, inserções sequenciais mistas interior/exterior) confirmaram
propriedade de Delaunay + volume idêntico ao `scipy` em praticamente todos
os casos.

**Um bug pré-existente e não relacionado foi descoberto no processo**, e
depois **corrigido** - mas o caminho até a correção vale ser registrado,
porque o primeiro diagnóstico estava *errado*.

Sintoma: `tetrahedralize()` raramente (~1 em 40 configurações aleatórias)
produzia uma tetraedralização com um tetraedro fino, válido e adjacente ao
casco convexo faltando de verdade.

**Primeiro diagnóstico (errado)**: culpei a precisão de ponto flutuante no
predicado de circumsfera, já que o caso que falhava envolvia um tetraedro
quase degenerado ("sliver") com volume de só `~2.7%` do que suas arestas
sugeririam. Era uma hipótese razoável - é exatamente a razão pela qual
CGAL, TetGen e Qhull usam aritmética exata/adaptativa nesses predicados.
Isso levou à construção de predicados de aritmética exata (`BigInt` +
`RobustPredicates`, descritos abaixo).

**Causa real**: reimplementei o algoritmo inteiro em aritmética racional
exata (`Fraction` do Python - zero ponto flutuante em qualquer lugar) e o
mesmo buraco persistiu, com diferença de volume exata e não-nula. Isso
descartou precisão numérica de vez. Rastreando passo a passo, cada
inserção individual era perfeitamente consistente em volume - o problema
estava no **tetraedro-super ser pequeno demais**. Seus quatro vértices são
um substituto *finito* para "pontos no infinito"; quando não estão longe o
bastante, a circumsfera finita de um tetraedro que os inclui pode testar
erradamente como contendo um ponto real próximo, deletando um tetraedro
que deveria sobreviver - e um sliver junto ao casco é exatamente onde esse
erro aparece. Medido: no caso original, `R = 20*deltaMax` dava 6 tetraedros
e volume `0.00096` abaixo do casco convexo real, enquanto `R >= 200*deltaMax`
dava os 7 corretos batendo com o qhull exatamente.

**Correção aplicada**: aumentar o tetraedro-super (com margem generosa
sobre a maior falha observada - custa nada em runtime, são quatro pontos
descartados no final). Validado com **200 execuções aleatórias contra o
qhull: 200/200** batendo volume total e propriedade de Delaunay (antes:
falhava ~1 em 40). Mais 100 execuções combinando `tetrahedralize()` +
`insertSteinerPoint()`: 100/100. Adicionado teste C++ de regressão com os
7 pontos exatos que expuseram o bug, verificando contra o volume do casco
convexo calculado independentemente (via `ConvexHull` do qhull, que não usa
Delaunay nenhum) - vale notar que `satisfiesDelaunayProperty()` **não**
pegava esse bug: a saída defeituosa o satisfazia perfeitamente, porque
circumsfera-vazia não diz nada sobre cobertura completa do casco.

**Lição registrada no código**: "tem um sliver envolvido, logo é problema
de precisão" era plausível e estava errado. Refazer a conta em aritmética
exata é um teste barato e decisivo que deveria ter vindo *antes* de
construir infraestrutura baseada na teoria da precisão estar certa.

## Predicados geométricos exatos (BigInt + RobustPredicates)

Construídos durante a investigação acima. **Não foram o que corrigiu aquele
bug** (a correção foi o tetraedro-super), mas são mantidos como melhoria de
robustez legítima por si só - protegem contra a classe de erro que a
hipótese original descrevia, que é real mesmo não sendo a causa daquele
caso.

`aether::core::BigInt` é um inteiro de largura fixa (1280 bits) feito
sob medida para esses predicados, não uma biblioteca de precisão
arbitrária de uso geral. Limbs de `uint32_t` (não `uint64_t`)
especificamente para que a multiplicação nunca precise de um tipo
intermediário de 128 bits - o MSVC não tem `__int128` portátil, e todo
produto `uint32 * uint32` cabe exatamente num `uint64_t`.

`aether::mesh::RobustPredicates` usa isso para `orientation3D()` e
`inSphere3D()`. A exatidão vem de um fato simples: todo `double` IEEE754 é
*exatamente* `significando * 2^expoente` (um racional diádico) - nenhuma
aproximação é necessária para converter um, só para computar com ele em
ponto flutuante depois. Cada coordenada vira esse par exato, e toda
operação seguinte (subtração, multiplicação, adição) roda sobre essa
representação, então nenhum arredondamento reentra na conta. Só o **sinal**
final é extraído, que é tudo que os dois predicados precisam.
`DelaunayTetrahedralization3D` agora usa ambos em todos os seus testes de
orientação e circumsfera - inclusive permitindo trocar uma checagem de
degeneração baseada em tolerância (`|volume| < 1e-12`) por uma comparação
com zero exato, sem limiar arbitrário para adivinhar.

**Estendido também para 2D**: `orientation2D()` e `inCircle2D()` (mesma
abordagem, uma dimensão abaixo), agora usados por `DelaunayTriangulation2D`
e `PolygonTriangulation2D` - este último tinha *três* testes de sinal em
ponto flutuante (direção de curva no ear clipping, ponto-em-triângulo, e
in-circumcircle nos flips), todos são testes de **sinal** puro, então
exatidão é precisamente o que importa e nenhuma magnitude é necessária.

**O mesmo bug do super-triângulo estava vivo em 2D** - e desta vez foi
medido antes de afirmar, não presumido. Fiz um A/B direto: com a escala
antiga (`20*deltaMax`), 200 conjuntos de pontos aleatórios comparados
contra a área do casco convexo (calculada independentemente via
`ConvexHull` do qhull - uma triangulação precisa preencher exatamente seu
casco) falharam **10 vezes (190/200)**; com a escala corrigida, **200/200**.
Adicionado teste C++ de regressão com um dos casos que falhavam,
verificando cobertura de área, não só validade - e vale repetir que
`satisfiesDelaunayProperty()` **passava** na saída defeituosa, porque
circumcírculo-vazio não diz nada sobre cobertura.

Todos os 5 suites de teste C++ passam (`ctest -C Debug`, 5/5).

O solver físico oficial segue o método dos Volumes Finitos (FVM), conforme
o roadmap do projeto; os protótipos em `python/research/lbm_prototypes/`
foram uma etapa de validação/experimentação e não são estendidos.

`TaylorGreenVortexSolver2D` é deliberadamente escopado a domínio duplamente
periódico (sem paredes sólidas) para poder validar contra uma solução exata.
`LidDrivenCavitySolver2D` adiciona paredes sólidas (a cavidade com tampa
deslizante clássica) reaproveitando o mesmo método - ambos usam grade
colocalizada (não staggered), o que é mais simples de implementar mas pode
gerar oscilação tipo "checkerboard" sem interpolação de Rhie-Chow (não
implementada); a divergência pós-projeção fica em torno de ~4e-3 (Taylor-
Green) a ~0.05-0.19 (cavidade, que tem gradientes de parede mais fortes) em
vez de zero de máquina - um artefato de discretização conhecido e
documentado, não um bug. A validação da cavidade evita deliberadamente
comparar com tabelas de benchmark da literatura (ex. Ghia et al. 1982), que
teriam que ser citadas de memória - em vez disso valida propriedades
deriváveis: repouso exato com tampa parada, conservação de massa, e a
reversão de escoamento que a conservação de massa força numa caixa fechada.

`StaggeredNavierStokesSolver3D` estende Navier-Stokes para 3D **e** resolve
de vez o problema do checkerboard acima: em vez de grade colocalizada, usa
o arranjo staggered clássico (MAC - Marker-and-Cell), com pressão no centro
da célula e cada componente de velocidade na sua própria face (u nas faces
x, v nas faces y, w nas faces z). Isso torna o gradiente de pressão e o
divergente operadores adjuntos discretos exatos (cada um usa só as duas
células/faces diretamente vizinhas, não uma diferença central de 2h) - a
composição dos dois vira exatamente o Laplaciano compacto de 7 pontos
resolvido para a pressão, sem o gap de checkerboard que a grade colocalizada
tem. Medido diretamente: divergência pós-projeção ~9.4e-13 (praticamente
precisão de máquina), não os ~4e-3 a ~0.05-0.19 das versões colocalizadas.
O termo de convecção precisa de mais interpolações (os componentes cruzados,
ex. v e w interpolados até a posição de uma face de u, para montar os termos
de fluxo uv/uw) - mais contabilidade que o esquema colocalizado, mas não
conceitualmente novo.

Ainda restrito a domínio triplamente periódico (sem paredes sólidas), para
poder validar contra uma solução exata genuinamente 3D: um escoamento de
Beltrami (`rot(u) = u`), o modo ABC único (`u0=sin(z)+cos(y)`,
`v0=sin(x)+cos(z)`, `w0=sin(y)+cos(x)`). Como `rot(u0)=u0`, a identidade
vetorial `(u.grad)u = grad(|u|^2/2) - u x rot(u)` faz o termo advectivo não-
linear colapsar num gradiente puro (`u x u = 0`), absorvível numa pressão
modificada - então `u(t) = u0 * exp(-nu*t)` é solução exata das equações de
Navier-Stokes incompressíveis completas (cada componente cartesiano também
satisfaz `laplaciano(u0) = -u0` por diferenciação direta, batendo com a
mesma taxa de decaimento no lado difusivo). Essa derivação é autocontida
(identidade vetorial padrão + verificação direta por substituição), não um
dado de benchmark citado de memória, e genuinamente 3D - diferente de embutir
um escoamento 2D invariante em z numa grade 3D, que nunca exercitaria os
novos termos cruzados de fato. Um bug real de indexação foi encontrado e
corrigido durante a validação: os índices lineares já achatados (`im`, `jm`,
`km`, usados para acessar vizinhos do mesmo campo ao longo do próprio eixo)
estavam sendo reusados como se fossem coordenadas escalares individuais ao
montar índices para os termos cruzados entre campos (ex. `v` interpolado
numa posição de `u`) - um erro sutil que só se manifestava como segfault
dependente do layout de memória, não de forma consistente, o que tornou o
diagnóstico mais difícil que os bugs de sinal encontrados nos solvers
anteriores.

`StaggeredLidDrivenCavitySolver3D` estende isso pra paredes sólidas nas seis
faces + tampa deslizante no topo (a cavidade 3D), seguindo o mesmo padrão já
usado em 2D (`TaylorGreenVortexSolver2D` periódico primeiro,
`LidDrivenCavitySolver2D` com paredes depois). Diferença estrutural chave em
relação à versão periódica: como uma parede sólida tem uma face física
genuína (diferente do wraparound periódico, onde a face `nx` coincide com a
face `0`), cada componente de velocidade agora tem um tamanho de array
diferente - `u` cobre `nx+1` posições ao longo do próprio eixo (as duas
pontas são os valores de contorno diretamente, sempre 0, nunca resolvidas),
o mesmo para `v` e `w` nos seus eixos. As direções *tangenciais* de cada
componente (onde ele não tem face própria) usam mirror-Dirichlet fantasma
(`ghost = 2*valor_parede - interior`, a mesma técnica já validada em
`LidDrivenCavitySolver2D`) - só o mirror de `u` na parede do topo carrega o
valor não-nulo da tampa (`lidVelocity * sin²(pi*x/Lx) * sin²(pi*y/Ly)`,
afunilando a zero nas quatro arestas do topo, a generalização 3D do
afunilamento nos dois cantos da versão 2D); todo outro mirror em toda a
classe é homogêneo (zero). Cada uma das três equações de momento foi
deduzida reaproveitando exatamente a mesma derivação já validada da versão
periódica (mesmo padrão cíclico x→y→z, u→v→w), só trocando os acessos
`wrap()`-periódicos por acessores que decidem entre índice direto (eixo
próprio) e mirror fantasma (eixos tangenciais) - o que deu confiança
adicional na correção, já que a álgebra em si já tinha sido validada antes.
Funcionou corretamente checado via Python antes mesmo de escrever o teste
C++: repouso exato com tampa parada, e com a tampa em movimento (Re=10),
divergência pós-projeção em ~2.5e-14 - a mesma precisão de máquina da versão
periódica, confirmando que a propriedade de adjunto discreto exato da malha
staggered carrega corretamente pro caso de parede sólida - e a mesma
topologia de vórtice primário já validada em 2D (camada superior arrastada
na direção da tampa, camada inferior revertida pela conservação de massa).

`MixingLengthChannelFlowSolver1D` (Módulo 6) usa o mesmo problema 1D de
canal já validado no Poiseuille laminar, agora com viscosidade turbulenta
`nu_t = l_m^2 |du/dy|` (comprimento de mistura de Prandtl, `l_m = kappa*y`
perto da parede, limitado a `0.09*(altura/2)` longe dela). A velocidade de
atrito `u_tau` vem do balanço de momento exato da equação (`sqrt(fonte *
altura/2)`), não de uma estimativa por diferença finita perto da parede -
uma versão inicial usava essa estimativa com viscosidade só molecular e
divergia do valor exato em ~38%, porque o solver de fato usa uma
viscosidade efetiva (molecular + turbulenta da célula vizinha) nessa mesma
face. A validação cobre só a inclinação da lei da parede logarítmica
(1/kappa - praticamente garantida pela própria construção do modelo), não a
constante aditiva B nem qualquer tabela de benchmark de canal turbulento -
ambas exigiriam valores citados de memória.

`KEpsilonChannelFlowSolver1D` (Módulo 6, segunda etapa) - k-epsilon padrão
(C_mu=0.09, C_eps1=1.44, C_eps2=1.92) com funções de parede de equilíbrio
(malha não resolve até a parede; primeira célula fica na região log-law,
y+~30-100). Três bugs reais apareceram e foram corrigidos durante o
desenvolvimento: (1) partida fria - com velocidade inicial zero a produção
de turbulência começa em zero e k colapsa antes do campo de velocidade se
desenvolver, corrigido com "warm-start" a partir do próprio
`MixingLengthChannelFlowSolver1D`; (2) o termo de destruição quadrático de
epsilon (`-C_eps2*epsilon^2/k`) tratado explicitamente é instável (um
problema clássico e documentado desse tipo de termo), corrigido com
linearização semi-implícita (regra de Patankar: termos de sumidouro sempre
no denominador); (3) condições de fluxo (Neumann) nas duas paredes da
equação de momento deixam a solução sem referência absoluta (mesmo
problema de espaço nulo já visto na pressão do Taylor-Green/cavidade, só
que aqui a variável tem valor físico real, não arbitrário) - corrigido
fixando a célula adjacente à parede pela própria lei logarítmica (usando
`u_tau` exato + a constante B=5.0 padrão, o único valor deste projeto de
turbulência citado de memória em vez de derivado/medido - só ancora o
nível absoluto do perfil, a inclinação validada não depende de B estar
exatamente certo).

`KOmegaSSTChannelFlowSolver1D` (Módulo 6, terceira etapa) - modelo k-omega
SST (Shear Stress Transport) de Menter, o fechamento de turbulência mais
sofisticado do projeto até agora: mistura dois conjuntos de constantes via
uma função de blending F1 dependente da distância à parede (conjunto 1,
k-omega de Wilcox, perto da parede; conjunto 2, derivado de k-epsilon,
longe dela), mais uma segunda função F2 que ativa um limitador de
tensão de cisalhamento (Bradshaw) na viscosidade turbulenta
(`nu_t = a1*k / max(a1*omega, |dU/dy|*F2)`, a1=0.31) e um termo de
difusão cruzada (`2*(1-F1)*sigma_w2/omega * dk/dy*domega/dy`) na equação
de omega. Essas constantes específicas (Menter 1994/2003) são valores
padrão da literatura citados de memória - como a constante B=5.0 da lei
logarítmica já usada pelo k-epsilon, isso não é algo derivável a partir de
primeiros princípios, e é o principal ponto desta classe que não pode ser
verificado de forma independente como a maior parte dos solvers deste
projeto. Por isso a validação segue a mesma filosofia: checa a inclinação
da lei da parede logarítmica e a simetria do perfil (auto-deriváveis), não
uma tabela de benchmark.

Mesmas condições de parede por função de equilíbrio do k-epsilon (`k_wall
= u_tau^2/sqrt(beta*)`, `omega_wall = u_tau/(kappa*y_wall*sqrt(beta*))` -
esta última derivada da consistência com `epsilon = beta**k*omega` e a
própria fórmula de `eps_wall` do k-epsilon, não um valor memorizado à
parte) e mesmo ancoramento de `u` pela lei log (mesmo problema de espaço
nulo do k-epsilon). Um bug real foi encontrado durante o desenvolvimento:
a viscosidade turbulenta nas duas células de parede nunca era atualizada
(a função que calcula o blending só cobre células interiores, já que
precisa do gradiente de velocidade), ficando presa no valor zero inicial
e usando só viscosidade molecular no fluxo da célula adjacente à parede -
medido diretamente: isso deixava a inclinação medida da lei log em 2.99
(22.5% de erro) em vez do valor esperado (1/kappa=2.439). Corrigido
fixando essas duas células ao valor de equilíbrio
`nu_t = kappa*u_tau*y_wall` (a mesma fórmula que o modelo de comprimento
de mistura e o próprio k-epsilon produzem nesse regime, sob a mesma
suposição de equilíbrio já usada para `k_wall`/`omega_wall`) em vez de
avaliar a fórmula do SST ali - o que reduziu o erro pra ~8.6% (2.65 medido),
comparável ao próprio k-epsilon (~6-9% no mesmo canal). Mesma prática de
linearização de Patankar para os dois termos de destruição bilineares
(`beta**k*omega` na equação de k, `beta*omega^2` na de omega) e mesmo
warm-start de `u` a partir do comprimento de mistura, pelas mesmas razões
já documentadas para o k-epsilon.

`MixingLengthLidDrivenCavitySolver2D` (Módulo 6, quarta etapa) - o primeiro
fechamento de turbulência acoplado a um escoamento 2D de verdade: todos os
três fechamentos acima (mistura, k-epsilon, k-omega SST) só viram o
problema 1D de canal totalmente desenvolvido, onde a convecção desaparece
por construção. Estende `LidDrivenCavitySolver2D` (mesmo tratamento de
célula-fantasma, mesma tampa regularizada, mesma projeção de Chorin) com
`nu_t(x,y) = l_m(x,y)^2 * S(x,y)`, onde `S` é a taxa de deformação 2D
completa (`sqrt(2*(du/dx)^2 + 2*(dv/dy)^2 + (du/dy+dv/dx)^2)`, não só
`|du/dy|` como no canal 1D - a cavidade tem deformação normal real além do
cisalhamento) e `l_m` usa a distância à parede mais próxima entre as
quatro (em vez de uma única coordenada normal à parede). `nu_t` entra no
termo de difusão do preditor como `(nu+nu_t)` médio em cada face, aplicado
independentemente por componente de velocidade - uma simplificação
("forma quase-laminar") que omite os termos cruzados que uma divergência
rigorosa do tensor de tensões com viscosidade variável incluiria; comum em
códigos RANS introdutórios/de engenharia e adequada como primeiro passo,
não uma implementação completa. `nu_t` é exatamente zero em todas as
quatro paredes por construção (comprimento de mistura zero lá).

Validado sem tabela de benchmark (mesma filosofia da cavidade laminar):
com a tampa parada, `nu_t` deve ser exatamente zero em toda a malha (não
há deformação nenhuma) - checado de forma exata, não aproximada. Com a
tampa em movimento: divergência limitada (mesmo estilo de checagem de
regressão da cavidade laminar, mesma limitação de grade colocalizada sem
Rhie-Chow), a mesma topologia de vórtice primário já validada na cavidade
laminar, e uma checagem estrutural de `nu_t` (não-negativo em toda parte,
maior no centro da cavidade que numa célula adjacente à parede). Não
resolve um regime realmente turbulento (exigiria malha muito mais fina e
Reynolds muito mais alto do que o validado aqui) - isso é uma limitação
reconhecida, não uma alegação de precisão nesse regime.

`KEpsilonLidDrivenCavitySolver2D` (Módulo 6, quinta etapa) estende o
fechamento de duas equações (antes só acoplado ao canal 1D) pra essa mesma
cavidade real - um passo materialmente mais difícil que o do comprimento
de mistura: `k` e `epsilon` agora são campos transportados (suas próprias
equações de advecção-difusão-produção-destruição, avançadas explicitamente
a cada passo junto com o momento, não resolvidas até convergência por um
laço de Picard externo como no canal 1D). **A escolha de projeto central
foi o tratamento de parede**, deliberadamente não funções de parede de
equilíbrio (o que o `KEpsilonChannelFlowSolver1D` usa): funções de parede
precisam da velocidade de atrito local `u_tau`, que no canal 1D vinha de um
balanço de forças *exato* específico de escoamento totalmente desenvolvido
- não existe uma relação exata assim pra um escoamento 2D recirculante em
geral, e estimar `u_tau` localmente a partir de um gradiente perto da
parede numa malha não refinada especificamente pra isso seria exatamente o
tipo de aproximação não validada que este projeto evita. Em vez disso, usa
a relação assintótica de baixo-Reynolds pra epsilon, derivada do fato de
que `k ~ y²` perto de qualquer parede sólida (verdade geral, não só pra
escoamentos em equilíbrio/lei-log): `epsilon_parede = 2*nu*k_celula/y_meio²`,
onde `y_meio = h/2` é a distância exata e conhecida do centro da célula até
a parede sob a convenção de mirror-fantasma desta classe. `k` é fixado
exatamente em 0 em toda parede (energia cinética turbulenta se anula numa
parede sem escorregamento por definição); `Cmu` não é amortecido perto da
parede (sem correção `f_mu` de baixo-Re tipo Launder-Sharma) - mais uma
simplificação deliberada, já que isso puxaria mais funções empíricas
citadas de memória.

O termo de destruição quadrático de epsilon (`-Ceps2*epsilon²/k`) usa a
mesma linearização de Patankar já necessária no k-epsilon 1D (denominador
`1 + dt*Ceps2*epsilon/k`, usando o `epsilon`/`k` da iteração anterior como
fator congelado), aplicada preventivamente desta vez em vez de descoberta
via colapso numérico - junto com o warm-start de velocidade a partir do
`MixingLengthLidDrivenCavitySolver2D` (mesmo problema de partida fria do
k-epsilon 1D: `u=0` dá produção zero, `k`/`epsilon` colapsariam antes do
campo de velocidade se desenvolver), **isso evitou reproduzir qualquer um
dos três bugs reais que tornaram o k-epsilon 1D "o arco de depuração mais
difícil da sessão"** - testado incrementalmente via Python em várias
resoluções e números de Reynolds (100 e 1000) sem NaN, sem explosão
numérica, com divergência limitada durante toda a execução.

Validado com a mesma filosofia sem-benchmark: com a tampa parada, `u`/`v`
ficam exatamente zero (produção e difusão do momento se anulam
independentemente de `k`/`epsilon`, já que multiplicam um campo de
velocidade uniformemente nulo) - `k`/`epsilon` NÃO são checados como
permanecendo no valor inicial (ao contrário do `nu_t` do comprimento de
mistura): com `k` fixado em 0 na parede mas inicializado num valor
uniforme não-nulo no interior, a difusão sozinha já evolui `k`/`epsilon`
em direção à parede mesmo em repouso - comportamento real e esperado de um
campo transportado, não uma regressão pra checar bit a bit. Com a tampa em
movimento (Re=100): conservação de massa, mesma topologia de vórtice, e
checagem estrutural de `k`/`epsilon`/`nu_t` (não-negativos em toda parte,
menores numa célula adjacente à parede que no centro). Medido diretamente
(não suposto): tanto `k` quanto `nu_t` **decaem** ao longo da execução
nesses números de Reynolds em vez de atingir um nível turbulento
estatisticamente estacionário - o comportamento fisicamente correto pra
k-epsilon nesse regime (a cavidade com tampa deslizante só transiciona pra
turbulência de verdade em Re da ordem de 10⁴, muito além do que qualquer
solver deste projeto já rodou) - não é um bug, e este teste não alega o
contrário.

`KOmegaSSTLidDrivenCavitySolver2D` estende a mesma cavidade real pro
k-omega SST, arquitetura quase idêntica ao `KEpsilonLidDrivenCavitySolver2D`
(mesmas paredes-fantasma, mesmo warm-start, mesma linearização de Patankar
pros termos de destruição bilineares, agora `beta*·k·omega` e `beta·omega²`).
A condição de contorno de `omega` perto da parede também é autoderivada: em
vez de recorrer à fórmula de baixo-Re de Wilcox (citada de memória),
substitui a mesma relação `epsilon_parede = 2·nu·k/y_meio²` já usada no
k-epsilon na identidade exata `epsilon = beta*·k·omega` e toma o limite
`y→0` (onde `k ~ a·y²` cancela a dependência em `k`), dando
`omega_parede = 2·nu/(beta*·y_meio²)` - uma forma autoconsistente com a
derivação do epsilon já validada, mas **não** alegada como batendo
exatamente com a fórmula mais cuidadosa de Wilcox (comumente citada com um
fator 6 em vez de 2, vinda de uma análise assintótica direta da própria
equação de omega) - uma simplificação documentada, não escondida.

**Um bug real foi encontrado e corrigido durante o desenvolvimento**: a
condição inicial de `omega` usava `omega0 = eps0/(beta*·k0)`, que quando
`lidVelocity=0` (o teste de repouso) dá `k0=0` exatamente, logo `0/0 = NaN`
- toda a malha de omega nascia `NaN` na própria construção, e o teste de
repouso falhava com `u`/`v` virando `NaN` no primeiro passo (o `nu_t`
corrompido por `NaN` se propaga pro termo de difusão do momento mesmo
multiplicando uma diferença de velocidade identicamente zero: `0*NaN=NaN`,
não zero). Corrigido simplificando algebricamente a mesma fórmula pra
`omega0 = sqrt(k0)/(L·beta*^0.25)` (substituindo `eps0 = beta*^0.75·k0^1.5/L`
na razão e cancelando o `k0` do denominador antes de calcular, em vez de
computar `eps0` e depois dividir por `k0` em tempo de execução) - matematicamente
idêntico quando `k0>0`, mas sem a divisão por zero quando `k0=0`. Achado
rapidamente ao rodar o mesmo teste de repouso já usado no k-epsilon.

Fora esse bug (já corrigido), testado incrementalmente em Re=100 na mesma
grade do k-epsilon: sem explosão numérica, e os números batem de forma
consistente entre os dois fechamentos independentes (topologia de vórtice
top/bottom quase idêntica, mesma ordem de grandeza de divergência) - uma
checagem cruzada tranquilizadora, já que ambos são implementações
independentes de modelos diferentes acoplados ao mesmo solver de momento
subjacente.

## Mixing-length acoplado à cavidade 3D real (primeira turbulência genuinamente 3D)

`MixingLengthLidDrivenCavitySolver3D` fecha o thread que ficou explicitamente
em aberto desde os fechamentos de turbulência 2D acima: até aqui, todo
fechamento que já viu convecção real (`MixingLengthLidDrivenCavitySolver2D`,
`KEpsilonLidDrivenCavitySolver2D`, `KOmegaSSTLidDrivenCavitySolver2D`) estava
confinado a uma malha 2D colocalizada plana. Esta classe acopla o mesmo
fechamento de comprimento de mistura de Prandtl à cavidade 3D em malha
staggered/MAC (`StaggeredLidDrivenCavitySolver3D`), a primeira turbulência
deste projeto num escoamento genuinamente 3D.

`nu_t = l_m² · |S|`, com `l_m = min(kappa·d_parede, 0.09·min(Lx,Ly,Lz)/2)`
(mesma fórmula de von-Kármán-perto-da-parede/cap de Escudier de todo
fechamento de comprimento de mistura deste projeto, agora com `d_parede` =
distância à mais próxima das *seis* paredes) e `|S| = sqrt(2·Sij·Sij)` o
tensor de taxa de deformação 3D completo (6 componentes independentes, a
generalização direta da fórmula 2D que só tinha 3). Cada componente da taxa
de deformação é calculado em pontos cell-centered (onde a pressão e agora
`nu_t` também vivem) combinando a derivada ao longo do próprio eixo (diferença
exata entre as duas faces que flanqueiam a célula, ex. `dudx` a partir das
duas faces-x de `u`) com as derivadas cruzadas (ex. `dudy`: média de duas
diferenças centrais em y, uma avaliada em cada face-x da célula) - reaproveita
exatamente o mesmo padrão de média em arestas que os próprios termos
convectivos cruzados do `step()` (`duvdy`, `dvudx`, ...) já usam nessa malha
staggered, em vez de inventar um esquema de interpolação novo.

**Simplificação deliberada, documentada**: a viscosidade efetiva usada na
difusão de cada componente é *exata* (varia por face, usando `nu_t` das duas
células de pressão vizinhas) ao longo do eixo próprio do componente, mas usa
um único valor amostrado localmente (média das mesmas duas células) nos dois
eixos transversais - uma simplificação adicional sobre a já simplificada
"forma quase-laminar" que `MixingLengthLidDrivenCavitySolver2D` documenta
(que por sua vez já omite os termos cruzados completos do tensor de
viscosidade variável). Um tratamento totalmente rigoroso precisaria resolver
`nu_t` em cada aresta da malha staggered (12 posições distintas por célula em
3D, cada uma com sua própria média de 4 células) - complexidade real por um
ganho de precisão que a filosofia de validação deste projeto (checagens
estruturais/de primeiros princípios, não batimento contra tabela de
literatura) não reivindica nem precisa.

**Validação em Python antes de qualquer teste C++** (prática estabelecida):
com a tampa parada, velocidade fica em zero exato (mesma checagem forte de
todo teste de repouso deste projeto) - e aqui a garantia é ainda mais forte
que nos fechamentos de duas equações (k-epsilon/k-omega SST), já que `nu_t`
é puramente algébrico (função da taxa de deformação local, sem equação de
transporte própria): com velocidade em zero exato em toda parte, a taxa de
deformação também é zero exata em toda parte, então `nu_t` fica em zero exato
- não apenas aproximadamente, medido: `max|u|=max|v|=max|w|=max(nu_t)=0.0`
depois de 20 passos. Com a tampa em movimento (Re=100, malha 12³, 400
passos): divergência máxima medida em `~2.9e-13` - mesma ordem de grandeza
de precisão de máquina que `StaggeredLidDrivenCavitySolver3D` puro já
alcançava (confirma que a vantagem da malha staggered sobre os fechamentos
2D colocalizados, que ficam em `~0.14-0.19`, sobrevive à adição do termo de
difusão turbulenta); topologia do vórtice primário preservada (camada
superior arrastada a `~0.151` na direção da tampa, camada inferior revertida
a `~-0.0049` por conservação de massa); `nu_t` não-negativo em toda parte e
maior no centro da cavidade (`~9.1e-4`) que numa célula de canto adjacente a
três paredes (`~1.6e-6`) - cresce com a distância à parede, como esperado.

Todos os 5 suites de teste C++ passam (`ctest -C Debug`, 5/5) na primeira
tentativa - nenhum bug novo neste incremento. Vinculado ao Python como
`aether.MixingLengthLidDrivenCavitySolver3D`.

**Ainda em aberto**: k-epsilon e k-omega SST acoplados à cavidade 3D
(análogo às duas etapas seguintes que já existem em 2D) ficam como próximo
passo natural, seguindo o mesmo padrão estabelecido de provar o fechamento
mais simples primeiro.

## k-epsilon acoplado à cavidade 3D real

`KEpsilonLidDrivenCavitySolver3D` dá o próximo passo imediato: o fechamento
de duas equações (antes só acoplado ao canal 1D e à cavidade 2D real)
estendido pra mesma cavidade 3D em malha staggered que acabou de receber o
comprimento de mistura. Arquitetura quase idêntica à versão 2D: mesmo
tratamento de parede (`epsilon_parede = 2·nu·k/y_meio²`, a mesma relação
assintótica auto-derivada, agora aplicada em qualquer um dos três eixos
conforme qual parede é cruzada), `k` fixado em 0 nas seis paredes, mesmo
warm-start (agora a partir de `MixingLengthLidDrivenCavitySolver3D`) e mesma
linearização de Patankar pro termo de destruição quadrático de epsilon.

**Uma simplificação real que a extensão pra 3D expôs**: diferente das
equações de momento (staggered, precisando do tratamento "eixo
próprio exato / eixos transversais aproximados" que
`MixingLengthLidDrivenCavitySolver3D` documenta), `k` e `epsilon` vivem nos
mesmos pontos cell-centered que a pressão e `nu_t` - então o transporte
deles usa um estêncil colocalizado comum de 6 vizinhos (2 por eixo), a
generalização direta e sem complicação extra do estêncil de 4 vizinhos da
versão 2D. Nenhuma simplificação adicional foi necessária aqui.

**Validado em Python antes de qualquer teste C++**: com a tampa parada,
velocidade em zero exato (mesma checagem forte de sempre); `k`/`epsilon`
não são checados contra o valor inicial (mesmo raciocínio da versão 2D -
com `k` fixado em 0 na parede mas não-nulo no interior, a difusão sozinha já
evolui os campos mesmo em repouso). Com a tampa em movimento (Re=100, malha
12³, 600 passos): divergência máxima medida em `~4.1e-11` - ainda muitas
ordens de grandeza abaixo do `~0.14-0.19` das versões 2D colocalizadas,
embora um pouco acima do `~2.9e-13` do comprimento de mistura puro (`k`/
`epsilon` realimentam `nu_t` a cada passo, adicionando uma fonte extra de
variação); topologia do vórtice preservada (camada superior a `~0.148`,
inferior a `~-0.0045`); `k`, `epsilon` e `nu_t` não-negativos em toda parte
e maiores no centro que num canto adjacente a três paredes. Como na versão
2D, `k` e `nu_t` **decaem** ao longo da execução nesse Reynolds (mesmo
comportamento fisicamente correto documentado lá - a cavidade só transiciona
pra turbulência de verdade em Re ~ 10⁴).

Nenhum bug novo neste incremento. Todos os 5 suites de teste C++ passam
(`ctest -C Debug`, 5/5) na primeira tentativa. Vinculado ao Python como
`aether.KEpsilonLidDrivenCavitySolver3D`.

## k-omega SST acoplado à cavidade 3D real (fecha a progressão de turbulência 3D)

`KOmegaSSTLidDrivenCavitySolver3D` fecha a última etapa natural: SST
estendido da cavidade 2D real pra mesma cavidade 3D staggered, arquitetura
quase idêntica à `KEpsilonLidDrivenCavitySolver3D` (mesmo tratamento de
parede pra `omega` via `omega_parede = 2·nu/(beta*·y_meio²)`, mesmo
warm-start a partir do `MixingLengthLidDrivenCavitySolver3D`, mesma
linearização de Patankar).

**A única física genuinamente nova que o 3D expôs**: vorticidade é um vetor,
não um escalar. O limitador de Bradshaw da versão 2D usa `|dv/dx - du/dy|`,
o componente z do rotacional - o único componente não-nulo possível num
escoamento genuinamente 2D. Em 3D, `curl(u) = (dw/dy-dv/dz, du/dz-dw/dx,
dv/dx-du/dy)` tem os três componentes em geral, então o limitador usa a
magnitude completa `sqrt(wx²+wy²+wz²)` - calculada a partir das mesmas seis
derivadas cruzadas (`dudy`, `dvdx`, `dudz`, `dwdx`, `dvdz`, `dwdy`) já
necessárias pro tensor de deformação, reaproveitando o mesmo padrão de média
em arestas mais uma vez, agora combinando os mesmos seis números de forma
antissimétrica (vorticidade) em vez de só simétrica (deformação).

**Validado em Python antes de qualquer teste C++, incluindo o cenário que
causou o bug na versão 2D**: diferente do desenvolvimento original do SST
2D (que teve um `0/0 = NaN` real inicializando `omega0 = eps0/(beta*·k0)`
com `k0=0`), esta versão 3D já nasceu com a fórmula algebricamente
simplificada `omega0 = sqrt(k0)/(L·beta*^0.25)` desde o início - o teste de
repouso passou limpo na primeira tentativa, sem NaN em lugar nenhum. Com a
tampa em movimento (Re=100, malha 12³, 600 passos): divergência máxima
medida em `~9.9e-13` - ainda mais perto da precisão de máquina que o próprio
`KEpsilonLidDrivenCavitySolver3D` (`~4.1e-11`); topologia de vórtice
praticamente idêntica ao k-epsilon 3D na mesma configuração (`~0.148`/
`~-0.0045` pra ambos) - a mesma checagem cruzada tranquilizadora entre
fechamentos independentes já documentada em 2D; `k`, `omega`, `nu_t`
não-negativos em toda parte e maiores no centro que num canto.

**Nenhum bug encontrado neste incremento** - segunda implementação limpa
seguida (depois do k-epsilon 3D), reforçando que aplicar preventivamente as
lições já aprendidas (aqui, a simplificação algébrica do `omega0` desde o
primeiro rascunho) evita os mesmos problemas que apareceram nas versões
originais 1D/2D.

Todos os 5 suites de teste C++ passam (`ctest -C Debug`, 5/5) na primeira
tentativa. Vinculado ao Python como `aether.KOmegaSSTLidDrivenCavitySolver3D`.

**Isso fecha por completo a progressão de turbulência 3D**: mixing-length,
k-epsilon e k-omega SST agora todos acoplados à cavidade 3D real, espelhando
exatamente a mesma progressão de três etapas que já existia em 2D. A
turbulência-convecção 3D deixa de ser um item em aberto do roadmap.

## Módulo 8: viewer 3D para os solvers da cavidade staggered

Com quatro solvers 3D novos construídos nesta sessão (`StaggeredLidDrivenCavitySolver3D`
e suas três extensões de turbulência) e nenhum viewer capaz de mostrá-los -
só diagnóstico de console -, `aether_unified_viewer` ganhou um quinto modo,
`cavity3d`: renderiza `StaggeredLidDrivenCavitySolver3D` (a base
compartilhada por todo fechamento 3D) como um campo de vetores de
velocidade 3D, reaproveitando o esqueleto de câmera orbital/perspectiva já
existente no modo `mesh`.

**Decisão deliberada: um modo só, não um por fechamento.** Os quatro
solvers 3D (staggered puro, mixing-length, k-epsilon, k-omega SST)
resolvem o mesmo campo de velocidade com a mesma topologia de vórtice -
renderizá-los seria visualmente quase idêntico. Os campos realmente
específicos de cada fechamento (`nu_t`, `k`, `epsilon`/`omega`) são volumes
escalares, mais bem servidos por um renderizador de iso-superfícies/cortes
do que por setas - o que exigiria marching cubes 3D, que ainda não existe
(`engine/postprocessing` só tem o marching squares 2D, já documentado como
"o primeiro passo rumo ao marching cubes 3D"). Construir esse renderizador
pela metade agora seria pior que não construir - fica como trabalho futuro
bem definido em vez de forçado nesta etapa.

**Simplificação documentada: setas sem ponta.** Uma seta 3D de verdade
precisa de um referencial perpendicular ao eixo do vetor pra desenhar a
ponta, e "perpendicular a uma linha 3D" não tem escolha única (ao contrário
do caso 2D do `cavity_mode`, onde só existe uma perpendicular no plano) -
resolver isso corretamente exigiria escolher/normalizar um referencial
arbitrário por seta. Como a cor (por velocidade) e a direção da própria
linha já comunicam o campo, as setas aqui são só o eixo (sem ponta) -
simplificação real, documentada no código, não uma omissão escondida.

Validado rodando o modo (malha 10³, Re=10, 400 passos) e capturando uma
screenshot real da janela (técnica `PrintWindow`, a mesma usada em toda a
sessão pra validar os viewers): o cubo de arame aparece corretamente, e as
setas mostram claramente velocidades maiores (vermelho/branco) perto da
tampa deslizante no topo e menores (azul) no resto da cavidade - exatamente
o padrão físico esperado, confirmando visualmente o que os testes numéricos
já validavam sem imagem nenhuma.

## LES (Smagorinsky): sai de RANS, entra em simulação de grandes escalas

`SmagorinskyLesLidDrivenCavitySolver3D` é a etapa seguinte do Módulo 6
depois do k-ω SST - e a primeira que **não é RANS**. A fórmula parece
enganosamente parecida com a do comprimento de mistura, então vale dizer
claramente o que muda:

- **RANS** (mistura, k-epsilon, k-ω SST) modela *toda* a turbulência. O
  comprimento característico vem da **geometria do escoamento** (distância
  à parede, ou uma escala transportada). Refinar a malha resolve melhor a
  mesma física modelada, mas não encolhe `nu_t` sistematicamente.
- **LES** (esta classe) resolve os grandes turbilhões diretamente e modela
  só o que cai abaixo do filtro da malha. O comprimento característico é a
  **largura de filtro** `Delta = cbrt(dx*dy*dz)`, atada à **malha**. Refinar
  encolhe `Delta`, encolhe `nu_sgs`, e no limite `Delta→0` tende a DNS.

`nu_sgs = l_s² · |S|`, com `l_s = min(Cs·Delta, kappa·d_parede)` e o mesmo
tensor de deformação 3D completo já derivado e validado no
`MixingLengthLidDrivenCavitySolver3D` (reaproveitado, não reinventado).

**Sobre o termo de parede, e por que não é van Driest**: Smagorinsky puro
superestima `nu_sgs` perto da parede. O remédio de manual é o amortecimento
de van Driest, `1 - exp(-y+/A+)`, que precisa de `y+` e portanto de `u_tau`
- que este projeto recusou repetidamente estimar em geometria genérica
(mesma decisão e mesmo raciocínio do `KEpsilonLidDrivenCavitySolver2D`).
Em vez disso reaproveita o cap `kappa·d_parede`: usa só a distância
exatamente conhecida às seis paredes, anula-se na parede como exigido, e é
o mesmo limite geométrico auto-justificável já validado aqui. É uma
substituição deliberada e documentada - **não** é van Driest e não alega
reproduzi-lo.

**Constante recitada, sinalizada honestamente**: `Cs = 0.17` é valor
empírico de literatura (citado entre 0.1 e 0.2, e genuinamente dependente
do escoamento - uma fraqueza conhecida do Smagorinsky puro, e a motivação
do procedimento dinâmico de Germano). É recitada, não derivada, e está
exposta como parâmetro do construtor em vez de enterrada no código.

**A validação que de fato distingue LES de RANS** - e que é *medível*, não
só afirmável. Mesmo problema físico (Re=100), mesmo tempo final, três
resoluções, rodando as duas classes:

| n | LES `nu_sgs` | RANS `nu_t` |
|---|---|---|
| 8 | 2.65e-3 | 3.81e-3 |
| 12 | 1.62e-3 | 5.66e-3 |
| 16 | 1.07e-3 | 6.26e-3 |

O LES **cai 2.5x** ao refinar (rumo a DNS); o RANS **sobe 1.6x** (a malha
mais fina resolve gradientes mais íngremes, aumentando `|S|`, enquanto seu
comprimento fixado pela geometria não encolhe). O teste C++ fixa as
*direções* (monotonicidade estrita no LES, e LES caindo enquanto RANS
sobe), não os números específicos - as direções são o que a física exige,
as magnitudes dependem de resolução e duração.

Demais validações: repouso em zero exato (inclusive `nu_sgs`, garantia
forte que só os modelos algébricos oferecem - não há equação de transporte
própria), divergência máxima em `~4.2e-13`, topologia de vórtice
preservada, `nu_sgs` não-negativo e menor num canto adjacente a três
paredes que no centro.

Nenhum bug neste incremento; os 5 suites C++ passam na primeira tentativa.
Vinculado ao Python como `aether.SmagorinskyLesLidDrivenCavitySolver3D`.

**Ainda em aberto no Módulo 6**: DES (precisa de LES + RANS combinados -
agora ambos existem em 3D, então destravou) e DNS.

## DES (Detached-Eddy Simulation): híbrido RANS/LES, fecha a lacuna que o LES deixou aberta

`DesSstLidDrivenCavitySolver3D` implementa a formulação clássica SST-DES
(Strelets 2001, a variante do DES97 original de Spalart baseada em k-ω SST)
- destravada assim que LES e RANS passaram a existir na mesma malha
staggered 3D. A ideia cabe numa frase: pega o `KOmegaSSTLidDrivenCavitySolver3D`
inteiro (mesmas equações de transporte, blending, tratamento de parede,
limitador de Bradshaw) e troca a escala de comprimento RANS implícita no
termo de destruição da equação de k, `L_RANS = sqrt(k)/(beta*·omega)`, por
`L_DES = min(L_RANS, C_DES·Delta)`, com `Delta = cbrt(dx·dy·dz)` - a mesma
largura de filtro do `SmagorinskyLesLidDrivenCavitySolver3D`.

Concretamente: destruição_k = beta*·k·omega = k^1.5/L_RANS normalmente.
Substituir `L_DES` por `L_RANS` equivale a multiplicar essa destruição por
`F_DES = max(L_RANS/(C_DES·Delta), 1.0)`. `F_DES=1` deixa o SST inalterado
(malha grossa demais pra valer a pena tentar resolver turbilhões
diretamente); `F_DES>1` empurra k (e portanto `nu_t`, pela mesma fórmula de
viscosidade turbulenta do SST) para valores estilo LES - o mesmo mecanismo
de "escala atada à malha" que define o LES, agora condicionado a qual das
duas escalas é menor, em vez de se aplicar sempre.

**`C_DES = 0.61` é constante de literatura recitada, sinalizada
honestamente** (valor de Strelets para SST-DES; o DES97 original de
Spalart-Allmaras usa 0.65 para uma escala de comprimento diferente, não
diretamente comparável). Exposta como parâmetro do construtor, mesma prática
do `Cs` do LES.

**Limitação conhecida, de manual, não corrigida aqui**: isto é DES97
clássico, não DDES (Delayed DES). Numa malha refinada na direção paralela à
parede mas que ainda não resolve a camada limite, `C_DES·Delta` pode ficar
menor que `L_RANS` *dentro* da camada limite, trocando prematuramente pro
modo LES ali - "modeled-stress depletion"/separação induzida pela malha, a
motivação documentada da função de blindagem do DDES. Não implementado
aqui: as malhas 3D deste projeto até agora são cartesianas uniformes (sem
estiramento na direção normal à parede), então a falha existe em princípio
mas não é exercitada pela validação usada (comparação de refino de malha
dominada pelo interior, não um caso de camada limite resolvida).

**Validação, toda medida, nada suposto**:
1. Repouso exato com a tampa parada (mesmo padrão de todo solver 3D aqui).
2. **Redução exata ao SST puro**: forçando `C_DES` gigante (1e6),
   `C_DES·Delta` nunca cai abaixo de `L_RANS` realista, então `F_DES=1.0`
   durante toda a trajetória (instrumentado, não suposto) - e com isso o DES
   reproduz o SST **bit a bit** (diferença máxima medida: exatamente `0.0`
   em `nu_t` e em `u`). Comparar com o `C_DES=0.61` padrão seria mais fraco:
   `F_DES` pode passar de 1.0 transitoriamente mesmo em casos que acabam
   assentando em `F_DES=1`, divergindo um pouco as duas trajetórias antes de
   reconvergir - efeito real, não bug, mas mais ruidoso de afirmar.
3. **A validação que de fato distingue DES de SST puro**: a Re=1000 (escoamento
   genuinamente turbulento, mesma escolha do teste de refino do LES),
   `F_DES` médio no núcleo do domínio cresce com o refino - medido, não
   suposto: n=4→2.01, n=6→2.89, n=8→3.63, estritamente crescente. Mesmo
   mecanismo "escala atada à malha" do LES, checado um nível acima, no fator
   de blending que o aciona dentro do fechamento RANS do DES.
4. Conservação de massa (`~5.8e-13`), topologia do vórtice primário, e
   não-negatividade de k/omega/`nu_t`/`F_DES` - mesmo estilo estrutural do
   teste do `KOmegaSSTLidDrivenCavitySolver3D`.

Nenhum bug neste incremento; os 5 suites C++ passam na primeira tentativa.
Vinculado ao Python como `aether.DesSstLidDrivenCavitySolver3D`
(`des_factor(i,j,k)`, `filter_width()`, `c_des` no construtor).

**Módulo 6 agora só tem DNS em aberto** - mistura → k-epsilon → k-ω SST →
LES → DES, todos em 3D real. DNS é o limite de resolução total (sem nenhum
modelo de turbulência), fora de alcance computacional prático nas malhas
usadas aqui - não está sendo perseguido como próximo passo natural, ao
contrário dos anteriores.

## CDT (tetraedralização restrita por contorno): fatia tratável da task #74

A generalidade completa - contorno não-convexo arbitrário, com garantia de
recuperação mesmo além da obstrução clássica do poliedro de Schönhardt - é
genuinamente nível de pesquisa, como o próprio cabeçalho da classe já
documentava. Em vez disso, `DelaunayTetrahedralization3D` ganhou a fatia
prática que geradores de malha de produção (TetGen etc.) também usam para
entrada bem-condicionada: a maioria dos facetas de contorno de um domínio
"razoável" já aparece na tetraedralização Delaunay *irrestrita* assim que
todo vértice de contorno entra como ponto comum; as poucas que não aparecem
costumam ser recuperáveis com um número limitado de pontos de Steiner.

Três métodos novos:
- `missingFacets(facets)`: quais facetas (dadas como três índices de ponto)
  ainda não são face de nenhum tetraedro - diagnóstico/teste, O(facetas ×
  tetraedros).
- `recoverFacets(facets, maxRounds=4)`: para cada faceta ausente, insere o
  centroide como ponto de Steiner e recursa nos três sub-triângulos
  resultantes, até `maxRounds` níveis. Retorna `recoveredFacets` (as
  faces elementares garantidamente presentes - a entrada original, sem
  necessidade de ajuda, mais qualquer subdivisão que precisou) e
  `unrecovered` (o que nem a recursão resolveu).
- `removeRegion(seed, walls)`: a partir do tetraedro que contém `seed`,
  remove por adjacência de face tudo que for alcançável sem atravessar
  nenhuma faceta em `walls` - a técnica clássica de "ponto no buraco" para
  esculpir uma cavidade interna, usando `recoverFacets()`'s own
  `recoveredFacets` como as paredes.

**Validação de ponta a ponta com forma escolhida deliberadamente**: octaedro
externo (R=2) com um buraco em forma de octaedro menor (r=0.5) esculpido no
meio, ambos centrados no mesmo ponto. Faces de octaedro são triângulos
genuínos, sem escolha de diagonal de quad - por quê isso importa, ver o
próximo parágrafo. Medido diretamente: as 16 facetas pedidas já vêm
presentes na tetraedralização irrestrita (0 ausentes), então
`recoverFacets()` é um no-op aqui (mas ainda testado como parte do
pipeline completo, não só o caso de sorte). `removeRegion()` esculpe o
buraco, `satisfiesDelaunayProperty()` continua válido depois (remoção só
apaga tetraedros, nunca adiciona pontos), e o volume somado bate
**exatamente** com a fórmula fechada do octaedro (bola L1),
`(4/3)·(R³-r³) = (4/3)·(8-0.125) = 10.5` - `diff` medido: `0.0`.

**Uma limitação real, honesta, encontrada testando - não escondida.** Um
cubo simples (8 vértices), pedindo os 12 triângulos de suas 6 faces
quadradas (cada face dividida por uma diagonal escolhida à mão): a
tetraedralização irrestrita é livre pra escolher *qualquer* diagonal de
cada face quadrada (ambas dão uma triangulação local igualmente válida e de
circunsfera vazia de um conjunto de 4 pontos exatamente coplanares) - medido
diretamente, faltaram exatamente 6 das 12 facetas pedidas (a outra
diagonal foi escolhida em 3 faces). A heurística de `recoverFacets()` **não
recupera essas**: o centroide de uma faceta ausente cai exatamente no mesmo
plano dos 4 cantos da face (por construção, já que a faceta é coplanar com
sua parceira de diagonal), então a inserção ou é rejeitada de imediato pela
própria guarda de degenerescência de `insertSteinerPoint()`, ou (observado
num caso relacionado com buraco durante o desenvolvimento) recursa por
vários sub-triângulos igualmente coplanares antes de ainda assim falhar.
Em vez de mascarar isso, o teste verifica exatamente o comportamento
honesto: `missingFacets` bate com os 6 medidos, e `recoverFacets` particiona
os 12 facetas de entrada inteiramente entre `recoveredFacets` (6) e
`unrecovered` (6) - nada perdido ou corrompido silenciosamente, mesmo
espírito do caminho de rejeição honesta que `insertSteinerPoint()` já
tinha para pontos exteriores degenerados.

Nenhum bug novo neste incremento (a limitação acima é comportamento
esperado da heurística escolhida, não um bug); os 5 suites C++ passam.
Vinculado ao Python (`missing_facets`, `recover_facets` retornando
`FacetRecoveryResult` com `recovered_facets`/`unrecovered`, `remove_region`).

**Ainda em aberto, pela razão certa agora**: recuperação combinatória de
faceta (busca de troca de diagonal antes de recorrer a pontos de Steiner,
a técnica que geradores de malha de produção usam primeiro) e o caso
genuinamente Schönhardt-difícil onde nenhuma subdivisão finita resolve sem
uma escolha de ponto mais cuidadosa; GMRES/BiCGSTAB; módulos 9-14 do
roadmap.

## BiCGSTAB: o primeiro sistema linear genuinamente não-simétrico do projeto

Todo sistema linear resolvido neste projeto até aqui era simétrico (Poisson
de pressão, difusão/Poiseuille, as varreduras Patankar-linearizadas de
k/omega) ou resolvido por Gauss-Seidel puro (que não se importa com
simetria) - então Gradiente Conjugado sempre foi a ferramenta certa e
suficiente, e nenhuma matriz genuinamente não-simétrica precisou ser
montada. O gatilho concreto identificado num checkpoint de planejamento
anterior, e agora executado: convecção-difusão 1D estacionária com
upwind de 1ª ordem é o caso mais simples que realmente precisa de um
método não-simétrico - o upwind pesa os dois vizinhos de forma desigual
conforme a direção do escoamento, quebrando a simetria que todo stencil
puramente difusivo aqui tem.

`ImplicitConvectionDiffusionSolver1D`: `u·dphi/dx = nu·d²phi/dx² + source`,
Dirichlet nas duas pontas, discretização por volumes finitos - diferença
central pra difusão (igual a todo solver de difusão aqui), upwind de 1ª
ordem pra convecção (a escolha padrão que garante limitação/estabilidade em
qualquer número de Peclet de célula, ao custo de difusão numérica/falsa
medida abaixo). Mesma convenção de parede espelhada Dirichlet de
`MixingLengthChannelFlowSolver1D`.

**Solução exata, derivada aqui, não recitada**: `u·phi' = nu·phi''` tem
solução geral `phi = A·exp(u·x/nu) + B`; aplicando as duas condições de
Dirichlet dá `phi(x) = esquerda + (direita-esquerda)·(exp(u·x/nu)-1)/(exp(u·L/nu)-1)`
- solução de EDO linear comum, não uma tabela de benchmark.

`solveGaussSeidel()` é a referência (relaxação diagonal/resíduo genérica -
o coeficiente diagonal exato de cada célula é conhecido em forma fechada
pra qualquer sinal de `u`, então não precisou derivar dois casos à mão).
`solveBiCGStab()` (Van der Vorst 1992) é o ponto real: matrix-free, no
primeiro operador genuinamente não-simétrico montado aqui.

**Validação, tudo medido antes de escolher tolerância**: a Pe_cell=0.1
(nx=40), Gauss-Seidel precisou de 2804 iterações pra convergir a 1e-13;
BiCGSTAB precisou de 52 - uma redução de **~54x**, o mesmo tipo de ganho que
CG já deu sobre Gauss-Seidel nos sistemas simétricos deste projeto, agora
demonstrado num não-simétrico. Os dois métodos concordam entre si a ~7.8e-11
(checa consistência de implementação, não verdade); contra a solução exata
da EDO, o erro máximo foi ~0.0151 - atribuído à difusão numérica conhecida
do upwind de 1ª ordem, não um bug, confirmado por um segundo teste que
mostra o erro caindo por refino de malha na razão certa (~0.52x ao dobrar
nx, consistente com O(h)). Um terceiro teste cobre termo de fonte e
escoamento reverso (velocidade negativa, o outro ramo do upwind),
checado só contra Gauss-Seidel.

Nenhum bug encontrado; os 5 suites C++ passam na primeira tentativa.
Vinculado ao Python como `aether.ImplicitConvectionDiffusionSolver1D`
(`solve_gauss_seidel`, `solve_bicgstab`, `value`, `cell_center_x`).

## GMRES(m): fecha o item "GMRES/BiCGSTAB" do Módulo 5

`solveGmres()` na mesma classe, sobre o mesmo operador não-simétrico:
GMRES reiniciado (Saad & Schultz 1986), Arnoldi com Gram-Schmidt
**modificado** (deliberadamente não o clássico, que perde ortogonalidade
catastroficamente em ponto flutuante) e rotações de Givens aplicadas
incrementalmente à matriz de Hessenberg - o que mantém uma fatoração QR
sempre atualizada e faz o resíduo ser conhecido exatamente a cada passo
interno, sem resolver o problema de mínimos quadrados.

**O que GMRES tem que BiCGSTAB não tem**: ele minimiza o resíduo sobre todo
o subespaço de Krylov a cada passo, então a norma do resíduo é
monotonicamente não-crescente **por construção** - inclusive entre
reinícios, já que cada ciclo parte do iterado anterior e a correção nula
sempre está no novo subespaço. BiCGSTAB não garante nada disso.

**A validação que de fato distingue os dois** - medida, não afirmada.
Violações de monotonicidade no histórico completo de resíduo de cada método:

| caso | GMRES | BiCGSTAB |
|---|---|---|
| Pe moderado | 0 | 4 |
| velocidade reversa, m=10 | 0 | 6 |
| com termo de fonte | 0 | 14 |

A afirmação que sustenta o teste é o **zero** do GMRES (é a garantia
matemática real). O número não-nulo do BiCGSTAB é checado também, mas só
como prova de que o teste não é vazio nessas entradas - **não** é uma
alegação de que BiCGSTAB sempre oscila: num quarto caso medido (Pe mais
alto) ele saiu monotônico também, e é exatamente por isso que esse caso
ficou de fora.

Segundo teste estrutural: sem reinício, num sistema de dimensão n, o
subespaço de Krylov não pode crescer além de n antes de se tornar
invariante, então GMRES converge em no máximo n iterações - garantia de
terminação finita que nenhum método estacionário (Gauss-Seidel) tem.
Medido: n=25 convergiu em exatamente 25.

**Um achado honesto que vale registrar: aqui GMRES é mais lento que
BiCGSTAB, não mais rápido.** No caso de Pe moderado, GMRES(30) precisou de
174 iterações contra 51 do BiCGSTAB - e mesmo contando produtos
matriz-vetor (GMRES faz 1 por iteração, BiCGSTAB faz 2), são 174 contra
102. A vantagem do GMRES não é velocidade neste problema, é robustez: a
monotonicidade acima. Vale ter os dois pelo mesmo motivo que se tem
Gauss-Seidel e CG - são trade-offs diferentes, não um substituindo o outro.

**Sobre precondicionadores, e por que nenhum foi implementado aqui**: o
candidato óbvio seria Jacobi (diagonal). Mas neste sistema o coeficiente
diagonal é `|u| + 2·nu/h` - que **não depende da célula**: é constante em
todo o domínio (medido: 42.0 para todas as células no caso padrão). Então
`D⁻¹A` é apenas `A` multiplicado por um escalar, e métodos de Krylov são
invariantes a reescalonamento escalar do sistema (o subespaço de Krylov é
idêntico, o minimizador é idêntico). Precondicionamento de Jacobi aqui
seria um **placebo demonstrável**, não uma otimização - então não foi
implementado. Um precondicionador que valha a pena (ILU, ou Jacobi num
sistema de coeficientes variáveis) precisa de um sistema com diagonal
não-constante, que o projeto ainda não tem em forma não-simétrica.

Nenhum bug; os 5 suites C++ passam na primeira tentativa. Vinculado ao
Python como `solve_gmres(restart, max_iterations, tolerance)` e
`residual_history()`.

## Coeficientes variáveis + precondicionadores: fecha o Módulo 5

O item anterior terminou dizendo que precondicionadores estavam travados
por falta de um sistema com diagonal variável. A forma honesta de seguir
era **remover o bloqueio**, não contorná-lo - e de um jeito fisicamente
motivado: todo fechamento de turbulência deste projeto já produz uma
difusividade efetiva `nu + nu_t` que varia célula a célula, então um
operador de convecção-difusão de coeficientes variáveis é exatamente a
forma que um solve implícito de transporte turbulento teria.

`setVelocityField()` / `setDiffusivityField()` (valores por célula, média
para as faces - o mesmo idioma que os solvers 3D já usam pra `nu_t`), mais
`setPreconditioner()` com `None` / `Jacobi` / `IncompleteLU`, aplicados à
esquerda em GMRES e BiCGSTAB.

**Extração da matriz por sondagem, não por transcrição à mão.** Os três
diagonais são recuperados aplicando o operador a vetores unitários, uma
coluna por vez - deliberadamente **não** copiando os coeficientes da
fórmula do stencil. Dois motivos: as linhas de contorno diferem das
interiores (o espelho fantasma dobra um coeficiente de vizinho de volta na
diagonal, com sinal trocado), e este projeto **já foi mordido exatamente
por esse modo de falha** - o `frictionVelocity()` do
`MixingLengthChannelFlowSolver1D` assumia um tratamento de parede que o
próprio `solve()` não usava, erro de ~38%. A sondagem custa nx produtos
matriz-vetor uma vez por solve e é correta por construção. De quebra, ela
**verifica** a tridiagonalidade em vez de assumi-la.

**ILU(0) é exato aqui, e isso é uma propriedade do 1D, não uma alegação
geral.** ILU(0) mantém só o padrão de esparsidade de A; a fatoração LU de
uma matriz tridiagonal não gera preenchimento fora desse padrão, então nada
é descartado e ILU(0) **é** a fatoração LU exata (a recorrência clássica de
Thomas). Logo `M⁻¹A` é a identidade e o Krylov converge em **uma
iteração** - medido: exatamente 1, para GMRES e BiCGSTAB, com e sem
coeficientes variáveis. Em 2D/3D o stencil não é tridiagonal, ILU(0)
descarta preenchimento real e vira aproximação genuína - registrado pra que
o resultado de 1 iteração não seja mal lido depois.

**Jacobi: placebo com coeficiente constante, ganho real com coeficiente
variável** - e o teste mede as duas coisas. Diagonal constante ⟹ `M = D` é
múltiplo escalar da identidade ⟹ Krylov é invariante (o próprio critério de
parada também, já que `‖M⁻¹r‖/‖M⁻¹b‖ = ‖r‖/‖b‖`). Medido: GMRES 251 sem vs
276 com Jacobi - nenhum ganho. Com coeficientes variáveis: 360 vs 287, ~20%
melhor.

**E o limite honesto do resultado do Jacobi: o ganho não é grande nem
monotônico no contraste da difusividade.** Medido variando o contraste
pico-parede na mesma tolerância: 19x → GMRES 360/287 (ajuda 20%); 100x →
304/356 (Jacobi **atrapalha**); 400x e 2500x → ambos batem o teto de 500
iterações (GMRES(30) reiniciado estagna de qualquer jeito). O caso de 19x é
o usado no teste por ser onde o efeito é real e estável, não por ser
representativo - Jacobi é um precondicionador fraco e isso está medido
assim em vez de escolher o enquadramento lisonjeiro.

### O bug mais instrutivo desta etapa: 1 ULP

O teste de Jacobi falhou de um jeito que custou tempo e três hipóteses
erradas. Sintoma: o executável C++ reportava *breakdown* do BiCGSTAB num
caso que as bindings Python resolviam em 107 iterações - **mesma
biblioteca, mesmos parâmetros**.

Hipóteses erradas, na ordem: (1) o limiar absoluto de breakdown estava mal
escalado - troquei por relativo, e **piorei**, porque usar `epsilon` como
fator dispara espuriamente (rho decai naturalmente à escala de epsilon
perto da convergência sem nada ter dado errado); corrigido pra `epsilon²`,
mas o sintoma persistiu idêntico; (2) build obsoleto - timestamps
desmentiram; (3) vazamento de estado entre objetos - não havia.

A causa real só apareceu ao **comparar os 60 valores do campo de entrada um
a um**: eles diferiam em **1 ULP em 16 das 60 entradas**. O MSVC vetoriza um
laço de `std::sin` sob `/O2` usando um seno vetorial, que pode divergir do
seno escalar no último bit - e diverge; o Python chama o escalar. Essa
diferença de 1 ULP na entrada bastava pra mandar o BiCGSTAB por um caminho
totalmente diferente: convergia em 107 a partir de um conjunto de entradas e
**quebrava na iteração 99** a partir do outro.

Duas correções, ambas reais: o perfil de difusividade passou de `sin²` para
a parábola `4t(1-t)` - mesma forma e mesma faixa, mas construída só com
`+ - * /`, operações que o IEEE-754 exige serem corretamente arredondadas,
logo idênticas bit a bit entre compiladores, níveis de otimização e
linguagens. E as asserções passaram a ser sobre **GMRES**, cuja contagem de
iterações é bit-estável (verificado: inalterada sob uma perturbação
deliberada de 1 ULP no campo inteiro), em vez de BiCGSTAB, cuja contagem
oscila algumas unidades sob a mesma perturbação.

Essa última medição é, por acaso, a demonstração mais concreta desta etapa
de **por que a monotonicidade do GMRES importa na prática**: sob a mesma
perturbação de 1 ULP, GMRES dá exatamente a mesma resposta e o BiCGSTAB não.

Mantidas também as melhorias de robustez do BiCGSTAB descobertas no
caminho: limiares de breakdown relativos (`epsilon²`, não `epsilon`) e -
genuinamente correto - um breakdown alcançado com resíduo já dentro da
tolerância é reportado como **convergência**, não falha. Um caso medido
pedia 1e-12, quebrava, e segurava resíduo relativo de 1.688e-12; a primeira
versão jogava essa resposta fora como fracasso.

Um bug latente também foi corrigido: trocar os campos de coeficientes
invalidava a extração da matriz mas **não** a fatoração ILU construída a
partir dela, deixando fatores obsoletos em uso.

Os 5 suites C++ passam. Vinculado ao Python (`set_velocity_field`,
`set_diffusivity_field`, `set_preconditioner` com o enum `Preconditioner`).

**Módulo 5 está completo como escopado**: FVM → Multigrid → CG → BiCGSTAB
→ GMRES → precondicionadores.

**Ainda em aberto**: estender a um sistema convecção-difusão 2D/3D real -
que é também o que daria sentido a um ILU(0) *aproximado* (ainda não é
necessidade concreta, já que os solvers de Navier-Stokes 2D/3D usam projeção
explícita, não um solve implícito); recuperação combinatória de faceta no
CDT; módulos 9-14 do roadmap.

## Módulo 9: UI com painéis (primeira etapa)

Primeiro módulo de produto do roadmap. A decisão de partida foi de escopo,
não de código: **construir a UI do zero em vez de adotar Dear ImGui**. Não
por preferência técnica - ImGui resolveria isso melhor e mais rápido - mas
porque o projeto tem uma decisão permanente de não ter dependências
externas de runtime (o `apps/common` existe justamente porque o bootstrap
WGL/GL-3.3 foi escrito à mão em vez de puxar GLFW/GLAD). Reverter isso não
é uma escolha a tomar de passagem. O custo é real e está registrado: são
algumas centenas de linhas que o ImGui daria de graça, com muito menos
widgets.

**`apps/common/Ui`** - GUI em modo imediato sobre o mesmo pipeline GL 3.3
dos viewers:

- **Fonte bitmap 5x7 embarcada**, cobrindo os 95 caracteres ASCII
  imprimíveis. Sem arquivo de fonte pra distribuir ou carregar. A tabela
  foi autorada e depois **verificada glifo a glifo**, renderizando o
  conjunto inteiro como arte ASCII, e só então convertida pra C++ -
  **programaticamente, não transcrita à mão**. Transcrever 665 números é
  exatamente o passo que planta um pixel errado que ninguém nota até uma
  letra específica sair estranha.
- **Um shader, uma textura, uma chamada de desenho por quadro.** Toda a
  geometria - fundo de painel, faces de botão, trilhos de slider e cada
  glifo de texto - acumula num único buffer de vértices. Retângulos sólidos
  amostram uma célula deliberadamente toda branca no atlas da fonte, então
  não precisam de shader nem caminho de desenho separado.
- Widgets: painel, label, separador, botão, slider e checkbox. Modo
  imediato: sem árvore retida, sem callbacks, sem invalidação - o estado
  mora nos objetos de simulação, não nos widgets. O único estado retido é
  qual slider está sendo arrastado, porque um arrasto atravessa quadros por
  definição.

**Modo `sim`** - o primeiro modo interativo no sentido de CFD, não no de
câmera: todos os anteriores resolviam o problema de uma vez e exibiam um
resultado fixo. Aqui o solver avança dentro do laço de render e o painel
escreve de volta nele. Dá pra trocar o fechamento de turbulência (laminar /
comprimento de mistura / k-epsilon / k-ω SST), editar resolução,
viscosidade e velocidade da tampa, controlar run/pause/passo/reiniciar e
ler diagnósticos ao vivo (Re, malha, passos, tempo, divergência máxima) -
tudo sem recompilar. Os quatro fechamentos 2D têm construtor e superfície
`u`/`v`/`time`/`maxDivergence` idênticos, então um ponteiro por fechamento
mais um visitante pequeno bastam; não foi preciso hierarquia de wrapper.

Duas simplificações deliberadas, ambas pra manter esta primeira etapa sobre
a UI e não sobre encanamento de renderer: o campo é desenhado como heatmap
por célula através do próprio `drawRect()` da UI (que já agrupa retângulos
em uma chamada só - mais barato e muito menos código que levantar um
segundo shader, e exercita a camada nova o bastante pra servir de teste de
fumaça); e é 2D, não 3D - o painel é agnóstico de dimensão, mas os
fechamentos 3D custam o suficiente por passo pra que passo interativo exija
uma thread de trabalho, o que é tarefa própria.

**Validado com captura real da janela** (técnica `PrintWindow`, a mesma de
todos os viewers): a fonte sai legível, o painel e os widgets desenham
corretos, e o campo mostra o resultado fisicamente certo - banda vermelha
rápida no topo **afinando nos dois cantos** (o taper regularizado da tampa,
que existe justamente pra remover a singularidade de pressão de uma tampa
descontínua) e o núcleo do vórtice primário como mancha de baixa velocidade
acima do centro, deslocado no sentido do movimento da tampa, que é onde ele
deve estar a Re=100.

### Campos numéricos editáveis por teclado

Complemento ao slider: `Ui::textField()` - clique pra focar, `WM_CHAR` pra
digitar (dígitos, `.`, `-`), Backspace edita, Enter confirma (analisa,
limita ao intervalo, escreve `*value`), Escape descarta, e **clique fora
também confirma** - a convenção usual de campo de texto, não só a de
arrasto do slider. Um buffer malformado no commit (vazio, ou só `"-"`)
deixa `*value` inalterado em vez de escrever lixo.

Reaproveita o mesmo padrão "enfileira no WndProc, drena por quadro" que
`mousePressed`/`mouseReleased` já usavam: `WM_CHAR` empilha em
`g_pendingText`, drenado pra `UiInput::textInput` uma vez por quadro. Um
detalhe que evitou um caminho de código a mais: `TranslateMessage()` (já
chamado no laço de mensagens) gera `WM_CHAR` também pra Backspace/Enter/Esc,
não só pra caracteres imprimíveis - então um único handler cobre digitação
e os três controles sem precisar de `WM_KEYDOWN`/`VK_*` separado.

Trocados por `textField()` os parâmetros que costumam precisar de um valor
preciso (um Reynolds específico) que uma trilha de slider de 260px não
alcança bem: viscosidade e velocidade da tampa. Resolução continua slider -
precisa cair num inteiro, e arrastar é o jeito natural de escolher um.

**Validado com captura real da janela, os dois caminhos de confirmação**:
clique + digitação + Enter aplicou o valor, disparou reconstrução (passos
voltam a 0) e reformatou a caixa; clique + digitação + clique-fora (sem
Enter) também aplicou - confirmando que o commit por clique externo, um
caminho de código distinto do de Enter, funciona de verdade e não só em
teoria.

### Modo `sim3d`: o painel sobre a cavidade 3D real, com thread de trabalho

**Este é o primeiro código multi-thread do projeto**, e a razão é concreta:
um passo em 3D custa O(n³) em vez de O(n²), e os construtores de
k-epsilon/k-omega SST/DES rodam um primer de 400 passos de mistura por
dentro antes do painel poder mostrar qualquer coisa - nenhum dos dois é
barato o bastante pra caber no laço de render a 60fps como no modo `sim`
2D. Então o solver passou a viver na própria thread.

**Contrato de threading, registrado explicitamente porque é fácil errar
sutilmente**: `rebuildIfRequested()`/`stepOnce()` só rodam na thread de
trabalho; `requestRebuild()`/`requestStepOnce()`/`setRunning()` só rodam na
thread de render/UI; `trySnapshot()` só roda na thread de render e **nunca
bloqueia** - é o único lugar onde um `lock_guard` comum seria errado, já
que uma leitura bloqueante ali congelaria o painel pelo tempo que uma
reconstrução levar.

A thread de render usa `std::mutex::try_lock()` uma vez por quadro pra
copiar o que precisa desenhar; se falhar (worker no meio de um passo ou
reconstrução), o quadro simplesmente redesenha o retrato anterior. Trocar
um parâmetro *pode* bloquear brevemente por até a duração de um passo (a
troca precisa da mesma trava que o passo usa), mas isso é uma ação rara do
usuário, não o laço de desenho - a resposividade que importa (painel,
câmera orbital) nunca espera por nada.

**Validado com captura real da janela em três cenários**: (1) rodando por
~6s o laminar acumulou **112.812 passos** com divergência em `1.63e-12`,
confirmando que o worker avança livre do teto de quadros do render; (2)
trocar de fechamento durante a execução manteve o painel inteiramente
interativo durante e depois da reconstrução, sem travar; (3) arrastar fora
do painel orbita a câmera corretamente (o cubo girou, o painel ficou
intacto), confirmando que `ui.wantsMouse()` separa corretamente cliques de
UI de arrasto de câmera.

Reaproveita a maior parte do `cavity3d_mode` já existente (câmera orbital,
setas coloridas por velocidade, caixa de arame) - a diferença é que aqui o
buffer de vértices é reconstruído a cada quadro a partir do retrato mais
recente em vez de uma vez só no início, e os seis fechamentos 3D (não só o
laminar) ficam disponíveis por um seletor no painel, com o mesmo padrão
"um ponteiro por variante mais um pequeno visitante" das outras telas.

**Ainda em aberto no Módulo 9**: janelas móveis/ancoráveis; árvore de
cena.

## Faxina antes dos módulos 9-14

Antes de entrar nos módulos de produto (UI, GPU, persistência, IA, API,
plugins), uma consolidação em três fases. O motivo do timing é concreto:
essas camadas todas se acoplam à API dos solvers, e o Módulo 10 (GPU)
exigiria portar o **mesmo** kernel duplicado seis vezes. Pagar essa dívida
depois custa muito mais.

**Fase 1 - versionamento.** O repositório tinha *um commit e um arquivo
rastreado* (o README): ~22 mil linhas de engine viviam fora do controle de
versão. Commitado tudo (121 arquivos), com `.gitignore` ampliado (artefatos
MSVC, `settings.local.json` que continha caminhos de uma máquina antiga) e
um `.gitattributes` novo normalizando fim de linha - relevante porque o
projeto é desenvolvido em mais de um dispositivo, e sem isso todo arquivo
tocado em outra máquina aparece como reescrito por inteiro no diff.

**Fase 2 - código morto.** Os quatro visualizadores standalone (`viewer`,
`field_viewer`, `cavity_viewer`, `turbulence_viewer`) já tinham virado
modos do `unified_viewer`, mas continuavam na árvore: uma segunda cópia do
mesmo código de render, condenada a divergir. Removidos (−1.092 linhas)
após confirmar que o unificado exercita tudo que eles faziam.

**Fase 3 - `StaggeredCavityBase3D`.** Seis classes de cavidade 3D
carregavam cada uma sua cópia literal da maquinaria de malha staggered.
Extraída para uma base comum: **−2.513 linhas** no módulo solver.

| classe | antes | depois |
|---|---:|---:|
| `StaggeredLidDrivenCavitySolver3D` | 391 | **19** |
| `MixingLengthLidDrivenCavitySolver3D` | 475 | 92 |
| `SmagorinskyLesLidDrivenCavitySolver3D` | 485 | 107 |
| `KEpsilonLidDrivenCavitySolver3D` | 648 | 280 |
| `KOmegaSSTLidDrivenCavitySolver3D` | 750 | 375 |
| `DesSstLidDrivenCavitySolver3D` | 754 | 381 |
| `StaggeredCavityBase3D` (novo) | — | 382 |

**A extração foi verificada como preservadora de comportamento antes de ser
feita, não depois.** Cada função foi comparada nas seis classes
normalizando espaço/comentários e comparando hash: `uAt`, `vAt`, `wAt`,
`maxDivergence`, `applyLaplacian`, `projectToDivergenceFree`,
`lidVelocityAt`, `wallDistanceAt` e `dot` eram **idênticas byte a byte** nas
seis; o preditor de momento era idêntico nas cinco turbulentas. Só duas
coisas divergiam, ambas compreendidas: `pAt` (a classe LES escrevia o clamp
com `std::clamp`, as outras com `if/else` - mesma semântica) e
`stableTimeStep` (a laminar não tem `nu_t` a somar).

**Sem funções virtuais, deliberadamente.** O desenho óbvio - um hook virtual
`eddyViscosityAt()` - poria despacho virtual no laço mais interno do
momento, seis vezes por célula por componente. Em vez disso cada fechamento
registra um ponteiro pro seu próprio `nut_` via `setEddyViscosityField()`;
`nutAt()` é função comum, inlinável, que retorna 0.0 quando o ponteiro é
nulo (o caso laminar). Um desvio previsível em vez de uma consulta de
vtable, e a classe continua não-polimórfica (daí o destrutor protegido e
não-virtual).

**O caso laminar mereceu medição, não suposição.** Unificar seu momento
significa trocar `nu·laplaciano` pela forma ponderada por face - com
`nu_t=0` as duas são algebricamente idênticas, mas **não** bit a bit (é
`nu*(a-b) - nu*(b-c)` contra `nu*(a-2b+c)`). Medido diretamente contra
valores de referência capturados antes da mudança, após 200 passos: 18.8%
das amostras idênticas bit a bit e diferença máxima de **9.09e-16** em
valores até 0.5 - cerca de 9 ULPs acumulados, exatamente a mudança de ordem
de arredondamento esperada, não de comportamento (divergência máxima
3.293e-13 contra 3.304e-13 antes). Com isso a classe laminar caiu de 391
para **19 linhas**.

**Dois bugs foram introduzidos pela refatoração e pegos pelos testes** -
vale registrar porque ambos vieram de automatizar demais. Ao reescrever os
construtores por regex, a substituição da lista de inicialização engoliu
`cDes_(cDes)` no DES, deixando o membro não inicializado (lixo → `NaN` no
`F_DES` já no primeiro passo); e o `setEddyViscosityField` não foi inserido
no k-epsilon porque o corpo do construtor dele começa com um comentário, não
com a linha que o padrão procurava. Os dois apareceram no teste de repouso
exato (`u == 0.0` bit a bit) - o tipo de asserção que parece pedante até
salvar uma refatoração. Também houve duas tentativas descartadas de
script antes disso (um splitter que cortava função por "linha igual a `}`",
engolindo funções de uma linha; e um filtro de header que comia o corpo da
classe), revertidas via git - que é exatamente por que a Fase 1 veio antes.

## Marching cubes 3D: fecha o pré-requisito do renderizador de iso-superfícies

`marchingCubes3D` fecha o item explicitamente deixado em aberto na task
anterior: extração de iso-superfícies 3D, o "primeiro passo rumo a" que o
`marchingSquares2D` já documentava desde que foi criado.

**Implementado via marching *tetraedros*, não a tabela clássica de 256
casos do cubo.** Recitar essa tabela de memória (ou os ~15 casos canônicos
+ rotações/reflexões/complementos em que costuma ser comprimida) seria
exatamente o tipo de risco de erro de transcrição que este projeto evita -
o mesmo raciocínio que já levou o `marchingSquares2D` a derivar seus casos
da contagem de arestas cruzadas em vez de uma tabela de 16 linhas. Marching
tetraedros evita o problema por completo: cada cubo é dividido em 6
tetraedros compartilhando a diagonal principal (decomposição de Freudenthal
padrão), e um tetraedro só tem 2⁴=16 combinações de sinal nos cantos -
poucas o bastante, e simples o bastante (nenhum caso ambíguo de face existe
num tetraedro, ao contrário de um cubo), pra derivar por contagem de cantos
acima do iso-valor em vez de recitar qualquer coisa: 0 ou 4 cantos acima -
sem triângulo; 1 ou 3 cantos acima - as 3 arestas do canto isolado até os
outros três são exatamente as arestas cruzadas, um triângulo; 2 cantos
acima e 2 abaixo - as 4 arestas conectando o par "acima" ao par "abaixo"
formam um quadrilátero (percorrido em ordem cíclica pelo emparelhamento
bipartido), dividido em 2 triângulos.

**Validado a decomposição em si antes de qualquer código**: os 6 tetraedros
de Freudenthal de um cubo unitário somam exatamente 1.0 de volume (`numpy`,
verificação independente) - a mesma disciplina de "medir um invariante
exato" já usada na validação por volume da bipirâmide do
`DelaunayTetrahedralization3D`.

**Orientação das normais resolvida de forma genérica, não por caso**: em
vez de raciocinar manualmente sobre a ordem de vértices de cada um dos 3
padrões de caso, cada triângulo gerado é reorientado (trocando dois
vértices se necessário) checando o produto escalar entre sua normal e a
direção até o lado "acima do iso-valor" conhecido - funciona igual pros
casos de 1, 2 ou 3 cantos acima, sem lógica duplicada.

**Validado em Python antes do teste C++, mesmo padrão do `marchingSquares2D`
com seu círculo exato**: `f(x,y,z) = distância até um centro` tem a esfera
de raio `r` como sua iso-superfície `f=r`, pra qualquer `r`. Medido
diretamente (malha 30³, r=0.6): erro máximo de raio de `~2.6e-3` (espaçamento
de malha `~0.067`, mais de uma ordem de grandeza maior), área total de
`~4.51` contra `4·pi·r²≈4.524` esperado, e - a checagem específica da
lógica de orientação por caso - **as 9216 normais geradas apontam pra fora
da esfera, sem uma única exceção**.

Nenhum bug encontrado - a decisão de derivar por contagem de cantos em vez
de tentar reconstruir a tabela de 256 casos evitou exatamente o tipo de
erro sutil que essa tabela costuma introduzir. Todos os 5 suites de teste
C++ passam (`ctest -C Debug`, 5/5) na primeira tentativa. Vinculado ao
Python como `aether.marching_cubes_3d` / `aether.Triangle3D`.

**Ainda em aberto**: um modo de viewer que renderize de fato as
iso-superfícies extraídas (ex.: `nu_t`/`k`/`epsilon`/`omega` dos
fechamentos de turbulência 3D) - o algoritmo em si está pronto e validado,
falta só a integração com `apps/unified_viewer`.

## Viewer de iso-superfície 3D: fecha a lacuna do Módulo 8

`apps/unified_viewer` ganhou um sexto modo, `isosurface`, fechando o item
deixado em aberto na task do marching cubes: renderiza a iso-superfície de
`nu_t` (viscosidade turbulenta) de `KOmegaSSTLidDrivenCavitySolver3D` - o
fechamento de turbulência mais completo deste projeto - como uma malha
sombreada de verdade, não setas. Reaproveita o pipeline de shading do modo
`mesh` (posição+normal, lambertiano de dois lados, sombreamento plano por
face) em vez do pipeline de linhas coloridas do `cavity3d`, já que um volume
escalar é exatamente pra isso que o marching cubes foi construído.

**Iso-valor relativo, não absoluto**: em vez de um valor fixo de `nu_t`,
usa uma fração (`0.3`) do máximo observado naquela execução - a magnitude
de `nu_t` depende de viscosidade/Re/resolução, então um valor absoluto
recitado de antemão poderia facilmente não cruzar o campo em lugar nenhum
(superfície vazia) ou cruzar tudo. Uma fração do máximo garante uma
superfície não-vazia e de tamanho razoável independente dos números exatos
que a execução produzir.

**Resolução/passos escolhidos por medição direta, não palpite**: a
primeira tentativa (malha 16³, 500 passos) levou `~32s` numa build Debug -
tempo de espera longo demais pra um viewer que deveria abrir rápido. Medido
o tempo em função do tamanho via Python antes de decidir: malha 12³ com 300
passos (a mesma resolução que o próprio teste C++ de
`KOmegaSSTLidDrivenCavitySolver3D` já usa) leva `~8s` - ainda perceptível
(esse é, de longe, o solver mais caro que qualquer modo deste viewer roda,
por causa do warm-start de 400 passos + transporte de k/omega), mas
razoável.

Validado com uma captura real de tela: um blob laranja liso e corretamente
sombreado aparece dentro do cubo de arame, com variação de brilho por face
condizente com a iluminação lambertiana de dois lados - confirmando
visualmente que a extração + sombreamento funcionam de ponta a ponta, não
só "compilou e não travou". `1132` triângulos gerados nesta execução.

Compilação limpa de primeira. Isso fecha por completo a sequência iniciada
com o gap de visualização do Módulo 8: os quatro solvers 3D agora têm
representação visual (vetores via `cavity3d`, e o campo escalar mais
interessante de cada fechamento de turbulência via `isosurface`).

`MultigridPoissonSolver2D` - a fatia mais valiosa do próximo item do roadmap
(Multigrid/GMRES/BiCGSTAB/precondicionadores): multigrid geométrico V-cycle
(esquema de correção) para a mesma equação de Poisson 2D com contorno de
Dirichlet que `SteadyDiffusionSolver` já resolve, mas com contagem de
iterações independente da resolução da malha - a vantagem clássica do
multigrid. Implementado como classe autônoma (não encaixado em
`DiffusionProblem`) porque o coarsening exige `nx`/`ny` potências de 2, algo
que a hierarquia genérica de 6 faces/anisotrópica nunca foi projetada para
exigir. Usa mirror-Dirichlet fantasma (`ghost = 2*valor_parede - interior`,
a mesma técnica já validada em `LidDrivenCavitySolver2D`) em vez da
convenção de "célula de contorno fixada" do `DiffusionProblem` - escolha
deliberada, porque assim toda célula do array é uma incógnita genuína, o
que faz o coarsening parear blocos 2x2 sem casos especiais na borda do
domínio. Restrição por média simples do bloco 2x2 (natural para células
uniformes tipo volume finito); prolongação por injeção constante por
partes (mais simples que interpolação bilinear, simplificação documentada
desta primeira versão). Validado contra a mesma solução de série de
Fourier já usada pelo `SteadyDiffusionSolver` (erro ~2.19 numa malha
64x64, o mesmo tipo de erro O(h) já documentado ali) e medido diretamente
contra CG/Gauss-Seidel no mesmo problema: **11 V-cycles** contra **244
iterações de CG** e **7425 de Gauss-Seidel** - cada V-cycle custa
aproximadamente o equivalente a um punhado de varreduras completas de
Gauss-Seidel (2 antes + 2 depois na malha mais fina, mais barato em cada
nível mais grosso), então isso é um ganho real, não só uma contagem menor
com custo por iteração escondido.

Também foi adicionado `SteadyDiffusionSolver::solvePreconditionedConjugateGradient()`
(CG pré-condicionado por Jacobi, `M^-1 = 1/diag(A)`) - medido diretamente e
inicialmente contraintuitivo até se pensar sobre o motivo: converge no
*mesmo* número de iterações que o CG simples (276 vs 276 numa malha
81x41 anisotrópica; 619 vs 618 numa 200x20) em vez de menos. Isso é um
resultado nulo real e explicável, não um bug: toda célula livre interior
nas malhas Cartesianas uniformes deste projeto compartilha exatamente o
mesmo stencil discreto, então `diag(A)` é um único valor constante repetido
em toda parte - pré-condicionamento de Jacobi vira só um reescalonamento
uniforme do sistema inteiro, que não muda o comportamento de convergência
do CG (que depende do espalhamento *relativo* dos autovalores, algo que um
fator constante global não altera). Jacobi ajudaria numa malha com
diagonal genuinamente variável por célula (malha graduada, difusividade
espacialmente variável) - um caso que os testes deste projeto ainda não
exercitam.

**GMRES e BiCGSTAB foram deliberadamente adiados**, não esquecidos: ambos
são métodos para sistemas lineares *não-simétricos* em geral, mas todo
sistema linear que este projeto resolve até agora (Poisson/difusão, todas
as equações de pressão dos solvers de Navier-Stokes) é simétrico positivo-
definido, onde CG já é o método ótimo - implementá-los agora significaria
construir e "validar" contra um sistema que não existe ainda no código,
exatamente o tipo de funcionalidade especulativa e não demonstrada que este
projeto evita. Fazem mais sentido quando uma discretização genuinamente
não-simétrica existir (ex. um esquema de convecção upwind), o que ainda não
é o caso.

`engine/postprocessing` (Módulo 7) - `Streamline2D` integra linhas de
corrente reais via Runge-Kutta de 4ª ordem sobre um campo de velocidade
cell-centered (mesma convenção usada por todos os solvers 2D colocalizados
deste projeto), interpolado bilinearmente entre os 4 centros de célula mais
próximos pra amostrar em pontos arbitrários (não alinhados à malha). Suporta
domínio periódico (wrap físico, mesma convenção do `TaylorGreenVortexSolver2D`)
ou não-periódico (a trajetória para se sair do domínio, mesma convenção do
`LidDrivenCavitySolver2D`). Validado contra o vórtice de Taylor-Green: como
`u=U0*cos(x)*sin(y)`, `v=-U0*sin(x)*cos(y)` tem função de corrente exata
`psi=-U0*cos(x)*cos(y)` (verificável por diferenciação direta,
`d(psi)/dy=u`, `-d(psi)/dx=v`), uma linha de corrente por definição preserva
`psi` constante - checagem forte e autocontida (medido diretamente: desvio
de `psi` de ~7.7e-7 ao longo de 3000 passos), não uma comparação com
software externo.

`marchingSquares2D` extrai iso-contornos de um campo escalar 2D
cell-centered - o primeiro passo rumo ao marching cubes 3D, escopado
deliberadamente pra 2D primeiro (mesma prática incremental usada em toda a
malha não-estruturada deste projeto: 2D antes de 3D). Em vez de uma tabela
de 16 casos decorada (risco real de erro de transcrição), os segmentos são
derivados diretamente da contagem de arestas cruzadas do quadrilátero
(sempre 0, 2 ou 4 - uma certeza topológica de um contorno de 4 vértices),
com o caso ambíguo de 4 cruzamentos ("sela") resolvido por um critério
simples de valor médio dos 4 cantos. Validado contra um campo com contorno
exato conhecido: `f(x,y)=x²+y²` tem o círculo de raio `r` como seu
iso-contorno `f=r²`, pra qualquer `r` - medido diretamente, erro máximo de
raio de ~1.2e-4 (a espaçamento de malha ~0.031, duas ordens de grandeza
maior) e comprimento total do contorno a ~8e-4 de `2*pi*r`.

`apps/cavity_viewer` foi atualizado pra desenhar streamlines reais (curvas
amarelas, RK4 via `Streamline2D`) sobrepostas ao campo de setas já
existente - confirmado visualmente via captura real de janela que elas
seguem o mesmo vórtice único de recirculação já validado numericamente.

`apps/viewer` (Módulo 8) foi modernizado de OpenGL fixo/immediate-mode
(`glBegin`/`glVertex`, GLU) para OpenGL 3.3 core profile de verdade: shaders
GLSL compilados em tempo de execução, VBO/VAO pra geometria, e uma matriz
MVP calculada à mão (`Matrix4x4`, novo em `engine/core` - `perspective()` e
`lookAt()` substituem `gluPerspective`/`gluLookAt`, que não existem num
contexto core profile). Sem depender de GLFW/GLAD/Vulkan (mantendo a
filosofia "nada além do Windows SDK" do projeto): os poucos ponteiros de
função GL 3.3/WGL necessários (`glGenVertexArrays`, `glCreateShader`,
`wglCreateContextAttribsARB`, etc.) são declarados e carregados à mão via
`wglGetProcAddress`, o padrão usual pra evitar uma biblioteca de loader só
por causa de uma dúzia de símbolos. Criar um contexto core profile exige o
truque clássico de "contexto descartável": criar uma janela+contexto legado
temporária só pra obter `wglCreateContextAttribsARB` (que só existe depois
que *algum* contexto já está ativo), destruí-la, e então criar a janela
real com um contexto 3.3 de verdade.

**Um bug real de janela foi encontrado e corrigido durante o
desenvolvimento**: a primeira versão usava o mesmo `WndProc` tanto pra
janela temporária quanto pra real - como `WndProc` chama
`PostQuitMessage(0)` em `WM_DESTROY`, destruir a janela temporária colocava
uma mensagem `WM_QUIT` na fila da *thread* (não da janela específica), que
a janela real capturava na primeiríssima iteração do loop de mensagens e
saía imediatamente, antes de renderizar um único frame - o processo
terminava com código de saída 0, sem nenhuma mensagem de erro, então o
`printf` inicial aparecia mas a janela nunca ficava visível. Encontrado ao
notar que o processo desaparecia da lista de tarefas segundos depois do
lançamento, mesmo sem erro impresso. Corrigido dando à janela temporária um
`WndProc` mínimo (`DefWindowProc` puro), separado do `WndProc` real.

Validado com uma captura real de janela (mesma técnica `PrintWindow` usada
em toda a sessão): o tetraedro de teste renderiza corretamente, com
sombreamento lambertiano de dois lados e teste de profundidade funcionando
(só a face voltada pra câmera aparece, as outras três ficam ocultas atrás,
como esperado para um sólido fechado).

**Depois de provar o pipeline no `viewer`, o toolkit de bootstrap (contexto
GL 3.3 core + WGL, compilação/link de shaders) foi extraído para
`apps/common` (`aether_app_common`)** - em vez de copiar e colar as mesmas
~150 linhas de boilerplate (e arriscar reintroduzir o mesmo bug de janela)
nos outros três apps, cada um agora só declara seus próprios shaders,
monta seu próprio VBO e usa o toolkit compartilhado pra bootstrap/compilação.
Os quatro visualizadores legados (`viewer`, `field_viewer`, `cavity_viewer`,
`turbulence_viewer`) foram todos migrados pra OpenGL 3.3 core profile:
`field_viewer` (heatmap por quads coloridos), `cavity_viewer` (setas +
streamlines por `GL_LINES`/`GL_LINE_STRIP`) e `turbulence_viewer` (curvas
log-law por `GL_LINE_STRIP` + marcadores `GL_POINTS`) - cada um validado
com uma captura real de janela confirmando renderização pixel-a-pixel
equivalente à versão anterior em OpenGL fixo. Adicionado `Matrix4x4::ortho()`
(projeção ortográfica, faltava desde a primeira versão do `viewer`, que só
precisava de `perspective()`/`lookAt()`), validado com um teste autoderivado
(o centro exato da caixa `ortho(0,2,0,1,-1,1)` deve mapear pra origem do
NDC, verdade por construção pra qualquer projeção ortográfica).

**Etapa seguinte: unificação num único executável.** `apps/unified_viewer`
reúne os quatro modos (`mesh`, `heatmap`, `cavity`, `turbulence`) num só
`main.cpp`, selecionados por `argv[1]`. Cada modo manteve seu código quase
verbatim (shaders, globais, `WndProc`, geração de VBO) dentro do seu próprio
namespace (`mesh_mode`, `heatmap_mode`, `cavity_mode`, `turbulence_mode`) -
decisão deliberada de *não* forçar uma abstração de render-loop única sobre
quatro necessidades genuinamente diferentes (uma cena 3D em perspectiva com
câmera orbital interativa vs. três overlays 2D ortográficos estáticos), já
que apenas um modo roda por invocação do processo mesmo. Os namespaces
existem só pra evitar colisão de símbolos (`WndProc`, `g_width`, `kVertexShaderSource`
etc. repetidos) agora que os quatro `.cpp` viraram uma única unidade de
tradução.

Os quatro executáveis originais (`aether_viewer`, `aether_field_viewer`,
`aether_cavity_viewer`, `aether_turbulence_viewer`) foram mantidos na árvore
por um tempo como pontos de entrada adicionais, e **removidos depois**, na
faxina anterior aos módulos 9-14: como cada um virou um modo do executável
unificado, mantê-los significava carregar uma segunda cópia do mesmo código
condenada a divergir da primeira. Ver a seção "Faxina antes dos módulos
9-14" mais abaixo.

Validado compilando e rodando os quatro modos (`mesh` com um tetraedro STL
de teste, `heatmap`, `cavity`, `turbulence` sem argumentos) e conferindo que
os diagnósticos impressos no console batem exatamente com os das versões
standalone (mesmas iterações até convergência, mesmo `u_tau` nos três
fechamentos de turbulência, mesma divergência máxima da cavidade) - a
migração de código é estrutural (namespaces + dispatch), não numérica,
então a validação relevante é "mesmo comportamento observável", não uma
nova derivação analítica.

## Build (C++ core + bindings)

Requer CMake >= 3.20, um compilador C++20 (MSVC, GCC ou Clang) e, para os
bindings Python, `pip install pybind11` no ambiente Python usado para configurar.

```
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Se o pybind11 não for encontrado, o CMake pula os bindings automaticamente
(o core e os testes C++ continuam buildando normalmente). Se o pybind11
estiver instalado num venv não-padrão, aponte o CMake pra ele com
`-Dpybind11_DIR=<venv>/Lib/site-packages/pybind11/share/cmake/pybind11`.

Para usar o pacote `aether` a partir do Python, coloque a pasta com o
`.pyd`/`.so` compilado no `PYTHONPATH` (ex.:
`build/bindings/python/Release` no Windows com MSVC):

```
$env:PYTHONPATH = "build/bindings/python/Release"
python -c "import aether; print(aether.Vector3(1,1,1).norm())"
```

## Visualizador

Um único executável com seleção de modo por argumento:

```
build\apps\unified_viewer\Release\aether_unified_viewer.exe <modo>
```

| modo | o que mostra |
|---|---|
| `mesh <arquivo.stl>` | malha STL/OBJ importada, câmera orbital em perspectiva |
| `heatmap` | placa com 3 lados frios e 1 quente (`SteadyDiffusionSolver`), heatmap azul-branco-vermelho |
| `cavity` | cavidade 2D com tampa deslizante (`LidDrivenCavitySolver2D`, Re=10): vetores + streamlines reais |
| `cavity3d` | campo de vetores 3D da cavidade staggered (`StaggeredLidDrivenCavitySolver3D`) |
| `isosurface` | iso-superfície real de `nu_t` do k-ω SST 3D, via marching cubes |
| `turbulence` | u+ vs ln(y+) dos três fechamentos de canal 1D sobre a lei da parede teórica |

No modo `mesh`/`cavity3d`/`isosurface`: arrastar com o botão esquerdo gira a
vista, a roda aproxima/afasta. Pipeline OpenGL 3.3 core profile (shaders
GLSL + VBOs/VAOs) em todos os modos.

Sobre o modo `turbulence`, vale a explicação que o app também imprime no
console: a curva laranja (k-epsilon) gruda na reta teórica desde o primeiro
ponto porque a célula adjacente à parede é *literalmente fixada* por essa
fórmula - é a própria condição de contorno. A azul (comprimento de mistura)
fica sistematicamente abaixo com a mesma inclinação: o modelo acerta o
expoente 1/kappa (validado em `solver_tests.cpp`) mas não foi ajustado pra
reproduzir a constante aditiva B=5.0, que não existe explicitamente nesse
fechamento mais simples (sem amortecimento de van Driest). Isso é esperado,
não um bug.
