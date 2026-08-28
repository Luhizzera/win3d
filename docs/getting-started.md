# Começando

Este guia leva de um clone novo até um resultado aberto no ParaView.

O `README.md` na raiz é um diário de desenvolvimento: registra o que foi
medido, o que falhou e por quê, em ordem cronológica. É útil para continuar
o desenvolvimento e péssimo para começar a usar. Este documento é o
contrário.

---

## 1. Instalar

Duas formas, dependendo do que você quer.

**Só usar o pacote Python**: `pip install .` (ou, para gerar o arquivo e
guardar/distribuir, `pip wheel . -w dist`) na raiz do repositório. Precisa
de CMake ≥ 3.20 e um compilador C++20 (MSVC, GCC ou Clang) — `pip` cuida de
buscar `scikit-build-core` e `pybind11` sozinho, num ambiente de build
isolado. Compila o motor, empacota as extensões dentro do próprio pacote
`aether`, e o resultado importa em qualquer venv, sem o repositório
clonado:

```
pip install .
python -c "import aether; print(aether.Vector3(1, 1, 1).norm())"
```

Verificado de ponta a ponta (não só "deveria funcionar"): wheel construída,
instalada numa venv nova sem nenhuma relação com este repositório, e um
solver rodado a partir dela sem erro. `requires-python = ">=3.12"` em
`pyproject.toml` é literal — é a única versão com a qual este projeto
realmente testou, não uma suposição de compatibilidade mais ampla.

**Desenvolver o motor** (rodar a suíte C++, editar `engine/`): `pip
install .` recompila do zero a cada chamada, o que é lento para iteração.
Use o fluxo de sempre em vez disso —

```
python build.py
```

— que configura, compila em Release e roda as 13 suítes de teste
incrementalmente. `python/` no `PYTHONPATH` encontra os módulos compilados
sozinho dentro do `build/` deste repositório, sem precisar instalar nada:

```python
import sys; sys.path.insert(0, "python")
import aether
print(aether.Vector3(1, 1, 1).norm())
```

Os dois caminhos coexistem sem conflito: `python/aether/__init__.py`
procura as extensões primeiro dentro de si mesmo (o layout que uma
instalação via `pip` produz), depois em `build/` (o layout de
desenvolvimento) — o que for encontrado primeiro vence.

---

## 2. A primeira simulação

```
python examples/flow_around_object.py
```

Roda o pipeline inteiro sobre um icosaedro construído em memória. Passe um
arquivo para usar sua própria geometria:

```
python examples/flow_around_object.py meu_objeto.stl
```

O exemplo está comentado passo a passo e é o melhor ponto de partida para
ler. O resto deste guia explica *por que* cada etapa existe.

---

## 3. O pipeline, etapa por etapa

### Geometria

```python
surface = aether.load_stl("objeto.stl")   # ou load_obj
```

`load_stl` solda vértices coincidentes automaticamente, e precisa: o
formato STL grava três vértices independentes por triângulo, sem
conectividade compartilhada, então sem soldar nenhum par de triângulos
divide uma aresta e a superfície inteira parece cheia de furos.
`load_obj` não solda porque OBJ já guarda uma lista indexada.

Para geometria vinda de um CAD real, há `load_step`:

```python
result = aether.load_step("peca.step")
if result.unsupported_features:
    print(result.unsupported_features)  # o que não pôde ser interpretado
surface = result.mesh
```

Só o subconjunto **BREP facetado** de STEP é entendido — toda face plana,
todo contorno um `POLY_LOOP` simples, sem furos. É o formato que a maioria
dos exportadores de CAD produz quando pedem uma malha "facetada"/"low
poly" em vez de superfícies exatas. Uma face curva (cilíndrica, B-spline)
ou uma representação tesselada não trava nem falha em silêncio: aparece
listada em `unsupported_features`, e `result.mesh` traz o que pôde ser
extraído do resto do arquivo. Não há suporte a IGES nem a CAD curvo de
verdade — isso exigiria um kernel CAD (OpenCASCADE), ainda não incluído.

**A superfície precisa ser fechada.** O mesher esculpe o *interior* do
objeto de dentro de uma caixa de fluido, e uma superfície com furos não tem
interior:

```python
assert surface.is_watertight()
```

### Malha de volume

```python
domain = aether.mesh_flow_around_object(surface, margin=2.5)
```

Constrói uma caixa em volta do objeto, tetraedraliza, recupera as facetas
da superfície e esculpe o sólido — sobra o fluido.

- `margin` é o quanto a caixa se estende além do objeto, em múltiplos do
  tamanho do próprio objeto. Pequeno demais e as paredes interferem no
  escoamento; grande demais e você gasta células em espaço vazio.
- `background_coarsening` controla a densidade da malha de fundo, em
  múltiplos do tamanho médio de triângulo da superfície. Menor = mais
  células.

Levanta `MeshGenerationError` se alguma faceta não pôde ser recuperada, se
o ponto-semente não caiu dentro do sólido, ou se o volume esculpido não bate
com caixa-menos-objeto dentro de 1%. Prefere falhar a devolver uma malha
que não passa na própria checagem.

**Escopo**: escoamento externo, em torno de um objeto. Escoamento interno
em duto (o interior de uma superfície fechada) não está coberto — é uma
capacidade diferente, não um ajuste desta.

### Condições de contorno

```python
lid = aether.driving_wall_velocity(domain, face="z_max",
                                   direction=(1.0, 0.0, 0.0), speed=1.0)
```

Há **duas formas** de dirigir a mesma malha, e a escolha é o que separa os
dois tipos de caso.

**Caixa fechada, agitada por uma tampa** (`driving_wall_velocity`): uma
face da caixa desliza no próprio plano; todo o resto — as outras cinco
faces e o objeto — é parede sem escorregamento.

**A direção precisa ser tangencial à face**, e a função recusa se não for.
Numa caixa selada, uma parede empurrando ao longo da própria normal injeta
massa que não tem para onde ir: não é um caso difícil, é um caso insolúvel.
Vale confirmar antes de gastar tempo:

```python
report = aether.check_closed_domain_conservation(domain.mesh, lid)
assert report.is_conservative
```

Geometria pura, sem marchar um passo.

**Escoamento passante** (`freestream_boundary`): entra por uma face e sai
pela oposta, passando pelo objeto — aerodinâmica externa comum.

```python
wall_velocity, is_outlet = aether.freestream_boundary(
    domain, inlet_face="x_min", outlet_face="x_max", velocity=(1.0, 0.0, 0.0))
solver = aether.UnstructuredCavitySolver3D(domain.mesh, 0.5, wall_velocity, is_outlet, 0.0)
```

Aqui a velocidade **precisa** ter componente ao longo da normal da entrada —
o oposto do que a caixa selada exige, e bem-posto justamente porque existe
uma saída para a mesma massa sair. A checagem a observar não é
`check_closed_domain_conservation` e sim o `net_boundary_flux()` do solver
depois de rodar: a afirmação aqui é "o que entra sai", não "o que entra é
zero".

### Vai rodar?

```python
stability = aether.measure_mesh_stability(domain.mesh, viscosity=0.5)
print(stability.spectral_radius, stability.is_stable)
```

**Não existe critério a-priori honesto de qualidade de malha neste
domínio**, e isso foi medido, não suposto: uma malha com não-ortogonalidade
2,24 roda bem e uma com 2,07 diverge. Então esta sonda não inventa um limiar
sobre métrica-proxy. Ela mede as duas coisas que de fato decidem:

- **repouso** — com todas as paredes paradas o campo tem que ficar
  *exatamente* em zero; qualquer deriva é defeito estrutural;
- **raio espectral** — o fator de crescimento de um passo. Acima de 1 a
  marcha diverge por menor que seja o `dt`.

Custa dezenas de passos contra os milhares de uma execução real. As métricas
de qualidade vêm no relatório como contexto para quem lê, e
deliberadamente **não** entram no veredito.

### Resolver

```python
solver = aether.UnstructuredCavitySolver3D(domain.mesh, 0.5, lid)
run = aether.run_to_steady_state(solver, domain.mesh, max_steps=600)
print(run.summary())
```

Marcha até o campo parar de mudar. A tolerância é *relativa* porque uma
absoluta é frouxa demais ou apertada demais conforme a escala do campo, que
o chamador não deveria precisar saber de antemão.

O relatório traz `converged`, `diverged`, divergência por faces, fluxo
líquido de contorno e velocidade máxima. Três leituras que valem o hábito:

| sinal | o que significa |
|---|---|
| `net_boundary_flux` ≉ 0 numa caixa selada | as condições de contorno não conservam massa |
| `max_velocity` > velocidade da tampa | não convergiu, ou está divergindo |
| `converged=False` com `last_change` grande | precisa de mais passos |

### Exportar

```python
aether.export_result_vtk("resultado.vtk", solver, domain.mesh)
```

Velocidade, pressão e `speed` derivado, em VTK legado que ParaView e VisIt
abrem direto.

**Exportar é a resposta, não um viewer embutido.** Corte, limiar, semeadura
de linhas de corrente e iso-superfície sobre malha irregular é exatamente o
conjunto de operações que um pós-processador geral existe para oferecer, e
que ParaView já faz melhor do que este repositório faria.

---

## 4. Além do básico

### Turbulência

```python
M = aether.UnstructuredCavitySolver3D.TurbulenceModel
solver = aether.UnstructuredCavitySolver3D(domain.mesh, 0.05, lid,
                                            turbulence=M.MIXING_LENGTH)
print(solver.eddy_viscosity(0), solver.wall_distance(0))
```

Comprimento de mistura de Prandtl: algébrico, sem equações de transporte
próprias. É o único fechamento disponível no solver não-estruturado;
k-ε e k-ω SST existem apenas nos solvers estruturados.

### Temperatura

```python
E = aether.UnstructuredCavitySolver3D.EnergyModel
solver.enable_energy(E.PASSIVE, thermal_diffusivity=0.01)
solver.set_wall_temperature(lambda p: p.z > 1.0 - 1e-9, 100.0)
```

- `PASSIVE` — o escoamento carrega o calor, o calor não move o escoamento
  (convecção forçada: resfriar um componente).
- `BOUSSINESQ` — adiciona o empuxo, com `reference_temperature`,
  `thermal_expansion` e `gravity` (convecção natural).

Face de contorno sem temperatura prescrita fica **adiabática**.

### Salvar e recarregar

```python
archive = aether.FieldArchive()
aether.save_tetrahedral_mesh(archive, domain.mesh)
archive.set_field("pressure", [solver.pressure(c) for c in range(domain.mesh.cell_count())])
archive.save("caso.aecf")

reloaded = aether.FieldArchive.load("caso.aecf")
mesh = aether.load_tetrahedral_mesh(reloaded)
```

A malha e os campos viajam no mesmo arquivo de propósito: um checkpoint em
que malha e campos podem se separar é um checkpoint que pode ser
silenciosamente descasado.

**Condições de contorno e configuração do solver ainda não são
serializadas** — reconstrua-as no código ao recarregar.

---

## 5. O que este motor ainda não faz

Registrado aqui para não ser descoberto no meio de um projeto:

- **Escoamento compressível** — só incompressível.
- **Multifásico** — só uma fase.
- **Multi-região / multi-material** — um domínio de fluido por vez.
- **STEP/IGES** — `aether.load_step` importa STEP (ISO 10303-21), mas só o
  subconjunto **BREP facetado** (toda face plana, todo contorno um
  `POLY_LOOP`, sem furos): o caso comum de CAD já exportado como facetas.
  Uma face curva ou a representação tesselada do AP242 é detectada e
  reportada em `StepLoadResult.unsupported_features`, não adivinhada. CAD
  curvo de verdade e IGES continuam de fora — exigem um kernel CAD completo
  (OpenCASCADE), decisão de dependência ainda em aberto.
- **Escoamento interno em duto** — o mesher esculpe o objeto de dentro de
  uma caixa (escoamento externo); o caso oposto, tetraedralizar o interior
  de uma superfície, é outra capacidade.
- **AMR** (refinamento adaptativo) — não existe.
- **GPU** — há um kernel CUDA validado, mas não um solver residente na GPU;
  não acelera nada de verdade ainda.
- **Wheel publicada (PyPI)** — `pip install .`/`pip wheel .` funciona (veja
  a seção 1) e produz uma wheel real e instalável, mas ela só existe depois
  de você mesmo construí-la; publicá-la (e gerar wheels por plataforma via
  CI) é trabalho separado, ainda não feito.

E uma limitação de escala que vale saber de antemão: a tetraedralização é
**O(N²)** no número de pontos. Alguns milhares de pontos rodam em segundos;
10⁵ ainda não é prático.

---

## 6. Onde procurar mais

| arquivo | o que é |
|---|---|
| `examples/flow_around_object.py` | o pipeline inteiro, comentado |
| `README.md` | diário de desenvolvimento, por módulo |
| `ROADMAP.md` | o que vem a seguir e por que nessa ordem |
| `DIVIDA_TECNICA.md` | defeitos conhecidos, com as medições que os fecharam |

Os cabeçalhos das classes em `engine/` são a referência de API real: cada um
explica não só o que a classe faz, mas por que foi feita assim e o que já
foi medido a respeito.
