# Changelog

Formato baseado em [Keep a Changelog](https://keepachangelog.com/pt-BR/1.0.0/).
Este projeto ainda não segue versionamento semântico estrito (pré-1.0): um
`0.x` significa "a API pode mudar", não "incompleto" — a física e a suíte de
testes por trás dela são o que este projeto mede com rigor, não o número da
versão.

Para o estado de desenvolvimento dia a dia (o que foi medido, o que falhou e
por quê), veja `README.md`. Para o que vem a seguir, `ROADMAP.md`. Para
defeitos conhecidos com a medição que os fechou, `DIVIDA_TECNICA.md`. Este
arquivo é só o resumo por release.

## [0.1.0-beta] - 2026-08-28

**Primeiro release marcado com tag.** Não é o começo do projeto — o motor já
tinha solvers estruturados e não-estruturados validados, turbulência (RANS
completo, LES, DES), um pipeline de escoamento externo geometria→malha→solver
e uma suíte de mais de 13 conjuntos de teste antes deste release existir.
O que esta tag marca é o motor deixar de ser "algo que só compila neste
checkout" e passar a ser instalável, com uma prova de integração ponta a
ponta em geometria importada e um primeiro ponto de comparação externo —
os quatro itens que fecham a linha de trabalho para uma beta "palpável".

### Added

- **Empacotamento real**: `pip install .` / `pip wheel .` funcionam.
  `pyproject.toml` novo com `scikit-build-core` orquestrando o
  `CMakeLists.txt` já existente; as 11 extensões compiladas passam a se
  instalar dentro do próprio pacote `aether`. Verificado numa venv nova,
  fora deste repositório, rodando um solver de verdade a partir da wheel
  instalada — não só "a wheel foi criada". Ver `docs/getting-started.md`
  seção 1.
- **Recuperação combinatória de facetas** na tetraedralização restrita
  (`DelaunayTetrahedralization3D::tryFlipCoplanarQuadDiagonal`): troca a
  diagonal de um quadrilátero coplano no lugar, sem ponto de Steiner. Fecha
  o caso concretamente medido (cubo de 8 vértices) parcialmente, e um
  cilindro real (prisma poligonal a 10/24/60 lados) por completo — 100% das
  facetas recuperadas nos três tamanhos.
- **Prova de integração ponta a ponta com geometria importada real**: o
  pipeline completo (`mesh_flow_around_object` → `freestream_boundary` →
  `UnstructuredCavitySolver3D`) rodado sobre um cilindro **importado via
  round-trip STL binário de verdade** (não mantido em memória) — 0 facetas
  perdidas, volume exato, desbalanço de massa de 0,0014% da entrada.
  Registrado como teste permanente em `python/tests/test_pipeline.py`.
- **Primeira comparação contra um benchmark publicado** (Ghia, Ghia & Shin,
  1982), com dados obtidos de fontes reais e cruzados, não de memória. O
  resultado é reportado com honestidade, não maquiado: o solver não
  reproduz a tabela diretamente, e a causa provável (o lid regularizado por
  `sin²`, uma decisão deliberada já existente, não um bug) é quantificada,
  não só sugerida. Ver `python/research/ghia_1982_validation.py`.
- **CI corrigido e verificado passando pela primeira vez na história do
  projeto** — a causa raiz era um `.gitignore` que mantinha
  `engine/testing/` inteiramente fora do controle de versão desde o
  primeiro commit deste repositório.
- Importação de STEP (ISO 10303-21) para o subconjunto BREP facetado.
- Esquema de convecção limitado (van Leer) exposto nos seis solvers 3D
  staggered, com `Central` mantido como default por decisão medida (não
  trocado sem uma comparação que justificasse).

### Known limitations

Lista completa e honesta em `docs/getting-started.md` seção 5. As mais
relevantes pra decidir se este release serve pro seu caso:

- Só escoamento externo incompressível, laminar ou RANS/LES/DES — sem
  compressível, multifásico ou calor conjugado.
- STEP só no subconjunto facetado; sem IGES; sem CAD curvo (exigiria
  OpenCASCADE, decisão de dependência ainda em aberto).
- Tetraedralização é O(N²) — milhares de pontos rodam em segundos, 10⁵
  ainda não é prático.
- GPU tem um kernel validado, mas nenhum solver residente na GPU ainda.
- Wheel funciona localmente; não há publicação no PyPI nem wheels por
  plataforma via CI ainda.
- Sem `LICENSE` neste repositório — decisão do dono do projeto, não tomada
  aqui.
