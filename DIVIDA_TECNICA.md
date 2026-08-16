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

### 1.1 Fluxo de saída inconsistente entre projeção e diagnóstico [criada nesta sessão]

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

**O que fazer**: a projeção precisa corrigir o fluxo de contorno do mesmo jeito
que corrige os internos — aplicar `−dt·∇p·A` às faces de saída dentro da
projeção, e medir o mesmo fluxo corrigido.

### 1.2 Poisson da pressão sem correção não-ortogonal no solver de NS [criada nesta sessão]

`UnstructuredDiffusionSolver` tem correção não-ortogonal por correção
diferida; `UnstructuredCavitySolver3D::applyPoissonOperator()` **não tem**.
São o mesmo operador matemático, com tratamentos diferentes.

**Por que importa**: a Fase 2.2 mediu que, sem essa correção, o erro **estagna**
(ordem 0,10 em vez de convergir). O solver de NS está hoje com a versão que
sabidamente não converge, na equação mais importante que ele resolve.

**O que fazer**: extrair o operador validado da Fase 2.2 e usá-lo, em vez de
manter duas implementações divergentes do mesmo Laplaciano.

### 1.3 Gradiente zero silencioso em célula de estêncil deficiente [criada nesta sessão]

`UnstructuredCavitySolver3D::scalarGradients()` devolve gradiente **zero**
quando a matriz de mínimos quadrados é singular. `UnstructuredDiffusionSolver`
cai para Green-Gauss no mesmo caso.

**Por que não dá para contornar**: gradiente zero não é "sem resposta", é uma
resposta errada e plausível. Numa malha de geometria real, células de estêncil
deficiente aparecem justamente perto de contornos complicados — onde o
gradiente importa mais.

**O que fazer**: usar o mesmo recuo para Green-Gauss, ou recusar a malha
explicitamente. Nunca devolver zero como se fosse resultado.

---

## 2. Duplicação que vai divergir

### 2.1 Dois solvers não-estruturados com a mesma geometria copiada [criada nesta sessão]

`UnstructuredDiffusionSolver` e `UnstructuredCavitySolver3D` têm cópias
independentes de: coeficientes de face (`|A|²/(A·d)`), decomposição
over-relaxed, estêncil de mínimos quadrados, inversa 3×3 simétrica com guarda
de posto, e o laço de CG.

**Por que não dá para contornar**: é exatamente a situação que motivou extrair
`StaggeredCavityBase3D` das seis cavidades 3D antes dos Módulos 9-14 — e a
razão registrada lá vale igual aqui: **portar uma correção seis vezes é a
versão cara do problema**. O item 1.2 acima já é a primeira divergência entre
as duas cópias, com menos de um dia de vida.

**O que fazer**: extrair uma base compartilhada (`UnstructuredFvmBase`) com
geometria de face, gradientes e CG, do mesmo jeito e pelo mesmo motivo que a
base staggered foi extraída. Fazer isso **antes** de escrever um terceiro
solver não-estruturado.

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

1. **1.1** (fluxo de saída) — sozinho, pequeno, causa já identificada
2. **2.1** (extrair base compartilhada) — antes que 1.2 e 1.3 virem três cópias
3. **1.2 e 1.3** — de graça depois do 2.1, já que passam a ter um lugar só
4. **5.2** (bindings) — destrava todo o resto da investigação
5. **3.2** (solução manufaturada) — mede o que hoje é inferido
6. **3.1** (esquema de alta ordem) — só depois que 3.2 der uma régua confiável
7. **4.1** (margem nos estruturados) — independente, pode entrar quando quiser
8. **5.1** (push) — um comando

O item **6** (tetraedralização restrita) fica por último não por prioridade,
mas porque é o único que não se resolve com disciplina — se resolve com
pesquisa.
