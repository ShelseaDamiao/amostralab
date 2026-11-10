# AmostraLab

Calculadora de tamanho de amostra para estudos em Moçambique, uma iniciativa da **ComuniDados Moz**.

O AmostraLab calcula a amostra, distribui-a pelas províncias e distritos, sorteia unidades de forma reprodutível e escreve a justificação para o protocolo, em português e inglês.

## O que faz

| Ferramenta | Para quê |
| --- | --- |
| Prevalência e média | Inquéritos: margem de erro, população finita, efeito de desenho, não resposta |
| Duas proporções e duas médias | Comparar grupos, com ajuste pela linha de base (ANCOVA) |
| RCT por conglomerados | Escolas, aldeias ou unidades sanitárias aleatorizadas (ICC, R², correcção t) |
| Efeito mínimo detectável | Modo inverso: o que um número fixo de conglomerados consegue detectar |
| Estratificação | 11 províncias e 161 distritos e cidades, 4 métodos de alocação, margem por estrato, orçamento e dias de campo, exportação CSV |
| Pesquisa rápida de dados | População do Censo 2017 por distrito, peso na província e no país |
| Aleatorização reprodutível | Sorteio com semente, estratificado, para 2 ou 3 braços, exportação CSV |

## Inspiração e créditos

As fórmulas inspiram-se em obras clássicas de amostragem, ensaios clínicos e avaliação de impacto (Cochran, Kish, Lohr, Neyman, Fleiss, Chow, Cohen, Bloom, Hayes e Moulton, Donner e Klar, Dong e Maynard, Duflo, Glennerster, Kraft, Bruhn e McKenzie, entre outros), citadas na própria plataforma.

O AmostraLab é uma ferramenta de apoio ao planeamento e ao ensino, escrita de raiz. Não reproduz nem substitui essas obras e não tem qualquer ligação aos seus autores ou editoras. As implementações foram simplificadas e adaptadas por nós; qualquer erro é da nossa responsabilidade. Para estudos reais, consulte as obras originais, confirme os parâmetros com dados piloto e compare com ferramentas validadas (PowerUp!, G*Power, `pwr` em R).

## Dados

`data/mz-censo-2017.json` e `data/mz-censo-2017-distritos.csv`: população do IV Recenseamento Geral da População e Habitação 2017 (INE Moçambique, resultados definitivos), compilada por [citypopulation.de](https://www.citypopulation.de/en/mozambique/admin/). Os valores não estão ajustados pela omissão censitária (cerca de 3,7%). Distritos criados depois de 2017 não aparecem em separado.

## Correr localmente

Site estático, sem dependências nem compilação:

```bash
python3 -m http.server 8000
# abrir http://localhost:8000
```

## Publicar no Vercel

Importe o repositório em [vercel.com/new](https://vercel.com/new). Framework preset: **Other**. Sem comando de build; o directório de saída é a raiz. Cada `git push` para `main` publica uma nova versão.

## Licença

MIT © 2026 ComuniDados Moz
