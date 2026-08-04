# Pipeline 16S de mel de abelhas sem ferrão

Este repositório reúne os scripts usados no processamento de amplicons 16S rRNA V3–V4, na classificação taxonômica, na construção de objetos `phyloseq` e na produção de análises e figuras. 

Os comentários metodológicos e os parâmetros necessários à reprodução foram mantidos. 


## Estrutura

| Caminho | Finalidade |
|---|---|
| `R/pipeline_bootstrap.R` | Define os caminhos, as pastas de saída e os contratos entre etapas. |
| `R/funcoes_estatisticas_exatas.R` | Reúne funções de testes exatos por enumeração das alocações únicas. |
| `scripts/` | Contém as etapas principais disponibilizadas. |
| `config/environment.example` | Exemplo de variáveis de ambiente para uma execução local. |

## Descrição dos scripts

| Script | Descrição básica |
|---|---|
| `01_conferencia_primers.R` | Confere pareamento dos FASTQ, integridade dos identificadores e ocorrência/orientação dos primers. |
| `02_dada2.R` | Remove primers, filtra reads, estima erros por corrida, infere ASVs, une tabelas, remove quimeras e exporta matrizes e auditorias. |
| `03_auditoria_quimeras_recuperacao_ASVs.R` | Audita ASVs removidas e produz cenários de sensibilidade sem substituir o resultado principal. |
| `04a_silva138.R` | Classifica as ASVs com SILVA 138.2 usando DADA2. |
| `04b_rdp19.R` | Classifica as ASVs com RDP 19 usando DADA2. |
| `04c_Greengenes2.R` | Classifica as ASVs com Greengenes2 e preserva a nomenclatura derivada do GTDB. |
| `04e_GSRDB_V3V4_QIIME2.R` | Executa a classificação com GSR-DB V3–V4 por QIIME 2 e salva a análise principal e a sensibilidade. |
| `05_rblast.R` | Alinha ASVs contra um banco local NCBI 16S e organiza evidências por ASV. |
| `06_comparacao_banco_de_dados.R` | Integra as classificações por prioridade hierárquica e registra conflitos entre bancos. |
| `07_phyloseq.R` | Constrói e valida os objetos `phyloseq`, exporta componentes e calcula resultados descritivos. |
| `10_graficos.R` | Lê resultados oficiais das etapas anteriores e gera as figuras; não usa resultados numéricos fixos nas legendas. |
| `11_faprotax.R` | Realiza a inferência funcional por FAPROTAX. |


## Entradas esperadas

O pipeline pressupõe:

- FASTQ paired-end em `data/raw/` ou no caminho indicado por `PIPELINE_RAW_DIR`;
- `metadados.tsv` na mesma pasta dos dados brutos;
- bancos taxonômicos locais em `bancodados/`;
- Cutadapt disponível no `PATH` ou indicado por `CUTADAPT_BIN`;
- R e os pacotes carregados por cada script;
- QIIME 2 para a etapa GSR-DB;
- BLAST+ e um banco local NCBI 16S para a etapa de alinhamento.

Os bancos e os FASTQ não devem ser enviados ao GitHub. Para cada banco, documente nome, versão, fonte, licença, nome esperado do arquivo e checksum SHA-256.

## Configuração

Execute a partir da raiz do repositório. As variáveis abaixo podem ser ajustadas diretamente no terminal ou carregadas a partir do exemplo em `config/environment.example`.

```bash
export PIPELINE_PROJECT_DIR="$PWD"
export PIPELINE_LIB_DIR="$PWD/R"
export PIPELINE_RAW_DIR="$PWD/data/raw"
export PIPELINE_OUTPUT_DIR="$PWD/results/output_V1"
export PIPELINE_VERSION="V1"
export CUTADAPT_BIN="$(command -v cutadapt)"
```

## Ordem de execução

```bash
Rscript scripts/01_conferencia_primers.R
Rscript scripts/02_dada2.R
Rscript scripts/03_auditoria_quimeras_recuperacao_ASVs.R
Rscript scripts/04a_silva138.R
Rscript scripts/04b_rdp19.R
Rscript scripts/04c_Greengenes2.R
Rscript scripts/04d_BEExact.R
Rscript scripts/04e_GSRDB_V3V4_QIIME2.R
Rscript scripts/05_rblast.R
Rscript scripts/06_comparacao_banco_de_dados.R
Rscript scripts/07_phyloseq.R
Rscript scripts/08_inext.R
Rscript scripts/08a_analises_ecologicas.R
Rscript scripts/09_ancombc2.R
Rscript scripts/10_graficos.R
Rscript scripts/11_faprotax.R
```

O Script 01 contém um checkpoint para confirmar o cenário de filtragem. A decisão deve ser predefinida e registrada antes de prosseguir para as análises a jusante.

## Reprodutibilidade e publicação

Antes de tornar o repositório público:

1. inclua os quatro módulos ausentes listados acima;
2. acrescente um `renv.lock` gerado no ambiente definitivo;
3. publique um metadado anonimizado e seu dicionário de campos;
4. valide a sintaxe de todos os scripts com `parse(file = ...)`;
5. execute o fluxo em uma clonagem limpa;
6. compare os resultados finais com os números relatados no artigo;
7. escolha uma licença e crie uma versão marcada para o artigo.


