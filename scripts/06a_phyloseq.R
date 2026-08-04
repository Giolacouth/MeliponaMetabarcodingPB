# -----------------------------------------------------------------------------
# Infraestrutura compartilhada: caminhos, I/O, auditoria e metadados.
# -----------------------------------------------------------------------------
.file_args <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grepl("^--file=", commandArgs(trailingOnly = FALSE))
])
.script_dir <- if (length(.file_args)) {
  dirname(normalizePath(.file_args[[1L]], mustWork = FALSE))
} else {
  getwd()
}
.bootstrap_candidates <- unique(Filter(nzchar, c(
  Sys.getenv("PIPELINE_LIB_DIR", unset = ""),
  file.path(.script_dir, "..", "R"),
  file.path(getwd(), "R"),
  file.path(getwd(), "..", "R")
)))
.bootstrap_files <- file.path(.bootstrap_candidates, "pipeline_bootstrap.R")
.bootstrap_files <- .bootstrap_files[file.exists(.bootstrap_files)]
if (!length(.bootstrap_files)) {
  stop("pipeline_bootstrap.R nao localizado; defina PIPELINE_LIB_DIR.", call. = FALSE)
}
source(.bootstrap_files[[1L]], local = .GlobalEnv)
rm(.file_args, .script_dir, .bootstrap_candidates, .bootstrap_files)
run_pipeline_script("06a_phyloseq.R", "phyloseq", function(ctx) {
###############################################################################
# SCRIPT 06 — PHYLOSEQ 16S / ASV
#
#
# Estrutura analitica:
#   1. Construir objeto phyloseq completo com 10 amostras (ps_all10)
#   2. Separar dois conjuntos:
#        ps_core9  — 9 amostras da primeira run (analise principal)
#        ps_plus10 — ps_core9 + fasciculata_auxiliar (analise de sensibilidade)
#   3. Descrever e exportar ps_core9, sem duplicar inferencia do Script 08.
#   4. Descrever e exportar ps_plus10 separadamente. A inferencia ecologica
#      nao filogenetica pertence ao Script 08; a filogenetica, ao Script 07.
#
#
###############################################################################

options(encoding = "UTF-8", stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)     # Objeto central e funcoes de microbioma
  library(Biostrings)   # DNAStringSet / writeXStringSet
  library(ggplot2)      # Visualizacao
  library(dplyr)        # Manipulacao de data.frames
  library(vegan)        # adonis2 / betadisper / permutest
  library(ggrepel)
  # dada2 NAO e necessario no Script 6; o seqtab e lido como RDS pronto.
  # grid e pacote base; grid::unit() e acessivel sem library(grid).
})

VERSAO        <- "4.2_core9_plus10_premissas_sem_rarefacao"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Semente compartilhada para reprodutibilidade das permutacoes.
SEED_GLOBAL <- 1234L

log_msg <- function(msg, tipo = "INFO")
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))

# Funcoes compartilhadas para testes exatos por enumeracao. O arquivo e
# localizado pelo diretorio informado pelo orquestrador ou pelo proprio script.
arq_funcoes_exatas <- file.path(.pipeline_lib_dir, "funcoes_estatisticas_exatas.R")
if (!file.exists(arq_funcoes_exatas)) {
  stop("funcoes_estatisticas_exatas.R nao encontrado em: ", arq_funcoes_exatas, call. = FALSE)
}
sys.source(arq_funcoes_exatas, envir = .GlobalEnv)

cat("=============================================================\n")
cat("PHYLOSEQ 16S / ASV — v", VERSAO, "\n", sep = "")
cat("Data:", DATA_EXECUCAO, "\n")
cat("=============================================================\n\n")

###############################################################################
# 0. PARAMETROS GLOBAIS
###############################################################################

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$stage$root
plot_path <- ctx$stage$figures
raw_path <- ctx$raw_path

# Helper compartilhado para status nativo/introduzido.
helper_metadata <- file.path(.pipeline_lib_dir, "utils_metadata.R")
if (file.exists(helper_metadata)) source(helper_metadata)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path,   recursive = TRUE, showWarnings = FALSE)

plot_core9  <- file.path(plot_path, "core9_primeira_run")
plot_plus10 <- file.path(plot_path, "plus10_sensibilidade_auxiliar")
dir.create(plot_core9,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_plus10, recursive = TRUE, showWarnings = FALSE)

arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
arq_taxa   <- ctx$contracts[["taxa_consenso"]]
arq_contaminantes <- ctx$contracts[["contaminantes"]]
arq_asvmap <- ctx$contracts[["asv_sequences"]]
arq_meta_taxa <- ctx$contracts[["metadata_taxonomia"]]
arq_taxa_gsr07 <- ctx$contracts[["taxa_gsr07_alinhada"]]

arq_meta <- ctx$contracts[["metadata_final"]]
if (!file.exists(arq_meta)) {
  stop("metadata_final.tsv nao encontrado em ", output_path, ". Execute o Script 01 antes do Script 06.")
}

# Contrato do desenho experimental.
RUN_PRINCIPAL <- "run_main"
RUN_AUXILIAR  <- "run_aux"
PADRAO_AUXILIAR <- "S10"

# Ranks taxonomicos esperados (ordem canonica DADA2 / SILVA)
RANKS_ESPERADOS <- c("Kingdom", "Phylum", "Class", "Order",
                     "Family", "Genus", "Species")

# Contrato estrito: Script 5 e Script 6 devem representar exatamente o mesmo
# universo apos remover os contaminantes. Nenhuma perda silenciosa e aceita.
LIMIAR_TAXA_MISS <- 0

# Top N para graficos
TOP_N_PHYLUM <- 10L
TOP_N_GENUS  <- 15L
TOP_N_HEAT   <- 30L
N_PERM <- 9999L
###############################################################################
# 1. FUNCOES AUXILIARES
# Funcoes puras (sem side-effects) e funcoes de I/O bem delimitadas.
###############################################################################

# Predicado: sequencia DNA 16S valida (>= 50 nt, alfabeto IUPAC)
eh_dna <- function(x)
  grepl("^[ACGTURYKMSWBDHVNacgturykmswbdhvn.-]+$", x) & nchar(x) >= 50

# Matriz de ASVs com táxons nas linhas (independente da orientacao do objeto)
matriz_taxa_linhas <- function(ps_obj) {
  mat <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) mat <- t(mat)
  mat
}

# Remover contaminantes definidos exclusivamente pelo Script 5.
# O Script 6 nao decide o que e contaminante; apenas aplica a lista exportada.
remover_contaminantes_seqtab <- function(seqtab, arq_contaminantes, output_path) {
  if (!file.exists(arq_contaminantes)) {
    stop("asvs_contaminantes_excluir.csv ausente. Execute o Script 5 antes do Script 6.")
  }

  cont <- read.csv(arq_contaminantes, stringsAsFactors = FALSE)
  if (!all(c("ASV_ID", "ASV_seq", "Motivo") %in% colnames(cont))) {
    stop("asvs_contaminantes_excluir.csv deve conter ASV_ID, ASV_seq e Motivo.")
  }

  reads_antes <- sum(seqtab)

  if (nrow(cont) == 0) {
    write.csv(data.frame(Reads_antes = reads_antes,
                         Reads_removidas = 0,
                         Reads_depois = reads_antes,
                         Pct_reads_removidas = 0),
              file.path(output_path, "resumo_remocao_contaminantes_script6.csv"),
              row.names = FALSE)
    log_msg("Arquivo de contaminantes vazio; nenhum contaminante removido.", "OK")
    return(seqtab)
  }

  seqs_remover <- unique(cont$ASV_seq)
  seqs_presentes <- intersect(colnames(seqtab), seqs_remover)
  reads_removidas <- sum(seqtab[, seqs_presentes, drop = FALSE])

  seqtab_filtrado <- seqtab[, !(colnames(seqtab) %in% seqs_remover), drop = FALSE]
  reads_depois <- sum(seqtab_filtrado)

  write.csv(data.frame(ASV_ID = cont$ASV_ID[match(seqs_presentes, cont$ASV_seq)],
                       ASV_seq = seqs_presentes,
                       Motivo = "removida_antes_do_phyloseq_por_lista_do_script5",
                       stringsAsFactors = FALSE),
            file.path(output_path, "asvs_contaminantes_removidas_script6.csv"),
            row.names = FALSE)

  write.csv(data.frame(Reads_antes = reads_antes,
                       Reads_removidas = reads_removidas,
                       Reads_depois = reads_depois,
                       Pct_reads_removidas = round(100 * reads_removidas / reads_antes, 4)),
            file.path(output_path, "resumo_remocao_contaminantes_script6.csv"),
            row.names = FALSE)

  log_msg(sprintf("%d ASVs contaminantes removidas; %.4f%% dos reads removidos.",
                  length(seqs_presentes), 100 * reads_removidas / reads_antes), "OK")

  seqtab_filtrado
}
# Transformacao para abundancia relativa por amostra
transformar_relativo <- function(ps_obj)
  transform_sample_counts(ps_obj, function(x) if (sum(x) > 0) x / sum(x) else x)

# Salvar ggplot em PDF (ponto unico de chamada a ggsave)
salvar_plot <- function(p, arquivo, dir_saida, width = 8, height = 6)
  ggsave(file.path(dir_saida, arquivo), plot = p, width = width, height = height)

# Exportar distancia como CSV e RDS
exportar_distancia <- function(dist_obj, nome, prefixo) {
  write.csv(
    as.matrix(dist_obj),
    file.path(output_path, paste0(prefixo, "_beta_", nome, ".csv"))
  )
  saveRDS(dist_obj,
          file.path(output_path, paste0(prefixo, "_dist_", nome, ".rds")))
}

# Exportar matriz de contagens de ASVs, taxonomia, metadados e FASTA de um objeto phyloseq
exportar_componentes_phyloseq <- function(ps_obj, prefixo) {

  otu_df <- as.data.frame(matriz_taxa_linhas(ps_obj))
  otu_df <- cbind(ASV_ID = rownames(otu_df), otu_df)
  write.csv(otu_df,
            file.path(output_path, paste0(prefixo, "_asv_count_table.csv")),
            row.names = FALSE)

  taxa_df <- cbind(
    ASV_ID = taxa_names(ps_obj),
    as.data.frame(tax_table(ps_obj), stringsAsFactors = FALSE)
  )
  write.csv(taxa_df,
            file.path(output_path, paste0(prefixo, "_taxa_table.csv")),
            row.names = FALSE)

  write.csv(
    as(sample_data(ps_obj), "data.frame"),
    file.path(output_path, paste0(prefixo, "_sample_data.csv")),
    row.names = TRUE)

  Biostrings::writeXStringSet(
    refseq(ps_obj),
    filepath = file.path(output_path, paste0(prefixo, "_refseq_ASV.fasta"))
  )

  log_msg(paste("Componentes exportados:", prefixo), "SAVE")
}

# Validacao minima de um objeto phyloseq.
# ASVs zeradas devem ser removidas e auditadas ANTES desta funcao. Aqui elas
# interrompem a execucao para impedir uma poda silenciosa em objetos inesperados.
validar_phyloseq_basico <- function(ps_obj, nome_obj) {
  stopifnot(inherits(ps_obj, "phyloseq"))
  if (nsamples(ps_obj) < 2)
    stop(nome_obj, ": menos de 2 amostras.")
  if (ntaxa(ps_obj) < 2)
    stop(nome_obj, ": menos de 2 ASVs.")
  if (any(sample_sums(ps_obj) == 0))
    stop(nome_obj, ": amostra(s) com zero reads.")
  if (!all(grepl("^ASV_", taxa_names(ps_obj))))
    log_msg(paste(nome_obj, ": taxa_names fora do padrao ASV_"), "WARN")
  if (!"BeeSpecies" %in% sample_variables(ps_obj))
    stop(nome_obj, ": variavel BeeSpecies ausente.")
  zero_taxa <- taxa_names(ps_obj)[taxa_sums(ps_obj) == 0]
  if (length(zero_taxa) > 0L) {
    stop(
      nome_obj, ": ", length(zero_taxa),
      " ASV(s) com soma zero. Remova-as explicitamente e registre a decisao antes da validacao: ",
      paste(head(zero_taxa, 20L), collapse = ";")
    )
  }
  ps_obj
}

###############################################################################
# 2. VALIDACOES DE ARQUIVOS
###############################################################################

cat("=== VALIDACOES DE ARQUIVOS ===\n\n")

for (arq in c(arq_seqtab, arq_taxa, arq_meta, arq_contaminantes,
                arq_asvmap, arq_meta_taxa)) {
  if (!file.exists(arq)) {
    stop(
      "Arquivo ausente: ", arq, "\n",
      "Execute o Script 5 antes do Script 6. O Script 5 deve gerar ",
      "asvs_contaminantes_excluir.csv, mesmo que vazio."
    )
  }
}

log_msg("seqtab_global_nochim.rds OK", "OK")
log_msg("taxa_consenso_final.rds OK",  "OK")
log_msg("metadata_consenso_taxonomico.csv OK", "OK")
log_msg(paste("metadata_final.tsv OK:", arq_meta), "OK")
cat("\n")

###############################################################################
# 2B. METADADOS DA TAXONOMIA
###############################################################################

meta_taxa_exec <- read.csv(
  arq_meta_taxa,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(meta_taxa_exec) < 1L ||
    !"Descricao_hierarquia" %in% colnames(meta_taxa_exec)) {
  stop(
    "metadata_consenso_taxonomico.csv invalido. Execute o Script 05 revisado."
  )
}

descricao_hierarquia_taxa <- as.character(
  meta_taxa_exec$Descricao_hierarquia[1L]
)

log_msg(
  paste("Hierarquia taxonomica:", descricao_hierarquia_taxa),
  "INFO"
)

###############################################################################
# 3. SEQTAB
###############################################################################

cat("=== SEQTAB ===\n\n")

seqtab <- tryCatch(
  as.matrix(readRDS(arq_seqtab)),
  error = function(e) {
    stop(
      "Falha ao ler seqtab_global_nochim.rds: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)

# Validacoes estruturais
if (nrow(seqtab) == 0L || ncol(seqtab) == 0L) {
  stop("seqtab: matriz vazia.", call. = FALSE)
}
if (!is.numeric(seqtab)) {
  stop("seqtab: as contagens devem ser numericas.", call. = FALSE)
}

if (is.null(rownames(seqtab)) || any(rownames(seqtab) == ""))
  stop("seqtab: nomes de amostras ausentes ou vazios.")

if (is.null(colnames(seqtab)) || any(colnames(seqtab) == ""))
  stop("seqtab: colnames nao contem sequencias validas.")

if (anyDuplicated(rownames(seqtab)) > 0)
  stop("seqtab: SampleIDs duplicados.")

if (anyDuplicated(colnames(seqtab)) > 0)
  stop("seqtab: sequencias ASV duplicadas.")

if (!all(eh_dna(colnames(seqtab))))
  stop("seqtab: colnames contem strings que nao parecem sequencias 16S IUPAC.")

if (any(is.na(seqtab)) || any(!is.finite(seqtab)))
  stop("seqtab: contem valores NA ou nao finitos.")

if (any(seqtab < 0))
  stop("seqtab: contem contagens negativas.")

if (max(abs(seqtab - round(seqtab))) > 1e-8)
  stop("seqtab: contagens nao inteiras; o objeto canonico deve conter reads brutas.")

if (any(rowSums(seqtab) == 0))
  stop("seqtab: amostra(s) com zero reads.")

n_zero_col <- sum(colSums(seqtab) == 0)

if (n_zero_col > 0) {
  stop(
    "seqtab canonico contem ", n_zero_col,
    " ASV(s) com soma zero. Isso indica contrato quebrado no Script 01; nenhuma ASV sera removida silenciosamente."
  )
}

cat(sprintf(
  "Seqtab antes da remocao de contaminantes: %d amostras | %d ASVs\n",
  nrow(seqtab),
  ncol(seqtab)
))

# Remover contaminantes identificados no Script 5.
# Esta etapa precisa ocorrer ANTES de criar ASV_IDs.
seqtab <- remover_contaminantes_seqtab(seqtab, arq_contaminantes, output_path)

if (ncol(seqtab) == 0) {
  stop("Todas as ASVs foram removidas apos exclusao de contaminantes.")
}

if (any(rowSums(seqtab) == 0)) {
  amostras_zero <- rownames(seqtab)[rowSums(seqtab) == 0]
  stop(
    "Apos remover contaminantes, as seguintes amostras ficaram com zero reads:\n",
    paste(amostras_zero, collapse = "\n")
  )
}

cat(sprintf(
  "Seqtab final para phyloseq: %d amostras | %d ASVs\n",
  nrow(seqtab),
  ncol(seqtab)
))

# Preservar IDs ASV canonicos do Script 01 apos a remocao de contaminantes.
asv_seqs <- colnames(seqtab)
asv_map <- tryCatch(
  read.delim(
    arq_asvmap,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  ),
  error = function(e) {
    stop(
      "Falha ao ler ASV_sequences.tsv: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)
if (!all(c("ASV_ID", "Sequence", "Origem") %in% colnames(asv_map))) {
  stop("ASV_sequences.tsv deve conter ASV_ID, Sequence e Origem.", call. = FALSE)
}
asv_map$ASV_ID <- trimws(as.character(asv_map$ASV_ID))
asv_map$Sequence <- toupper(trimws(as.character(asv_map$Sequence)))
asv_map$Origem <- trimws(as.character(asv_map$Origem))
if (anyNA(asv_map$ASV_ID) || any(asv_map$ASV_ID == "") ||
    anyNA(asv_map$Sequence) || any(asv_map$Sequence == "") ||
    anyNA(asv_map$Origem) || any(asv_map$Origem == "")) {
  stop("ASV_sequences.tsv contem valor obrigatorio ausente ou vazio.", call. = FALSE)
}
if (anyDuplicated(asv_map$ASV_ID) > 0 || anyDuplicated(asv_map$Sequence) > 0) {
  stop("ASV_sequences.tsv contem IDs ou sequencias duplicadas.", call. = FALSE)
}
if (!all(eh_dna(asv_map$Sequence))) {
  stop("ASV_sequences.tsv contem sequencia 16S/IUPAC invalida.", call. = FALSE)
}
idx_map <- match(asv_seqs, asv_map$Sequence)
if (any(is.na(idx_map)))
  stop("Ha ASVs do seqtab sem correspondencia no mapa canonico ASV_sequences.tsv.")
asv_ids <- asv_map$ASV_ID[idx_map]
seq2id  <- setNames(asv_ids, asv_seqs)
id2seq  <- setNames(asv_seqs, asv_ids)

otu_tab <- otu_table(seqtab, taxa_are_rows = FALSE)
colnames(otu_tab) <- asv_ids

write.csv(data.frame(ASV_ID = asv_ids,
                     ASV_seq = asv_seqs,
                     Origem = asv_map$Origem[idx_map],
                     stringsAsFactors = FALSE),
          file.path(output_path, "ASV_ID_map_script6_preservado.csv"),
          row.names = FALSE)

cat(sprintf("IDs ASV preservados do Script 01: %s ... %s\n\n", asv_ids[1], tail(asv_ids, 1)))

###############################################################################
# 4. TAXONOMIA
###############################################################################

cat("=== TAXONOMIA ===\n\n")

taxa_mat <- tryCatch(
  as.matrix(readRDS(arq_taxa)),
  error = function(e) {
    stop(
      "Falha ao ler taxa_consenso_final.rds: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)
if (nrow(taxa_mat) == 0L || ncol(taxa_mat) == 0L) {
  stop("taxa_consenso_final.rds esta vazio.", call. = FALSE)
}

if (is.null(rownames(taxa_mat)) ||
    anyNA(rownames(taxa_mat)) ||
    any(rownames(taxa_mat) == "")) {
  stop("taxa_consenso_final.rds: rownames ausentes ou vazios.", call. = FALSE)
}
if (anyDuplicated(rownames(taxa_mat)) > 0) {
  stop("taxa_consenso_final.rds: rownames duplicados.", call. = FALSE)
}
if (is.null(colnames(taxa_mat)) ||
    anyNA(colnames(taxa_mat)) ||
    any(colnames(taxa_mat) == "") ||
    anyDuplicated(colnames(taxa_mat)) > 0L) {
  stop("taxa_consenso_final.rds: ranks ausentes, vazios ou duplicados.", call. = FALSE)
}

# Garantir que todos os ranks canonicos existam
ranks_ausentes <- setdiff(RANKS_ESPERADOS, colnames(taxa_mat))
if (length(ranks_ausentes) > 0) {
  add_mat <- matrix(NA_character_, nrow(taxa_mat), length(ranks_ausentes),
                    dimnames = list(rownames(taxa_mat), ranks_ausentes))
  taxa_mat <- cbind(taxa_mat, add_mat)
  log_msg(sprintf("Ranks adicionados como NA: %s",
                  paste(ranks_ausentes, collapse = ", ")), "WARN")
}
taxa_mat <- taxa_mat[, RANKS_ESPERADOS, drop = FALSE]

# Contrato estrito: o Script 05 salva rownames como sequencias nucleotidicas.
# Nao aceitar fallback por ASV_ID, pois isso poderia mascarar mistura de execucoes.
idx_match <- match(asv_seqs, rownames(taxa_mat))
tipo_match <- "sequencia"

n_miss   <- sum(is.na(idx_match))
pct_miss <- n_miss / length(idx_match)

log_msg(sprintf("Correspondencia por %s: %d/%d ASVs com match; %d sem match (%.1f%%)",
                tipo_match, sum(!is.na(idx_match)), length(idx_match),
                n_miss, 100 * pct_miss),
        if (n_miss == 0) "OK" else "WARN")

if (pct_miss > LIMIAR_TAXA_MISS)
  stop(sprintf(
    "%.1f%% das ASVs sem correspondencia taxonomica (contrato esperado: 0%%).\n",
    100 * pct_miss),
    "Verificar se seqtab_global_nochim.rds, taxa_consenso_final.rds e ",
    "asvs_contaminantes_excluir.csv provem da mesma execucao do Script 05.")

chaves_taxa_validas <- if (tipo_match == "sequencia") asv_seqs else asv_ids
extras_taxa <- setdiff(rownames(taxa_mat), chaves_taxa_validas)
if (length(extras_taxa) > 0L) {
  stop(
    "taxa_consenso_final.rds contem ", length(extras_taxa),
    " linha(s) extra(s) apos a remocao dos contaminantes."
  )
}


taxa_ord <- matrix(NA_character_, length(asv_seqs), length(RANKS_ESPERADOS),
                   dimnames = list(asv_ids, RANKS_ESPERADOS))
validos <- !is.na(idx_match)
if (any(validos))
  taxa_ord[validos, ] <- taxa_mat[idx_match[validos], , drop = FALSE]

# Auditoria de ASVs sem correspondencia taxonomica
if (n_miss > 0) {
  asvs_miss <- data.frame(
    ASV_ID   = asv_ids[!validos],
    Sequence = asv_seqs[!validos],
    Motivo   = "Presente em seqtab; ausente em taxa_consenso_final",
    stringsAsFactors = FALSE)
  write.csv(asvs_miss,
            file.path(output_path, "asvs_sem_taxonomia.csv"), row.names = FALSE)
  log_msg(sprintf("%d ASVs auditadas em asvs_sem_taxonomia.csv", n_miss), "WARN")
}

taxa_ps <- tax_table(taxa_ord)

cat(sprintf("Ranks: %s\n", paste(rank_names(taxa_ps), collapse = ", ")))
cat(sprintf("ASVs com Genus   classificado: %d / %d\n",
            sum(!is.na(taxa_ps[, "Genus"])),   nrow(taxa_ps)))
cat(sprintf("ASVs com Species classificada: %d / %d\n\n",
            sum(!is.na(taxa_ps[, "Species"])), nrow(taxa_ps)))

###############################################################################
# 5. METADADOS
###############################################################################

cat("=== METADADOS ===\n\n")

meta <- read.table(arq_meta, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE,
                   na.strings = c("", "NA", "NaN"))

# Validacoes estruturais obrigatorias
if (!"SampleID" %in% colnames(meta))
  stop("metadata_final.tsv: coluna SampleID ausente.")
meta$SampleID <- trimws(as.character(meta$SampleID))
if (anyNA(meta$SampleID) || any(meta$SampleID == ""))
  stop("metadata_final.tsv: SampleID ausente ou vazio.")
if (anyDuplicated(meta$SampleID) > 0)
  stop("SampleID duplicado: ",
       paste(unique(meta$SampleID[duplicated(meta$SampleID)]), collapse = ", "))
if (!"BeeSpecies" %in% colnames(meta) && !"Species" %in% colnames(meta))
  stop("metadata_final.tsv: coluna BeeSpecies (ou alias legado Species) ausente.")
if (!"Run" %in% colnames(meta))
  stop("metadata_final.tsv: coluna Run ausente.")
if (!"SampleLabel" %in% colnames(meta)) {
  stop(
    "metadata_final.tsv: coluna SampleLabel ausente. ",
    "O SampleLabel e necessario para rotular os graficos."
  )
}

meta$SampleLabel <- trimws(as.character(meta$SampleLabel))

if (any(is.na(meta$SampleLabel)) || any(meta$SampleLabel == "")) {
  stop("metadata_final.tsv: existem valores ausentes ou vazios em SampleLabel.")
}

if (anyDuplicated(meta$SampleLabel) > 0) {
  duplicados <- unique(meta$SampleLabel[duplicated(meta$SampleLabel)])

  stop(
    "metadata_final.tsv: SampleLabel duplicado: ",
    paste(duplicados, collapse = ", "),
    ". Cada amostra deve possuir um SampleLabel unico."
  )
}


# colidir com Species bacteriana da tax_table em psmelt().
if (!"BeeSpecies" %in% colnames(meta)) meta$BeeSpecies <- meta$Species
if ("Species" %in% colnames(meta)) meta$Species <- NULL

# Verificar cobertura entre seqtab e metadados
amostras_seqtab <- rownames(seqtab)
faltando_meta   <- setdiff(amostras_seqtab, meta$SampleID)
if (length(faltando_meta) > 0)
  stop("Amostras no seqtab sem metadados:\n", paste(faltando_meta, collapse = "\n"))

extra_meta <- setdiff(meta$SampleID, amostras_seqtab)
if (length(extra_meta) > 0) {
  stop(
    "metadata_final.tsv contem amostra(s) ausente(s) do seqtab; nenhuma sera ignorada silenciosamente: ",
    paste(extra_meta, collapse = ", ")
  )
}

rownames(meta) <- meta$SampleID
meta <- meta[amostras_seqtab, , drop = FALSE]

if (!identical(rownames(meta), rownames(seqtab)))
  stop("Ordem das amostras diverge entre metadados e seqtab apos reindexacao.")
if (nrow(meta) != 10L ||
    anyNA(meta$Run) ||
    any(trimws(as.character(meta$Run)) == "") ||
    sum(as.character(meta$Run) == RUN_PRINCIPAL) != 9L ||
    sum(as.character(meta$Run) == RUN_AUXILIAR) != 1L) {
  stop("metadata_final.tsv nao respeita o desenho de duas corridas (9 + 1).")
}

# Validar BeeSpecies
if (any(is.na(meta$BeeSpecies)))
  stop("BeeSpecies contem NA.")
meta$BeeSpecies <- factor(meta$BeeSpecies, levels = c(
  "Melipona fasciculata", "Melipona scutellaris", "Melipona subnitida"))
if (any(is.na(meta$BeeSpecies)))
  stop("BeeSpecies contem nomes divergentes dos niveis esperados.")

# Variavel derivada para responder a pergunta nativo/introduzido.
if (exists("adicionar_status_origem")) {
  meta <- adicionar_status_origem(meta, "BeeSpecies")
} else {
  meta$Nativo_Introduzido <- factor(ifelse(as.character(meta$BeeSpecies) == "Melipona fasciculata",
                                           "Introduzida", "Nativa"),
                                    levels = c("Introduzida", "Nativa"))
}

meta$Run <- factor(meta$Run)

# Identificar amostra da segunda run
is_auxiliar <- grepl(PADRAO_AUXILIAR, meta$SampleID)
if (sum(is_auxiliar) != 1)
  stop("Esperava exatamente 1 amostra com padrao '", PADRAO_AUXILIAR,
       "'; encontrado: ", sum(is_auxiliar), ".")

# Variavel de bookkeeping para auditoria (nao usada em testes)
meta$AnalysisSet <- factor(
  ifelse(is_auxiliar, "plus10_auxiliar", "core9_primeira_run"),
  levels = c("core9_primeira_run", "plus10_auxiliar"))

# Verificacao de coordenadas sem ponto decimal (ex: -7027003 vs -7.027003)
for (col_coord in c("latitude","longitude","lat","long","lon",
                    "Latitude","Longitude")) {
  if (!col_coord %in% colnames(meta)) next
  vals <- suppressWarnings(as.numeric(meta[[col_coord]]))
  eh_la <- col_coord %in% c("latitude","lat","Latitude")
  prob  <- !is.na(vals) & ((eh_la & abs(vals) > 90) | (!eh_la & abs(vals) > 180))
  if (any(prob))
    log_msg(paste0("'", col_coord,
                   "' parece ter coordenadas sem ponto decimal."), "WARN")
}

meta_ps <- sample_data(meta)

cat(sprintf("Metadados: %d amostras x %d variaveis\n", nrow(meta), ncol(meta)))
cat("BeeSpecies por Run:\n");         print(table(meta$BeeSpecies, meta$Run))
cat("\nNativo_Introduzido:\n");       print(table(meta$Nativo_Introduzido))
cat("\nConjuntos analiticos:\n");     print(table(meta$AnalysisSet))
cat("\nAmostra auxiliar:", meta$SampleID[is_auxiliar], "\n\n")

###############################################################################
# 6. REFSEQS
###############################################################################

cat("=== REFSEQS ===\n\n")

ref_seqs <- Biostrings::DNAStringSet(asv_seqs)
names(ref_seqs) <- asv_ids

cat(sprintf("refseqs: %d sequencias | comprimento: %d -- %d pb\n\n",
            length(ref_seqs), min(width(ref_seqs)), max(width(ref_seqs))))

###############################################################################
# 7. CONSTRUIR PS_ALL10
###############################################################################

cat("=== CONSTRUINDO PS_ALL10 ===\n\n")

ps_all10 <- phyloseq(otu_tab, taxa_ps, meta_ps, ref_seqs)
if (any(taxa_sums(ps_all10) == 0)) {
  stop("ps_all10 contem ASV(s) com soma zero; isso indica falha de integracao.")
}
if (any(sample_sums(ps_all10) == 0)) {
  stop("ps_all10 contem amostra(s) com zero reads; nenhuma amostra sera removida silenciosamente.")
}
ps_all10 <- validar_phyloseq_basico(ps_all10, "ps_all10")

saveRDS(ps_all10, file.path(output_path, "phyloseq_meliponini_completo_10.rds"))
log_msg("phyloseq_meliponini_completo_10.rds salvo (nao rarefado)", "SAVE")
cat(sprintf("ps_all10: %d amostras | %d ASVs\n\n",
            nsamples(ps_all10), ntaxa(ps_all10)))

###############################################################################
# 8. DEFINIR CORE9 E PLUS10
###############################################################################

cat("=== DEFININDO CONJUNTOS ANALITICOS ===\n\n")

sd_all <- as(sample_data(ps_all10), "data.frame")

sample_auxiliar <- sd_all$SampleID[grepl(PADRAO_AUXILIAR, sd_all$SampleID)]

if (length(sample_auxiliar) != 1) {
  stop("Falha ao identificar exatamente uma amostra auxiliar no objeto phyloseq.")
}

run_auxiliar <- as.character(sd_all$Run[sd_all$SampleID == sample_auxiliar])

if (length(run_auxiliar) != 1L || is.na(run_auxiliar) || run_auxiliar != RUN_AUXILIAR) {
  stop(
    "A amostra auxiliar deve pertencer exclusivamente a ", RUN_AUXILIAR,
    "; valor encontrado: ", paste(run_auxiliar, collapse = ", ")
  )
}

runs_disponiveis <- sort(unique(as.character(sd_all$Run)))
runs_esperadas <- sort(c(RUN_PRINCIPAL, RUN_AUXILIAR))
if (!identical(runs_disponiveis, runs_esperadas)) {
  stop(
    "Runs do phyloseq divergem do desenho experimental. Esperadas: ",
    paste(runs_esperadas, collapse = ", "), "; encontradas: ",
    paste(runs_disponiveis, collapse = ", ")
  )
}

samples_core9 <- rownames(sd_all)[as.character(sd_all$Run) == RUN_PRINCIPAL]

if (length(samples_core9) != 9) {
  stop(
    "A analise principal deveria conter 9 amostras da primeira run, mas encontrou ",
    length(samples_core9), "."
  )
}

ps_core9 <- prune_samples(samples_core9, ps_all10)
asvs_zero_core9 <- taxa_names(ps_core9)[taxa_sums(ps_core9) == 0]
write.csv(
  data.frame(ASV_ID = asvs_zero_core9, Motivo = "ausente_na_corrida_principal",
             stringsAsFactors = FALSE),
  file.path(output_path, "core9_asvs_exclusivas_corrida_auxiliar.csv"),
  row.names = FALSE
)
if (length(asvs_zero_core9) > 0L) {
  ps_core9 <- prune_taxa(taxa_sums(ps_core9) > 0, ps_core9)
  log_msg(sprintf("core9: %d ASV(s) exclusivas da corrida auxiliar removidas explicitamente.",
                  length(asvs_zero_core9)), "INFO")
}
ps_core9 <- validar_phyloseq_basico(ps_core9, "ps_core9")

ps_plus10 <- validar_phyloseq_basico(ps_all10, "ps_plus10")

saveRDS(ps_core9,  file.path(output_path, "phyloseq_core9_primeira_run.rds"))
saveRDS(ps_plus10, file.path(output_path, "phyloseq_plus10_com_auxiliar.rds"))

# ---------------------------------------------------------------------------
# 8B. OBJETOS PARA SENSIBILIDADE TAXONOMICA GSR 0.7
#
# A matriz de contagens, metadados, refseqs e lista de contaminantes permanecem
# identicos ao consenso principal. 
# ---------------------------------------------------------------------------

gsr_taxonomia_disponivel <- file.exists(arq_taxa_gsr07)
ps_gsr_all10 <- NULL
ps_gsr_core9 <- NULL
ps_gsr_plus10 <- NULL

if (gsr_taxonomia_disponivel) {
  taxa_gsr <- as.matrix(readRDS(arq_taxa_gsr07))

  if (is.null(rownames(taxa_gsr)) || anyDuplicated(rownames(taxa_gsr)) > 0L) {
    stop("taxa_gsr07_alinhada_analise.rds sem rownames validos.")
  }

  ranks_gsr_ausentes <- setdiff(RANKS_ESPERADOS, colnames(taxa_gsr))
  if (length(ranks_gsr_ausentes) > 0L) {
    add_gsr <- matrix(
      NA_character_,
      nrow = nrow(taxa_gsr),
      ncol = length(ranks_gsr_ausentes),
      dimnames = list(rownames(taxa_gsr), ranks_gsr_ausentes)
    )
    taxa_gsr <- cbind(taxa_gsr, add_gsr)
  }
  taxa_gsr <- taxa_gsr[, RANKS_ESPERADOS, drop = FALSE]

  # O Script 05 alinha a GSR e preserva rownames como sequencias. Nao aceitar
  # fallback por ASV_ID, pois isso poderia ocultar mistura de execucoes.
  idx_gsr <- match(asv_seqs, rownames(taxa_gsr))

  if (anyNA(idx_gsr)) {
    stop(
      "GSR 0.7 nao corresponde ao universo apos remocao de contaminantes: ",
      sum(is.na(idx_gsr)), " ASV(s) sem match."
    )
  }

  taxa_gsr_ord <- taxa_gsr[idx_gsr, RANKS_ESPERADOS, drop = FALSE]
  rownames(taxa_gsr_ord) <- asv_ids

  ps_gsr_all10 <- phyloseq(
    otu_tab,
    tax_table(taxa_gsr_ord),
    meta_ps,
    ref_seqs
  )
  if (any(taxa_sums(ps_gsr_all10) == 0) || any(sample_sums(ps_gsr_all10) == 0)) {
    stop("ps_gsr_all10 contem taxa ou amostra com soma zero; falha de integracao GSR.")
  }
  ps_gsr_all10 <- validar_phyloseq_basico(ps_gsr_all10, "ps_gsr_all10")

  ps_gsr_core9 <- prune_samples(samples_core9, ps_gsr_all10)
  zero_gsr_core9 <- taxa_names(ps_gsr_core9)[taxa_sums(ps_gsr_core9) == 0]
  if (!setequal(zero_gsr_core9, asvs_zero_core9)) {
    stop(
      "A taxonomia GSR gerou conjunto inesperado de ASVs zeradas no core9. ",
      "Consenso=", length(asvs_zero_core9), "; GSR=", length(zero_gsr_core9)
    )
  }
  if (length(zero_gsr_core9) > 0L) {
    ps_gsr_core9 <- prune_taxa(taxa_sums(ps_gsr_core9) > 0, ps_gsr_core9)
  }
  ps_gsr_core9 <- validar_phyloseq_basico(ps_gsr_core9, "ps_gsr_core9")

  ps_gsr_plus10 <- ps_gsr_all10

  saveRDS(
    ps_gsr_all10,
    file.path(output_path, "phyloseq_all10_gsr07_sensibilidade.rds")
  )
  saveRDS(
    ps_gsr_core9,
    file.path(output_path, "phyloseq_core9_gsr07_sensibilidade.rds")
  )
  saveRDS(
    ps_gsr_plus10,
    file.path(output_path, "phyloseq_plus10_gsr07_sensibilidade.rds")
  )

  log_msg(
    "Objetos phyloseq GSR 0.7 salvos para sensibilidade taxonomica.",
    "SAVE"
  )
} else {
  log_msg(
    paste(
      "Taxonomia GSR 0.7 ausente; objetos phyloseq de sensibilidade",
      "taxonomica nao foram criados."
    ),
    "WARN"
  )
}

cat("Amostra auxiliar: unidade amostral distinta, pertencente a segunda run:\n")
cat(sample_auxiliar, "\n\n")

cat("Run principal usada na analise estatistica:\n")
cat(RUN_PRINCIPAL, "\n\n")

cat("Run da amostra auxiliar:\n")
cat(run_auxiliar, "\n\n")

cat(sprintf("ps_core9 : %d amostras | %d ASVs (analise principal)\n",
            nsamples(ps_core9), ntaxa(ps_core9)))
cat(sprintf("ps_plus10: %d amostras | %d ASVs (sensibilidade)\n\n",
            nsamples(ps_plus10), ntaxa(ps_plus10)))

cat("Distribuicao core9:\n")
print(table(sample_data(ps_core9)$BeeSpecies, sample_data(ps_core9)$Run))

cat("\nDistribuicao plus10:\n")
print(table(sample_data(ps_plus10)$BeeSpecies, sample_data(ps_plus10)$Run))
cat("\n")
###############################################################################
# 9. FUNCAO PRINCIPAL: ANALISE EM FUNIL INVERTIDO
#
###############################################################################

obter_mapa_samplelabel <- function(ps_obj) {

  meta_df <- as(
    phyloseq::sample_data(ps_obj),
    "data.frame"
  )

  if (!"SampleLabel" %in% colnames(meta_df)) {
    stop("SampleLabel ausente no sample_data do objeto phyloseq.")
  }

  ids <- rownames(meta_df)
  labels <- trimws(as.character(meta_df$SampleLabel))

  if (any(is.na(labels)) || any(labels == "")) {
    stop("SampleLabel possui valores ausentes ou vazios.")
  }

  if (anyDuplicated(labels) > 0) {
    stop("SampleLabel deve ser unico para cada amostra.")
  }

  setNames(labels, ids)
}


aplicar_samplelabel <- function(x, mapa) {

  x <- as.character(x)
  rotulos <- unname(mapa[x])

  # Fallback seguro: se algum ID nao estiver no mapa,
  # conservar o proprio SampleID.
  ausentes <- is.na(rotulos) | rotulos == ""
  rotulos[ausentes] <- x[ausentes]

  rotulos
}


# Auditoria explícita das premissas dos testes por BeeSpecies. 

auditar_premissas_testes <- function(ps_obj, prefixo, output_path) {
  meta_p <- as(phyloseq::sample_data(ps_obj), "data.frame")
  n_sp <- table(droplevels(as.factor(meta_p$BeeSpecies)))
  n_run <- table(droplevels(as.factor(meta_p$Run)))
  mel_info <- resumir_estrutura_meliponario(ps_obj)
  tem_mel_compartilhado <- nrow(mel_info$compartilhados) > 0L
  uma_run <- length(n_run) == 1L
  run_ajustavel <- uma_run

  linhas <- data.frame(
    Conjunto = prefixo,
    Analise = c(
      "Teste exato de postos", "Teste exato de postos",
      "PERMANOVA", "PERMANOVA", "BETADISPER",
      "DESeq2/analises diferenciais posteriores", "Todas"
    ),
    Premissa = c(
      "Grupos com pelo menos duas unidades",
      "Permutabilidade/independencia sob a hipotese nula",
      "Grupos com pelo menos duas unidades",
      "Permutabilidade sem efeito de corrida nao controlado",
      "Pelo menos tres unidades por grupo para diagnostico estavel",
      "Corrida estimavel separadamente de BeeSpecies",
      "Independencia completa em relacao ao meliponario"
    ),
    Atendida_diretamente = c(
      all(n_sp >= 2L),
      uma_run && !tem_mel_compartilhado,
      all(n_sp >= 2L),
      uma_run,
      all(n_sp >= 3L),
      run_ajustavel,
      !tem_mel_compartilhado
    ),
    Status = c(
      if (all(n_sp >= 2L)) "ATENDIDA" else "NAO_ATENDIDA",
      if (uma_run && !tem_mel_compartilhado) "ATENDIDA" else "LIMITADA",
      if (all(n_sp >= 2L)) "ATENDIDA" else "NAO_ATENDIDA",
      if (uma_run) "ATENDIDA" else "NAO_ATENDIDA",
      if (all(n_sp >= 3L)) "ATENDIDA_COM_BAIXO_PODER" else "NAO_ESTIMAVEL",
      if (run_ajustavel) "NAO_NECESSARIA_UMA_RUN" else "NAO_ATENDIDA",
      if (!tem_mel_compartilhado) "ATENDIDA" else "LIMITADA"
    ),
    Consequencia = c(
      paste(names(n_sp), n_sp, sep = "=", collapse = "; "),
      if (uma_run && !tem_mel_compartilhado) {
        "Permutacoes livres sao compatíveis com o desenho observado."
      } else {
        "P-valores permanecem exploratorios; enumeracao exata nao corrige dependencia, meliponario ou batch."
      },
      paste(names(n_sp), n_sp, sep = "=", collapse = "; "),
      if (uma_run) {
        "Nao ha variacao de corrida dentro do conjunto."
      } else {
        "A unica amostra da corrida auxiliar pertence a M. fasciculata; especie e corrida nao sao separaveis."
      },
      if (all(n_sp >= 3L)) {
        "Calculavel, mas n=3 por grupo ainda oferece baixo poder para detectar heterogeneidade."
      } else {
        "O diagnostico e omitido quando algum grupo tem n<3."
      },
      if (run_ajustavel) {
        "Todos os dados pertencem a mesma corrida."
      } else {
        "Modelos plus10 sem Run estimam uma combinacao inseparavel de efeito biologico e tecnico."
      },
      mel_info$resumo
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    linhas,
    file.path(output_path, paste0(prefixo, "_auditoria_premissas_testes.csv")),
    row.names = FALSE
  )
  linhas
}

# Resume a independencia das unidades amostrais sem confundi-la com
# independencia em relacao ao local de coleta.
resumir_estrutura_meliponario <- function(ps_obj) {
  meta_design <- as(phyloseq::sample_data(ps_obj), "data.frame")
  meta_design$SampleID <- rownames(meta_design)

  if (!"Meliponary" %in% colnames(meta_design)) {
    return(list(
      resumo = paste(
        "As unidades amostrais correspondem a amostras de mel distintas.",
        "A coluna Meliponary nao esta disponivel neste objeto; portanto,",
        "a independencia em relacao ao meliponario nao pode ser verificada."
      ),
      compartilhados = data.frame()
    ))
  }

  mel <- trimws(as.character(meta_design$Meliponary))
  validos <- !is.na(mel) & nzchar(mel)
  if (!any(validos)) {
    return(list(
      resumo = paste(
        "As unidades amostrais correspondem a amostras de mel distintas.",
        "Meliponary esta vazio; a independencia em relacao ao meliponario",
        "nao pode ser verificada."
      ),
      compartilhados = data.frame()
    ))
  }

  meta_validos <- meta_design[validos,
    c("SampleID", "BeeSpecies", "Meliponary"), drop = FALSE]
  meta_validos$Meliponary <- trimws(as.character(meta_validos$Meliponary))
  grupos <- split(meta_validos, meta_validos$Meliponary, drop = TRUE)

  linhas <- lapply(names(grupos), function(nome_meliponario) {
    grupo <- grupos[[nome_meliponario]]
    especies <- sort(unique(as.character(grupo$BeeSpecies)))
    if (nrow(grupo) < 2L || length(especies) < 2L) return(NULL)
    data.frame(
      Meliponary = nome_meliponario,
      N_amostras = nrow(grupo),
      BeeSpecies = paste(especies, collapse = "; "),
      SampleID = paste(grupo$SampleID, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  linhas <- linhas[!vapply(linhas, is.null, logical(1))]
  compartilhados <- if (length(linhas)) {
    do.call(rbind, linhas)
  } else {
    data.frame()
  }

  if (nrow(compartilhados) == 0L) {
    resumo <- paste(
      
    )
  } else {
    detalhes <- paste(
      paste0(compartilhados$Meliponary, " [", compartilhados$BeeSpecies, "]"),
      collapse = "; "
    )
    resumo <- paste0(
      "As unidades amostrais correspondem a amostras de mel distintas; ",
      "contudo, a independencia nao e completa em relacao ao meliponario. ",
      "Os seguintes meliponarios contribuem com amostras de mais de uma ",
      "BeeSpecies: ", detalhes, ". Amostras do mesmo meliponario podem ",
      "compartilhar efeitos ambientais e de manejo; por isso, os testes por ",
      "BeeSpecies devem ser interpretados como exploratorios."
    )
  }

  list(resumo = resumo, compartilhados = compartilhados)
}

analisar_funil_invertido <- function(ps_obj,
                                     prefixo,
                                     dir_fig,
                                     executar_testes     = TRUE,
                                     incluir_covariaveis = TRUE) {

  cat("=============================================================\n")
  cat("ANALISE:", prefixo, "\n")
  cat("=============================================================\n\n")

  ps_obj <- validar_phyloseq_basico(ps_obj, prefixo)
  ps_rel <- transformar_relativo(ps_obj)
  mapa_samplelabel <- obter_mapa_samplelabel(ps_obj)
  n_por_especie <- table(sample_data(ps_obj)$BeeSpecies)
  estrutura_meliponario <- resumir_estrutura_meliponario(ps_obj)

  cat("Numero de unidades amostrais distintas por BeeSpecies:\n")
  print(n_por_especie)
  cat("\n")

  if (nrow(estrutura_meliponario$compartilhados) > 0L) {
    write.csv(
      estrutura_meliponario$compartilhados,
      file.path(
        output_path,
        paste0(prefixo, "_meliponarios_compartilhados_entre_especies.csv")
      ),
      row.names = FALSE
    )
  }

  if (executar_testes) {
    if (identical(prefixo, "plus10")) {
      log_msg(
        paste0(
          "plus10 inclui uma amostra da corrida run_aux. Os testes sao ",
          "sensibilidades exploratorias separadas; o efeito de Run nao pode ser ",
          "estimado porque a segunda corrida possui somente uma amostra e esta ",
          "amostra pertence a Melipona fasciculata."
        ),
        "WARN"
      )
    }
    log_msg(estrutura_meliponario$resumo, "WARN")
    log_msg(
      paste0(
        "A interpretacao inferencial tambem e limitada pelo baixo n por ",
        "BeeSpecies: ",
        paste(names(n_por_especie), n_por_especie, sep = "=", collapse = "; "),
        "."
      ),
      "WARN"
    )
  }

  if (executar_testes && any(n_por_especie < 2)) {
    stop(
      "Ha BeeSpecies com menos de 2 amostras. ",
      "Comparacoes por grupo, PERMANOVA e BETADISPER nao devem ser executadas."
    )
  }
  saveRDS(ps_rel,
          file.path(output_path, paste0(prefixo, "_phyloseq_relativo.rds")))

  # -------------------------------------------------------------------------
  # 9.1 NIVEL 1 — AMOSTRAS INDIVIDUAIS
  # -------------------------------------------------------------------------

  cat("--- Nivel 1: amostras individuais ---\n\n")

  qc_amostras <- data.frame(
    SampleID = sample_names(ps_obj),

    SampleLabel = aplicar_samplelabel(
      sample_names(ps_obj),
      mapa_samplelabel
    ),

    Reads = sample_sums(ps_obj),

    BeeSpecies = as.character(
      sample_data(ps_obj)$BeeSpecies
    ),
    Nativo_Introduzido = if ("Nativo_Introduzido" %in% sample_variables(ps_obj))
      as.character(sample_data(ps_obj)$Nativo_Introduzido) else NA_character_,
    Run        = as.character(sample_data(ps_obj)$Run),
    stringsAsFactors = FALSE)
  write.csv(qc_amostras,
            file.path(output_path, paste0(prefixo, "_nivel1_amostras_qc.csv")),
            row.names = FALSE)

  otu_taxa <- matriz_taxa_linhas(ps_obj)
  qc_asvs <- data.frame(
    ASV_ID     = taxa_names(ps_obj),
    TotalReads = taxa_sums(ps_obj),
    Prevalence = rowSums(otu_taxa > 0),
    stringsAsFactors = FALSE)
  write.csv(qc_asvs,
            file.path(output_path, paste0(prefixo, "_nivel1_asvs_prevalencia.csv")),
            row.names = FALSE)

  p_depth <- ggplot(
    qc_amostras,
    aes(
      x = reorder(SampleID, Reads),
      y = Reads,
      fill = BeeSpecies
    )
  ) +
    geom_col(width = 0.75) +
    coord_flip() +

    scale_x_discrete(
      labels = function(x) {
        aplicar_samplelabel(x, mapa_samplelabel)
      }
    ) +

    scale_fill_brewer(palette = "Set2") +
    theme_classic(base_size = 11) +
    labs(title = paste0(prefixo, " — profundidade por amostra"),x = NULL,
      y = "Reads", fill = "Especie de abelha"
    )
  salvar_plot(p_depth, paste0(prefixo, "_nivel1_sample_depth.pdf"),
              dir_fig, width = 9, height = 6)
  print(qc_amostras); cat("\n")

  # -------------------------------------------------------------------------
  # 9.2 DIVERSIDADE ALFA
  # -------------------------------------------------------------------------

  cat("--- Diversidade alfa ---\n\n")

  medidas_alfa <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
  # Chao1 e ACE sao exportados apenas como diagnostico suplementar. Como o
  # Script 01 aplica filtro de frequencia, singletons/doubletons podem ter sido
  # removidos e esses estimadores nao devem sustentar a inferencia principal.
  medidas_alfa_teste <- c("Observed", "Shannon", "Simpson", "InvSimpson")

  alfa <- estimate_richness(ps_obj, measures = medidas_alfa)
  alfa$SampleID <- rownames(alfa)

  # Unir metadados ao data.frame de alfa para acesso a BeeSpecies e Run
  meta_df <- as(sample_data(ps_obj), "data.frame")
  meta_df$SampleID <- rownames(meta_df)
  alfa <- dplyr::left_join(alfa, meta_df, by = "SampleID")

  write.csv(alfa,
            file.path(output_path, paste0(prefixo, "_nivel2_alpha_por_amostra.csv")),
            row.names = FALSE)

  alfa_por_sp <- alfa |>
    dplyr::group_by(BeeSpecies) |>
    dplyr::summarise(
      n           = dplyr::n(),
      Observed_m  = round(mean(Observed,  na.rm = TRUE), 1),
      Observed_sd = round(sd(Observed,    na.rm = TRUE), 1),
      Chao1_m     = round(mean(Chao1,     na.rm = TRUE), 1),
      Chao1_sd    = round(sd(Chao1,       na.rm = TRUE), 1),
      Shannon_m   = round(mean(Shannon,   na.rm = TRUE), 3),
      Shannon_sd  = round(sd(Shannon,     na.rm = TRUE), 3),
      Simpson_m   = round(mean(Simpson,   na.rm = TRUE), 3),
      Simpson_sd  = round(sd(Simpson,     na.rm = TRUE), 3),
      .groups = "drop")
  write.csv(alfa_por_sp,
            file.path(output_path, paste0(prefixo, "_nivel2_alpha_por_BeeSpecies.csv")),
            row.names = FALSE)
  print(alfa_por_sp); cat("\n")

  p_alfa <- plot_richness(ps_obj, x = "BeeSpecies", color = "BeeSpecies",
                          measures = c("Observed","Chao1","ACE","Shannon","Simpson")) +
    geom_boxplot(aes(fill = BeeSpecies), alpha = 0.35, width = 0.5,
                 outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.5, alpha = 0.85) +
    scale_fill_brewer(palette  = "Set2") +
    scale_color_brewer(palette = "Set2") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(title    = paste0(prefixo, " — diversidade alfa"),
         subtitle = "Contagens nao rarefadas; Chao1/ACE apenas suplementares apos filtro de frequencia",
         x = NULL, y = "Estimativa")
  salvar_plot(p_alfa, paste0(prefixo, "_nivel2_alpha.pdf"),
              dir_fig, width = 12, height = 7)

  # -------------------------------------------------------------------------
  # 9.3 TESTES ALFA: KRUSKAL-WALLIS (apenas se executar_testes)
  # AVISO ESTATISTICO:
  #   -
  # -------------------------------------------------------------------------

  if (executar_testes) {
    cat("--- Testes alfa: teste de postos exato por BeeSpecies ---\n\n")
    log_msg(
      "AVISO ESTATISTICO: unidades amostrais distintas, com possivel compartilhamento de ambiente e manejo dentro de meliponario; baixo n por BeeSpecies; interpretar o teste exato como exploratorio; a enumeracao nao corrige dependencia ou confundimento.", "WARN")

    kw_alfa <- function(medida) {
      z <- tryCatch(
        teste_postos_exato_grupos_fixos(
          alfa[[medida]], alfa$BeeSpecies,
          tamanhos_esperados = if (identical(prefixo, "core9")) c(2L, 3L, 4L) else c(3L, 3L, 4L),
          max_alocacoes = 1000000L
        ),
        error = function(e) e
      )
      if (inherits(z, "error")) {
        log_msg(paste("Teste exato falhou para", medida, ":", conditionMessage(z)), "WARN")
        return(data.frame(
          Medida = medida, Statistic = NA_real_, P_value = NA_real_,
          P_value_assintotico_diagnostico = NA_real_, df = NA_integer_,
          N_permutacoes_exatas = NA_real_, p_min_teorico = NA_real_,
          Epsilon2 = NA_real_, Metodo = "Teste exato por enumeracao",
          Observacao = paste("Falhou:", conditionMessage(z)),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        Medida = medida,
        Statistic = z$H,
        P_value = z$p_exato,
        P_value_assintotico_diagnostico = z$p_assintotico_diagnostico,
        df = z$df,
        N_permutacoes_exatas = z$N_permutacoes_exatas,
        p_min_teorico = z$p_min_teorico,
        Epsilon2 = z$Epsilon2,
        Metodo = z$Metodo,
        Observacao = paste0(
          "Exploratorio; ", z$Nota, " n por BeeSpecies: ",
          paste(names(n_por_especie), n_por_especie, sep = "=", collapse = "; ")
        ),
        stringsAsFactors = FALSE
      )
    }

    testes_alfa <- dplyr::bind_rows(lapply(medidas_alfa_teste, kw_alfa))
    testes_alfa$P_value_ajustado_BH <- p.adjust(testes_alfa$P_value, method = "BH")
    testes_alfa$Significativo_BH <- !is.na(testes_alfa$P_value_ajustado_BH) &
      testes_alfa$P_value_ajustado_BH < 0.05
    write.csv(testes_alfa,
              file.path(output_path, paste0(prefixo, "_teste_alpha_kruskal.csv")),
              row.names = FALSE)
    print(testes_alfa); cat("\n")

  } else {
    log_msg("Testes alfa ignorados (analise de sensibilidade).", "INFO")
  }

  # -------------------------------------------------------------------------
  # 9.4 COMPOSICAO POR BeeSpecies — BARPLOTS
  #
  # -------------------------------------------------------------------------

  cat("--- Composicao por BeeSpecies: barplots ---\n\n")

  gerar_barplot_tax <- function(nivel, top_n, arquivo) {
    if (!nivel %in% rank_names(ps_obj)) {
      log_msg(sprintf("Nivel '%s' ausente. Barplot ignorado.", nivel), "WARN")
      return(invisible(NULL))
    }
    ps_agr     <- tax_glom(ps_obj, taxrank = nivel, NArm = FALSE)
    ps_agr_rel <- transformar_relativo(ps_agr)
    melt_df    <- psmelt(ps_agr_rel)

    # NA substituido por "Unclassified" somente na camada de visualizacao
    melt_df$TaxonPlot <- as.character(melt_df[[nivel]])
    melt_df$TaxonPlot[is.na(melt_df$TaxonPlot) | melt_df$TaxonPlot == ""] <-
      "Unclassified"

    top_taxa <- melt_df |>
      dplyr::group_by(TaxonPlot) |>
      dplyr::summarise(Media = mean(Abundance), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(Media)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::pull(TaxonPlot)
    melt_df$TaxonPlot[!melt_df$TaxonPlot %in% top_taxa] <- "Outros"

    p <- ggplot(
      melt_df,
      aes(x = Sample, y = Abundance, fill = TaxonPlot)
    ) +
      geom_bar(stat = "identity", width = 0.9) +

      facet_wrap(
        ~ BeeSpecies,
        scales = "free_x"
      ) +

      scale_x_discrete(
        labels = function(x) {
          aplicar_samplelabel(x, mapa_samplelabel)
        }
      ) +

      scale_fill_viridis_d(option = "turbo") +
      theme_classic(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = grid::unit(0.35, "cm")
      ) +
      labs(
        title = paste0(prefixo, " — composicao relativa: ", nivel),
        subtitle = paste("Top", top_n, "taxa; demais agrupados como Outros"),
        x = NULL,
        y = "Abundancia relativa",
        fill = nivel
      )
    salvar_plot(p, arquivo, dir_fig, width = 13, height = 6)
    invisible(p)
  }

  gerar_barplot_tax("Phylum", TOP_N_PHYLUM,
                    paste0(prefixo, "_nivel2_barplot_phylum.pdf"))
  gerar_barplot_tax("Genus",  TOP_N_GENUS,
                    paste0(prefixo, "_nivel2_barplot_genus.pdf"))

  # -------------------------------------------------------------------------
  # 9.5 BETA-DIVERSIDADE
  # Bray-Curtis sobre abundancias relativas (ps_rel).
  # Jaccard binario sobre presenca/ausencia (ps_obj, binary = TRUE).
  # -------------------------------------------------------------------------

  cat("--- Beta-diversidade ---\n\n")

  dist_bray    <- phyloseq::distance(ps_rel, method = "bray")
  dist_jaccard <- phyloseq::distance(ps_obj, method = "jaccard", binary = TRUE)

  exportar_distancia(dist_bray,    "bray_rel",       prefixo)
  exportar_distancia(dist_jaccard, "jaccard_binary",  prefixo)

  gerar_ordination <- function(dist_obj, metodo, nome_dist, arquivo) {
    set.seed(SEED_GLOBAL)
    ord <- tryCatch(
      {
        if (metodo == "NMDS") {
          ordinate(
            ps_obj,
            method = metodo,
            distance = dist_obj,
            trymax = 200,
            trace = FALSE
          )
        } else {
          ordinate(
            ps_obj,
            method = metodo,
            distance = dist_obj
          )
        }
      },
      error = function(e) {
        log_msg(
          paste("Falha:", metodo, nome_dist, conditionMessage(e)),
          "WARN"
        )
        NULL
      }
    )
    if (is.null(ord)) return(invisible(NULL))

    saveRDS(ord, file.path(output_path,
                           paste0(prefixo, "_ordination_", metodo, "_", nome_dist, ".rds")))

    meta_plot <- as(sample_data(ps_obj), "data.frame")
    usar_shape_run <- "Run" %in% colnames(meta_plot) &&
      length(unique(meta_plot$Run)) > 1

    p <- if (usar_shape_run) {
      plot_ordination(
        ps_obj,
        ord,
        color = "BeeSpecies",
        shape = "Run"
      )
    } else {
      plot_ordination(
        ps_obj,
        ord,
        color = "BeeSpecies"
      )
    }

    p <- p +
      geom_point(size = 3.5, alpha = 0.85) +
      scale_color_brewer(palette = "Set2") +
      theme_classic(base_size = 11) +
      labs(
        title = paste0(
          prefixo,
          " — ",
          metodo,
          " / ",
          nome_dist
        ),
        color = "Especie de abelha",
        shape = "Corrida"
      )

    if ("SampleLabel" %in% colnames(p$data)) {
      p <- p +
        ggrepel::geom_text_repel(
          aes(label = SampleLabel),
          size = 3,
          show.legend = FALSE,
          max.overlaps = Inf,
          box.padding = 0.35,
          point.padding = 0.25,
          min.segment.length = 0
        )
    }
    # Elipse de confianca (t-distribution, 95%): apenas quando todos os grupos
    # tiverem n >= 4. Com n = 3, a elipse nao e estimavel de forma confiavel.
    n_grupo <- table(as(sample_data(ps_obj), "data.frame")$BeeSpecies)
    if (all(n_grupo >= 4)) {
      p <- p + stat_ellipse(aes(group = BeeSpecies), type = "t",
                            linetype = 2, level = 0.95, alpha = 0.5,
                            show.legend = FALSE)
    } else {
      log_msg(paste0(prefixo, ": elipses omitidas — grupos com n < 4."), "WARN")
    }
    salvar_plot(p, arquivo, dir_fig, width = 8, height = 5)

    # Scree plot somente para PCoA (variancia explicada por eixo)
    if (metodo == "PCoA") {
      p_scree <- tryCatch(
        plot_scree(ord, paste("Scree —", prefixo, nome_dist)) +
          theme_classic(base_size = 11),
        error = function(e) {
          log_msg(
            paste0(prefixo, "/", nome_dist,
                   ": scree plot omitido por erro: ", conditionMessage(e)),
            "WARN"
          )
          NULL
        })
      if (!is.null(p_scree))
        salvar_plot(p_scree, sub("\\.pdf$", "_scree.pdf", arquivo),
                    dir_fig, width = 8, height = 5)
    }
    invisible(ord)
  }

  gerar_ordination(dist_bray,    "PCoA", "bray_rel",
                   paste0(prefixo, "_nivel2_pcoa_bray.pdf"))
  gerar_ordination(dist_bray,    "NMDS", "bray_rel",
                   paste0(prefixo, "_nivel2_nmds_bray.pdf"))
  gerar_ordination(dist_jaccard, "PCoA", "jaccard_binary",
                   paste0(prefixo, "_nivel2_pcoa_jaccard.pdf"))
  gerar_ordination(dist_jaccard, "NMDS", "jaccard_binary",
                   paste0(prefixo, "_nivel2_nmds_jaccard.pdf"))

  # -------------------------------------------------------------------------
  # 9.6 PERMANOVA E BETADISPER (apenas se executar_testes)
  #
  # PERMANOVA:
  #   testa diferencas entre centroides das comunidades por BeeSpecies.
  #
  # BETADISPER:
  #   testa se os grupos possuem dispersoes semelhantes.
  #   Deve ser interpretado junto com a PERMANOVA.
  #
  # Interpretacao:
  #   PERMANOVA significativa + BETADISPER nao significativo:
  #     evidencia exploratoria de diferenca de composicao entre BeeSpecies.
  #
  #   PERMANOVA significativa + BETADISPER significativo:
  #     resultado pode refletir diferenca de dispersao entre grupos, nao apenas
  #     diferenca entre centroides.
  #
  #   PERMANOVA nao significativa:
  #     ausencia de evidencia estatistica, mas nao ausencia de diferenca biologica,
  #     pois o n amostral e pequeno.
  # -------------------------------------------------------------------------

  if (executar_testes) {

    cat("--- PERMANOVA e BETADISPER ---\n\n")

    log_msg(
      paste0(
        "AVISO ESTATISTICO: unidades amostrais distintas, mas a independencia nao e completa em relacao ao meliponario. ",
        if (identical(prefixo, "core9")) {
          "O core9 pertence a uma unica corrida. "
        } else {
          "O plus10 inclui uma unica amostra da corrida auxiliar; Run e BeeSpecies nao podem ser separados. "
        },
        "Baixo n por BeeSpecies; PERMANOVA e BETADISPER devem ser interpretados como exploratorios. ",
        "n por grupo: ",
        paste(names(n_por_especie), n_por_especie, sep = "=", collapse = "; ")
      ),
      "WARN"
    )

    rodar_permanova_betadisper <- function(dist_obj,
                                           nome_dist,
                                           grupo = "BeeSpecies",
                                           n_perm = N_PERM) {

      if (!inherits(dist_obj, "dist")) {
        stop(nome_dist, ": dist_obj precisa ser da classe 'dist'.")
      }

      meta_ord <- as(sample_data(ps_obj), "data.frame")

      if (!grupo %in% colnames(meta_ord)) {
        stop("Variavel de grupo ausente em sample_data: ", grupo)
      }

      # Garantir que a ordem dos metadados seja exatamente a ordem da matriz de distancia.
      meta_ord <- meta_ord[labels(dist_obj), , drop = FALSE]

      if (!identical(rownames(meta_ord), labels(dist_obj))) {
        stop("Ordem entre matriz de distancia e metadados diverge em: ", nome_dist)
      }

      grupo_vec <- droplevels(as.factor(meta_ord[[grupo]]))
      n_grupo_local <- table(grupo_vec)

      if (length(n_grupo_local) < 2) {
        stop(nome_dist, ": menos de 2 grupos em ", grupo, ".")
      }

      if (any(n_grupo_local < 2)) {
        stop(
          nome_dist,
          ": ha grupo com menos de 2 amostras. ",
          "PERMANOVA/BETADISPER nao devem ser executados."
        )
      }

      # ---------------------------------------------------------------------
      # PERMANOVA
      # ---------------------------------------------------------------------

      set.seed(SEED_GLOBAL)

      formula_perm <- as.formula(paste("dist_obj ~", grupo))

      perm <- vegan::adonis2(
        formula_perm,
        data = meta_ord,
        permutations = n_perm,
        by = "margin"
      )

      perm_df <- as.data.frame(perm)
      perm_df$Termo <- rownames(perm_df)
      rownames(perm_df) <- NULL
      perm_df$Distancia <- nome_dist
      perm_df$Grupo <- grupo
      perm_df$Permutacoes <- n_perm
      perm_df$Interpretacao <- "PERMANOVA exploratoria; interpretar junto com BETADISPER"

      write.csv(
        perm_df,
        file.path(output_path, paste0(prefixo, "_permanova_", nome_dist, ".csv")),
        row.names = FALSE
      )

      saveRDS(
        perm,
        file.path(output_path, paste0(prefixo, "_permanova_", nome_dist, ".rds"))
      )

      cat("\nPERMANOVA —", nome_dist, "\n")
      print(perm)

      # ---------------------------------------------------------------------
      # BETADISPER
      # GUARDA: betadisper estima a dispersao de cada grupo em torno do
      # centroide. Com n < 3 a estimativa tem 1 g.l. (instavel) e
      # bias.adjust = TRUE pode produzir valores degenerados ou erro.
      # Quando algum grupo tem n < 3, o teste e omitido e registrado como
      # nao realizavel por desenho amostral.
      # ---------------------------------------------------------------------

      betadisper_realizavel <- all(n_grupo_local >= 3)

      if (!betadisper_realizavel) {
        log_msg(
          sprintf("%s: BETADISPER omitido — grupo(s) com n < 3 (%s). Nao realizavel por desenho.",
                  nome_dist,
                  paste(names(n_grupo_local), n_grupo_local, sep = "=", collapse = "; ")),
          "WARN"
        )
        bd          <- NULL
        bd_anova    <- NULL
        bd_perm     <- NULL
        bd_anova_df <- data.frame(
          Termo = "Groups", Distancia = nome_dist, Grupo = grupo,
          Metodo = "BETADISPER nao realizavel (grupo com n < 3)",
          Bias_adjust = NA, `Pr(>F)` = NA_real_,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        write.csv(
          bd_anova_df,
          file.path(output_path, paste0(prefixo, "_betadisper_anova_", nome_dist, ".csv")),
          row.names = FALSE
        )
      } else {

      bd <- vegan::betadisper(
        dist_obj,
        group = grupo_vec,
        type = "centroid",
        bias.adjust = TRUE
      )

      bd_anova <- anova(bd)

      set.seed(SEED_GLOBAL)

      bd_perm <- vegan::permutest(
        bd,
        permutations = n_perm,
        pairwise = TRUE
      )

      saveRDS(
        bd,
        file.path(output_path, paste0(prefixo, "_betadisper_", nome_dist, ".rds"))
      )

      saveRDS(
        bd_perm,
        file.path(output_path, paste0(prefixo, "_betadisper_permutest_", nome_dist, ".rds"))
      )

      bd_anova_df <- as.data.frame(bd_anova)
      bd_anova_df$Termo <- rownames(bd_anova_df)
      rownames(bd_anova_df) <- NULL
      bd_anova_df$Distancia <- nome_dist
      bd_anova_df$Grupo <- grupo
      bd_anova_df$Metodo <- "ANOVA sobre distancias ao centroide"
      bd_anova_df$Bias_adjust <- TRUE

      write.csv(
        bd_anova_df,
        file.path(output_path, paste0(prefixo, "_betadisper_anova_", nome_dist, ".csv")),
        row.names = FALSE
      )

      capture.output(
        bd_perm,
        file = file.path(output_path, paste0(prefixo, "_betadisper_permutest_", nome_dist, ".txt"))
      )

      cat("\nBETADISPER — ANOVA —", nome_dist, "\n")
      print(bd_anova)

      cat("\nBETADISPER — PERMUTEST —", nome_dist, "\n")
      print(bd_perm)

      # ---------------------------------------------------------------------
      # Distancia de cada amostra ao centroide do grupo
      # ---------------------------------------------------------------------

      bd_distancias <- data.frame(
        SampleID = names(bd$distances),
        Grupo = as.character(bd$group),
        Distancia_ao_centroide = as.numeric(bd$distances),
        Distancia = nome_dist,
        stringsAsFactors = FALSE
      )

      meta_export <- meta_ord[bd_distancias$SampleID, , drop = FALSE]
      meta_export$SampleID <- rownames(meta_export)

      bd_distancias <- dplyr::left_join(
        bd_distancias,
        meta_export,
        by = "SampleID"
      )

      write.csv(
        bd_distancias,
        file.path(output_path, paste0(prefixo, "_betadisper_distancias_amostras_", nome_dist, ".csv")),
        row.names = FALSE
      )

      bd_resumo_grupo <- bd_distancias |>
        dplyr::group_by(Grupo) |>
        dplyr::summarise(
          n = dplyr::n(),
          Media_distancia = mean(Distancia_ao_centroide, na.rm = TRUE),
          DP_distancia = sd(Distancia_ao_centroide, na.rm = TRUE),
          Mediana_distancia = median(Distancia_ao_centroide, na.rm = TRUE),
          Min_distancia = min(Distancia_ao_centroide, na.rm = TRUE),
          Max_distancia = max(Distancia_ao_centroide, na.rm = TRUE),
          .groups = "drop"
        )

      write.csv(
        bd_resumo_grupo,
        file.path(output_path, paste0(prefixo, "_betadisper_resumo_grupos_", nome_dist, ".csv")),
        row.names = FALSE
      )

      # ---------------------------------------------------------------------
      # Grafico diagnostico do BETADISPER
      # ---------------------------------------------------------------------

      plot_bd <- ggplot(bd_distancias, aes(x = Grupo, y = Distancia_ao_centroide, fill = Grupo)) +
        geom_boxplot(alpha = 0.35, width = 0.55, outlier.shape = NA) +
        geom_jitter(width = 0.10, size = 2.8, alpha = 0.85) +
        scale_fill_brewer(palette = "Set2") +
        theme_classic(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none"
        ) +
        labs(
          title = paste0(prefixo, " — BETADISPER / ", nome_dist),
          subtitle = "Distancia de cada amostra ao centroide do grupo",
          x = NULL,
          y = "Distancia ao centroide"
        )

      salvar_plot(
        plot_bd,
        paste0(prefixo, "_betadisper_distancias_", nome_dist, ".pdf"),
        dir_fig,
        width = 8,
        height = 5
      )

      }  # fim do else (betadisper realizavel)

      # ---------------------------------------------------------------------
      # Tabela interpretativa integrada
      # ---------------------------------------------------------------------

      p_perm <- NA_real_
      r2_perm <- NA_real_

      if (grupo %in% perm_df$Termo) {
        p_perm <- perm_df$`Pr(>F)`[perm_df$Termo == grupo]
        r2_perm <- perm_df$R2[perm_df$Termo == grupo]
      }

      p_betadisper <- NA_real_

      if ("Pr(>F)" %in% colnames(bd_anova_df)) {
        p_betadisper <- bd_anova_df$`Pr(>F)`[bd_anova_df$Termo == "Groups"]
      }

      interpretacao_integrada <- dplyr::case_when(
        is.na(p_perm) | is.na(p_betadisper) ~
          "Nao foi possivel gerar interpretacao automatica.",
        p_perm < 0.05 & p_betadisper >= 0.05 ~
          "PERMANOVA significativa e BETADISPER nao significativo: evidencia exploratoria de diferenca entre centroides das comunidades.",
        p_perm < 0.05 & p_betadisper < 0.05 ~
          "PERMANOVA significativa e BETADISPER significativo: diferenca pode refletir heterogeneidade de dispersao entre grupos.",
        p_perm >= 0.05 & p_betadisper < 0.05 ~
          "PERMANOVA nao significativa, mas BETADISPER significativo: grupos diferem em dispersao, sem evidencia de diferenca entre centroides.",
        p_perm >= 0.05 & p_betadisper >= 0.05 ~
          "PERMANOVA e BETADISPER nao significativos: sem evidencia estatistica de diferenca entre centroides ou dispersoes; interpretar considerando baixo n.",
        TRUE ~
          "Interpretacao nao classificada."
      )

      interpretacao_df <- data.frame(
        Prefixo = prefixo,
        Distancia = nome_dist,
        Grupo = grupo,
        N_amostras = nsamples(ps_obj),
        N_por_grupo = paste(names(n_grupo_local), n_grupo_local, sep = "=", collapse = "; "),
        PERMANOVA_R2 = r2_perm,
        PERMANOVA_p = p_perm,
        BETADISPER_p = p_betadisper,
        Permutacoes = n_perm,
        Inclui_segunda_run = !identical(prefixo, "core9"),
        Run_separavel_de_BeeSpecies = identical(prefixo, "core9"),
        Permutabilidade_plena = identical(prefixo, "core9") &&
          nrow(resumir_estrutura_meliponario(ps_obj)$compartilhados) == 0L,
        Status_inferencia = if (identical(prefixo, "core9")) {
          "principal_exploratoria_por_baixo_n_e_meliponario"
        } else {
          "sensibilidade_exploratoria_com_batch_confundido"
        },
        Interpretacao = if (identical(prefixo, "core9")) {
          interpretacao_integrada
        } else {
          paste0(
            interpretacao_integrada,
            " Plus10 inclui uma unica amostra de outra corrida; o resultado nao separa efeito de especie de efeito tecnico."
          )
        },
        stringsAsFactors = FALSE
      )

      write.csv(
        interpretacao_df,
        file.path(output_path, paste0(prefixo, "_interpretacao_permanova_betadisper_", nome_dist, ".csv")),
        row.names = FALSE
      )

      invisible(list(
        permanova = perm,
        betadisper = bd,
        betadisper_anova = bd_anova,
        betadisper_permutest = bd_perm,
        distancias_centroide = if (betadisper_realizavel) bd_distancias else NULL,
        interpretacao = interpretacao_df
      ))
    }

    resultado_bray <- tryCatch(
      rodar_permanova_betadisper(dist_bray, "bray_rel"),
      error = function(e) stop("PERMANOVA/BETADISPER bray_rel falhou: ",
                               conditionMessage(e), call. = FALSE)
    )
    resultado_jaccard <- tryCatch(
      rodar_permanova_betadisper(dist_jaccard, "jaccard_binary"),
      error = function(e) stop("PERMANOVA/BETADISPER jaccard_binary falhou: ",
                               conditionMessage(e), call. = FALSE)
    )

    cat("\n")

  } else {
    log_msg("PERMANOVA/BETADISPER nao executados para este conjunto por configuracao.", "INFO")
  }

  # -------------------------------------------------------------------------
  # 9.7 NIVEL 3 — COVARIAVEIS
  # Avalia cada covariavel candidata quanto a: completude, variacao,
  # confundimento com BeeSpecies e elegibilidade para modelo inferencial.
  # A inclusao em PERMANOVA nao e automatica para evitar superajuste
  # dado o n pequeno.
  # -------------------------------------------------------------------------

  cat("--- Nivel 3: covariaveis ---\n\n")

  COVARIAVEIS_CANDIDATAS <- c("Coexistence_Level", "Environment",
                              "Municipality", "Meliponary")
  meta_cov <- as(sample_data(ps_obj), "data.frame")

  # avaliar_covariavel: retorna 1 linha do resumo para uma covariavel
  avaliar_covariavel <- function(cv) {

    if (!cv %in% colnames(meta_cov)) {
      return(data.frame(
        Covariavel = cv,
        N_na = NA_integer_,
        N_niveis = NA_integer_,
        Usavel_descritivo = FALSE,
        Usavel_modelo = FALSE,
        Observacao = "Ausente nos metadados",
        stringsAsFactors = FALSE
      ))
    }

    x <- meta_cov[[cv]]

    n_na <- sum(is.na(x))
    n_niveis <- length(unique(x[!is.na(x)]))

    quase_id_amostra <- n_niveis >= (0.7 * nsamples(ps_obj))

    confundida <- FALSE

    if (n_niveis >= 2) {
      tab <- table(meta_cov$BeeSpecies, x, useNA = "no")

      cov_determinada_por_especie <- all(rowSums(tab > 0) <= 1L)
      especie_determinada_por_cov <- all(colSums(tab > 0) <= 1L)
      confundida <- cov_determinada_por_especie || especie_determinada_por_cov
    }

    usavel_descritivo <- n_niveis >= 2

    usavel_modelo <- incluir_covariaveis &&
      executar_testes &&
      n_na == 0 &&
      n_niveis >= 2 &&
      !confundida &&
      !quase_id_amostra &&
      nsamples(ps_obj) >= 9

    obs <- dplyr::case_when(
      !usavel_descritivo ~
        "Sem variacao suficiente",

      n_na > 0 ~
        "Possui NA; somente uso descritivo",

      quase_id_amostra ~
        "Muitos niveis para o n amostral; usar apenas descritivamente",

      confundida ~
        "Confundida com BeeSpecies; excluir de modelo inferencial",

      !usavel_modelo ~
        "Uso descritivo",

      TRUE ~
        "Pode ser testada com cautela"
    )

    capture.output(
      table(meta_cov$BeeSpecies, x, useNA = "ifany"),
      file = file.path(
        output_path,
        paste0(prefixo, "_covariavel_", cv, ".txt")
      )
    )

    data.frame(
      Covariavel = cv,
      N_na = n_na,
      N_niveis = n_niveis,
      Usavel_descritivo = usavel_descritivo,
      Usavel_modelo = usavel_modelo,
      Observacao = obs,
      stringsAsFactors = FALSE
    )
  }
  resumo_cov <- dplyr::bind_rows(lapply(COVARIAVEIS_CANDIDATAS, avaliar_covariavel))
  write.csv(resumo_cov,
            file.path(output_path, paste0(prefixo, "_nivel3_covariaveis.csv")),
            row.names = FALSE)
  print(resumo_cov)

  if (executar_testes && any(resumo_cov$Usavel_modelo))
    log_msg(
      "Covariaveis testaveis detectadas. Adicionar manualmente a PERMANOVA se justificado.",
      "WARN")
  cat("\n")

  # -------------------------------------------------------------------------
  # 9.8 HEATMAP EXPLORATORIO (nivel Genus)
  # -------------------------------------------------------------------------

  cat("--- Heatmap exploratorio ---\n\n")

  gerar_heatmap <- function(nivel, top_n, arquivo) {
    if (!nivel %in% rank_names(ps_obj)) {
      log_msg(sprintf("Nivel '%s' ausente. Heatmap ignorado.", nivel), "WARN")
      return(invisible(NULL))
    }
    ps_agr     <- tax_glom(ps_obj, taxrank = nivel, NArm = FALSE)
    ps_agr_rel <- transformar_relativo(ps_agr)
    top_taxa   <- names(sort(taxa_sums(ps_agr_rel), decreasing = TRUE))[
      seq_len(min(top_n, ntaxa(ps_agr_rel)))]
    ps_top <- prune_taxa(top_taxa, ps_agr_rel)

    p <- tryCatch(
      plot_heatmap(ps_top, method = "NMDS", distance = "bray",
                   sample.label = "BeeSpecies", taxa.label = nivel,
                   title = paste(prefixo, "— heatmap top", top_n, nivel)) +
        theme_classic(base_size = 9),
      error = function(e) {
        log_msg(paste("Heatmap falhou:", conditionMessage(e)), "WARN"); NULL })
    if (!is.null(p))
      salvar_plot(p, arquivo, dir_fig, width = 10, height = 8)
    invisible(p)
  }

  gerar_heatmap("Genus", TOP_N_HEAT,
                paste0(prefixo, "_nivel2_heatmap_genus.pdf"))

  # -------------------------------------------------------------------------
  # 9.9 CLUSTERING HIERARQUICO (metodo UPGMA / average)
  # -------------------------------------------------------------------------

  cat("--- Clustering hierarquico ---\n\n")

  gerar_hclust <- function(dist_obj, nome_dist) {
    hc <- hclust(dist_obj, method = "average")

    # Preservar SampleID no objeto analitico salvo. A substituicao por
    # SampleLabel ocorre somente em uma copia destinada a visualizacao.
    saveRDS(hc, file.path(output_path,
                          paste0(prefixo, "_hclust_", nome_dist, ".rds")))

    hc_plot <- hc
    hc_plot$labels <- aplicar_samplelabel(hc$labels, mapa_samplelabel)

    if (any(is.na(hc_plot$labels)) || any(hc_plot$labels == "")) {
      stop(nome_dist, ": SampleLabel ausente no dendrograma.")
    }
    if (anyDuplicated(hc_plot$labels) > 0L) {
      stop(nome_dist, ": SampleLabel duplicado no dendrograma.")
    }

    pdf(file.path(dir_fig, paste0(prefixo, "_hclust_", nome_dist, ".pdf")),
        width = 9, height = 6)
    plot(hc_plot, main = paste(prefixo, "—", nome_dist),
         xlab = NULL,
         sub = "Metodo: UPGMA (average) | rotulos: SampleLabel",
         cex = 0.8)
    dev.off()
    invisible(hc)
  }

  gerar_hclust(dist_bray,    "bray_rel")
  gerar_hclust(dist_jaccard, "jaccard_binary")

  # -------------------------------------------------------------------------
  # 9.10 AUDITORIA DAS PREMISSAS DOS TESTES
  # -------------------------------------------------------------------------
  # As curvas de rarefacao/iNEXT foram externalizadas para o script
  # diagnostico_completude_amostral_inext.R. O Script 06 permanece responsavel
  # apenas pela construcao e validacao dos objetos phyloseq.

  premissas_testes <- auditar_premissas_testes(
    ps_obj = ps_obj,
    prefixo = prefixo,
    output_path = output_path
  )
  cat("\n--- Auditoria das premissas ---\n")
  print(premissas_testes)
  cat("\n")

  log_msg(paste("Analise concluida:", prefixo), "FINAL")
  cat("\n")

  invisible(list(
    ps           = ps_obj,
    ps_rel       = ps_rel,
    alfa         = alfa,
    dist_bray    = dist_bray,
    dist_jaccard = dist_jaccard,
    covariaveis  = resumo_cov,
    premissas    = premissas_testes
  ))
}

###############################################################################
# 10. ANALISE PRINCIPAL — CORE9
###############################################################################

res_core9 <- analisar_funil_invertido(
  ps_obj              = ps_core9,
  prefixo             = "core9",
  dir_fig             = plot_core9,
  executar_testes     = FALSE,
  incluir_covariaveis = TRUE)

###############################################################################
# 11. ANALISE DE SENSIBILIDADE — PLUS10
###############################################################################

res_plus10 <- analisar_funil_invertido(
  ps_obj              = ps_plus10,
  prefixo             = "plus10",
  dir_fig             = plot_plus10,
  executar_testes     = FALSE,
  incluir_covariaveis = TRUE)

###############################################################################
# 12. COMPARACAO CORE9 VS PLUS10
###############################################################################

cat("=== COMPARACAO CORE9 VS PLUS10 ===\n\n")

comp_resumo <- data.frame(
  Conjunto     = c("core9", "plus10"),
  Amostras     = c(nsamples(ps_core9),  nsamples(ps_plus10)),
  ASVs         = c(ntaxa(ps_core9),     ntaxa(ps_plus10)),
  Reads_totais = c(sum(sample_sums(ps_core9)), sum(sample_sums(ps_plus10))),
  Prof_min     = c(min(sample_sums(ps_core9)), min(sample_sums(ps_plus10))),
  Prof_mediana = c(median(sample_sums(ps_core9)), median(sample_sums(ps_plus10))),
  Prof_max     = c(max(sample_sums(ps_core9)), max(sample_sums(ps_plus10))),
  stringsAsFactors = FALSE)
write.csv(comp_resumo, file.path(output_path, "comparacao_core9_vs_plus10.csv"),
          row.names = FALSE)
print(comp_resumo)

# Sensibilidade: destacar a amostra fasciculata_auxiliar dentro do plus10
alfa_plus <- res_plus10$alfa
alfa_plus$Eh_auxiliar <- grepl(PADRAO_AUXILIAR, alfa_plus$SampleID)

write.csv(
  alfa_plus,
  file.path(output_path, "sensibilidade_plus10_alpha_auxiliar_destacada.csv"),
  row.names = FALSE
)

p_comp <- ggplot(
  alfa_plus,
  aes(
    x = BeeSpecies,
    y = Shannon,
    color = Eh_auxiliar,
    shape = Run
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9,
    position = position_jitter(width = 0.08)
  ) +

  ggrepel::geom_text_repel(
    aes(label = SampleLabel),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  )  +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(
    title = "Analise de sensibilidade: inclusao da amostra Auxiliar",
    subtitle = "Amostra da segunda corrida destacada; testes plus10 sao exploratorios e separados",
    x = NULL,
    y = "Shannon",
    color = "Amostra Auxiliar",
    shape = "Corrida"
  )

salvar_plot(
  p_comp,
  "sensibilidade_plus10_shannon_auxiliar_destacada.pdf",
  plot_path,
  width = 9,
  height = 6
)
cat("\n")
###############################################################################
# 13. EXPORTACOES FINAIS
###############################################################################

cat("=== EXPORTACOES FINAIS ===\n\n")

exportar_componentes_phyloseq(ps_all10,  "all10")
exportar_componentes_phyloseq(ps_core9,  "core9")
exportar_componentes_phyloseq(ps_plus10, "plus10")

###############################################################################
# 14. METADADOS DE EXECUCAO
###############################################################################

run_metadata <- data.frame(
  Script                 = "6_phyloseq_V",
  Versao                 = VERSAO,
  Data_execucao          = DATA_EXECUCAO,
  Amostras_all10         = nsamples(ps_all10),
  Amostras_core9         = nsamples(ps_core9),
  Amostras_plus10        = nsamples(ps_plus10),
  ASVs_all10             = ntaxa(ps_all10),
  ASVs_core9             = ntaxa(ps_core9),
  ASVs_plus10            = ntaxa(ps_plus10),
  Amostra_auxiliar         = sample_auxiliar,
  Analise_principal      = "core9 (9 amostras, primeira run)",
  Analise_sensibilidade  = "plus10 (10 amostras; testes exploratorios separados; Run nao estimavel)",
  Taxa_source            = "taxa_consenso_final.rds",
  Hierarquia_taxonomica  = descricao_hierarquia_taxa,
  GSR07_phyloseq_sensibilidade = gsr_taxonomia_disponivel,
  GSR07_uso = if (gsr_taxonomia_disponivel) {
    "taxonomia alternativa; mesmas contagens, metadados e ASVs do consenso principal"
  } else {
    "nao disponivel"
  },
  Testes_plus10          = "Executados separadamente como sensibilidade exploratoria; sem ajuste para Run",
  NA_em_taxa_preservado  = "TRUE — conversao para Unclassified somente em plots",
  stringsAsFactors = FALSE)
write.csv(run_metadata,
          file.path(output_path, "metadata_execucao_script6.csv"),
          row.names = FALSE)
log_msg("metadata_execucao_script6.csv salvo", "SAVE")

cat("\n=============================================================\n")
cat("PHYLOSEQ — CONCLUIDO\n")
cat(sprintf("ps_all10 : %d amostras | %d ASVs\n",
            nsamples(ps_all10), ntaxa(ps_all10)))
cat(sprintf("ps_core9 : %d amostras | %d ASVs (analise principal)\n",
            nsamples(ps_core9), ntaxa(ps_core9)))
cat(sprintf("ps_plus10: %d amostras | %d ASVs (sensibilidade)\n",
            nsamples(ps_plus10), ntaxa(ps_plus10)))
cat("Species bacteriana: NA nos objetos; Unclassified somente em plots.\n")
cat("Phyloseq GSR 0.7 de sensibilidade:", gsr_taxonomia_disponivel, "\n")
cat("=============================================================\n\n")

rm(seqtab, taxa_mat)
gc(verbose = FALSE)
log_msg("Finalizado", "FINAL")
})
