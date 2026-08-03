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
run_pipeline_script("11a_faprotax_graficos.R", "faprotax", function(ctx) {
###############################################################################
# SCRIPT 11 — FAPROTAX / INFERENCIA FUNCIONAL ECOLOGICA 16S
# Versao revisada: microeco::trans_func
#
# Projeto: Microbioma de mel de Melipona spp. (abelhas sem ferrao)
#
# Entradas principais, geradas pelos scripts anteriores:
#   output_<V>/06_phyloseq/phyloseq_core9_primeira_run.rds
#   output_<V>/06_phyloseq/phyloseq_plus10_com_auxiliar.rds
#   output_<V>/01_dada2/ASV_sequences.tsv
#
# Objetivo:
#   Inferir potencial funcional ecologico/metabolico por FAPROTAX a partir da
#   taxonomia consenso integrada aos objetos phyloseq.
#
# Motor FAPROTAX:
#   microeco::trans_func$new(dataset = microtable)
#   trans_func$cal_func(prok_database = "FAPROTAX")
#
# Decisao metodologica:
#   - Este script NAO usa faprotax_function_taxon_map.tsv.
#   - A anotacao FAPROTAX e feita pelo motor interno do microeco.
#   - As funcoes FAPROTAX nao sao mutuamente exclusivas: uma mesma ASV pode
#     contribuir para mais de uma funcao. Portanto, a soma das abundancias
#     relativas entre funcoes pode exceder 100%.
#   - Os testes inferenciais sao exploratorios devido ao baixo n por grupo.
#
# Estrutura analitica:
#   core9  — 9 amostras principais: descritivo + inferencia exploratoria.
#   plus10 — 10 amostras: sensibilidade exploratoria com testes separados; Run nao estimavel.
#
# Saidas:
#   output_<V>/11_faprotax/core9_*
#   output_<V>/11_faprotax/plus10_*
#   output_<V>/11_faprotax/figuras/*
#
# Observacao:
#   Este script nao altera arquivos dos Scripts 01-10. Ele apenas le entradas
#   existentes e escreve novas saidas em subpastas proprias.
###############################################################################

options(encoding = "UTF-8", stringsAsFactors = FALSE, warn = 1)

###############################################################################
# 0. PACOTES
###############################################################################

pacotes_obrigatorios <- c(
  "phyloseq", "microeco", "dplyr", "tidyr", "ggplot2", "vegan", "scales", "ggrepel"
)

for (pkg in pacotes_obrigatorios) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Pacote ausente: ", pkg, "\n",
      "Instale antes de rodar o Script 11. Exemplo: install.packages('", pkg, "')",
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(microeco)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(scales)
  library(ggrepel)
})

###############################################################################
# 1. PARAMETROS GLOBAIS
###############################################################################

VERSAO <- "4.1_core9_plus10_testes_relatorio"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Caminho principal sobrescrevivel por variavel de ambiente.
base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$output_root
plot_root <- ctx$stage$figures
faprotax_out <- ctx$stage$root
faprotax_plot <- ctx$stage$figures
for (d in c(faprotax_out, faprotax_plot)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
}
arq_core9 <- ctx$contracts[["phyloseq_core9"]]
arq_plus10 <- ctx$contracts[["phyloseq_plus10"]]
arq_asvmap <- ctx$contracts[["asv_sequences"]]
arq_meta_taxa <- ctx$contracts[["metadata_taxonomia"]]

TOP_N_FUNCOES <- 40L
TOP_N_TAXA_LEGENDA <- 2L
USAR_LEGENDA_AUXILIAR_TAXON <- TRUE
N_PERM <- 9999L
ALPHA  <- 0.05
SEED   <- 1234L

RANKS_CANONICOS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Paletas mantidas coerentes com os graficos principais do pipeline.
CORES_ESP <- c(
  "Melipona fasciculata" = "#FF6361",
  "Melipona scutellaris" = "#000489",
  "Melipona subnitida"   = "#21EBC9"
)

LABELS_ESP <- c(
  "Melipona fasciculata" = "M. fasciculata",
  "Melipona scutellaris" = "M. scutellaris",
  "Melipona subnitida"   = "M. subnitida"
)

CORES_STATUS <- c(
  "Introduzida" = "#FF6361",
  "Nativa"      = "#000489"
)

LABELS_STATUS <- c(
  "Introduzida" = "Introduzida",
  "Nativa"      = "Nativa"
)

CORES_MULTI <- c(
  "#d45087", "#6A3D9A", "#1F78B4", "#33A02C", "#FF7F00",
  "#00BCD4", "#F06292", "#FFD600", "#00897B", "#FB8C00",
  "#8E24AA", "#1565C0", "#D81B60", "#43A047", "#2f4b7c",
  "#546E7A", "#665191", "#0288D1", "#558B2F", "#EF6C00",
  "#4527A0", "#00695C", "#AD1457", "#827717", "#37474F",
  "#1A237E", "#004D40", "#F57F17", "#880E4F"
)

TEMA_PUB <- theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "grey60"),
    strip.text       = element_text(face = "bold.italic", size = 10),
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(face = "bold"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey30")
  )

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", round(p, 3), ifelse(p < ALPHA, "*", " n.s."))
}

# Rótulo canônico das amostras. SampleID permanece como chave interna;
# SampleLabels é aceito somente como alias e padronizado para SampleLabel.
validar_samplelabel_meta <- function(meta_df) {
  if (!"SampleID" %in% colnames(meta_df)) meta_df$SampleID <- rownames(meta_df)
  if (!"SampleLabel" %in% colnames(meta_df) && "SampleLabels" %in% colnames(meta_df)) {
    colnames(meta_df)[colnames(meta_df) == "SampleLabels"] <- "SampleLabel"
  }
  if (!"SampleLabel" %in% colnames(meta_df))
    stop("SampleLabel ausente no sample_data; execute novamente o Script 06 atualizado.", call. = FALSE)
  meta_df$SampleID <- trimws(as.character(meta_df$SampleID))
  if (any(is.na(meta_df$SampleID)) || any(meta_df$SampleID == ""))
    stop("SampleID ausente ou vazio no sample_data.", call. = FALSE)
  if (anyDuplicated(meta_df$SampleID) > 0)
    stop("SampleID duplicado no sample_data.", call. = FALSE)
  meta_df$SampleLabel <- trimws(as.character(meta_df$SampleLabel))
  if (any(is.na(meta_df$SampleLabel)) || any(meta_df$SampleLabel == ""))
    stop("SampleLabel ausente ou vazio no sample_data.", call. = FALSE)
  if (anyDuplicated(meta_df$SampleLabel) > 0)
    stop("SampleLabel duplicado no sample_data; os rotulos dos graficos devem ser unicos.", call. = FALSE)
  meta_df
}

gerar_rotulos_amostras <- function(meta_df) {
  meta_df <- validar_samplelabel_meta(meta_df)
  as.character(meta_df$SampleLabel)
}

rotular_sampleid <- function(sample_id, meta_df) {
  meta_df <- validar_samplelabel_meta(meta_df)
  mapa <- setNames(meta_df$SampleLabel, meta_df$SampleID)
  ids <- as.character(sample_id)
  out <- unname(mapa[ids])
  out[is.na(out) | out == ""] <- ids[is.na(out) | out == ""]
  out
}

###############################################################################
# 2. FUNCOES UTILITARIAS — I/O, VALIDACAO, LOG E CHECKPOINTS
###############################################################################

log_msg <- function(msg, tipo = "INFO") {
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))
}

# Funcoes compartilhadas para testes exatos por enumeracao. O arquivo e
# localizado pelo diretorio informado pelo orquestrador ou pelo proprio script.
arq_funcoes_exatas <- file.path(.pipeline_lib_dir, "funcoes_estatisticas_exatas.R")
if (!file.exists(arq_funcoes_exatas)) {
  stop("funcoes_estatisticas_exatas.R nao encontrado em: ", arq_funcoes_exatas, call. = FALSE)
}
sys.source(arq_funcoes_exatas, envir = .GlobalEnv)

valor_escalar_seguro <- function(x, padrao = NA_character_) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(padrao)
  as.character(x[[1L]])
}

OBS_CONFUNDIMENTO <- paste(
  "Nativo_Introduzido e derivado de BeeSpecies;",
  "qualquer comparacao por origem e descritiva/exploratoria e nao estima",
  "um efeito independente da especie."
)

validar_arquivo <- function(arq, desc = basename(arq)) {
  if (!file.exists(arq)) {
    stop(desc, " nao encontrado: ", arq, call. = FALSE)
  }
  tamanho <- file.size(arq)
  if (is.na(tamanho) || tamanho == 0L) {
    stop(desc, " existe, mas esta vazio ou inacessivel: ", arq, call. = FALSE)
  }
  log_msg(sprintf("%s OK (%.2f MB)", desc, file.size(arq) / 1024 / 1024), "OK")
  invisible(TRUE)
}

salvar_csv <- function(df, arq) {
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, arq, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
  if (!file.exists(arq) || file.size(arq) == 0) {
    stop("Falha ao salvar CSV: ", arq, call. = FALSE)
  }
  invisible(arq)
}

salvar_tsv <- function(df, arq) {
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  write.table(df, arq, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  if (!file.exists(arq) || file.size(arq) == 0) {
    stop("Falha ao salvar TSV: ", arq, call. = FALSE)
  }
  invisible(arq)
}

salvar_rds <- function(obj, arq) {
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, arq)
  if (!file.exists(arq) || file.size(arq) == 0) {
    stop("Falha ao salvar RDS: ", arq, call. = FALSE)
  }
  invisible(arq)
}

salvar_matriz <- function(mat, arq) {
  df <- cbind(Feature = rownames(mat), as.data.frame(mat, check.names = FALSE))
  salvar_csv(df, arq)
}

salvar_plot_duplo <- function(p, nome_base, plot_dir, width = 10, height = 7) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  arq_pdf <- file.path(plot_dir, paste0(nome_base, ".pdf"))
  arq_png <- file.path(plot_dir, paste0(nome_base, ".png"))

  ggplot2::ggsave(arq_pdf, plot = p, width = width, height = height, bg = "white")
  ggplot2::ggsave(arq_png, plot = p, width = width, height = height, dpi = 300, bg = "white")

  if (!file.exists(arq_pdf) || !file.exists(arq_png)) {
    stop("Falha ao salvar grafico: ", nome_base, call. = FALSE)
  }

  invisible(c(pdf = arq_pdf, png = arq_png))
}

salvar_erro_etapa <- function(prefixo, etapa, erro, out_dir = faprotax_out) {
  df <- data.frame(
    Prefixo = prefixo,
    Etapa = etapa,
    Erro = conditionMessage(erro),
    Data_execucao = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  arq <- file.path(out_dir, paste0(prefixo, "_ERRO_", etapa, ".csv"))
  salvar_csv(df, arq)
  invisible(arq)
}

executar_plot_seguro <- function(expr, prefixo, etapa) {
  tryCatch(
    expr,
    error = function(e) {
      log_msg(paste(prefixo, etapa, "falhou; ver arquivo de erro:", conditionMessage(e)), "WARN")
      salvar_erro_etapa(prefixo, paste0("plot_", etapa), e)
      invisible(NULL)
    }
  )
}

matriz_taxa_linhas <- function(ps_obj) {
  mat <- as(phyloseq::otu_table(ps_obj), "matrix")
  if (!phyloseq::taxa_are_rows(ps_obj)) mat <- t(mat)
  storage.mode(mat) <- "numeric"

  if (is.null(rownames(mat)) || any(rownames(mat) == "")) {
    stop("otu_table sem rownames de ASV/taxon.", call. = FALSE)
  }
  if (is.null(colnames(mat)) || any(colnames(mat) == "")) {
    stop("otu_table sem colnames de amostras.", call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop("otu_table contem valores nao finitos.", call. = FALSE)
  }
  if (any(mat < 0, na.rm = TRUE)) {
    stop("otu_table contem contagens negativas.", call. = FALSE)
  }

  mat
}

adicionar_status_origem <- function(meta_df) {
  if (!"BeeSpecies" %in% colnames(meta_df)) {
    stop("BeeSpecies ausente no sample_data; nao e possivel derivar Nativo_Introduzido.", call. = FALSE)
  }

  meta_df$Nativo_Introduzido <- ifelse(
    as.character(meta_df$BeeSpecies) == "Melipona fasciculata",
    "Introduzida",
    "Nativa"
  )

  meta_df$Nativo_Introduzido <- factor(
    meta_df$Nativo_Introduzido,
    levels = c("Introduzida", "Nativa")
  )

  meta_df
}

validar_phyloseq_entrada <- function(ps_obj, nome_obj) {
  if (!inherits(ps_obj, "phyloseq")) {
    stop(nome_obj, ": objeto nao e phyloseq.", call. = FALSE)
  }
  if (phyloseq::nsamples(ps_obj) < 2) {
    stop(nome_obj, ": menos de 2 amostras.", call. = FALSE)
  }
  if (phyloseq::ntaxa(ps_obj) < 2) {
    stop(nome_obj, ": menos de 2 ASVs.", call. = FALSE)
  }
  if (is.null(phyloseq::otu_table(ps_obj, errorIfNULL = FALSE)) ||
      is.null(phyloseq::sample_data(ps_obj, errorIfNULL = FALSE)) ||
      is.null(phyloseq::tax_table(ps_obj, errorIfNULL = FALSE))) {
    stop(nome_obj, ": otu_table, sample_data ou tax_table ausente.", call. = FALSE)
  }
  mat <- matriz_taxa_linhas(ps_obj)
  if (anyNA(mat) || any(!is.finite(mat)) || any(mat < 0)) {
    stop(nome_obj, ": otu_table contem NA, valor nao finito ou negativo.", call. = FALSE)
  }
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    stop(nome_obj, ": ASV ou amostra com soma zero.", call. = FALSE)
  }
  if (max(abs(mat - round(mat))) > 1e-8) {
    stop(nome_obj, ": FAPROTAX deve receber contagens brutas inteiras.", call. = FALSE)
  }
  req_meta <- c("Run", "BeeSpecies")
  faltam_meta <- setdiff(req_meta, phyloseq::sample_variables(ps_obj))
  if (length(faltam_meta) > 0L) {
    stop(nome_obj, ": sample_data sem ", paste(faltam_meta, collapse = ", "), call. = FALSE)
  }
  ranks_faltantes <- setdiff(RANKS_CANONICOS, phyloseq::rank_names(ps_obj))
  if (length(ranks_faltantes) > 0L) {
    stop(nome_obj, ": tax_table sem rank(s): ", paste(ranks_faltantes, collapse = ", "), call. = FALSE)
  }
  if (is.null(phyloseq::taxa_names(ps_obj)) || any(phyloseq::taxa_names(ps_obj) == "")) {
    stop(nome_obj, ": taxa_names ausentes ou vazios.", call. = FALSE)
  }
  if (anyDuplicated(phyloseq::taxa_names(ps_obj)) > 0 ||
      anyDuplicated(phyloseq::sample_names(ps_obj)) > 0) {
    stop(nome_obj, ": taxa_names ou sample_names duplicados.", call. = FALSE)
  }
  invisible(TRUE)
}

validar_mapa_asv <- function(arq_asvmap, taxa_ids) {
  validar_arquivo(arq_asvmap, "ASV_sequences.tsv")

  asv_map <- read.delim(
    arq_asvmap,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  req <- c("ASV_ID", "Sequence", "Origem")
  faltam <- setdiff(req, colnames(asv_map))
  if (length(faltam) > 0) {
    stop("ASV_sequences.tsv sem coluna(s): ", paste(faltam, collapse = ", "), call. = FALSE)
  }

  if (anyDuplicated(asv_map$ASV_ID) > 0) {
    stop("ASV_sequences.tsv contem ASV_ID duplicados.", call. = FALSE)
  }
  if (anyDuplicated(asv_map$Sequence) > 0) {
    stop("ASV_sequences.tsv contem sequencias duplicadas.", call. = FALSE)
  }

  ids_ok <- all(taxa_ids %in% asv_map$ASV_ID)
  seqs_ok <- all(taxa_ids %in% asv_map$Sequence)
  if (!ids_ok && !seqs_ok) {
    stop(
      "taxa_names do phyloseq nao coincidem integralmente com ASV_ID nem com Sequence. ",
      "A rastreabilidade ASV foi interrompida.",
      call. = FALSE
    )
  }

  asv_map
}

###############################################################################
# 3. CONVERSAO PHYLOSEQ -> MICROECO E PREDICAO FAPROTAX
###############################################################################

limpar_taxonomia_para_microeco <- function(tax_df, prefixo) {
  tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE, check.names = FALSE)

  for (rk in RANKS_CANONICOS) {
    if (!rk %in% colnames(tax_df)) {
      tax_df[[rk]] <- ""
    }
  }

  tax_df <- tax_df[, unique(c(RANKS_CANONICOS, setdiff(colnames(tax_df), RANKS_CANONICOS))), drop = FALSE]

  for (cc in colnames(tax_df)) {
    tax_df[[cc]] <- as.character(tax_df[[cc]])
    tax_df[[cc]][is.na(tax_df[[cc]])] <- ""
    tax_df[[cc]] <- trimws(tax_df[[cc]])
    tax_df[[cc]][tax_df[[cc]] %in% c("NA", "Na", "na", "NULL", "None")] <- ""
  }

  if (all(tax_df$Kingdom == "")) {
    tax_df$Kingdom <- "Bacteria"
    log_msg(paste(prefixo, ": Kingdom ausente/vazio; definido como Bacteria para FAPROTAX 16S."), "WARN")
  }

  tax_df
}

phyloseq_para_microeco <- function(ps_obj, prefixo) {
  # microeco::microtable$new() exige data.frame tradicional, não matrix/tibble.
  # Algumas versões do microeco aplicam droplevels() internamente em colunas
  # categóricas; por isso, metadados e taxonomia são convertidos para factor
  # antes da criação do microtable. A otu_table permanece estritamente numérica.

  otu_mat <- matriz_taxa_linhas(ps_obj)
  otu_df <- as.data.frame(otu_mat, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(otu_df) <- rownames(otu_mat)

  # Garante que todas as colunas da matriz de contagens de ASVs sejam numéricas.
  otu_df[] <- lapply(otu_df, function(x) as.numeric(x))
  rownames(otu_df) <- rownames(otu_mat)

  tax_df <- as.data.frame(
    phyloseq::tax_table(ps_obj),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  tax_df <- limpar_taxonomia_para_microeco(tax_df, prefixo)

  sample_df <- as.data.frame(
    phyloseq::sample_data(ps_obj),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!identical(rownames(otu_df), rownames(tax_df))) {
    stop(prefixo, ": rownames da otu_table e tax_table nao coincidem.", call. = FALSE)
  }

  faltam_meta <- setdiff(colnames(otu_df), rownames(sample_df))
  if (length(faltam_meta) > 0) {
    stop(prefixo, ": amostras da otu_table ausentes no sample_data: ",
         paste(faltam_meta, collapse = ", "), call. = FALSE)
  }

  sample_df <- sample_df[colnames(otu_df), , drop = FALSE]

  # Checkpoints explícitos antes de chamar microeco.
  if (!is.data.frame(otu_df)) {
    stop(prefixo, ": otu_table nao esta em formato data.frame antes do microeco.", call. = FALSE)
  }
  if (!all(vapply(otu_df, is.numeric, logical(1)))) {
    stop(prefixo, ": otu_table contem coluna nao numerica antes do microeco.", call. = FALSE)
  }
  if (!is.data.frame(sample_df)) {
    stop(prefixo, ": sample_table nao esta em formato data.frame antes do microeco.", call. = FALSE)
  }
  if (!is.data.frame(tax_df)) {
    stop(prefixo, ": tax_table nao esta em formato data.frame antes do microeco.", call. = FALSE)
  }

  # Compatibilidade com versões do microeco que chamam droplevels() em colunas
  # individualmente. Nessas versões, colunas character E numeric em sample_table
  # podem disparar erro. Como o microeco só precisa dos metadados como grupos
  # aqui, convertemos TODAS as colunas de sample_table e tax_table para factor.
  # A otu_table permanece estritamente numérica.
  sample_df[] <- lapply(sample_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    factor(x)
  })

  tax_df[] <- lapply(tax_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    factor(x)
  })

  # Mantém rownames depois das conversões de coluna.
  rownames(sample_df) <- colnames(otu_df)
  rownames(tax_df) <- rownames(otu_df)

  mt <- microeco::microtable$new(
    otu_table = otu_df,
    sample_table = sample_df,
    tax_table = tax_df,
    auto_tidy = FALSE
  )

  # check_microtable() chamado por trans_func$new() também executa tidy_dataset().
  # Rodamos aqui para falhar cedo e com mensagem controlada.
  mt$tidy_dataset()

  if (!identical(rownames(mt$otu_table), rownames(mt$tax_table))) {
    stop(prefixo, ": microeco alterou desalinhamento entre otu_table e tax_table.", call. = FALSE)
  }
  if (!identical(colnames(mt$otu_table), rownames(mt$sample_table))) {
    stop(prefixo, ": microeco alterou desalinhamento entre otu_table e sample_table.", call. = FALSE)
  }

  mt
}

rodar_microeco_faprotax <- function(ps_obj, prefixo) {
  mt <- phyloseq_para_microeco(ps_obj, prefixo)

  tf <- microeco::trans_func$new(dataset = mt)

  for_what_atual <- tf$for_what
  for_what_txt <- if (
    is.null(for_what_atual) || length(for_what_atual) == 0L ||
    all(is.na(for_what_atual))
  ) {
    "NA"
  } else {
    as.character(for_what_atual[[1L]])
  }

  if (for_what_txt != "prok") {
    log_msg(
      paste0(prefixo, ": microeco detectou for_what = '", for_what_txt,
             "'; definido manualmente como 'prok' para FAPROTAX 16S."),
      "WARN"
    )
    tf$for_what <- "prok"
  }

  tf$cal_func(prok_database = "FAPROTAX")

  if (is.null(tf$res_func)) {
    stop(prefixo, ": microeco nao gerou tf$res_func.", call. = FALSE)
  }

  asv_func_bin <- as.data.frame(tf$res_func, check.names = FALSE)
  if (nrow(asv_func_bin) == 0 || ncol(asv_func_bin) == 0) {
    stop(prefixo, ": tf$res_func esta vazio.", call. = FALSE)
  }

  ids_asv <- rownames(mt$otu_table)
  n_match_linhas <- if (is.null(rownames(asv_func_bin))) 0L else sum(rownames(asv_func_bin) %in% ids_asv)
  n_match_colunas <- if (is.null(colnames(asv_func_bin))) 0L else sum(colnames(asv_func_bin) %in% ids_asv)
  if (n_match_colunas > n_match_linhas) {
    asv_func_bin <- as.data.frame(t(as.matrix(asv_func_bin)), check.names = FALSE)
    log_msg(paste(prefixo, ": tf$res_func transposto para ASVs nas linhas."), "INFO")
    n_match_linhas <- sum(rownames(asv_func_bin) %in% ids_asv)
  }
  if (n_match_linhas < 1L) {
    stop(
      prefixo, ": nao foi possivel identificar ASVs nas linhas ou colunas de tf$res_func.",
      call. = FALSE
    )
  }
  if (is.null(rownames(asv_func_bin)) || any(rownames(asv_func_bin) == "") ||
      anyDuplicated(rownames(asv_func_bin)) > 0L) {
    stop(prefixo, ": tf$res_func sem rownames de ASV validos e unicos.", call. = FALSE)
  }
  if (is.null(colnames(asv_func_bin)) || any(colnames(asv_func_bin) == "") ||
      anyDuplicated(colnames(asv_func_bin)) > 0L) {
    stop(prefixo, ": tf$res_func sem nomes de funcao validos e unicos.", call. = FALSE)
  }

  for (cc in colnames(asv_func_bin)) {
    asv_func_bin[[cc]] <- suppressWarnings(as.numeric(as.character(asv_func_bin[[cc]])))
  }
  if (anyNA(asv_func_bin) || any(!is.finite(as.matrix(asv_func_bin)))) {
    stop(prefixo, ": tf$res_func contem valor nao numerico ou nao finito.", call. = FALSE)
  }
  if (any(as.matrix(asv_func_bin) < 0)) {
    stop(prefixo, ": tf$res_func contem valor negativo.", call. = FALSE)
  }
  valores_nao_binarios <- setdiff(unique(as.numeric(as.matrix(asv_func_bin))), c(0, 1))
  if (length(valores_nao_binarios) > 0L) {
    log_msg(
      paste(prefixo, ": tf$res_func continha valores nao binarios; convertido para presenca/ausencia."),
      "WARN"
    )
  }
  asv_func_bin[] <- lapply(asv_func_bin, function(x) as.integer(x > 0))

  # Remove funcoes totalmente ausentes, mantendo a tabela mais interpretavel.
  asv_func_bin <- asv_func_bin[, colSums(asv_func_bin != 0) > 0, drop = FALSE]

  if (ncol(asv_func_bin) == 0) {
    stop(prefixo, ": nenhuma funcao FAPROTAX presente apos remover colunas zero.", call. = FALSE)
  }

  n_asv_mapeadas <- sum(rowSums(asv_func_bin != 0) > 0)
  if (n_asv_mapeadas == 0) {
    stop(prefixo, ": nenhuma ASV foi mapeada para funcoes FAPROTAX.", call. = FALSE)
  }

  # Calcula tambem a tabela de redundancia funcional do proprio microeco.
  # O script usa a matriz exata ASV x funcao para contagens/relativas, mas salva
  # res_func_FR como checkpoint metodologico e para reprodutibilidade.
  microeco_fr <- tryCatch({
    tf$cal_func_FR(abundance_weighted = TRUE, adj_tax = FALSE, perc = FALSE, remove_zero = TRUE)
    as.data.frame(tf$res_func_FR, check.names = FALSE)
  }, error = function(e) {
    log_msg(paste(prefixo, ": cal_func_FR falhou; analise segue com matriz binaria res_func:", conditionMessage(e)), "WARN")
    NULL
  })

  list(
    microtable = mt,
    trans_func = tf,
    asv_func_bin = asv_func_bin,
    microeco_fr = microeco_fr
  )
}

###############################################################################
# 4. COLAPSO FUNCIONAL, COBERTURA E ESTATISTICA
###############################################################################

calcular_abundancia_funcional_microeco <- function(otu_mat, asv_func_bin, prefixo) {
  comuns <- intersect(rownames(otu_mat), rownames(asv_func_bin))
  if (length(comuns) == 0) {
    stop(prefixo, ": nenhuma ASV em comum entre otu_table e matriz FAPROTAX do microeco.", call. = FALSE)
  }
  extras_func <- setdiff(rownames(asv_func_bin), rownames(otu_mat))
  if (length(extras_func) > 0L) {
    stop(prefixo, ": tf$res_func contem ASVs fora do universo do phyloseq.", call. = FALSE)
  }

  otu_use <- otu_mat[comuns, , drop = FALSE]
  fun_use <- as.matrix(asv_func_bin[comuns, , drop = FALSE])
  storage.mode(fun_use) <- "numeric"

  fun_counts <- t(fun_use) %*% otu_use
  fun_counts <- as.matrix(fun_counts)
  storage.mode(fun_counts) <- "numeric"

  # Remove funcoes sem contagem total.
  fun_counts <- fun_counts[rowSums(fun_counts) > 0, , drop = FALSE]

  if (nrow(fun_counts) == 0) {
    stop(prefixo, ": nenhuma funcao manteve abundancia apos colapso.", call. = FALSE)
  }

  total_reads <- colSums(otu_mat)
  fun_rel <- sweep(fun_counts, 2, total_reads, "/")
  fun_rel[is.na(fun_rel)] <- 0

  mapped_asvs <- rownames(asv_func_bin)[rowSums(asv_func_bin != 0) > 0]
  mapped_asvs <- intersect(mapped_asvs, rownames(otu_mat))
  mapped_reads <- colSums(otu_mat[mapped_asvs, , drop = FALSE])

  cobertura <- data.frame(
    SampleID = names(total_reads),
    Reads_total = as.numeric(total_reads),
    Reads_ASVs_mapeadas_unicas = as.numeric(mapped_reads),
    Pct_reads_mapeados_unicos = round(100 * mapped_reads / total_reads, 4),
    stringsAsFactors = FALSE
  )

  list(counts = fun_counts, relative = fun_rel, cobertura = cobertura)
}

asv_funcoes_long <- function(asv_func_bin) {
  df <- as.data.frame(asv_func_bin, check.names = FALSE)
  df$ASV_ID <- rownames(df)

  long <- tidyr::pivot_longer(
    df,
    cols = -ASV_ID,
    names_to = "Function",
    values_to = "Presente"
  )

  long <- long[!is.na(long$Presente) & long$Presente > 0, , drop = FALSE]
  long <- long[, c("ASV_ID", "Function"), drop = FALSE]
  as.data.frame(long, stringsAsFactors = FALSE)
}

limpar_taxon_texto_saida <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- gsub("^[a-zA-Z]__", "", x)
  x[x %in% c("NA", "Na", "na", "NULL", "None", "none")] <- ""
  x
}

taxon_informativo <- function(x) {
  x <- limpar_taxon_texto_saida(x)
  nzchar(x) &
    !grepl(
      "^(uncultured|unclassified|unidentified|unknown|metagenome|bacterium|bacteria|sp\\.?|)$",
      x,
      ignore.case = TRUE
    )
}

extrair_taxonomia_asv <- function(ps_obj) {
  tax_df <- as.data.frame(
    phyloseq::tax_table(ps_obj),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (rk in RANKS_CANONICOS) {
    if (!rk %in% colnames(tax_df)) tax_df[[rk]] <- ""
    tax_df[[rk]] <- limpar_taxon_texto_saida(tax_df[[rk]])
  }

  tax_df$ASV_ID <- rownames(tax_df)

  especie_resolvida <- tax_df$Species
  usar_binomio <- taxon_informativo(tax_df$Species) &
    taxon_informativo(tax_df$Genus) &
    !grepl("\\s", tax_df$Species)

  especie_resolvida[usar_binomio] <- paste(tax_df$Genus[usar_binomio], tax_df$Species[usar_binomio])

  tax_df$Taxon_associado <- dplyr::case_when(
    taxon_informativo(especie_resolvida) ~ especie_resolvida,
    taxon_informativo(tax_df$Genus)      ~ paste0(tax_df$Genus, " sp."),
    taxon_informativo(tax_df$Family)     ~ paste0(tax_df$Family, " (fam.)"),
    taxon_informativo(tax_df$Order)      ~ paste0(tax_df$Order, " (ord.)"),
    taxon_informativo(tax_df$Class)      ~ paste0(tax_df$Class, " (classe)"),
    taxon_informativo(tax_df$Phylum)     ~ paste0(tax_df$Phylum, " (filo)"),
    TRUE                                 ~ "Taxonomia não resolvida"
  )

  tax_df$Nivel_taxonomico_usado <- dplyr::case_when(
    taxon_informativo(especie_resolvida) ~ "Species",
    taxon_informativo(tax_df$Genus)      ~ "Genus",
    taxon_informativo(tax_df$Family)     ~ "Family",
    taxon_informativo(tax_df$Order)      ~ "Order",
    taxon_informativo(tax_df$Class)      ~ "Class",
    taxon_informativo(tax_df$Phylum)     ~ "Phylum",
    TRUE                                 ~ "Unresolved"
  )

  tax_df$Species_resolvida <- especie_resolvida

  tax_df[, c(
    "ASV_ID", RANKS_CANONICOS, "Species_resolvida",
    "Taxon_associado", "Nivel_taxonomico_usado"
  ), drop = FALSE]
}

criar_tabela_asv_func_taxonomia <- function(asv_funcoes, ps_obj, otu_mat) {
  tax_df <- extrair_taxonomia_asv(ps_obj)

  asv_reads <- data.frame(
    ASV_ID = rownames(otu_mat),
    Total_reads_ASV = as.numeric(rowSums(otu_mat)),
    Pct_reads_total_ASV = round(100 * rowSums(otu_mat) / sum(otu_mat), 6),
    stringsAsFactors = FALSE
  )

  asv_funcoes$ASV_ID <- as.character(asv_funcoes$ASV_ID)
  tax_df$ASV_ID <- as.character(tax_df$ASV_ID)

  out <- asv_funcoes |>
    dplyr::left_join(tax_df, by = "ASV_ID") |>
    dplyr::left_join(asv_reads, by = "ASV_ID") |>
    dplyr::arrange(Function, dplyr::desc(Total_reads_ASV), ASV_ID)

  as.data.frame(out, stringsAsFactors = FALSE)
}

calcular_contribuicao_taxa_funcoes <- function(asv_funcoes_tax, otu_mat, meta_df = NULL) {
  otu_df <- as.data.frame(otu_mat, check.names = FALSE)
  otu_df$ASV_ID <- rownames(otu_mat)

  otu_long <- tidyr::pivot_longer(
    otu_df,
    cols = -ASV_ID,
    names_to = "SampleID",
    values_to = "Reads"
  )

  cols_tax <- intersect(
    c("ASV_ID", "Function", RANKS_CANONICOS, "Species_resolvida",
      "Taxon_associado", "Nivel_taxonomico_usado"),
    colnames(asv_funcoes_tax)
  )

  leituras <- asv_funcoes_tax[, cols_tax, drop = FALSE] |>
    dplyr::inner_join(otu_long, by = "ASV_ID") |>
    dplyr::filter(.data$Reads > 0)

  if (nrow(leituras) == 0) {
    vazio <- data.frame()
    return(list(global = vazio, por_amostra = vazio))
  }

  global <- leituras |>
    dplyr::group_by(
      Function, Taxon_associado, Nivel_taxonomico_usado,
      Kingdom, Phylum, Class, Order, Family, Genus, Species_resolvida
    ) |>
    dplyr::summarise(
      ASVs = dplyr::n_distinct(ASV_ID),
      Reads_funcao_taxon = sum(Reads, na.rm = TRUE),
      Amostras_com_reads = dplyr::n_distinct(SampleID),
      .groups = "drop"
    )

  total_funcao <- leituras |>
    dplyr::group_by(Function) |>
    dplyr::summarise(
      Total_reads_funcao = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    )

  global <- global |>
    dplyr::left_join(total_funcao, by = "Function") |>
    dplyr::mutate(
      Pct_dentro_funcao = round(100 * Reads_funcao_taxon / Total_reads_funcao, 4)
    ) |>
    dplyr::arrange(Function, dplyr::desc(Pct_dentro_funcao), dplyr::desc(Reads_funcao_taxon))

  por_amostra <- leituras |>
    dplyr::group_by(
      SampleID, Function, Taxon_associado, Nivel_taxonomico_usado,
      Kingdom, Phylum, Class, Order, Family, Genus, Species_resolvida
    ) |>
    dplyr::summarise(
      ASVs = dplyr::n_distinct(ASV_ID),
      Reads_funcao_taxon_sample = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    )

  total_funcao_sample <- leituras |>
    dplyr::group_by(SampleID, Function) |>
    dplyr::summarise(
      Total_reads_funcao_sample = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    )

  por_amostra <- por_amostra |>
    dplyr::left_join(total_funcao_sample, by = c("SampleID", "Function")) |>
    dplyr::mutate(
      Pct_dentro_funcao_sample = round(
        100 * Reads_funcao_taxon_sample / Total_reads_funcao_sample,
        4
      )
    )

  if (!is.null(meta_df)) {
    meta_cols <- intersect(
      c("SampleID", "SampleLabel", "BeeSpecies", "Nativo_Introduzido"),
      colnames(meta_df)
    )
    por_amostra <- por_amostra |>
      dplyr::left_join(meta_df[, meta_cols, drop = FALSE], by = "SampleID")
  }

  por_amostra <- por_amostra |>
    dplyr::arrange(SampleID, Function, dplyr::desc(Pct_dentro_funcao_sample))

  list(
    global = as.data.frame(global, stringsAsFactors = FALSE),
    por_amostra = as.data.frame(por_amostra, stringsAsFactors = FALSE)
  )
}


# Cria rótulos integrados para gráficos: Função metabólica + táxon associado.
# Por padrão, mostra os 2 táxons mais abundantes dentro de cada função.
quebrar_rotulo <- function(x, largura = 58) {
  vapply(
    as.character(x),
    function(z) paste(strwrap(z, width = largura), collapse = "\n"),
    character(1)
  )
}

criar_labels_funcao_taxon <- function(contrib_global,
                                      funcoes,
                                      top_n_taxa = TOP_N_TAXA_LEGENDA,
                                      incluir_pct = TRUE,
                                      largura = 58) {
  funcoes <- as.character(funcoes)
  labels_padrao <- setNames(funcoes, funcoes)

  if (is.null(contrib_global) || nrow(contrib_global) == 0) {
    return(labels_padrao)
  }

  req <- c("Function", "Taxon_associado", "Pct_dentro_funcao", "Reads_funcao_taxon")
  if (!all(req %in% colnames(contrib_global))) {
    warning("contrib_global sem colunas esperadas; legenda usará apenas funções.")
    return(labels_padrao)
  }

  cg <- contrib_global |>
    dplyr::filter(.data$Function %in% funcoes) |>
    dplyr::mutate(
      Taxon_associado = ifelse(
        is.na(.data$Taxon_associado) | !nzchar(as.character(.data$Taxon_associado)),
        "Taxonomia não resolvida",
        as.character(.data$Taxon_associado)
      )
    ) |>
    dplyr::group_by(.data$Function) |>
    dplyr::arrange(
      dplyr::desc(.data$Pct_dentro_funcao),
      dplyr::desc(.data$Reads_funcao_taxon),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n_taxa) |>
    dplyr::mutate(
      Taxon_legenda = if (isTRUE(incluir_pct)) {
        paste0(.data$Taxon_associado, " (", round(.data$Pct_dentro_funcao, 1), "%)")
      } else {
        .data$Taxon_associado
      }
    ) |>
    dplyr::summarise(
      Taxa_legenda = paste(.data$Taxon_legenda, collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Label_funcao_taxon = paste0(.data$Function, " — ", .data$Taxa_legenda)
    )

  labels <- labels_padrao
  idx <- match(cg$Function, names(labels))
  labels[idx[!is.na(idx)]] <- cg$Label_funcao_taxon[!is.na(idx)]
  labels <- quebrar_rotulo(labels, largura = largura)
  labels
}

# Legenda auxiliar: cria códigos curtos F01, F02, ... para os gráficos e
# uma tabela/figura separada associando cada código à função metabólica e aos
# principais microrganismos contribuintes. Esta abordagem evita que legends
# longas sejam cortadas em barplots e heatmaps.
criar_legenda_funcao_taxon_df <- function(contrib_global,
                                          funcoes,
                                          top_n_taxa = TOP_N_TAXA_LEGENDA,
                                          incluir_pct = TRUE) {
  funcoes <- unique(as.character(funcoes))

  base <- data.frame(
    Function = funcoes,
    Funcao_codigo = sprintf("F%02d", seq_along(funcoes)),
    stringsAsFactors = FALSE
  )

  if (is.null(contrib_global) || nrow(contrib_global) == 0) {
    base$Taxa_associados <- "Taxonomia não resolvida ou contribuição indisponível"
    base$Taxon_principal <- "Taxonomia não resolvida"
    base$Nivel_taxonomico_principal <- "Unresolved"
    base$Pct_taxon_principal_dentro_funcao <- NA_real_
    base$Legenda_completa <- paste0(base$Funcao_codigo, " | ", base$Function, " — ", base$Taxa_associados)
    return(base)
  }

  req <- c("Function", "Taxon_associado", "Nivel_taxonomico_usado", "Pct_dentro_funcao", "Reads_funcao_taxon")
  faltam <- setdiff(req, colnames(contrib_global))
  if (length(faltam) > 0) {
    warning("contrib_global sem coluna(s): ", paste(faltam, collapse = ", "),
            "; legenda auxiliar será gerada sem táxons.")
    base$Taxa_associados <- "Taxonomia não resolvida ou contribuição indisponível"
    base$Taxon_principal <- "Taxonomia não resolvida"
    base$Nivel_taxonomico_principal <- "Unresolved"
    base$Pct_taxon_principal_dentro_funcao <- NA_real_
    base$Legenda_completa <- paste0(base$Funcao_codigo, " | ", base$Function, " — ", base$Taxa_associados)
    return(base)
  }

  cg <- contrib_global |>
    dplyr::filter(.data$Function %in% funcoes) |>
    dplyr::mutate(
      Taxon_associado = as.character(.data$Taxon_associado),
      Taxon_associado = ifelse(
        is.na(.data$Taxon_associado) | !nzchar(.data$Taxon_associado),
        "Taxonomia não resolvida",
        .data$Taxon_associado
      ),
      Nivel_taxonomico_usado = as.character(.data$Nivel_taxonomico_usado),
      Nivel_taxonomico_usado = ifelse(
        is.na(.data$Nivel_taxonomico_usado) | !nzchar(.data$Nivel_taxonomico_usado),
        "Unresolved",
        .data$Nivel_taxonomico_usado
      )
    ) |>
    dplyr::group_by(.data$Function) |>
    dplyr::arrange(
      dplyr::desc(.data$Pct_dentro_funcao),
      dplyr::desc(.data$Reads_funcao_taxon),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n_taxa) |>
    dplyr::mutate(
      Taxon_label = if (isTRUE(incluir_pct)) {
        paste0(.data$Taxon_associado, " [", .data$Nivel_taxonomico_usado, "; ",
               round(.data$Pct_dentro_funcao, 1), "%]")
      } else {
        paste0(.data$Taxon_associado, " [", .data$Nivel_taxonomico_usado, "]")
      }
    ) |>
    dplyr::summarise(
      Taxa_associados = paste(.data$Taxon_label, collapse = "; "),
      Taxon_principal = dplyr::first(.data$Taxon_associado),
      Nivel_taxonomico_principal = dplyr::first(.data$Nivel_taxonomico_usado),
      Pct_taxon_principal_dentro_funcao = dplyr::first(.data$Pct_dentro_funcao),
      Reads_taxon_principal = dplyr::first(.data$Reads_funcao_taxon),
      .groups = "drop"
    )

  out <- base |>
    dplyr::left_join(cg, by = "Function") |>
    dplyr::mutate(
      Taxa_associados = ifelse(
        is.na(.data$Taxa_associados) | !nzchar(.data$Taxa_associados),
        "Taxonomia não resolvida ou sem reads atribuídos",
        .data$Taxa_associados
      ),
      Taxon_principal = ifelse(
        is.na(.data$Taxon_principal) | !nzchar(.data$Taxon_principal),
        "Taxonomia não resolvida",
        .data$Taxon_principal
      ),
      Nivel_taxonomico_principal = ifelse(
        is.na(.data$Nivel_taxonomico_principal) | !nzchar(.data$Nivel_taxonomico_principal),
        "Unresolved",
        .data$Nivel_taxonomico_principal
      ),
      Legenda_completa = paste0(.data$Funcao_codigo, " | ", .data$Function, " — ", .data$Taxa_associados)
    )

  as.data.frame(out, stringsAsFactors = FALSE)
}

aplicar_codigo_funcao <- function(df, legenda_df) {
  if (is.null(legenda_df) || nrow(legenda_df) == 0) {
    df$Function_codigo <- as.character(df$Function)
    return(df)
  }

  mapa <- setNames(legenda_df$Funcao_codigo, legenda_df$Function)
  cod <- mapa[as.character(df$Function)]
  cod[is.na(cod)] <- as.character(df$Function)[is.na(cod)]
  df$Function_codigo <- cod
  df
}

plot_legenda_funcao_taxon_auxiliar <- function(legenda_df, prefixo, plot_dir) {
  if (is.null(legenda_df) || nrow(legenda_df) == 0) return(invisible(NULL))

  df <- legenda_df
  df$Linha <- rev(seq_len(nrow(df)))
  df$Function_plot <- quebrar_rotulo(df$Function, largura = 42)
  df$Taxa_plot <- quebrar_rotulo(df$Taxa_associados, largura = 72)

  altura <- max(5.5, 1.4 + 0.38 * nrow(df))

  p <- ggplot(df) +
    geom_rect(
      aes(xmin = -0.02, xmax = 1.64, ymin = Linha - 0.45, ymax = Linha + 0.45),
      fill = "white",
      color = NA,
      show.legend = FALSE
    ) +
    geom_text(aes(x = 0.00, y = Linha, label = Funcao_codigo),
              hjust = 0, size = 3.2, fontface = "bold") +
    geom_text(aes(x = 0.14, y = Linha, label = Function_plot),
              hjust = 0, size = 2.8, fontface = "bold") +
    geom_text(aes(x = 0.72, y = Linha, label = Taxa_plot),
              hjust = 0, size = 2.7) +
    annotate("text", x = 0.00, y = max(df$Linha) + 0.75,
             label = "Código", hjust = 0, size = 3.3, fontface = "bold") +
    annotate("text", x = 0.14, y = max(df$Linha) + 0.75,
             label = "Função metabólica", hjust = 0, size = 3.3, fontface = "bold") +
    annotate("text", x = 0.72, y = max(df$Linha) + 0.75,
             label = paste0("Microrganismos associados — top ", TOP_N_TAXA_LEGENDA,
                            " por contribuição dentro da função"),
             hjust = 0, size = 3.3, fontface = "bold") +
    coord_cartesian(xlim = c(-0.02, 1.64),
                    ylim = c(0.5, max(df$Linha) + 1.15),
                    clip = "off") +
    theme_void(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(size = 10, color = "grey30"),
      plot.margin      = margin(10, 20, 10, 10)
    ) +
    labs(
      title = paste0("Legenda auxiliar — função metabólica e microrganismo associado: ", prefixo),
      subtitle = "Os códigos F01, F02, ... são usados nos barplots e heatmaps para manter os gráficos legíveis. Percentuais entre colchetes indicam contribuição do táxon dentro da função."
    )

  salvar_plot_duplo(p, paste0(prefixo, "_G00_legenda_auxiliar_funcao_taxon"), plot_dir, width = 17, height = altura)
}

resumo_funcoes <- function(fun_counts, fun_rel) {
  data.frame(
    Function = rownames(fun_counts),
    Total_counts = rowSums(fun_counts),
    Prevalence_samples = rowSums(fun_counts > 0),
    Mean_relative = rowMeans(fun_rel),
    Max_relative = apply(fun_rel, 1, max),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(Total_counts), dplyr::desc(Prevalence_samples))
}

calcular_alpha_funcional <- function(fun_counts, meta_df) {
  mat <- t(fun_counts) # amostras x funcoes
  if (!setequal(rownames(mat), meta_df$SampleID)) {
    stop("Matriz funcional e metadados possuem universos de amostras diferentes.", call. = FALSE)
  }
  idx <- match(rownames(mat), meta_df$SampleID)
  if (anyNA(idx)) stop("Falha ao alinhar metadados da diversidade funcional.", call. = FALSE)

  alpha <- data.frame(
    SampleID = rownames(mat),
    Observed_functions = rowSums(mat > 0),
    # Soma de atribuicoes funcionais; pode exceder o total de reads porque uma
    # ASV pode contribuir para mais de uma funcao FAPROTAX.
    Total_function_assignments = rowSums(mat),
    meta_df[idx, setdiff(colnames(meta_df), "SampleID"), drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  alpha
}

kruskal_por_grupo <- function(alpha_df, grupo, prefixo, out_dir) {
  metricas <- c("Observed_functions")
  if (!grupo %in% colnames(alpha_df)) {
    log_msg(paste(prefixo, ": grupo ausente no alpha funcional:", grupo), "WARN")
    return(NULL)
  }

  res <- lapply(metricas, function(met) {
    df <- alpha_df[!is.na(alpha_df[[grupo]]) & !is.na(alpha_df[[met]]), , drop = FALSE]
    if (length(unique(df[[grupo]])) < 2) return(NULL)

    n_grupo <- table(df[[grupo]])
    if (any(n_grupo < 2)) {
      log_msg(
        sprintf("%s/%s/%s: teste exato omitido; grupo(s) com n < 2.", prefixo, grupo, met),
        "WARN"
      )
      return(NULL)
    }

    esperado <- if (identical(grupo, "BeeSpecies")) {
      if (identical(prefixo, "core9")) c(2L, 3L, 4L) else c(3L, 3L, 4L)
    } else NULL
    teste <- tryCatch(
      teste_postos_exato_grupos_fixos(
        df[[met]], df[[grupo]], tamanhos_esperados = esperado,
        max_alocacoes = 1000000L
      ),
      error = function(e) e
    )
    if (inherits(teste, "error")) {
      log_msg(
        sprintf("%s/%s/%s: teste exato falhou: %s",
                prefixo, grupo, met, conditionMessage(teste)),
        "WARN"
      )
      salvar_erro_etapa(prefixo, paste0("teste_exato_", grupo, "_", met), teste, out_dir)
      return(NULL)
    }

    data.frame(
      Prefixo = prefixo,
      Grupo = grupo,
      Metrica = met,
      Statistic = teste$H,
      Epsilon2 = round(teste$Epsilon2, 3),
      P_value = teste$p_exato,
      P_assintotico_diagnostico = teste$p_assintotico_diagnostico,
      N_permutacoes_exatas = teste$N_permutacoes_exatas,
      p_min_teorico = teste$p_min_teorico,
      N_por_grupo = paste(names(n_grupo), n_grupo, sep = "=", collapse = "; "),
      Metodo = teste$Metodo,
      Independencia_fator = if (grupo == "Nativo_Introduzido") {
        "nao_independente_de_BeeSpecies"
      } else {
        "nao_estimavel_separadamente_sem_modelo_balanceado"
      },
      stringsAsFactors = FALSE
    )
  })

  res <- dplyr::bind_rows(res)
  if (nrow(res) > 0) {
    res$P_adj_BH <- p.adjust(res$P_value, method = "BH")
    salvar_csv(res, file.path(out_dir, paste0(prefixo, "_alpha_kruskal_", grupo, ".csv")))
  }
  res
}

preparar_matriz_funcional_para_beta <- function(fun_rel, meta_df, prefixo, grupo) {
  mat <- t(fun_rel) # amostras x funcoes

  faltam <- setdiff(meta_df$SampleID, rownames(mat))
  if (length(faltam) > 0) {
    stop(prefixo, "/", grupo, ": amostras dos metadados ausentes na matriz funcional: ",
         paste(faltam, collapse = ", "), call. = FALSE)
  }

  mat <- mat[meta_df$SampleID, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  if (ncol(mat) < 2) {
    stop(prefixo, "/", grupo, ": menos de 2 funcoes com abundancia para beta diversidade.", call. = FALSE)
  }

  linhas_zero <- rowSums(mat) == 0
  if (any(linhas_zero)) {
    log_msg(
      paste0(prefixo, "/", grupo, ": removendo ", sum(linhas_zero),
             " amostra(s) sem funcao mapeada antes da beta diversidade."),
      "WARN"
    )
    mat <- mat[!linhas_zero, , drop = FALSE]
    meta_df <- meta_df[match(rownames(mat), meta_df$SampleID), , drop = FALSE]
  }

  list(mat = mat, meta = meta_df)
}

permanova_funcional <- function(fun_rel, meta_df, grupo, prefixo, out_dir) {
  if (!grupo %in% colnames(meta_df)) {
    log_msg(paste(prefixo, ": grupo ausente na PERMANOVA funcional:", grupo), "WARN")
    return(NULL)
  }

  prep <- tryCatch(
    preparar_matriz_funcional_para_beta(fun_rel, meta_df, prefixo, grupo),
    error = function(e) e
  )

  if (inherits(prep, "error")) {
    salvar_erro_etapa(prefixo, paste0("permanova_bray_", grupo), prep, out_dir)
    return(NULL)
  }

  mat <- prep$mat
  meta_ord <- prep$meta

  if (nrow(mat) < 3 || length(unique(meta_ord[[grupo]])) < 2) {
    log_msg(sprintf("PERMANOVA funcional omitida (%s/%s): amostras ou grupos insuficientes.", prefixo, grupo), "WARN")
    return(NULL)
  }

  dist_bray <- vegan::vegdist(mat, method = "bray")
  meta_ord <- meta_ord[match(labels(dist_bray), meta_ord$SampleID), , drop = FALSE]
  rownames(meta_ord) <- meta_ord$SampleID

  if (!identical(rownames(meta_ord), labels(dist_bray))) {
    stop(prefixo, "/", grupo, ": distancia e metadados fora de ordem.", call. = FALSE)
  }

  meta_ord$.grupo_perm <- droplevels(as.factor(meta_ord[[grupo]]))
  n_grp <- table(meta_ord$.grupo_perm)

  if (length(n_grp) < 2) {
    log_msg(sprintf("PERMANOVA funcional omitida (%s/%s): menos de 2 grupos apos filtros.", prefixo, grupo), "WARN")
    return(NULL)
  }

  grupo_nomeado <- setNames(as.character(meta_ord$.grupo_perm), rownames(meta_ord))
  esperado <- if (identical(grupo, "BeeSpecies")) {
    if (identical(prefixo, "core9")) c(2L, 3L, 4L) else c(3L, 3L, 4L)
  } else NULL
  perm <- permanova_exata_grupos_fixos(
    dist_bray, grupo_nomeado, tamanhos_esperados = esperado,
    max_alocacoes = 100000L
  )
  if (!is.na(perm$Erro)) {
    log_msg(paste(prefixo, grupo, "PERMANOVA exata nao executada:", perm$Erro), "WARN")
  }
  p_min_alcancavel <- perm$p_min_teorico

  bd_p <- NA_real_
  if (all(n_grp >= 3)) {
    bd <- tryCatch({
      bd_obj  <- vegan::betadisper(dist_bray, meta_ord$.grupo_perm)
      vegan::permutest(bd_obj, permutations = N_PERM)
    }, error = function(e) e)
    if (!inherits(bd, "error")) {
      bd_p <- bd$tab$`Pr(>F)`[1]
    } else {
      salvar_erro_etapa(prefixo, paste0("betadisper_bray_", grupo), bd, out_dir)
    }
  } else {
    log_msg(sprintf("BETADISPER funcional omitido (%s/%s): grupo(s) com n < 3.", prefixo, grupo), "WARN")
  }

  res <- data.frame(
    Prefixo = prefixo,
    Grupo = grupo,
    Distancia = "Bray-Curtis funcional microeco",
    N_amostras = nrow(meta_ord),
    N_por_grupo = paste(names(n_grp), n_grp, sep = "=", collapse = "; "),
    PERMANOVA_Df = perm$df_model,
    PERMANOVA_R2 = perm$R2,
    PERMANOVA_F = perm$F,
    PERMANOVA_p = perm$p_exato,
    N_alocacoes_exatas = perm$N_alocacoes_exatas,
    BETADISPER_p = bd_p,
    p_min_alcancavel = p_min_alcancavel,
    Metodo = paste(perm$Metodo, "Potencial funcional inferido por taxonomia; exploratorio."),
    Independencia_fator = if (grupo == "Nativo_Introduzido") {
      "nao_independente_de_BeeSpecies"
    } else {
      "nao_estimavel_separadamente_sem_modelo_balanceado"
    },
    stringsAsFactors = FALSE
  )

  salvar_csv(res, file.path(out_dir, paste0(prefixo, "_permanova_bray_", grupo, ".csv")))
  salvar_rds(dist_bray, file.path(out_dir, paste0(prefixo, "_dist_bray_funcional_", grupo, ".rds")))

  res
}

###############################################################################
# 5. GRAFICOS
###############################################################################

gerar_paleta_funcoes <- function(funcoes) {
  funcoes <- unique(as.character(funcoes))
  n <- length(funcoes)

  cores <- if (n <= length(CORES_MULTI)) {
    CORES_MULTI[seq_len(n)]
  } else {
    grDevices::colorRampPalette(CORES_MULTI)(n)
  }

  setNames(cores, funcoes)
}

preparar_top_funcoes_long <- function(fun_counts, fun_rel, meta_df, top_n = TOP_N_FUNCOES) {
  resumo <- resumo_funcoes(fun_counts, fun_rel)
  if (nrow(resumo) == 0) stop("Nenhuma funcao para plotar.", call. = FALSE)

  top_funcoes <- head(resumo$Function, min(top_n, nrow(resumo)))

  df <- as.data.frame(t(fun_rel[top_funcoes, , drop = FALSE]), check.names = FALSE)
  df$SampleID <- rownames(df)

  df <- tidyr::pivot_longer(
    df,
    cols = tidyselect::all_of(top_funcoes),
    names_to = "Function",
    values_to = "Relative_abundance"
  )

  idx_meta <- match(df$SampleID, meta_df$SampleID)
  if (anyNA(idx_meta)) stop("Amostra funcional sem metadado correspondente.", call. = FALSE)
  meta_cols <- setdiff(colnames(meta_df), "SampleID")
  df <- cbind(df, meta_df[idx_meta, meta_cols, drop = FALSE], stringsAsFactors = FALSE)
  df$Function <- factor(df$Function, levels = top_funcoes)
  df$SampleID <- factor(df$SampleID, levels = meta_df$SampleID)
  df$SampleLabel <- factor(df$SampleLabel, levels = meta_df$SampleLabel)

  df
}

agregar_funcoes_por_grupo <- function(fun_rel, meta_df, grupo, top_n = TOP_N_FUNCOES) {
  if (!grupo %in% colnames(meta_df)) {
    stop("Grupo ausente no metadata: ", grupo, call. = FALSE)
  }

  resumo <- data.frame(
    Function = rownames(fun_rel),
    Mean_relative = rowMeans(fun_rel),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(Mean_relative))

  top_funcoes <- head(resumo$Function, min(top_n, nrow(resumo)))
  mat <- t(fun_rel[top_funcoes, , drop = FALSE])
  mat <- mat[meta_df$SampleID, , drop = FALSE]

  df <- as.data.frame(mat, check.names = FALSE)
  df$SampleID <- rownames(df)
  idx_meta <- match(df$SampleID, meta_df$SampleID)
  if (anyNA(idx_meta)) stop("Amostra funcional sem metadado correspondente.", call. = FALSE)
  df[[grupo]] <- meta_df[[grupo]][idx_meta]

  long <- tidyr::pivot_longer(
    df,
    cols = tidyselect::all_of(top_funcoes),
    names_to = "Function",
    values_to = "Relative_abundance"
  )

  out <- long |>
    dplyr::group_by(.data[[grupo]], Function) |>
    dplyr::summarise(
      Mean_relative = mean(Relative_abundance, na.rm = TRUE),
      SD_relative = stats::sd(Relative_abundance, na.rm = TRUE),
      N = dplyr::n(),
      .groups = "drop"
    )

  colnames(out)[colnames(out) == grupo] <- "Grupo"
  out$Function <- factor(out$Function, levels = top_funcoes)
  out
}

plot_cobertura_faprotax <- function(cobertura_df, meta_df, prefixo, plot_dir) {
  if (anyDuplicated(cobertura_df$SampleID)) {
    stop("Tabela de cobertura FAPROTAX contem SampleID duplicado.", call. = FALSE)
  }
  idx_meta <- match(cobertura_df$SampleID, meta_df$SampleID)
  if (anyNA(idx_meta)) stop("Cobertura funcional sem metadado correspondente.", call. = FALSE)
  meta_cols <- setdiff(colnames(meta_df), "SampleID")
  df <- cbind(
    cobertura_df,
    meta_df[idx_meta, meta_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )
  df$SampleID <- factor(df$SampleID, levels = meta_df$SampleID)
  df$SampleLabel <- factor(df$SampleLabel, levels = meta_df$SampleLabel)

  p <- ggplot(df, aes(x = SampleLabel, y = Pct_reads_mapeados_unicos / 100, fill = BeeSpecies)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = CORES_ESP, labels = LABELS_ESP, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    TEMA_PUB +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(
      title = paste0("Microeco — Cobertura funcional por amostra: ", prefixo),
      subtitle = "Percentual de reads pertencentes a ASVs mapeadas a pelo menos uma função metabólica",
      x = "Amostra",
      y = "Reads mapeados a funcoes metabólicas",
      fill = "Espécie"
    )

  salvar_plot_duplo(p, paste0(prefixo, "_G01_cobertura_faprotax"), plot_dir, width = 10, height = 6)
}

plot_top_funcoes_amostras <- function(fun_counts, fun_rel, meta_df, prefixo, plot_dir,
                                      top_n = TOP_N_FUNCOES,
                                      contrib_global = NULL,
                                      legenda_funcao_taxon = NULL,
                                      top_n_taxa_legenda = TOP_N_TAXA_LEGENDA) {
  df <- preparar_top_funcoes_long(fun_counts, fun_rel, meta_df, top_n)
  fun_levels <- levels(df$Function)

  if (is.null(legenda_funcao_taxon)) {
    legenda_funcao_taxon <- criar_legenda_funcao_taxon_df(
      contrib_global = contrib_global,
      funcoes = fun_levels,
      top_n_taxa = top_n_taxa_legenda,
      incluir_pct = TRUE
    )
  } else {
    legenda_funcao_taxon <- legenda_funcao_taxon[legenda_funcao_taxon$Function %in% fun_levels, , drop = FALSE]
  }

  df <- aplicar_codigo_funcao(df, legenda_funcao_taxon)
  df$Function_codigo <- factor(df$Function_codigo, levels = legenda_funcao_taxon$Funcao_codigo)
  pal_fun <- gerar_paleta_funcoes(levels(df$Function_codigo))

  p <- ggplot(df, aes(x = SampleLabel, y = Relative_abundance, fill = Function_codigo)) +
    geom_col(position = "fill", width = 0.85) +
    facet_grid(. ~ BeeSpecies, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_fun, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)) +
    TEMA_PUB +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.text = element_text(size = 8)
    ) +
    labs(
      title = paste0("Microeco — Top ", min(top_n, length(unique(df$Function))), " funções por amostra: ", prefixo),
      subtitle = "Barras normalizadas; use a legenda auxiliar G00 para associar cada código à função e aos microrganismos",
      x = "Amostra",
      y = "Abundância relativa normalizada no gráfico",
      fill = "Código da função"
    )

  salvar_plot_duplo(p, paste0(prefixo, "_G02_top_funcoes_amostras"), plot_dir, width = 12, height = 7)
}

plot_top_funcoes_grupo <- function(fun_rel, meta_df, grupo, prefixo, plot_dir,
                                  top_n = TOP_N_FUNCOES,
                                  contrib_global = NULL,
                                  legenda_funcao_taxon = NULL,
                                  top_n_taxa_legenda = TOP_N_TAXA_LEGENDA) {
  df <- agregar_funcoes_por_grupo(fun_rel, meta_df, grupo, top_n)
  fun_levels <- levels(df$Function)

  if (is.null(legenda_funcao_taxon)) {
    legenda_funcao_taxon <- criar_legenda_funcao_taxon_df(
      contrib_global = contrib_global,
      funcoes = fun_levels,
      top_n_taxa = top_n_taxa_legenda,
      incluir_pct = TRUE
    )
  } else {
    legenda_funcao_taxon <- legenda_funcao_taxon[legenda_funcao_taxon$Function %in% fun_levels, , drop = FALSE]
  }

  df <- aplicar_codigo_funcao(df, legenda_funcao_taxon)
  df$Function_codigo <- factor(df$Function_codigo, levels = legenda_funcao_taxon$Funcao_codigo)
  pal_fun <- gerar_paleta_funcoes(levels(df$Function_codigo))

  p <- ggplot(df, aes(x = Grupo, y = Mean_relative, fill = Function_codigo)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.15, position = "fill") +
    scale_fill_manual(values = pal_fun, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)) +
    TEMA_PUB +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.text = element_text(size = 8)
    ) +
    labs(
      title = paste0("Microeco — Composição funcional média por ", grupo, ": ", prefixo),
      subtitle = "Média entre amostras; use a legenda auxiliar G00 para associar cada código à função e aos microrganismos",
      x = grupo,
      y = "Abundância relativa média normalizada no gráfico",
      fill = "Código da função"
    )

  nome <- paste0(prefixo, "_G03_top_funcoes_media_", grupo)
  salvar_plot_duplo(p, nome, plot_dir, width = 12, height = 7)
}


formatar_pct_heatmap <- function(x) {
  # x deve estar em proporção (0.001 = 0,1%).
  # Mantém zeros reais como "0%" e evita que valores positivos pequenos apareçam como 0%.
  x <- as.numeric(x)
  pct <- 100 * x
  out <- rep("", length(x))

  idx_na <- is.na(x) | !is.finite(x)
  out[idx_na] <- ""

  idx_zero <- !idx_na & x == 0
  out[idx_zero] <- "0%"

  idx_pos <- !idx_na & x > 0

  idx_muito_baixo <- idx_pos & pct < 0.001
  idx_baixo       <- idx_pos & pct >= 0.001 & pct < 0.01
  idx_medio_baixo <- idx_pos & pct >= 0.01  & pct < 0.1
  idx_medio       <- idx_pos & pct >= 0.1   & pct < 1
  idx_alto        <- idx_pos & pct >= 1

  out[idx_muito_baixo] <- "<0,001%"
  out[idx_baixo]       <- paste0(formatC(pct[idx_baixo],       format = "f", digits = 3, decimal.mark = ","), "%")
  out[idx_medio_baixo] <- paste0(formatC(pct[idx_medio_baixo], format = "f", digits = 2, decimal.mark = ","), "%")
  out[idx_medio]       <- paste0(formatC(pct[idx_medio],       format = "f", digits = 2, decimal.mark = ","), "%")
  out[idx_alto]        <- paste0(formatC(pct[idx_alto],        format = "f", digits = 1, decimal.mark = ","), "%")

  out
}

plot_heatmap_funcoes <- function(fun_counts, fun_rel, meta_df, prefixo, plot_dir,
                                  top_n = TOP_N_FUNCOES,
                                  contrib_global = NULL,
                                  legenda_funcao_taxon = NULL,
                                  top_n_taxa_legenda = TOP_N_TAXA_LEGENDA) {
  df <- preparar_top_funcoes_long(fun_counts, fun_rel, meta_df, top_n)
  fun_levels <- levels(df$Function)

  if (is.null(legenda_funcao_taxon)) {
    legenda_funcao_taxon <- criar_legenda_funcao_taxon_df(
      contrib_global = contrib_global,
      funcoes = fun_levels,
      top_n_taxa = top_n_taxa_legenda,
      incluir_pct = TRUE
    )
  } else {
    legenda_funcao_taxon <- legenda_funcao_taxon[
      legenda_funcao_taxon$Function %in% fun_levels,
      ,
      drop = FALSE
    ]
  }

  # Garante que a ordem dos códigos siga exatamente a ordem das funções selecionadas.
  legenda_funcao_taxon <- legenda_funcao_taxon[
    match(fun_levels, legenda_funcao_taxon$Function),
    ,
    drop = FALSE
  ]

  if (any(is.na(legenda_funcao_taxon$Function))) {
    legenda_funcao_taxon <- criar_legenda_funcao_taxon_df(
      contrib_global = contrib_global,
      funcoes = fun_levels,
      top_n_taxa = top_n_taxa_legenda,
      incluir_pct = TRUE
    )
  }

  df <- aplicar_codigo_funcao(df, legenda_funcao_taxon)

  eixo_y_df <- legenda_funcao_taxon |>
    dplyr::mutate(
      Function_y = paste0(.data$Funcao_codigo, " — ", .data$Function)
    )

  mapa_y <- setNames(eixo_y_df$Function_y, eixo_y_df$Function)
  df$Function_y <- mapa_y[as.character(df$Function)]
  df$Function_y[is.na(df$Function_y)] <- as.character(df$Function)[is.na(df$Function_y)]
  df$Function_y <- factor(df$Function_y, levels = rev(eixo_y_df$Function_y))

  # Rótulos com precisão adaptativa: evita que valores positivos pequenos apareçam como 0%.
  df$Label_pct <- formatar_pct_heatmap(df$Relative_abundance)

  # Cor do texto por contraste. Usa o percentil 75 dos valores positivos para evitar texto branco em células pouco intensas.
  positivos <- df$Relative_abundance[df$Relative_abundance > 0]
  limiar <- stats::quantile(positivos, probs = 0.75, na.rm = TRUE, names = FALSE)
  if (!is.finite(limiar)) limiar <- Inf
  df$Texto_cor <- ifelse(df$Relative_abundance >= limiar, "white", "black")

  p <- ggplot(df, aes(x = SampleLabel, y = Function_y, fill = Relative_abundance)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(
      aes(label = Label_pct, color = Texto_cor),
      size = 2.05,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_color_identity() +
    facet_grid(. ~ BeeSpecies, scales = "free_x", space = "free_x") +
    scale_y_discrete(labels = function(x) quebrar_rotulo(x, largura = 42)) +
    scale_fill_gradient(
      low = "#21EBC9",
      high = "#000489",
      trans = "sqrt",
      labels = scales::label_percent(accuracy = 0.001, decimal.mark = ","),
      na.value = "grey95"
    ) +
    TEMA_PUB +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = paste0("Microeco — Heatmap das funções mais abundantes: ", prefixo),
      subtitle = "Rótulos mostram o % real por amostra; Use G00 para microrganismos associados.",
      x = "Amostra",
      y = "Código — função metabólica",
      fill = "Abundância relativa"
    )

  salvar_plot_duplo(p, paste0(prefixo, "_G04_heatmap_funcoes"), plot_dir, width = 13, height = 9)
}

plot_alpha_funcional <- function(alpha_df, grupo, prefixo, plot_dir) {
  if (!grupo %in% colnames(alpha_df)) return(invisible(NULL))

  metricas <- c("Observed_functions")
  df <- alpha_df[, c("SampleID", grupo, metricas), drop = FALSE]

  long <- tidyr::pivot_longer(
    df,
    cols = tidyselect::all_of(metricas),
    names_to = "Metrica",
    values_to = "Valor"
  )

  long$Metrica <- factor(
    long$Metrica,
    levels = metricas,
    labels = c("Funções observadas")
  )

  fill_scale <- if (grupo == "BeeSpecies") {
    scale_fill_manual(values = CORES_ESP, labels = LABELS_ESP, drop = FALSE)
  } else {
    scale_fill_manual(values = CORES_STATUS, labels = LABELS_STATUS, drop = FALSE)
  }

  p <- ggplot(long, aes(x = .data[[grupo]], y = Valor, fill = .data[[grupo]])) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, size = 2.2, shape = 21, color = "black", alpha = 0.9) +
    facet_wrap(~ Metrica, scales = "free_y") +
    fill_scale +
    TEMA_PUB +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none") +
    labs(
      title = paste0("Microeco: Diversidade alfa funcional por ", grupo, ": ", prefixo),
      subtitle = "Riqueza funcional; funções metabólicas não são mutuamente exclusivas",
      x = grupo,
      y = "Funções observadas"
    )

  nome <- paste0(prefixo, "_G05_alpha_funcional_", grupo)
  salvar_plot_duplo(p, nome, plot_dir, width = 10, height = 6)
}

plot_pcoa_funcional <- function(fun_rel, meta_df, grupo, prefixo, plot_dir) {
  if (!grupo %in% colnames(meta_df)) return(invisible(NULL))

  prep <- tryCatch(
    preparar_matriz_funcional_para_beta(fun_rel, meta_df, prefixo, grupo),
    error = function(e) e
  )
  if (inherits(prep, "error")) stop(conditionMessage(prep), call. = FALSE)

  mat <- prep$mat
  meta_use <- prep$meta

  if (nrow(mat) < 3 || sum(rowSums(mat) > 0) < 3) {
    log_msg(paste("PCoA funcional omitida para", prefixo, grupo, "por amostras insuficientes."), "WARN")
    return(invisible(NULL))
  }

  meta_use <- validar_samplelabel_meta(meta_use)

  dist_bray <- vegan::vegdist(mat, method = "bray")
  k_use <- min(3L, nrow(mat) - 1L)
  ord <- stats::cmdscale(dist_bray, eig = TRUE, k = k_use)

  pontos <- as.data.frame(ord$points, stringsAsFactors = FALSE)
  colnames(pontos) <- paste0("PC", seq_len(ncol(pontos)))
  pontos$SampleID <- rownames(pontos)

  df <- dplyr::left_join(
    pontos,
    meta_use,
    by = "SampleID"
  )

  if (any(is.na(df[[grupo]]))) {
    stop(prefixo, " / ", grupo, ": PCoA gerou amostra(s) sem grupo no metadata.", call. = FALSE)
  }

  eig <- ord$eig
  soma_pos <- sum(eig[eig > 0])
  pct <- if (soma_pos > 0) round(100 * eig / soma_pos, 1) else rep(NA_real_, length(eig))

  eixo_lab <- function(i) {
    if (length(pct) < i || is.na(pct[i])) {
      paste0("PCoA", i)
    } else {
      sprintf("PCoA%d [%.1f%%]", i, pct[i])
    }
  }

  config_grupo <- function(grupo, valores) {
    valores <- unique(as.character(valores))

    if (grupo == "BeeSpecies") {
      cores <- CORES_ESP
      labels <- LABELS_ESP
      shapes <- c(
        "Melipona fasciculata" = 16,
        "Melipona scutellaris" = 17,
        "Melipona subnitida"   = 15
      )
      nome <- "Espécie"
    } else if (grupo == "Nativo_Introduzido") {
      cores <- CORES_STATUS
      labels <- LABELS_STATUS
      shapes <- c("Introduzida" = 17, "Nativa" = 16)
      nome <- "Origem"
    } else {
      cores <- setNames(grDevices::rainbow(length(valores)), valores)
      labels <- setNames(valores, valores)
      shapes <- setNames(rep(c(16, 17, 15, 18, 3), length.out = length(valores)), valores)
      nome <- grupo
    }

    extras <- setdiff(valores, names(cores))
    if (length(extras) > 0) {
      cores <- c(cores, setNames(grDevices::rainbow(length(extras)), extras))
      labels <- c(labels, setNames(extras, extras))
      shapes <- c(shapes, setNames(rep(19, length(extras)), extras))
    }

    list(cores = cores, labels = labels, shapes = shapes, nome = nome)
  }

  cfg <- config_grupo(grupo, df[[grupo]])

  calcular_hull <- function(df_plot, eixo_x, eixo_y) {
    df_plot |>
      dplyr::group_by(.data[[grupo]]) |>
      dplyr::filter(
        dplyr::n() >= 3,
        dplyr::n_distinct(.data[[eixo_x]]) >= 2,
        dplyr::n_distinct(.data[[eixo_y]]) >= 2
      ) |>
      dplyr::slice(grDevices::chull(.data[[eixo_x]], .data[[eixo_y]])) |>
      dplyr::ungroup()
  }

  calcular_perm_label <- function(dist_obj, meta_df, grupo, prefixo) {
    meta_ord <- meta_df[labels(dist_obj), , drop = FALSE]
    if (!identical(rownames(meta_ord), labels(dist_obj))) {
      return("PERMANOVA funcional indisponível: metadados fora de ordem.")
    }

    grupo_perm <- droplevels(as.factor(meta_ord[[grupo]]))
    n_grp <- table(grupo_perm)

    if (length(n_grp) < 2 || any(n_grp < 2)) {
      return("PERMANOVA funcional omitida: grupo(s) com n < 2.")
    }

    names(grupo_perm) <- rownames(meta_ord)
    esperado <- if (identical(grupo, "BeeSpecies")) {
      if (identical(prefixo, "core9")) c(2L, 3L, 4L) else c(3L, 3L, 4L)
    } else NULL
    perm <- permanova_exata_grupos_fixos(
      dist_obj, grupo_perm, tamanhos_esperados = esperado,
      max_alocacoes = 100000L
    )
    if (!is.na(perm$Erro)) {
      return(paste("PERMANOVA funcional indisponível:", perm$Erro))
    }
    sprintf(
      "PERMANOVA funcional exata: F = %.3f | R² = %.3f | %s%s",
      perm$F, perm$R2, fmt_p(perm$p_exato),
      ifelse(any(n_grp < 3), " | exploratória: grupo(s) com n < 3", "")
    )
  }

  perm_label <- calcular_perm_label(dist_bray, meta_use, grupo, prefixo)

  montar_pcoa_padrao_script10 <- function(df_plot, eixo_x, eixo_y, titulo, nome_base) {
    centroides <- df_plot |>
      dplyr::group_by(.data[[grupo]]) |>
      dplyr::summarise(
        cx = mean(.data[[eixo_x]], na.rm = TRUE),
        cy = mean(.data[[eixo_y]], na.rm = TRUE),
        .groups = "drop"
      )

    segmentos <- dplyr::left_join(df_plot, centroides, by = grupo)
    hulls <- calcular_hull(df_plot, eixo_x, eixo_y)

    p <- ggplot(df_plot, aes(x = .data[[eixo_x]], y = .data[[eixo_y]],
                             color = .data[[grupo]], shape = .data[[grupo]])) +
      geom_polygon(
        data = hulls,
        aes(x = .data[[eixo_x]], y = .data[[eixo_y]],
            fill = .data[[grupo]], group = .data[[grupo]]),
        inherit.aes = FALSE,
        alpha = 0.12, color = NA, show.legend = FALSE
      ) +
      geom_segment(
        data = segmentos,
        aes(x = cx, y = cy, xend = .data[[eixo_x]], yend = .data[[eixo_y]],
            color = .data[[grupo]]),
        linewidth = 0.6, alpha = 0.55, show.legend = FALSE
      ) +
      geom_point(size = 5, alpha = 0.95, stroke = 1.2) +
      geom_point(
        data = centroides,
        aes(x = cx, y = cy, color = .data[[grupo]]),
        inherit.aes = FALSE,
        size = 3.5, shape = 3, stroke = 1.8, show.legend = FALSE
      ) +
      ggrepel::geom_text_repel(
        aes(label = SampleLabel),
        size = 3, fontface = "bold",
        box.padding = 0.4, point.padding = 0.3,
        segment.color = "grey60", segment.size = 0.4,
        show.legend = FALSE
      ) +
      scale_color_manual(values = cfg$cores, labels = cfg$labels,
                         name = cfg$nome, drop = FALSE) +
      scale_fill_manual(values = cfg$cores, labels = cfg$labels,
                        name = cfg$nome, drop = FALSE) +
      scale_shape_manual(values = cfg$shapes, labels = cfg$labels,
                         name = cfg$nome, drop = FALSE) +
      TEMA_PUB +
      theme(legend.text = element_text(face = if (grupo == "BeeSpecies") "italic" else "plain")) +
      labs(
        title = titulo,
        subtitle = perm_label,
        x = eixo_lab(as.integer(sub("PC", "", eixo_x))),
        y = eixo_lab(as.integer(sub("PC", "", eixo_y)))
      )

    salvar_plot_duplo(p, nome_base, plot_dir, width = 10, height = 7)
    invisible(p)
  }

  nome_grupo <- ifelse(grupo == "BeeSpecies", "por espécie", "por origem")

  g12 <- montar_pcoa_padrao_script10(
    df,
    eixo_x = "PC1",
    eixo_y = "PC2",
    titulo = paste0("FAPROTAX/microeco — PCoA funcional Bray-Curtis ", nome_grupo, ": ", prefixo),
    nome_base = paste0(prefixo, "_G06_pcoa_funcional_", grupo)
  )

  g13 <- NULL
  if ("PC3" %in% colnames(df)) {
    g13 <- montar_pcoa_padrao_script10(
      df,
      eixo_x = "PC1",
      eixo_y = "PC3",
      titulo = paste0("FAPROTAX/microeco — PCoA funcional Bray-Curtis PC1 x PC3 ", nome_grupo, ": ", prefixo),
      nome_base = paste0(prefixo, "_G06b_pcoa_funcional_", grupo, "_PC1xPC3")
    )
  } else {
    log_msg(paste(prefixo, grupo, ": PC3 ausente; PCoA PC1 x PC3 omitida."), "WARN")
  }

  invisible(list(PC1xPC2 = g12, PC1xPC3 = g13))
}

gerar_graficos_faprotax <- function(fun_counts, fun_rel, cobertura_df, alpha_df, meta_df, prefixo, plot_dir,
                                    contrib_global = NULL,
                                    legenda_funcao_taxon = NULL) {
  if (isTRUE(USAR_LEGENDA_AUXILIAR_TAXON) && !is.null(legenda_funcao_taxon)) {
    executar_plot_seguro(
      plot_legenda_funcao_taxon_auxiliar(legenda_funcao_taxon, prefixo, plot_dir),
      prefixo,
      "G00_legenda_auxiliar_funcao_taxon"
    )
  }

  executar_plot_seguro(plot_cobertura_faprotax(cobertura_df, meta_df, prefixo, plot_dir), prefixo, "G01_cobertura")
  executar_plot_seguro(plot_top_funcoes_amostras(fun_counts, fun_rel, meta_df, prefixo, plot_dir,
                                                 contrib_global = contrib_global,
                                                 legenda_funcao_taxon = legenda_funcao_taxon), prefixo, "G02_top_funcoes_amostras")
  executar_plot_seguro(plot_top_funcoes_grupo(fun_rel, meta_df, "BeeSpecies", prefixo, plot_dir,
                                              contrib_global = contrib_global,
                                              legenda_funcao_taxon = legenda_funcao_taxon), prefixo, "G03_BeeSpecies")
  executar_plot_seguro(plot_top_funcoes_grupo(fun_rel, meta_df, "Nativo_Introduzido", prefixo, plot_dir,
                                              contrib_global = contrib_global,
                                              legenda_funcao_taxon = legenda_funcao_taxon), prefixo, "G03_Nativo_Introduzido")
  executar_plot_seguro(plot_heatmap_funcoes(fun_counts, fun_rel, meta_df, prefixo, plot_dir,
                                            contrib_global = contrib_global,
                                            legenda_funcao_taxon = legenda_funcao_taxon), prefixo, "G04_heatmap")
  executar_plot_seguro(plot_alpha_funcional(alpha_df, "BeeSpecies", prefixo, plot_dir), prefixo, "G05_alpha_BeeSpecies")
  executar_plot_seguro(plot_alpha_funcional(alpha_df, "Nativo_Introduzido", prefixo, plot_dir), prefixo, "G05_alpha_Nativo_Introduzido")
  executar_plot_seguro(plot_pcoa_funcional(fun_rel, meta_df, "BeeSpecies", prefixo, plot_dir), prefixo, "G06_pcoa_BeeSpecies")
  executar_plot_seguro(plot_pcoa_funcional(fun_rel, meta_df, "Nativo_Introduzido", prefixo, plot_dir), prefixo, "G06_pcoa_Nativo_Introduzido")

  invisible(TRUE)
}

###############################################################################
# 6. PIPELINE POR CONJUNTO
###############################################################################

rodar_faprotax_para_phyloseq <- function(ps_obj, prefixo, inferencia = FALSE) {
  log_msg(paste("Iniciando FAPROTAX/microeco para", prefixo), "RUN")

  validar_phyloseq_entrada(ps_obj, prefixo)

  otu_mat <- matriz_taxa_linhas(ps_obj)

  meta_df <- as(phyloseq::sample_data(ps_obj), "data.frame")
  faltam_meta <- setdiff(colnames(otu_mat), rownames(meta_df))
  extras_meta <- setdiff(rownames(meta_df), colnames(otu_mat))
  if (length(faltam_meta) > 0L || length(extras_meta) > 0L) {
    stop(
      prefixo, ": otu_table e sample_data divergem: faltantes=",
      length(faltam_meta), "; extras=", length(extras_meta),
      call. = FALSE
    )
  }
  meta_df <- meta_df[colnames(otu_mat), , drop = FALSE]
  if (!identical(rownames(meta_df), colnames(otu_mat))) {
    stop(prefixo, ": ordem do sample_data diverge da otu_table.", call. = FALSE)
  }
  meta_df$SampleID <- rownames(meta_df)
  meta_df <- validar_samplelabel_meta(meta_df)
  meta_df <- adicionar_status_origem(meta_df)

  res_microeco <- rodar_microeco_faprotax(ps_obj, prefixo)
  asv_func_bin <- res_microeco$asv_func_bin

  abund <- calcular_abundancia_funcional_microeco(otu_mat, asv_func_bin, prefixo)
  alpha <- calcular_alpha_funcional(abund$counts, meta_df)
  resumo <- resumo_funcoes(abund$counts, abund$relative)
  asv_funcoes <- asv_funcoes_long(asv_func_bin)
  asv_funcoes_tax <- criar_tabela_asv_func_taxonomia(asv_funcoes, ps_obj, otu_mat)
  contrib_taxa_funcoes <- calcular_contribuicao_taxa_funcoes(asv_funcoes_tax, otu_mat, meta_df)

  top_funcoes_legenda <- head(resumo$Function, min(TOP_N_FUNCOES, nrow(resumo)))
  legenda_funcao_taxon <- criar_legenda_funcao_taxon_df(
    contrib_global = contrib_taxa_funcoes$global,
    funcoes = top_funcoes_legenda,
    top_n_taxa = TOP_N_TAXA_LEGENDA,
    incluir_pct = TRUE
  )

  asv_mapeadas <- unique(asv_funcoes$ASV_ID)
  asv_nao_mapeadas <- setdiff(rownames(otu_mat), asv_mapeadas)
  nao_mapeadas_df <- data.frame(
    ASV_ID = asv_nao_mapeadas,
    Total_reads = rowSums(otu_mat[asv_nao_mapeadas, , drop = FALSE]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(Total_reads))

  checkpoints <- data.frame(
    Prefixo = prefixo,
    Amostras = ncol(otu_mat),
    ASVs_total = nrow(otu_mat),
    Funcoes_detectadas = nrow(abund$counts),
    ASVs_mapeadas = length(asv_mapeadas),
    ASVs_nao_mapeadas = length(asv_nao_mapeadas),
    Media_pct_reads_mapeados = round(mean(abund$cobertura$Pct_reads_mapeados_unicos), 4),
    Min_pct_reads_mapeados = round(min(abund$cobertura$Pct_reads_mapeados_unicos), 4),
    Max_pct_reads_mapeados = round(max(abund$cobertura$Pct_reads_mapeados_unicos), 4),
    microeco_for_what = valor_escalar_seguro(res_microeco$trans_func$for_what),
    microeco_database = valor_escalar_seguro(
      res_microeco$trans_func$database,
      padrao = "FAPROTAX"
    ),
    Taxonomia_origem = descricao_hierarquia_taxa,
    Interpretacao_origem = OBS_CONFUNDIMENTO,
    stringsAsFactors = FALSE
  )

  if (checkpoints$Media_pct_reads_mapeados < 10) {
    log_msg(
      sprintf(
        "%s: cobertura media FAPROTAX baixa (%.2f%%); interpretar perfis funcionais com cautela.",
        prefixo, checkpoints$Media_pct_reads_mapeados
      ),
      "WARN"
    )
  }

  # Saidas principais do conjunto analisado.
  salvar_matriz(abund$counts,   file.path(faprotax_out, paste0(prefixo, "_faprotax_counts.csv")))
  salvar_matriz(abund$relative, file.path(faprotax_out, paste0(prefixo, "_faprotax_relative.csv")))

  salvar_rds(abund$counts,   file.path(faprotax_out, paste0(prefixo, "_faprotax_counts.rds")))
  salvar_rds(abund$relative, file.path(faprotax_out, paste0(prefixo, "_faprotax_relative.rds")))

  salvar_csv(asv_funcoes,     file.path(faprotax_out, paste0(prefixo, "_asv_faprotax_funcoes.csv")))
  salvar_csv(asv_funcoes_tax, file.path(faprotax_out, paste0(prefixo, "_asv_faprotax_funcoes_taxonomia.csv")))
  salvar_csv(contrib_taxa_funcoes$global,
             file.path(faprotax_out, paste0(prefixo, "_funcao_taxa_contribuicao_global.csv")))
  salvar_csv(contrib_taxa_funcoes$por_amostra,
             file.path(faprotax_out, paste0(prefixo, "_funcao_taxa_contribuicao_por_amostra.csv")))
  salvar_csv(legenda_funcao_taxon,
             file.path(faprotax_out, paste0(prefixo, "_legenda_auxiliar_funcao_taxon.csv")))
  salvar_csv(nao_mapeadas_df, file.path(faprotax_out, paste0(prefixo, "_asvs_nao_mapeadas_faprotax.csv")))
  salvar_csv(abund$cobertura, file.path(faprotax_out, paste0(prefixo, "_cobertura_faprotax.csv")))
  salvar_csv(alpha,           file.path(faprotax_out, paste0(prefixo, "_alpha_funcional_faprotax.csv")))
  salvar_csv(resumo,          file.path(faprotax_out, paste0(prefixo, "_resumo_funcoes_faprotax.csv")))
  salvar_csv(checkpoints,     file.path(faprotax_out, paste0(prefixo, "_checkpoints_faprotax_microeco.csv")))

  # Saidas adicionais de rastreabilidade microeco.
  salvar_matriz(as.matrix(asv_func_bin), file.path(faprotax_out, paste0(prefixo, "_microeco_res_func_binary.csv")))
  salvar_rds(res_microeco$trans_func, file.path(faprotax_out, paste0(prefixo, "_microeco_trans_func.rds")))

  if (!is.null(res_microeco$microeco_fr)) {
    ids_fr <- rownames(res_microeco$microeco_fr)
    if (is.null(ids_fr) || any(ids_fr == "")) {
      ids_fr <- as.character(seq_len(nrow(res_microeco$microeco_fr)))
      log_msg(paste(prefixo, ": res_func_FR sem rownames; IDs sequenciais usados apenas na auditoria."), "WARN")
    }
    microeco_fr_out <- data.frame(
      Row_ID = ids_fr,
      as.data.frame(res_microeco$microeco_fr, check.names = FALSE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    salvar_csv(microeco_fr_out, file.path(faprotax_out, paste0(prefixo, "_microeco_res_func_FR.csv")))
  }

  gerar_graficos_faprotax(abund$counts, abund$relative, abund$cobertura, alpha, meta_df, prefixo, faprotax_plot,
                         contrib_global = contrib_taxa_funcoes$global,
                         legenda_funcao_taxon = legenda_funcao_taxon)

  if (inferencia) {
    if (identical(prefixo, "plus10")) {
      log_msg(
        paste0(
          "plus10 inclui uma unica amostra da corrida run_aux. ",
          "Os testes funcionais sao sensibilidades exploratorias; Run nao pode ",
          "ser estimada separadamente de BeeSpecies."
        ),
        "WARN"
      )
    }
    log_msg(
      paste(
        prefixo,
        "com inferencia exploratoria; Nativo_Introduzido e confundido com BeeSpecies"
      ),
      "STAT"
    )
    kruskal_por_grupo(alpha, "BeeSpecies", prefixo, faprotax_out)
    kruskal_por_grupo(alpha, "Nativo_Introduzido", prefixo, faprotax_out)
    permanova_funcional(abund$relative, meta_df, "BeeSpecies", prefixo, faprotax_out)
    permanova_funcional(abund$relative, meta_df, "Nativo_Introduzido", prefixo, faprotax_out)
  } else {
    log_msg(paste(prefixo, "executado sem inferencia por configuracao."), "INFO")
  }

  log_msg(
    sprintf("%s: %d funcoes | %d/%d ASVs mapeadas | cobertura media %.2f%%",
            prefixo, nrow(abund$counts), length(asv_mapeadas), nrow(otu_mat),
            mean(abund$cobertura$Pct_reads_mapeados_unicos)),
    "OK"
  )

  list(
    counts = abund$counts,
    relative = abund$relative,
    asv_funcoes = asv_funcoes,
    asv_funcoes_taxonomia = asv_funcoes_tax,
    funcao_taxa_global = contrib_taxa_funcoes$global,
    funcao_taxa_por_amostra = contrib_taxa_funcoes$por_amostra,
    legenda_funcao_taxon = legenda_funcao_taxon,
    cobertura = abund$cobertura,
    alpha = alpha,
    resumo = resumo,
    checkpoints = checkpoints,
    microeco_trans_func = res_microeco$trans_func
  )
}

###############################################################################
# 7. EXECUCAO
###############################################################################

cat("=============================================================\n")
cat("FAPROTAX 16S — Inferencia funcional ecologica v", VERSAO, "\n", sep = "")
cat("Data/Hora:", DATA_EXECUCAO, "\n")
cat("Motor: microeco::trans_func + prok_database='FAPROTAX'\n")
cat("=============================================================\n\n")

cat("=== VALIDACOES ===\n\n")
validar_arquivo(arq_core9,  "phyloseq_core9_primeira_run.rds")
validar_arquivo(arq_plus10, "phyloseq_plus10_com_auxiliar.rds")
validar_arquivo(arq_asvmap, "ASV_sequences.tsv")
validar_arquivo(arq_meta_taxa, "metadata_consenso_taxonomico.csv")

meta_taxa_exec <- read.csv(
  arq_meta_taxa,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (nrow(meta_taxa_exec) < 1L ||
    !"Descricao_hierarquia" %in% colnames(meta_taxa_exec)) {
  stop(
    "metadata_consenso_taxonomico.csv invalido; execute novamente o Script 05.",
    call. = FALSE
  )
}
descricao_hierarquia_taxa <- as.character(
  meta_taxa_exec$Descricao_hierarquia[[1L]]
)
if (is.na(descricao_hierarquia_taxa) || !nzchar(descricao_hierarquia_taxa)) {
  stop("Descricao_hierarquia vazia no metadado taxonomico.", call. = FALSE)
}

ps_core9 <- tryCatch(
  readRDS(arq_core9),
  error = function(e) {
    stop("Falha ao ler ps_core9: ", conditionMessage(e), call. = FALSE)
  }
)
ps_plus10 <- tryCatch(
  readRDS(arq_plus10),
  error = function(e) {
    stop("Falha ao ler ps_plus10: ", conditionMessage(e), call. = FALSE)
  }
)

validar_phyloseq_entrada(ps_core9, "ps_core9")
validar_phyloseq_entrada(ps_plus10, "ps_plus10")

if (phyloseq::nsamples(ps_core9) != 9L || phyloseq::nsamples(ps_plus10) != 10L) {
  stop("Esperados core9=9 e plus10=10 amostras.", call. = FALSE)
}
meta_core_check <- as(phyloseq::sample_data(ps_core9), "data.frame")
meta_plus_check <- as(phyloseq::sample_data(ps_plus10), "data.frame")
if (any(as.character(meta_core_check$Run) != "run_main")) {
  stop("ps_core9 deve conter apenas run_main.", call. = FALSE)
}
if (sum(as.character(meta_plus_check$Run) == "run_main") != 9L ||
    sum(as.character(meta_plus_check$Run) == "run_aux") != 1L) {
  stop("ps_plus10 deve conter 9 amostras run_main e 1 run_aux.", call. = FALSE)
}
if (!"S10" %in% rownames(meta_plus_check) ||
    as.character(meta_plus_check["S10", "Run"]) != "run_aux") {
  stop("A amostra auxiliar S10 nao esta corretamente identificada.", call. = FALSE)
}
if (!setequal(phyloseq::sample_names(ps_core9),
              setdiff(phyloseq::sample_names(ps_plus10), "S10"))) {
  stop("Universo de amostras core9 nao corresponde ao plus10 sem Auxiliar.", call. = FALSE)
}
if (!all(phyloseq::taxa_names(ps_core9) %in% phyloseq::taxa_names(ps_plus10))) {
  stop("ASVs do core9 nao sao subconjunto do plus10.", call. = FALSE)
}

asv_map <- validar_mapa_asv(arq_asvmap, taxa_names(ps_plus10))

cat("\n=== EXECUTANDO FAPROTAX VIA MICROECO ===\n\n")

res_core9 <- rodar_faprotax_para_phyloseq(
  ps_obj = ps_core9,
  prefixo = "core9",
  inferencia = TRUE
)

res_plus10 <- rodar_faprotax_para_phyloseq(
  ps_obj = ps_plus10,
  prefixo = "plus10",
  inferencia = TRUE
)

###############################################################################
# 8. METADADOS DE EXECUCAO E README LOCAL
###############################################################################

metadata_exec <- data.frame(
  Script = "11_faprotax",
  Versao = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Base_path = base_path,
  FAPROTAX_engine = "microeco::trans_func",
  FAPROTAX_database = "FAPROTAX",
  Res_func_orientacao = "detectada_por_match_de_ASV; transposta_quando_necessario",
  Res_func_tratamento = "presenca_ausencia_binaria",
  Abundancia_funcional_denominador = "reads_totais_por_amostra; funcoes_nao_mutuamente_exclusivas",
  Taxonomia_source = "taxa_consenso_final.rds",
  Taxonomia_hierarquia = descricao_hierarquia_taxa,
  microeco_version = as.character(utils::packageVersion("microeco")),
  Core9_amostras = phyloseq::nsamples(ps_core9),
  Core9_ASVs = phyloseq::ntaxa(ps_core9),
  Core9_funcoes = nrow(res_core9$counts),
  Core9_ASVs_mapeadas = res_core9$checkpoints$ASVs_mapeadas,
  Core9_media_pct_reads_mapeados = res_core9$checkpoints$Media_pct_reads_mapeados,
  Plus10_amostras = phyloseq::nsamples(ps_plus10),
  Plus10_ASVs = phyloseq::ntaxa(ps_plus10),
  Plus10_funcoes = nrow(res_plus10$counts),
  Plus10_ASVs_mapeadas = res_plus10$checkpoints$ASVs_mapeadas,
  Plus10_media_pct_reads_mapeados = res_plus10$checkpoints$Media_pct_reads_mapeados,
  Interpretacao = "Potencial funcional inferido por taxonomia; nao confirmacao direta de genes/metabolismo.",
  Confundimento = OBS_CONFUNDIMENTO,
  stringsAsFactors = FALSE
)

salvar_csv(metadata_exec, file.path(faprotax_out, "metadata_faprotax.csv"))

readme_txt <- c(
  "# FAPROTAX — Script 11",
  "",
  "Este diretorio contem a inferencia funcional ecologica/metabolica baseada em FAPROTAX via microeco::trans_func.",
  "",
  "## Entradas",
  "- phyloseq_core9_primeira_run.rds",
  "- phyloseq_plus10_com_auxiliar.rds",
  "- ASV_sequences.tsv",
  "",
  "## Motor analitico",
  "- microeco::microtable para organizar otu_table, tax_table e sample_table.",
  "- microeco::trans_func com prok_database = 'FAPROTAX'.",
  "- Nao usa faprotax_function_taxon_map.tsv.",
  "",
  "## Saidas principais",
  "- core9_faprotax_counts.csv/.rds: funcoes x amostras, contagens funcionais sobrepostas.",
  "- core9_faprotax_relative.csv/.rds: funcoes x amostras, abundancia relativa sobre reads totais.",
  "- plus10_faprotax_counts.csv/.rds: analise exploratoria de sensibilidade com 10 amostras, executada separadamente do core9.",
  "- *_asv_faprotax_funcoes.csv: rastreio ASV -> funcao.",
  "- *_asv_faprotax_funcoes_taxonomia.csv: rastreio ASV -> funcao com taxonomia microbiana associada.",
  "- *_funcao_taxa_contribuicao_global.csv: contribuição global de cada táxon/espécie dentro de cada função.",
  "- *_funcao_taxa_contribuicao_por_amostra.csv: contribuição de cada táxon/espécie dentro de cada função por amostra.",
  "- *_legenda_auxiliar_funcao_taxon.csv: códigos F01, F02, ... usados nos gráficos e associação função-táxon.",
  "- *_microeco_res_func_binary.csv: matriz ASV x funcao gerada pelo microeco.",
  "- *_microeco_trans_func.rds: objeto trans_func para auditoria/reuso.",
  "- *_cobertura_faprotax.csv: percentual de reads de ASVs mapeadas a pelo menos uma funcao.",
  "- *_checkpoints_faprotax_microeco.csv: checkpoints de cobertura, ASVs e funcoes detectadas.",
  "",
  "## Figuras",
  paste0("As figuras sao salvas em plots_", pipeline_version, "/faprotax como PDF e PNG:"),
  "- G00: legenda auxiliar associando código, função metabólica e microrganismos contribuintes.",
  "- G01: cobertura funcional por amostra.",
  "- G02: composicao das principais funcoes por amostra.",
  "- G03: composicao funcional media por BeeSpecies e Nativo_Introduzido.",
  "- G04: heatmap das principais funcoes, com percentual de abundância relativa em cada célula.",
  "- G05: riqueza funcional.",
  "- G06: PCoA Bray-Curtis funcional.",
  "",
  "## Interpretacao",
  "FAPROTAX infere potencial funcional a partir de taxonomia. O resultado nao confirma genes, vias ou metabolismo real.",
  "As funcoes nao sao mutuamente exclusivas; a soma das abundancias relativas entre funcoes pode exceder 100%.",
  "",
  "## Estatistica",
  "Core9 e plus10 recebem testes exploratorios separados. O plus10 inclui uma unica amostra da corrida auxiliar; seus resultados nao substituem o core9 e nao permitem separar efeito de Run de BeeSpecies.",
  "Com grupos pequenos, p-valores devem ser reportados como exploratorios e sempre acompanhados de n por grupo.",
  "Nativo_Introduzido e derivado de BeeSpecies e nao representa um efeito independente da especie."
)
writeLines(readme_txt, file.path(faprotax_out, "README_FAPROTAX.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(faprotax_out, "sessionInfo_faprotax.txt"))

cat("\n=============================================================\n")
cat("FAPROTAX — CONCLUIDO\n")
cat("Motor: microeco::trans_func\n")
cat("Saidas:", faprotax_out, "\n")
cat("Figuras:", faprotax_plot, "\n")
cat("Interpretacao: potencial funcional inferido; nao confirmacao direta.\n")
cat("=============================================================\n\n")

rm(ps_core9, ps_plus10, res_core9, res_plus10, asv_map)
gc(verbose = FALSE)
log_msg("Script 11 concluido.", "FINAL")
})
