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
  file.path(.script_dir, "..", "..", "R"),
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
run_pipeline_script("09_deseq2.R", "deseq2", function(ctx) {
###############################################################################
# SCRIPT 09 — DESEQ2: ABUNDÂNCIA DIFERENCIAL
###############################################################################

options(encoding = "UTF-8", stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
})

###############################################################################
# 0. PARÂMETROS GLOBAIS
###############################################################################

VERSAO        <- 
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

base_path <- ctx$base_path
pipeline_version <- ctx$version
out_path <- ctx$output_root
deseq_path <- ctx$stage$root
diag_path <- ctx$stage$figures
for (d in c(deseq_path, diag_path)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
}
ALPHA <- 0.05
LFC_THRESHOLD <- 1.0
PREV_MIN <- 2L
arq_core9 <- ctx$contracts[[]]
arq_plus10 <- ctx$contracts[[]]

###############################################################################
# 1. FUNÇÕES AUXILIARES
###############################################################################

log_msg <- function(msg, tipo = "INFO")
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))


  tryCatch(
    {
      out <- lfcShrink(dds_obj, res = res_raw, type = "ashr")
      attr(out, "shrink_type_used") <- "ashr"
      out
    },
    error = function(e) {
      log_msg(paste("ashr indisponivel ou falhou; usando type='normal':",
                    conditionMessage(e)), "WARN")
      out <- lfcShrink(dds_obj, contrast = contraste, res = res_raw, type = "normal")
      attr(out, "shrink_type_used") <- "normal"
      out
    }
  )
}



filter_threshold_seguro <- function(res_obj) {
  ft <- metadata(res_obj)$filterThreshold
  if (is.null(ft) || length(ft) == 0) return(NA_real_)
  round(as.numeric(ft), 2)
}

anotar_taxonomia <- function(df, taxa_tbl) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (is.null(rownames(df)) || any(rownames(df) == "")) {
    stop("Resultado DESeq2 sem rownames de ASV.", call. = FALSE)
  }
  if (!all(c("ASV_ID", ranks) %in% colnames(taxa_tbl))) {
    stop("tax_table sem ASV_ID ou ranks canonicos.", call. = FALSE)
  }
  if (anyDuplicated(taxa_tbl$ASV_ID) > 0L) {
    stop("tax_table contem ASV_ID duplicado.", call. = FALSE)
  }
  idx <- match(rownames(df), taxa_tbl$ASV_ID)
  if (anyNA(idx)) {
    stop(sum(is.na(idx)), " ASV(s) do DESeq2 sem taxonomia correspondente.", call. = FALSE)
  }
  out <- data.frame(
    ASV_ID = rownames(df),
    df,
    taxa_tbl[idx, ranks, drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out
}

salvar_csv <- function(x, arq) {
  write.csv(x, arq, row.names = FALSE, quote = TRUE, na = "")
  if (!file.exists(arq) || is.na(file.size(arq)) || file.size(arq) == 0L) {
    stop("Falha ao salvar CSV: ", arq, call. = FALSE)
  }
  invisible(arq)
}

###############################################################################
# 2. VALIDAÇÃO E CARREGAMENTO
###############################################################################

cat("=== VALIDAÇÃO E CARREGAMENTO ===\n\n")

if (!file.exists(arq_core9)) stop("Arquivo ausente: ", arq_core9, call. = FALSE)
if (is.na(file.size(arq_core9)) || file.size(arq_core9) == 0L) {
  stop("Arquivo core9 esta vazio: ", arq_core9, call. = FALSE)
}

ps_raw <- tryCatch(
  readRDS(arq_core9),
  error = function(e) stop("Falha ao ler core9: ", conditionMessage(e), call. = FALSE)
)
if (!inherits(ps_raw, "phyloseq")) stop("core9 nao e um objeto phyloseq.", call. = FALSE)
if (nsamples(ps_raw) != 9L) stop("core9 deve conter exatamente 9 amostras; encontrado: ", nsamples(ps_raw), call. = FALSE)
if (ntaxa(ps_raw) < 1L) stop("core9 nao possui ASVs.", call. = FALSE)
if (is.null(otu_table(ps_raw, errorIfNULL = FALSE))) stop("core9 sem otu_table.", call. = FALSE)
if (is.null(sample_data(ps_raw, errorIfNULL = FALSE))) stop("core9 sem sample_data.", call. = FALSE)
if (is.null(tax_table(ps_raw, errorIfNULL = FALSE))) stop("core9 sem tax_table.", call. = FALSE)
if (anyDuplicated(sample_names(ps_raw)) > 0L) stop("core9 possui SampleID duplicado.", call. = FALSE)
if (anyDuplicated(taxa_names(ps_raw)) > 0L) stop("core9 possui ASV_ID duplicado.", call. = FALSE)

meta_check <- as(sample_data(ps_raw), "data.frame")
req_meta <- c("Run", "BeeSpecies")
faltam_meta <- setdiff(req_meta, colnames(meta_check))
if (length(faltam_meta) > 0L) stop("sample_data sem coluna(s): ", paste(faltam_meta, collapse = ", "), call. = FALSE)
if (!identical(sort(rownames(meta_check)), sort(sample_names(ps_raw)))) {
  stop("sample_data e otu_table possuem universos de amostras diferentes.", call. = FALSE)
}
if (any(is.na(meta_check$Run)) || any(as.character(meta_check$Run) != "run_main")) {
  stop("core9 deve conter exclusivamente a corrida run_main.", call. = FALSE)
}

cat(sprintf("core9 carregado: %d amostras | %d ASVs\n",
            nsamples(ps_raw), ntaxa(ps_raw)))

# DESeq2 exige contagens inteiras brutas, finitas e nao negativas.
count_check <- as(otu_table(ps_raw), "matrix")
if (!taxa_are_rows(ps_raw)) count_check <- t(count_check)
storage.mode(count_check) <- "numeric"
if (anyNA(count_check) || any(!is.finite(count_check)) || any(count_check < 0)) {
  stop("otu_table possui NA, valores nao finitos ou contagens negativas.", call. = FALSE)
}
if (any(rowSums(count_check) == 0)) stop("core9 possui ASV(s) com soma zero.", call. = FALSE)
if (any(colSums(count_check) == 0)) stop("core9 possui amostra(s) com soma zero.", call. = FALSE)
desvio_max <- if (length(count_check) == 0L) NA_real_ else max(abs(count_check - round(count_check)))
if (!is.finite(desvio_max) || desvio_max > 1e-8) {
  stop("DESeq2 exige contagens brutas inteiras. Desvio maximo: ", desvio_max, call. = FALSE)
}
log_msg("Contagens brutas inteiras, finitas e nao negativas: OK", "OK")
cat("\n")

###############################################################################
# 3. FILTRO DE PREVALÊNCIA
###############################################################################

cat("=== FILTRO DE PREVALÊNCIA ===\n\n")
cat(sprintf("Critério: ASV presente em >= %d amostras (qualquer grupo)\n", PREV_MIN))

n_antes <- ntaxa(ps_raw)
ps_filt <- filter_taxa(ps_raw, function(x) sum(x > 0) >= PREV_MIN, TRUE)
n_depois <- ntaxa(ps_filt)
if (n_depois < 1L) stop("O filtro de prevalencia removeu todas as ASVs.", call. = FALSE)
cat(sprintf("ASVs: antes = %d | depois = %d | removidas = %d\n",
            n_antes, n_depois, n_antes - n_depois))

# Matriz de contagens (taxa nas linhas, amostras nas colunas)
count_mat <- as(otu_table(ps_filt), "matrix")
if (!taxa_are_rows(ps_filt)) count_mat <- t(count_mat)
count_mat <- round(count_mat)

# Metadados alinhados explicitamente com as colunas da matriz.
meta_all <- as(sample_data(ps_filt), "data.frame")
if (!setequal(rownames(meta_all), colnames(count_mat))) {
  stop("sample_data e count_mat possuem universos de amostras diferentes.", call. = FALSE)
}
meta_d <- meta_all[colnames(count_mat), , drop = FALSE]
if (!identical(rownames(meta_d), colnames(count_mat))) {
  stop("Falha ao alinhar sample_data e count_mat.", call. = FALSE)
}
# Nível de referência = Melipona fasciculata.
meta_d$BeeSpecies <- factor(
  meta_d$BeeSpecies,
  levels = c()
)
if (any(is.na(meta_d$BeeSpecies)))
  stop("BeeSpecies contem níveis inesperados no DESeq2.")
meta_d$Nativo_Introduzido <- factor(ifelse(as.character(meta_d$BeeSpecies) == ),
                                    levels = c())

cat("\nDistribuição por BeeSpecies:\n")
print(table(meta_d$BeeSpecies))
cat("Nível de referência (denominador):", levels(meta_d$BeeSpecies)[1], "\n\n")

# AVISO ESTATÍSTICO — réplicas por grupo:

n_por_grupo <- table(meta_d$BeeSpecies)
esperado_grupos <- c(
 
)
if (!identical(as.integer(n_por_grupo[names(esperado_grupos)]), as.integer(esperado_grupos))) {
  stop(
    "Distribuicao inesperada por BeeSpecies. Esperado 2/3/4; observado: ",
    paste(names(n_por_grupo), as.integer(n_por_grupo), sep = "=", collapse = "; "),
    call. = FALSE
  )
}
grupos_n_baixo <- names(n_por_grupo[n_por_grupo < 3])
if (length(grupos_n_baixo) > 0) {
  detalhe <- paste(sprintf("%s (n=%d)",
                           grupos_n_baixo, n_por_grupo[grupos_n_baixo]),
                   collapse = "; ")
  log_msg(
    paste0("AVISO ESTATISTICO: grupo(s) com n < 3 (minimo recomendado DESeq2): ",
           detalhe, ". Comparacoes que os envolvem sao exploratorias."),
    "WARN")
}

###############################################################################
# 4. CONSTRUIR E AJUSTAR DESEQ2
###############################################################################

cat("=== MODELO DESeq2 (design: ~ BeeSpecies) ===\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta_d,
  design    = ~ BeeSpecies
)

# Tenta ajuste paramétrico; cai para local se a convergência falhar
dds <- tryCatch(
  {
    log_msg("fitType = 'parametric'...", "INFO")
    DESeq(dds, test = "Wald", fitType = "parametric", sfType = "poscounts")
  },
  error = function(e) {
    log_msg(paste("Paramétrico falhou; usando 'local':",
                  conditionMessage(e)), "WARN")
    DESeq(dds, test = "Wald", fitType = "local", sfType = "poscounts")
  }
)
log_msg("DESeq2 ajustado", "OK")

###############################################################################
# 5. FATORES DE TAMANHO (diagnóstico)
###############################################################################

cat("\n=== FATORES DE TAMANHO ===\n\n")

sf <- sizeFactors(dds)
cat("Size factors por amostra:\n")
print(round(sf, 4))

if (anyNA(sf) || any(!is.finite(sf)) || any(sf <= 0)) {
  stop("DESeq2 produziu size factors ausentes, nao finitos ou nao positivos.", call. = FALSE)
}
razao <- max(sf) / min(sf)
cat(sprintf("Razão max/min: %.2f\n", razao))
if (razao > 5)
  log_msg(
    sprintf("Razão max/min = %.1f (> 5); interpretar resultados com cautela.",
            razao), "WARN")

sf_df <- data.frame(
  SampleID   = names(sf),
  SizeFactor = round(sf, 4),
  BeeSpecies = as.character(meta_d[names(sf), "BeeSpecies"]),
  stringsAsFactors = FALSE
)
salvar_csv(sf_df, file.path(deseq_path, "core9_size_factors.csv"))
log_msg("core9_size_factors.csv salvo", "SAVE")

###############################################################################
# 6. DIAGNÓSTICO: GRÁFICO DE DISPERSÃO
###############################################################################

cat("\n=== GRÁFICO DIAGNÓSTICO: DISPERSÃO ===\n\n")

grDevices::pdf(file.path(diag_path, "core9_dispersao.pdf"), width = 8, height = 6)
tryCatch(
  plotDispEsts(dds,
    main = paste0("Estimativas de Dispersão — DESeq2 / core9\n",
                  "Preto = por ASV | Curva = tendência | Azul = encolhida")),
  finally = grDevices::dev.off()
)
log_msg("core9_dispersao.pdf salvo", "SAVE")
cat("Muitos pontos distantes indicam ajuste ruim.\n")

###############################################################################
# 7. TABELA TAXONÔMICA PARA ANOTAÇÃO
###############################################################################

taxa_tbl <- as.data.frame(tax_table(ps_filt), stringsAsFactors = FALSE, check.names = FALSE)
ranks_canonicos <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
faltam_ranks <- setdiff(ranks_canonicos, colnames(taxa_tbl))
if (length(faltam_ranks) > 0L) {
  stop("tax_table sem rank(s) canonico(s): ", paste(faltam_ranks, collapse = ", "), call. = FALSE)
}
taxa_tbl$ASV_ID <- rownames(taxa_tbl)
if (!setequal(taxa_tbl$ASV_ID, rownames(count_mat))) {
  stop("tax_table e count_mat possuem universos de ASVs diferentes.", call. = FALSE)
}

###############################################################################
# 8. COMPARAÇÕES PAIRWISE 
###############################################################################

cat("\n=== COMPARAÇÕES PAIRWISE ===\n\n")

comparacoes <- list(
  list(
    nome      = "",
    contraste = c("")
  ),
  list(
    nome      = "",
    contraste = c("")
  ),
  list(
    nome      = "",
    contraste = c("")
  )
)

resumo_list <- vector("list", length(comparacoes))
resultado_completo_list <- vector("list", length(comparacoes))
p_global_list <- vector("list", length(comparacoes))
shrink_used_list <- character(length(comparacoes))

for (i in seq_along(comparacoes)) {
  comp <- comparacoes[[i]]
  nome <- comp$nome
  cat(sprintf("\n──────────────────────────────\n%s\n──────────────────────────────\n",
              nome))

  # Verificar se a comparação envolve grupo com n < 3
  grupos_comp <- comp$contraste[2:3]
  n_min_comp  <- min(n_por_grupo[grupos_comp])
  comp_exploratoria <- n_min_comp < 3
  if (comp_exploratoria)
    log_msg(sprintf("%s: envolve grupo com n = %d (< 3). Exploratória.",
                    nome, n_min_comp), "WARN")

  # Resultados brutos (sem shrinkage)
  res_bruto <- results(dds,
    contrast      = comp$contraste,
    alpha         = ALPHA,
    pAdjustMethod = "BH")

  n_filt_ind  <- sum(is.na(res_bruto$padj) & !is.na(res_bruto$pvalue))
  thresh_filt <- filter_threshold_seguro(res_bruto)
  cat(sprintf("ASVs filtradas (filtro independente DESeq2): %d | threshold: %s\n",
              n_filt_ind, ifelse(is.na(thresh_filt), "N/A", as.character(thresh_filt))))

  # lfcShrink: ashr se disponível; normal como fallback
  res_shrunk <- shrink_lfc(dds, comp$contraste, res_bruto)
  shrink_used <- attr(res_shrunk, "shrink_type_used")
  if (is.null(shrink_used) || !nzchar(shrink_used)) shrink_used <- "desconhecido"
  shrink_used_list[i] <- shrink_used
  cat("Sumário (LFC encolhido):\n")
  print(summary(res_shrunk, alpha = ALPHA))

  # MA plot; o dispositivo e sempre fechado, inclusive em erro.
  arq_ma <- file.path(diag_path, paste0("core9_MA_", nome, ".pdf"))
  grDevices::pdf(arq_ma, width = 8, height = 5)
  erro_ma <- NULL
  tryCatch(
    {
      DESeq2::plotMA(
        res_shrunk,
        alpha = ALPHA,
        main = paste0("MA — ", gsub("_", " ", nome),
                      "\nLFC encolhido | vermelho: padj < 0.05"),
        ylim = c(-6, 6)
      )
      graphics::abline(h = c(-LFC_THRESHOLD, LFC_THRESHOLD),
                       lty = 2, col = "grey50")
    },
    error = function(e) erro_ma <<- conditionMessage(e),
    finally = grDevices::dev.off()
  )
  if (!is.null(erro_ma)) {
    log_msg(paste("Falha no MA plot", nome, ":", erro_ma), "WARN")
  } else {
    log_msg(paste("MA plot salvo:", basename(arq_ma)), "SAVE")
  }

  # Converter para data.frame e anotar com taxonomia
  res_df          <- as.data.frame(res_shrunk)
  res_df$log2FC_bruto <- as.data.frame(res_bruto)$log2FoldChange
  res_ann <- anotar_taxonomia(res_df, taxa_tbl)

  # Filtrar significativos: padj < ALPHA e |LFC| > LFC_THRESHOLD
  res_sig <- res_ann[
    !is.na(res_ann$padj) &
      res_ann$padj < ALPHA &
      abs(res_ann$log2FoldChange) >= LFC_THRESHOLD, ]
  res_sig <- res_sig[order(res_sig$padj), ]

  n_up   <- sum(res_sig$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_sig$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("Significativas (padj < %.2f e |LFC| > %.1f): %d [up = %d | down = %d]\n",
              ALPHA, LFC_THRESHOLD, nrow(res_sig), n_up, n_down))

  # Armazenar para incluir a correcao BH global entre os tres contrastes.
  resultado_completo_list[[i]] <- res_ann
  p_global_list[[i]] <- data.frame(
    Comparacao = nome,
    ASV_ID = rownames(res_df),
    pvalue = res_df$pvalue,
    stringsAsFactors = FALSE
  )

  resumo_list[[i]] <- data.frame(
    Comparacao           = nome,
    N_min_grupo          = n_min_comp,
    Exploratoria_n_baixo = comp_exploratoria,
    ASVs_testadas        = sum(!is.na(res_df$padj)),
    ASVs_filtradas_indep = n_filt_ind,
    Threshold_filtro     = thresh_filt,
    N_sig                = nrow(res_sig),
    N_up                 = n_up,
    N_down               = n_down,
    LFC_threshold        = LFC_THRESHOLD,
    Shrinkage            = shrink_used,
    stringsAsFactors     = FALSE
  )
}

# Correcao adicional entre todas as hipoteses dos tres contrastes.

p_global <- dplyr::bind_rows(p_global_list)
p_global$padj_global_3_contrastes_BH <- NA_real_
idx_p <- !is.na(p_global$pvalue)
p_global$padj_global_3_contrastes_BH[idx_p] <- p.adjust(p_global$pvalue[idx_p], method = "BH")
salvar_csv(p_global, file.path(deseq_path, "core9_pvalores_tres_contrastes_BH_global.csv"))

for (i in seq_along(comparacoes)) {
  nome <- comparacoes[[i]]$nome
  res_ann <- resultado_completo_list[[i]]
  qg <- p_global[p_global$Comparacao == nome,
                 c("ASV_ID", "padj_global_3_contrastes_BH"), drop = FALSE]
  idx_qg <- match(res_ann$ASV_ID, qg$ASV_ID)
  if (anyNA(idx_qg)) stop("Falha ao alinhar padj global em ", nome, call. = FALSE)
  res_ann$padj_global_3_contrastes_BH <- qg$padj_global_3_contrastes_BH[idx_qg]
  res_sig <- res_ann[
    !is.na(res_ann$padj) & res_ann$padj < ALPHA &
      abs(res_ann$log2FoldChange) >= LFC_THRESHOLD,
    , drop = FALSE
  ]
  res_sig <- res_sig[order(res_sig$padj), , drop = FALSE]
  arq_sig  <- file.path(deseq_path, paste0("core9_sig_", nome, ".csv"))
  arq_full <- file.path(deseq_path, paste0("core9_completo_", nome, ".csv"))
  salvar_csv(res_sig, arq_sig)
  salvar_csv(res_ann, arq_full)
  resumo_list[[i]]$N_sig_BH_global <- sum(
    !is.na(res_ann$padj_global_3_contrastes_BH) &
      res_ann$padj_global_3_contrastes_BH < ALPHA &
      abs(res_ann$log2FoldChange) >= LFC_THRESHOLD
  )
  log_msg(paste("Salvos:", basename(arq_sig), "/", basename(arq_full)), "SAVE")
}

###############################################################################
# 8B. CONTRASTE DESCRITIVO
###############################################################################
cat("\n=== MODELO DESeq2 ===\n\n")

status_path <- file.path(deseq_path, "")
dir.create(status_path, recursive = TRUE, showWarnings = FALSE)

dds_status <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta_d,
  design    = ~ Nativo_Introduzido
)

dds_status <- tryCatch(
  {
    log_msg("Nativo/Introduzido: fitType = 'parametric'...", "INFO")
    DESeq(dds_status, test = "Wald", fitType = "parametric", sfType = "poscounts")
  },
  error = function(e) {
    log_msg(paste("Nativo/Introduzido paramétrico falhou; usando 'local':",
                  conditionMessage(e)), "WARN")
    DESeq(dds_status, test = "Wald", fitType = "local", sfType = "poscounts")
  }
)

contraste_status <- c("")
res_status_bruto <- results(dds_status,
  contrast = contraste_status,
  alpha = ALPHA,
  pAdjustMethod = "BH")
res_status_shrunk <- shrink_lfc(dds_status, contraste_status, res_status_bruto)
shrink_status <- attr(res_status_shrunk, "shrink_type_used")
if (is.null(shrink_status) || !nzchar(shrink_status)) shrink_status <- "desconhecido"

res_status_df <- as.data.frame(res_status_shrunk)
res_status_df$log2FC_bruto <- as.data.frame(res_status_bruto)$log2FoldChange
res_status_ann <- anotar_taxonomia(res_status_df, taxa_tbl)
res_status_sig <- res_status_ann[
  !is.na(res_status_ann$padj) &
    res_status_ann$padj < ALPHA &
    abs(res_status_ann$log2FoldChange) >= LFC_THRESHOLD, ]
res_status_sig <- res_status_sig[order(res_status_sig$padj), ]

salvar_csv(
  res_status_ann,
  file.path(status_path, "")
)
salvar_csv(
  res_status_sig,
  file.path(status_path, "c")
)

resumo_status <- data.frame(
  Comparacao = "a",
  Design = 
  PREV_MIN = PREV_MIN,
  N_introduzida = sum(meta_d$Nativo_Introduzido == "Introduzida"),
  N_nativa = sum(meta_d$Nativo_Introduzido == "Nativa"),
  ASVs_testadas = sum(!is.na(res_status_df$padj)),
  N_sig = nrow(res_status_sig),
  Exploratoria_n_baixo = TRUE,
  Confundida_com_BeeSpecies = TRUE,
  Shrinkage = shrink_status,
  stringsAsFactors = FALSE)
salvar_csv(
  resumo_status,
  file.path(status_path, "core9_resumo_nativa_vs_introduzida.csv")
)
log_msg("DESeq2 nativo/introduzido salvo em deseq2/nativo_introduzido_prev2", "SAVE")


###############################################################################
# 8C. SENSIBILIDADE EXPLORATORIA PLUS10 — DESEQ2
################################################################################
cat("\n=== PLUS10 — DESeq2 EXPLORATORIO SEPARADO ===\n\n")
log_msg(
  paste0(
    "sensibilidade exploratoria"
  ),
  "WARN"
)

if (!file.exists(arq_plus10) || is.na(file.size(arq_plus10)) || file.size(arq_plus10) == 0L) {
  stop("Arquivo plus10 ausente ou vazio: ", arq_plus10, call. = FALSE)
}
ps_p10_raw <- tryCatch(
  readRDS(arq_plus10),
  error = function(e) stop("Falha ao ler plus10: ", conditionMessage(e), call. = FALSE)
)
if (!inherits(ps_p10_raw, "phyloseq") || nsamples(ps_p10_raw) != 10L) {
  stop("plus10 deve ser phyloseq com exatamente 10 amostras.", call. = FALSE)
}
meta_p10_check <- as(sample_data(ps_p10_raw), "data.frame")
meta_p10_check <- meta_p10_check[sample_names(ps_p10_raw), , drop = FALSE]
if (!all(c("Run", "BeeSpecies") %in% colnames(meta_p10_check)) ||
    sum(as.character(meta_p10_check$Run) == "run_main") != 9L ||
    sum(as.character(meta_p10_check$Run) == "run_aux") != 1L ||
    !"S10" %in% rownames(meta_p10_check) ||
    as.character(meta_p10_check["S10", "Run"]) != "run_aux") {
  stop("plus10 nao respeita o desenho 9+1 com Auxiliar na segunda corrida.", call. = FALSE)
}

count_p10_all <- as(otu_table(ps_p10_raw), "matrix")
if (!taxa_are_rows(ps_p10_raw)) count_p10_all <- t(count_p10_all)
storage.mode(count_p10_all) <- "numeric"
if (anyNA(count_p10_all) || any(!is.finite(count_p10_all)) || any(count_p10_all < 0) ||
    max(abs(count_p10_all - round(count_p10_all))) > 1e-8) {
  stop("plus10: contagens invalidas para DESeq2.", call. = FALSE)
}

prev_p10 <- rowSums(count_p10_all > 0)
reads_p10_total <- rowSums(count_p10_all)
presenca_auxiliar_p10 <- count_p10_all[, "S10"] > 0
presenca_outras_p10 <- rowSums(
  count_p10_all[, setdiff(colnames(count_p10_all), "S10"), drop = FALSE] > 0
) > 0
auditoria_universo_p10 <- data.frame(
  ASV_ID = rownames(count_p10_all),
  Reads_totais_plus10 = as.numeric(reads_p10_total),
  Prevalencia_plus10 = as.integer(prev_p10),
  Presente_Auxiliar = as.logical(presenca_auxiliar_p10),
  Presente_outras_amostras = as.logical(presenca_outras_p10),
  Exclusiva_Auxiliar = as.logical(presenca_auxiliar_p10 & !presenca_outras_p10),
  Elegivel_DESeq2_PREV2 = prev_p10 >= PREV_MIN,
  Motivo_nao_elegivel = ifelse(
    prev_p10 >= PREV_MIN,
    NA_character_,
    paste(
      "Prevalencia inferior a", PREV_MIN,
      "; efeito biologico nao separavel de uma unica unidade amostral."
    )
  ),
  stringsAsFactors = FALSE
)
taxa_p10_all <- as.data.frame(
  tax_table(ps_p10_raw), stringsAsFactors = FALSE, check.names = FALSE
)
taxa_p10_all$ASV_ID <- rownames(taxa_p10_all)
auditoria_universo_p10 <- merge(
  auditoria_universo_p10, taxa_p10_all,
  by = "ASV_ID", all.x = TRUE, sort = FALSE
)
auditoria_universo_p10 <- auditoria_universo_p10[
  match(rownames(count_p10_all), auditoria_universo_p10$ASV_ID), , drop = FALSE
]
salvar_csv(
  auditoria_universo_p10,
  file.path(deseq_path, "plus10_universo_ASVs_elegibilidade_DESeq2.csv")
)

exclusivas_p10 <- names(prev_p10)[prev_p10 == 1L]
if (length(exclusivas_p10)) {
  amostra_unica <- vapply(exclusivas_p10, function(id) {
    colnames(count_p10_all)[which(count_p10_all[id, ] > 0)[1L]]
  }, character(1))
  audit_exclusivas <- data.frame(
    ASV_ID = exclusivas_p10,
    Amostra_unica = amostra_unica,
    Reads = vapply(seq_along(exclusivas_p10), function(ii) {
      count_p10_all[exclusivas_p10[ii], amostra_unica[ii]]
    }, numeric(1)),
    Incluida_DESeq2_PREV2 = FALSE,
    Motivo = paste(
      "Prevalencia=1; efeito biologico nao separavel de amostra/run;",
      "mantida em auditoria, nao forçada no teste diferencial."
    ),
    stringsAsFactors = FALSE
  )
  salvar_csv(audit_exclusivas,
             file.path(deseq_path, "plus10_asvs_prevalencia1_nao_testaveis.csv"))
}

ps_p10 <- filter_taxa(ps_p10_raw, function(x) sum(x > 0) >= PREV_MIN, TRUE)
if (ntaxa(ps_p10) < 1L) stop("plus10: filtro de prevalencia removeu todas as ASVs.")
count_p10 <- as(otu_table(ps_p10), "matrix")
if (!taxa_are_rows(ps_p10)) count_p10 <- t(count_p10)
count_p10 <- round(count_p10)
meta_p10 <- as(sample_data(ps_p10), "data.frame")
meta_p10 <- meta_p10[colnames(count_p10), , drop = FALSE]
meta_p10$BeeSpecies <- factor(
  meta_p10$BeeSpecies,
  levels = c("Melipona fasciculata", "Melipona scutellaris", "Melipona subnitida")
)
if (anyNA(meta_p10$BeeSpecies) ||
    !identical(as.integer(table(meta_p10$BeeSpecies)), c(3L, 3L, 4L))) {
  stop("plus10: distribuicao BeeSpecies deve ser 3/3/4.")
}
meta_p10$Nativo_Introduzido <- factor(
  ifelse(as.character(meta_p10$BeeSpecies) == "Melipona fasciculata",
         "Introduzida", "Nativa"),
  levels = c("Introduzida", "Nativa")
)

# Tabela de desenho e impossibilidade de ajuste de Run
mat_design_p10 <- model.matrix(~ Run + BeeSpecies, data = meta_p10)
rank_design_p10 <- qr(mat_design_p10)$rank
rank_full_p10 <- ncol(mat_design_p10)
design_audit_p10 <- data.frame(
  Conjunto = "plus10",
  N_amostras = nrow(meta_p10),
  Distribuicao_Run = paste(names(table(meta_p10$Run)), table(meta_p10$Run), sep = "=", collapse = ";"),
  Distribuicao_BeeSpecies = paste(names(table(meta_p10$BeeSpecies)), table(meta_p10$BeeSpecies), sep = "=", collapse = ";"),
  Formula_avaliada = "~ Run + BeeSpecies",
  Colunas_modelo = rank_full_p10,
  Posto_modelo = rank_design_p10,
  Modelo_com_Run_posto_completo = rank_design_p10 == rank_full_p10,
  Modelo_usado = "~ BeeSpecies",
  Run_ajustada = FALSE,
  Justificativa = paste(
    "A segunda corrida possui uma unica amostra e esta pertence a M. fasciculata;",
    "mesmo quando a matriz formal tem posto completo, nao ha replicacao para",
    "estimar variacao dentro de run_aux nem separar batch da amostra Auxiliar."
  ),
  stringsAsFactors = FALSE
)
salvar_csv(design_audit_p10, file.path(deseq_path, "plus10_auditoria_design_run.csv"))

dds_p10 <- DESeqDataSetFromMatrix(
  countData = count_p10,
  colData = meta_p10,
  design = ~ BeeSpecies
)
dds_p10 <- tryCatch(
  {
    DESeq(dds_p10, test = "Wald", fitType = "parametric", sfType = "poscounts")
  },
  error = function(e) {
    log_msg(paste("plus10: ajuste parametrico falhou; usando local:",
                  conditionMessage(e)), "WARN")
    DESeq(dds_p10, test = "Wald", fitType = "local", sfType = "poscounts")
  }
)

sf_p10 <- sizeFactors(dds_p10)
sf_p10_df <- data.frame(
  SampleID = names(sf_p10), SizeFactor = as.numeric(sf_p10),
  BeeSpecies = as.character(meta_p10[names(sf_p10), "BeeSpecies"]),
  Run = as.character(meta_p10[names(sf_p10), "Run"]),
  Exploratoria_segunda_run = TRUE,
  stringsAsFactors = FALSE
)
salvar_csv(sf_p10_df, file.path(deseq_path, "plus10_size_factors.csv"))

grDevices::pdf(file.path(diag_path, "plus10_dispersao.pdf"), width = 8, height = 6)
tryCatch(
  plotDispEsts(dds_p10,
               main = "Estimativas de dispersao — DESeq2 / plus10 exploratorio"),
  finally = grDevices::dev.off()
)

taxa_p10 <- as.data.frame(tax_table(ps_p10), stringsAsFactors = FALSE, check.names = FALSE)
taxa_p10$ASV_ID <- rownames(taxa_p10)

resumo_p10_list <- vector("list", length(comparacoes))
completos_p10 <- vector("list", length(comparacoes))
pglob_p10_list <- vector("list", length(comparacoes))
shrink_p10_used <- character(length(comparacoes))

for (i in seq_along(comparacoes)) {
  comp <- comparacoes[[i]]
  nome <- comp$nome
  res_raw <- results(
    dds_p10, contrast = comp$contraste, alpha = ALPHA, pAdjustMethod = "BH"
  )
  res_shr <- shrink_lfc(dds_p10, comp$contraste, res_raw)
  shrink_nm <- attr(res_shr, "shrink_type_used")
  if (is.null(shrink_nm) || !nzchar(shrink_nm)) shrink_nm <- "desconhecido"
  shrink_p10_used[i] <- shrink_nm

  arq_ma <- file.path(diag_path, paste0("plus10_MA_", nome, ".pdf"))
  grDevices::pdf(arq_ma, width = 8, height = 5)
  tryCatch(
    {
      DESeq2::plotMA(
        res_shr, alpha = ALPHA,
        main = paste0("MA — plus10 exploratorio — ", gsub("_", " ", nome)),
        ylim = c(-6, 6)
      )
      graphics::abline(h = c(-LFC_THRESHOLD, LFC_THRESHOLD), lty = 2, col = "grey50")
    },
    finally = grDevices::dev.off()
  )

  df <- as.data.frame(res_shr)
  df$log2FC_bruto <- as.data.frame(res_raw)$log2FoldChange
  ann <- anotar_taxonomia(df, taxa_p10)
  ann$Conjunto <- "plus10"
  ann$Exploratoria_segunda_run <- TRUE
  ann$Run_ajustada <- FALSE
  ann$Nota_batch <- paste(
    "Uma unica amostra pertence a run_aux; mudancas podem refletir batch;",
    "interpretar apenas pela concordancia com core9."
  )
  sig <- ann[
    !is.na(ann$padj) & ann$padj < ALPHA &
      abs(ann$log2FoldChange) >= LFC_THRESHOLD,
    , drop = FALSE
  ]
  sig <- sig[order(sig$padj), , drop = FALSE]
  salvar_csv(ann, file.path(deseq_path, paste0("plus10_completo_", nome, ".csv")))
  salvar_csv(sig, file.path(deseq_path, paste0("plus10_sig_", nome, ".csv")))
  completos_p10[[i]] <- ann
  pglob_p10_list[[i]] <- data.frame(
    Comparacao = nome, ASV_ID = rownames(df), pvalue = df$pvalue,
    stringsAsFactors = FALSE
  )
  resumo_p10_list[[i]] <- data.frame(
    Conjunto = "plus10", Comparacao = nome,
    N_min_grupo = min(table(meta_p10$BeeSpecies)[comp$contraste[2:3]]),
    ASVs_antes_filtro = ntaxa(ps_p10_raw),
    ASVs_apos_PREV2 = ntaxa(ps_p10),
    ASVs_testadas = sum(!is.na(df$padj)),
    N_sig = nrow(sig), Shrinkage = shrink_nm,
    Exploratoria_segunda_run = TRUE, Run_ajustada = FALSE,
    stringsAsFactors = FALSE
  )
}

pglob_p10 <- dplyr::bind_rows(pglob_p10_list)
pglob_p10$padj_global_3_contrastes_BH <- NA_real_
idx_pg10 <- !is.na(pglob_p10$pvalue)
pglob_p10$padj_global_3_contrastes_BH[idx_pg10] <-
  p.adjust(pglob_p10$pvalue[idx_pg10], method = "BH")
salvar_csv(pglob_p10,
           file.path(deseq_path, "plus10_pvalores_tres_contrastes_BH_global.csv"))

for (i in seq_along(comparacoes)) {
  nome <- comparacoes[[i]]$nome
  ann <- completos_p10[[i]]
  qg <- pglob_p10[pglob_p10$Comparacao == nome,
                  c("ASV_ID", "padj_global_3_contrastes_BH"), drop = FALSE]
  ann$padj_global_3_contrastes_BH <-
    qg$padj_global_3_contrastes_BH[match(ann$ASV_ID, qg$ASV_ID)]
  salvar_csv(ann, file.path(deseq_path, paste0("plus10_completo_", nome, ".csv")))

  core_file <- file.path(deseq_path, paste0("core9_completo_", nome, ".csv"))
  if (file.exists(core_file)) {
    core_tab <- read.csv(core_file, stringsAsFactors = FALSE, check.names = FALSE)
    cmp <- merge(
      core_tab[, c("ASV_ID", "log2FoldChange", "pvalue", "padj")],
      ann[, c("ASV_ID", "log2FoldChange", "pvalue", "padj")],
      by = "ASV_ID", all = TRUE, suffixes = c("_core9", "_plus10")
    )
    cmp$Direcao_core9 <- sign(cmp$log2FoldChange_core9)
    cmp$Direcao_plus10 <- sign(cmp$log2FoldChange_plus10)
    cmp$Direcao_concordante <- with(
      cmp,
      is.finite(log2FoldChange_core9) & is.finite(log2FoldChange_plus10) &
        Direcao_core9 == Direcao_plus10
    )
    cmp$Sig_core9 <- !is.na(cmp$padj_core9) & cmp$padj_core9 < ALPHA &
      abs(cmp$log2FoldChange_core9) >= LFC_THRESHOLD
    cmp$Sig_plus10 <- !is.na(cmp$padj_plus10) & cmp$padj_plus10 < ALPHA &
      abs(cmp$log2FoldChange_plus10) >= LFC_THRESHOLD
    cmp$Sig_concordante <- cmp$Sig_core9 == cmp$Sig_plus10
    salvar_csv(cmp,
               file.path(deseq_path, paste0("sensibilidade_core9_vs_plus10_", nome, ".csv")))
  }
}

# Contraste composto Nativa/Introduzida no plus10
dds_status_p10 <- DESeqDataSetFromMatrix(
  countData = count_p10, colData = meta_p10, design = ~ Nativo_Introduzido
)
dds_status_p10 <- tryCatch(
  DESeq(dds_status_p10, test = "Wald", fitType = "parametric", sfType = "poscounts"),
  error = function(e) DESeq(dds_status_p10, test = "Wald", fitType = "local", sfType = "poscounts")
)
contraste_status_p10 <- c("Nativo_Introduzido", "Nativa", "Introduzida")
res_status_p10_raw <- results(dds_status_p10, contrast = contraste_status_p10,
                              alpha = ALPHA, pAdjustMethod = "BH")
res_status_p10_shr <- shrink_lfc(dds_status_p10, contraste_status_p10, res_status_p10_raw)
status_p10_df <- as.data.frame(res_status_p10_shr)
status_p10_df$log2FC_bruto <- as.data.frame(res_status_p10_raw)$log2FoldChange
status_p10_ann <- anotar_taxonomia(status_p10_df, taxa_p10)
status_p10_ann$Conjunto <- "plus10"
status_p10_ann$Exploratoria_segunda_run <- TRUE
status_p10_ann$Run_ajustada <- FALSE
status_p10_sig <- status_p10_ann[
  !is.na(status_p10_ann$padj) & status_p10_ann$padj < ALPHA &
    abs(status_p10_ann$log2FoldChange) >= LFC_THRESHOLD,
  , drop = FALSE
]
salvar_csv(status_p10_ann,
           file.path(status_path, "plus10_completo_nativa_vs_introduzida.csv"))
salvar_csv(status_p10_sig,
           file.path(status_path, "plus10_sig_nativa_vs_introduzida.csv"))

resumo_p10_df <- dplyr::bind_rows(resumo_p10_list)
salvar_csv(resumo_p10_df, file.path(deseq_path, "plus10_resumo_comparacoes.csv"))

# Resumo de estabilidade dos tres contrastes
sens_resumo <- dplyr::bind_rows(lapply(comparacoes, function(comp) {
  arq_cmp <- file.path(deseq_path,
                       paste0("sensibilidade_core9_vs_plus10_", comp$nome, ".csv"))
  if (!file.exists(arq_cmp)) return(NULL)
  z <- read.csv(arq_cmp, stringsAsFactors = FALSE)
  ok_lfc <- is.finite(z$log2FoldChange_core9) & is.finite(z$log2FoldChange_plus10)
  data.frame(
    Comparacao = comp$nome,
    ASVs_comuns_com_LFC = sum(ok_lfc),
    Proporcao_direcao_concordante = if (sum(ok_lfc)) mean(z$Direcao_concordante[ok_lfc]) else NA_real_,
    Spearman_LFC = if (sum(ok_lfc) >= 3L)
      suppressWarnings(cor(z$log2FoldChange_core9[ok_lfc],
                           z$log2FoldChange_plus10[ok_lfc], method = "spearman"))
      else NA_real_,
    N_sig_core9 = sum(z$Sig_core9, na.rm = TRUE),
    N_sig_plus10 = sum(z$Sig_plus10, na.rm = TRUE),
    N_sig_em_ambos = sum(z$Sig_core9 & z$Sig_plus10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
salvar_csv(sens_resumo,
           file.path(deseq_path, "sensibilidade_deseq2_core9_vs_plus10_resumo.csv"))
log_msg("DESeq2 plus10 exploratorio e tabelas de concordancia salvos", "SAVE")

###############################################################################
# 9. TABELA RESUMO DAS COMPARAÇÕES
###############################################################################

cat("\n=== TABELA RESUMO ===\n\n")
resumo_df <- dplyr::bind_rows(resumo_list)
print(resumo_df)
salvar_csv(
  resumo_df,
  file.path(deseq_path, "core9_resumo_comparacoes.csv")
)
log_msg("core9_resumo_comparacoes.csv salvo", "SAVE")

###############################################################################
# 10. INFORMAÇÕES DA SESSÃO
###############################################################################

cat("\n=== SESSÃO R ===\n")
print(sessionInfo())

###############################################################################
# 11. METADADOS DE EXECUÇÃO
###############################################################################

run_meta <- data.frame(
  Script          = "9_deseq2",
  Versao          = VERSAO,
  Data_execucao   = DATA_EXECUCAO,
  Input           = "phyloseq_core9_primeira_run.rds; phyloseq_plus10_com_auxiliar.rds",
  Design          = "~ BeeSpecies; contraste composto ~ Nativo_Introduzido (confundido com especie)",
  Teste           = "Wald",
  Correcao_FDR    = paste(
    "BH por contraste (primaria); BH global adicional sobre",
    "todas as hipoteses dos 3 contrastes"
  ),
  LFC_threshold   = LFC_THRESHOLD,
  Alpha           = ALPHA,
  Filtro_prev_min = PREV_MIN,
  Shrinkage       = paste(unique(c(shrink_used_list, shrink_status)), collapse = ";"),
  Conjunto        = "core9 principal + plus10 sensibilidade exploratoria separada; Run nao ajustavel",
  stringsAsFactors = FALSE
)
salvar_csv(
  run_meta,
  file.path(deseq_path, "metadata_execucao_script9.csv")
)

cat("\n=============================================================\n")
cat("SCRIPT 09 (DESeq2) CONCLUÍDO\n")
cat("Saídas em:", deseq_path, "\n")
cat("=============================================================\n\n")

gc(verbose = FALSE)
log_msg("Finalizado", "FINAL")
})
