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
run_pipeline_script("10_graficos.R", "graficos", function(ctx) {
###############################################################################
# SCRIPT 10 — VISUALIZAÇÕES COMPLETAS (SEM RECALCULO ESTATISTICO)
#
# Gráficos gerados:
#   G01  Diversidade alfa        (jitter + medias + teste de postos exato)
#   G02  PCoA PC1 × PC2          (spider + convex hull + PERMANOVA)
#   G03  PCoA PC1 × PC3
#   G04  NMDS                    (elipses 90% + stress)
#   G05  Composição — Filo
#   G06  Composição — Gênero (todos identificados)
#   G07  Composição — Top 15 Gêneros
#   G08  iNEXT (rarefação/extrapolação)
#   G09  Whittaker (partição beta Sørensen)
#   G10  Heatmap top 20 Gêneros
#   G11  Diagrama de Venn (ASVs compartilhadas)
#   G12  Volcano DESeq2 (3 comparações)
#   G13  SIMPER top 10 ASVs por comparação
#   G14  IndVal ASVs indicadoras (FDR-BH < 0.05)
#
# LIMITAÇÃO DE DESENHO (relevante para a interpretação das figuras):
#   As perguntas de pesquisa 1 e 2 (quais microrganismos; a espécie influencia
#   o perfil) são abordadas por estas figuras. A pergunta 3 (a geografia indica
#   presença/ausência) NÃO é separável da espécie: espécie e município/ambiente
#   estão confundidos no desenho amostral (ver Script 08, seção 9). Cores e
#   agrupamentos por BeeSpecies nestas figuras carregam, inseparavelmente,
#   variação geográfica. Não atribuir diferenças exclusivamente à espécie.
#
# Dependências:
#   output_V1/phyloseq_core9_primeira_run.rds   (Script 06)
#   output_V1/phyloseq_plus10_com_auxiliar.rds     (Script 06)
#   output_V1/core9_dist_bray_rel.rds           (Script 06)
#   output_V1/core9_dist_jaccard_binary.rds     (Script 06)
#   output_V1/deseq2/core9_completo_*.csv       (Script 09)
#   output_V1/analises/core9_kruskal_wallis.csv (Script 08, opcional)
#   output_V1/analises/core9_simper_top10_*.csv (Script 08, opcional)
#   output_V1/analises/core9_indval_*.csv       (Script 08, opcional)
###############################################################################

options(encoding = "UTF-8", stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(iNEXT)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggVennDiagram)
  library(ggrepel)
  library(betapart)
  library(indicspecies)
})

###############################################################################
# 0. PARÂMETROS GLOBAIS
###############################################################################

VERSAO        <- "3.1_core9_plus10_inext_integrado"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

base_path <- ctx$base_path
pipeline_version <- ctx$version
out_path <- ctx$output_root
fig_path <- ctx$stage$figures
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(fig_path)) stop("Falha ao criar diretorio: ", fig_path, call. = FALSE)
N_PERM <- 9999L
SEED <- 1234L
ALPHA <- 0.05
LFC_THRESHOLD <- 1.0
LIMIAR_RAZAO_PROFUNDIDADE <- 10
arq_core9 <- ctx$contracts[["phyloseq_core9"]]
arq_dist_bray <- ctx$contracts[["core9_dist_bray"]]
arq_dist_jaccard <- ctx$contracts[["core9_dist_jaccard"]]
deseq_path <- ctx$layout$stages$deseq2$root
anal_path <- ctx$layout$stages$analises$root

###############################################################################
# 1. FUNÇÕES AUXILIARES
###############################################################################

log_msg <- function(msg, tipo = "INFO")
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))

# Funcoes compartilhadas para testes exatos por enumeracao. O arquivo e
# localizado pelo diretorio informado pelo orquestrador ou pelo proprio script.
arq_funcoes_exatas <- file.path(.pipeline_lib_dir, "funcoes_estatisticas_exatas.R")
if (!file.exists(arq_funcoes_exatas)) {
  stop("funcoes_estatisticas_exatas.R nao encontrado em: ", arq_funcoes_exatas, call. = FALSE)
}
sys.source(arq_funcoes_exatas, envir = .GlobalEnv)

salvar_plot <- function(p, arquivo, dir_saida, width = 10, height = 7)
  ggplot2::ggsave(file.path(dir_saida, arquivo), plot = p,
                  width = width, height = height, dpi = 300)

# Converte SampleID para o rótulo canônico de exibição SampleLabel.
# SampleID permanece como chave interna; SampleLabels é aceito como alias legado.
validar_samplelabel_meta <- function(meta) {
  if (!"SampleID" %in% colnames(meta)) meta$SampleID <- rownames(meta)
  if (!"SampleLabel" %in% colnames(meta) && "SampleLabels" %in% colnames(meta)) {
    colnames(meta)[colnames(meta) == "SampleLabels"] <- "SampleLabel"
  }
  if (!"SampleLabel" %in% colnames(meta))
    stop("SampleLabel ausente no sample_data; execute novamente o Script 06 atualizado.")
  meta$SampleLabel <- trimws(as.character(meta$SampleLabel))
  if (any(is.na(meta$SampleLabel)) || any(meta$SampleLabel == ""))
    stop("SampleLabel ausente ou vazio no sample_data.")
  if (anyDuplicated(meta$SampleLabel) > 0)
    stop("SampleLabel duplicado no sample_data; os rotulos dos graficos devem ser unicos.")
  meta
}

rotulo_amostra <- function(sample_id, meta = NULL) {
  if (is.null(meta)) meta <- get0("meta_c9", ifnotfound = NULL)
  if (is.null(meta)) return(as.character(sample_id))
  meta <- validar_samplelabel_meta(meta)
  mapa <- setNames(meta$SampleLabel, meta$SampleID)
  ids <- as.character(sample_id)
  out <- unname(mapa[ids])
  out[is.na(out) | out == ""] <- ids[is.na(out) | out == ""]
  out
}

# Matriz de ASVs com amostras nas linhas
matriz_amostras <- function(ps_obj) {
  mat <- as(otu_table(ps_obj), "matrix")
  if (taxa_are_rows(ps_obj)) mat <- t(mat)
  mat
}

relativo <- function(mat) {
  tot <- rowSums(mat)
  if (any(tot == 0))
    stop("Matriz contem amostra(s) com soma zero; abundancia relativa indefinida.")
  sweep(mat, 1, tot, "/")
}

gerar_rotulos_amostras <- function(meta) {
  meta <- validar_samplelabel_meta(meta)
  as.character(meta$SampleLabel)
}

alinhar_distancia <- function(dist_obj, sample_ids, nome) {
  if (!inherits(dist_obj, "dist")) stop(nome, " nao e objeto dist.", call. = FALSE)
  labs <- attr(dist_obj, "Labels")
  if (is.null(labs) || any(labs == "") || anyDuplicated(labs) > 0L) {
    stop(nome, " possui Labels ausentes ou duplicados.", call. = FALSE)
  }
  if (!setequal(labs, sample_ids)) {
    stop(nome, " e o conjunto informado possuem universos de amostras diferentes.", call. = FALSE)
  }
  if (!identical(labs, sample_ids)) {
    m <- as.matrix(dist_obj)[sample_ids, sample_ids, drop = FALSE]
    dist_obj <- as.dist(m)
  }
  dist_obj
}

###############################################################################
# 2. VALIDAÇÕES
###############################################################################

cat("=== VALIDAÇÕES ===\n\n")
for (arq in c(arq_core9)) {
  if (!file.exists(arq)) stop("Arquivo ausente: ", arq, call. = FALSE)
  if (is.na(file.size(arq)) || file.size(arq) == 0L) {
    stop("Arquivo vazio: ", arq, call. = FALSE)
  }
}
log_msg("Arquivos principais validados", "OK")
cat("\n")

###############################################################################
# 3. PALETA, TEMA E CONSTANTES VISUAIS
###############################################################################

# Cores por espécie (chave = nome completo, como em BeeSpecies)
CORES_ESP <- c(
  "Melipona fasciculata" = "#FF6361",
  "Melipona scutellaris" = "#000489",
  "Melipona subnitida"   = "#21EBC9"
)

# Rótulos curtos para legendas e eixos
LABELS_ESP <- c(
  "Melipona fasciculata" = "M. fasciculata",
  "Melipona scutellaris" = "M. scutellaris",
  "Melipona subnitida"   = "M. subnitida"
)

CORES_STATUS <- c("Introduzida" = "#FF6361", "Nativa" = "#2f4b7c")

TEMA_PUB <- theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey85", color = "grey60"),
    strip.text         = element_text(face = "bold.italic", size = 10),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(face = "bold"),
    legend.position    = "right",
    legend.title       = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey30")
  )

# Paleta de 29 cores base para filos/gêneros (interpolada se necessário)
CORES_MULTI <- c(
  "#d45087","#6A3D9A","#1F78B4","#33A02C","#FF7F00",
  "#00BCD4","#F06292","#FFD600","#00897B","#FB8C00",
  "#8E24AA","#1565C0","#D81B60","#43A047","#2f4b7c",
  "#546E7A","#665191","#0288D1","#558B2F","#EF6C00",
  "#4527A0","#00695C","#AD1457","#827717","#37474F",
  "#1A237E","#004D40","#F57F17","#880E4F"
)

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  paste0("p = ", round(p, 3), ifelse(p < ALPHA, "*", " n.s."))
}

###############################################################################
# 4. CARREGAR OBJETOS
###############################################################################

cat("=== CARREGAMENTO ===\n\n")

ps_c9 <- tryCatch(
  readRDS(arq_core9),
  error = function(e) stop("Falha ao ler core9: ", conditionMessage(e), call. = FALSE)
)
if (!inherits(ps_c9, "phyloseq")) stop("core9 nao e phyloseq.", call. = FALSE)
if (nsamples(ps_c9) != 9L) stop("core9 deve conter exatamente 9 amostras.", call. = FALSE)
if (ntaxa(ps_c9) < 1L) stop("core9 nao possui ASVs.", call. = FALSE)
if (is.null(sample_data(ps_c9, errorIfNULL = FALSE)) ||
    is.null(tax_table(ps_c9, errorIfNULL = FALSE)) ||
    is.null(otu_table(ps_c9, errorIfNULL = FALSE))) {
  stop("core9 sem sample_data, tax_table ou otu_table.", call. = FALSE)
}

meta_c9 <- as(sample_data(ps_c9), "data.frame")
meta_c9 <- validar_samplelabel_meta(meta_c9)
if (!setequal(rownames(meta_c9), sample_names(ps_c9))) {
  stop("sample_data e otu_table possuem universos de amostras diferentes.", call. = FALSE)
}
meta_c9 <- meta_c9[sample_names(ps_c9), , drop = FALSE]
if (!all(c("Run", "BeeSpecies") %in% colnames(meta_c9))) {
  stop("sample_data deve conter Run e BeeSpecies.", call. = FALSE)
}
if (any(is.na(meta_c9$Run)) || any(as.character(meta_c9$Run) != "run_main")) {
  stop("core9 deve conter exclusivamente run_main.", call. = FALSE)
}
meta_c9$BeeSpecies <- factor(
  meta_c9$BeeSpecies,
  levels = c("Melipona fasciculata", "Melipona scutellaris", "Melipona subnitida")
)
if (anyNA(meta_c9$BeeSpecies)) stop("BeeSpecies contem nivel inesperado.", call. = FALSE)
n_esp <- table(meta_c9$BeeSpecies)
if (!identical(as.integer(n_esp), c(2L, 3L, 4L))) {
  stop("Distribuicao por especie diferente de 2/3/4 no core9.", call. = FALSE)
}

asv_c9 <- matriz_amostras(ps_c9)
if (!identical(rownames(asv_c9), rownames(meta_c9))) {
  asv_c9 <- asv_c9[rownames(meta_c9), , drop = FALSE]
}
if (anyNA(asv_c9) || any(!is.finite(asv_c9)) || any(asv_c9 < 0) || any(rowSums(asv_c9) == 0)) {
  stop("otu_table do core9 contem valores invalidos ou amostra zerada.", call. = FALSE)
}
asv_pa_c9 <- as.data.frame((asv_c9 > 0) * 1L)

# Diagnostico da heterogeneidade da profundidade de sequenciamento. Esta etapa
# apenas registra a razao real maximo/minimo; nao rarefaz nem altera o objeto.
profundidade_core9 <- rowSums(asv_c9)
razao_profundidade_core9 <- max(profundidade_core9) / min(profundidade_core9)
diagnostico_profundidade <- data.frame(
  SampleID = names(profundidade_core9),
  SampleLabel = meta_c9[names(profundidade_core9), "SampleLabel"],
  BeeSpecies = as.character(meta_c9[names(profundidade_core9), "BeeSpecies"]),
  Reads = as.numeric(profundidade_core9),
  Profundidade_minima_core9 = min(profundidade_core9),
  Profundidade_maxima_core9 = max(profundidade_core9),
  Razao_max_min_core9 = razao_profundidade_core9,
  Limiar_diagnostico = LIMIAR_RAZAO_PROFUNDIDADE,
  Excede_limiar = razao_profundidade_core9 > LIMIAR_RAZAO_PROFUNDIDADE,
  stringsAsFactors = FALSE
)
write.csv(
  diagnostico_profundidade,
  file.path(fig_path, "diagnostico_profundidade_core9.csv"),
  row.names = FALSE
)
if (razao_profundidade_core9 > LIMIAR_RAZAO_PROFUNDIDADE) {
  log_msg(
    sprintf(
      paste0(
        "Profundidade core9: max/min = %.3f (> %.1f). ",
        "Riqueza observada e sensivel a profundidade; manter estimate_richness ",
        "apenas para continuidade e priorizar a estimativa iNEXT padronizada ",
        "por cobertura na interpretacao."
      ),
      razao_profundidade_core9, LIMIAR_RAZAO_PROFUNDIDADE
    ),
    "WARN"
  )
} else {
  log_msg(
    sprintf(
      "Profundidade core9: min=%d; max=%d; max/min=%.3f (<= %.1f).",
      min(profundidade_core9), max(profundidade_core9),
      razao_profundidade_core9, LIMIAR_RAZAO_PROFUNDIDADE
    ),
    "OK"
  )
}

taxa_tbl <- as.data.frame(tax_table(ps_c9), stringsAsFactors = FALSE, check.names = FALSE)
ranks_canonicos <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (!all(ranks_canonicos %in% colnames(taxa_tbl))) {
  stop("tax_table nao contem os sete ranks canonicos.", call. = FALSE)
}
if (!setequal(rownames(taxa_tbl), colnames(asv_c9))) {
  stop("tax_table e otu_table possuem universos de ASVs diferentes.", call. = FALSE)
}
taxa_tbl <- taxa_tbl[colnames(asv_c9), , drop = FALSE]
taxa_tbl$ASV_ID <- rownames(taxa_tbl)

# Variavel e rótulos derivados dos metadados
if (!"Nativo_Introduzido" %in% colnames(meta_c9)) {
  meta_c9$Nativo_Introduzido <- factor(ifelse(as.character(meta_c9$BeeSpecies) == "Melipona fasciculata",
                                              "Introduzida", "Nativa"),
                                       levels = c("Introduzida", "Nativa"))
}
meta_c9$SampleLabel <- gerar_rotulos_amostras(meta_c9)
mapa_samplelabel <- setNames(meta_c9$SampleLabel, meta_c9$SampleID)

cat(sprintf("ps_core9 : %d amostras | %d ASVs\n",
            nsamples(ps_c9), ntaxa(ps_c9)))

# Distâncias (carrega do Script 06 ou recomputa)
if (file.exists(arq_dist_bray)) {
  dist_bray <- tryCatch(
    readRDS(arq_dist_bray),
    error = function(e) {
      stop(
        "Falha ao ler core9_dist_bray_rel.rds: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  log_msg("dist_bray carregado (Script 06)", "OK")
} else {
  log_msg("Recomputando dist_bray (abundância relativa)", "WARN")
  ps_rel  <- transform_sample_counts(ps_c9, function(x) if (sum(x) > 0) x/sum(x) else x)
  dist_bray <- phyloseq::distance(ps_rel, method = "bray")
}
if (file.exists(arq_dist_jaccard)) {
  dist_jaccard <- tryCatch(
    readRDS(arq_dist_jaccard),
    error = function(e) {
      stop(
        "Falha ao ler core9_dist_jaccard_binary.rds: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  log_msg("dist_jaccard carregado (Script 06)", "OK")
} else {
  log_msg("Recomputando dist_jaccard (binário)", "WARN")
  dist_jaccard <- phyloseq::distance(ps_c9, method = "jaccard", binary = TRUE)
}
dist_bray <- alinhar_distancia(dist_bray, rownames(meta_c9), "dist_bray")
dist_jaccard <- alinhar_distancia(dist_jaccard, rownames(meta_c9), "dist_jaccard")
cat("\n")

###############################################################################
# 5. ESTATÍSTICAS CENTRAIS (calculadas uma vez; reutilizadas nos gráficos)
###############################################################################

cat("=== ESTATÍSTICAS CENTRAIS ===\n\n")

# ── 5.1 PERMANOVA (Bray-Curtis) ──────────────────────────────────────────────
# Fonte unica: Script 08. O Script 10 nao recalcula testes.
arq_perm_bray <- file.path(anal_path, "core9_permanova_global_exata.csv")
if (!file.exists(arq_perm_bray)) {
  stop("Resultado canonico ausente: ", arq_perm_bray,
       ". Execute o Script 08 antes do Script 10.", call. = FALSE)
}
perm_tab <- read.csv(arq_perm_bray, stringsAsFactors = FALSE, check.names = FALSE)
linha_bray <- perm_tab[perm_tab$Distancia == "bray_rel", , drop = FALSE]
if (nrow(linha_bray) != 1L || !all(c("R2", "F", "p_exato") %in% names(linha_bray))) {
  stop("core9_permanova_global_exata.csv invalido para Bray-Curtis.", call. = FALSE)
}
perm_R2 <- round(as.numeric(linha_bray$R2[[1L]]), 3)
perm_F <- round(as.numeric(linha_bray$F[[1L]]), 3)
perm_p <- as.numeric(linha_bray$p_exato[[1L]])
perm_label <- sprintf("PERMANOVA exata: F = %s | R² = %s | %s",
                      perm_F, perm_R2, fmt_p(perm_p))
cat(perm_label, "\n\n")

# ── 5.2 Diversidade alfa ─────────────────────────────────────────────────────
# Os indices sao calculados para desenho; os p-valores sao lidos do Script 08.
MEDIDAS_ALFA <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
alfa_df <- suppressWarnings(estimate_richness(ps_c9, measures = MEDIDAS_ALFA))
alfa_df$SampleID <- rownames(alfa_df)
alfa_df$SampleLabel <- meta_c9[rownames(alfa_df), "SampleLabel"]
alfa_df$BeeSpecies <- meta_c9[rownames(alfa_df), "BeeSpecies"]
alfa_df$Nativo_Introduzido <- meta_c9[rownames(alfa_df), "Nativo_Introduzido"]
alfa_df$Pielou <- with(alfa_df, ifelse(Observed > 1, Shannon / log(Observed), NA))

indices_kw <- c("Observed", "Shannon", "Simpson", "Pielou")
arq_kw_s08 <- file.path(anal_path, "core9_kruskal_wallis.csv")
if (!file.exists(arq_kw_s08)) {
  stop("Resultado canonico ausente: ", arq_kw_s08,
       ". Execute o Script 08 antes do Script 10.", call. = FALSE)
}
kw_s08 <- read.csv(arq_kw_s08, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Medida", "p_raw", "p_adj_BH", "Epsilon2") %in% names(kw_s08))) {
  stop("core9_kruskal_wallis.csv sem colunas canonicas.", call. = FALSE)
}
kw_p_vals <- setNames(rep(NA_real_, length(indices_kw)), indices_kw)
kw_p_adj_BH <- setNames(rep(NA_real_, length(indices_kw)), indices_kw)
comuns <- intersect(indices_kw, kw_s08$Medida)
kw_p_vals[comuns] <- kw_s08$p_raw[match(comuns, kw_s08$Medida)]
kw_p_adj_BH[comuns] <- kw_s08$p_adj_BH[match(comuns, kw_s08$Medida)]
log_msg("Estatisticas alfa reutilizadas do Script 08.", "OK")

# ── 5.3 iNEXT ────────────────────────────────────────────────────────────────
# Fonte unica: Script 06b. Ausencia ou incompatibilidade interrompe a execucao.
asv_inext <- as.matrix(otu_table(ps_c9))
if (!taxa_are_rows(ps_c9)) asv_inext <- t(asv_inext)
asv_inext <- asv_inext[rowSums(asv_inext) > 0, , drop = FALSE]
sample_ids_inext <- colnames(asv_inext)
rotulos_inext <- rotulo_amostra(sample_ids_inext, meta_c9)

# Paleta nomeada pelos SampleLabel usados pelo objeto iNEXT/ggiNEXT.
# Os nomes garantem correspondência determinística entre curvas e amostras.
if (anyNA(rotulos_inext) || any(!nzchar(rotulos_inext))) {
  stop(
    "SampleLabel ausente ou vazio para uma ou mais amostras do iNEXT.",
    call. = FALSE
  )
}

if (anyDuplicated(rotulos_inext)) {
  stop(
    "SampleLabel duplicado entre as amostras do iNEXT: ",
    paste(
      unique(rotulos_inext[duplicated(rotulos_inext)]),
      collapse = ", "
    ),
    call. = FALSE
  )
}

cores_inext <- setNames(
  rep_len(CORES_MULTI, length(rotulos_inext)),
  rotulos_inext
)

validar_resultado_inext_core9 <- function(obj) {
  if (!inherits(obj, "iNEXT")) return(FALSE)
  di <- tryCatch(as.data.frame(obj$DataInfo), error = function(e) NULL)
  if (is.null(di) || !all(c("n", "S.obs") %in% names(di))) return(FALSE)
  col_assem <- intersect(c("Assemblage", "assemblage", "site", "Site"), names(di))[1L]
  if (is.na(col_assem)) return(FALSE)
  idx <- match(rotulos_inext, as.character(di[[col_assem]]))
  if (anyNA(idx)) return(FALSE)
  all(as.numeric(di$n[idx]) == as.numeric(colSums(asv_inext))) &&
    all(as.numeric(di$S.obs[idx]) == as.numeric(colSums(asv_inext > 0)))
}
arq_inext_externo <- ctx$contracts[["inext_core9"]]

# Origem canônica do resultado iNEXT reutilizado pelo Script 10.
inext_fonte <- normalizePath(
  arq_inext_externo,
  winslash = "/",
  mustWork = FALSE
)

if (!file.exists(arq_inext_externo)) {
  stop("Resultado iNEXT ausente: ", arq_inext_externo,
       ". Execute o Script 06b antes do Script 10.", call. = FALSE)
}
rare_result <- readRDS(arq_inext_externo)
if (!validar_resultado_inext_core9(rare_result)) {
  stop("Resultado iNEXT externo incompatível com o phyloseq core9 atual.", call. = FALSE)
}
estimateD_cobertura <- NULL
cobertura_comum <- NA_real_
arq_estimateD_externo <- file.path(ctx$layout$stages$inext$root,
                                   "core9_iNEXT_estimateD_cobertura_comum.csv")
if (file.exists(arq_estimateD_externo)) {
  estimateD_cobertura <- read.csv(arq_estimateD_externo, check.names = FALSE)
  if ("SC" %in% names(estimateD_cobertura)) {
    z <- as.numeric(estimateD_cobertura$SC)
    z <- z[is.finite(z)]
    if (length(z)) cobertura_comum <- min(z)
  }
}
log_msg("iNEXT reutilizado do Script 06b; nenhum recalculo foi realizado.", "OK")

###############################################################################
# G01 — DIVERSIDADE ALFA
###############################################################################

cat("\nG01 - Alpha...\n")

alfa_long <- alfa_df |>
  tidyr::pivot_longer(cols = all_of(indices_kw),
                      names_to = "Indice", values_to = "Valor") |>
  dplyr::mutate(
    Indice_label = dplyr::case_when(
      Indice == "Observed" ~ paste0("Riqueza observada\n(KW exato: ", fmt_p(kw_p_adj_BH["Observed"]), ")"),
      Indice == "Shannon"  ~ paste0("Shannon (H')\n(KW exato: ", fmt_p(kw_p_adj_BH["Shannon"]), ")"),
      Indice == "Simpson"  ~ paste0("Simpson (1-D)\n(KW exato: ",  fmt_p(kw_p_adj_BH["Simpson"]),  ")"),
      Indice == "Pielou"   ~ paste0("Pielou (J)\n(KW exato: ",     fmt_p(kw_p_adj_BH["Pielou"]),   ")")
    ),
    Indice_label = factor(Indice_label, levels = unique(Indice_label))
  )

alfa_stats <- alfa_long |>
  dplyr::group_by(Indice_label, BeeSpecies) |>
  dplyr::summarise(media = mean(Valor, na.rm = TRUE),
                   ep    = sd(Valor, na.rm = TRUE) / sqrt(dplyr::n()),
                   .groups = "drop")

subtitle_alpha <- if (razao_profundidade_core9 > LIMIAR_RAZAO_PROFUNDIDADE) {
  sprintf(
    paste0(
      "Diamante = media | Barra = EP | p = KW exato com correcao BH | ",
      "profundidade max/min = %.2f; priorizar iNEXT por cobertura para riqueza"
    ),
    razao_profundidade_core9
  )
} else {
  sprintf(
    paste0(
      "Diamante = media | Barra = EP | p = KW exato com correcao BH | ",
      "profundidade max/min = %.2f"
    ),
    razao_profundidade_core9
  )
}

g01_alpha <- ggplot() +
  geom_linerange(data = alfa_stats,
                 aes(x = BeeSpecies, ymin = media - ep,
                     ymax = media + ep, color = BeeSpecies),
                 linewidth = 1.2, show.legend = FALSE) +
  geom_point(data = alfa_long,
             aes(x = BeeSpecies, y = Valor, fill = BeeSpecies),
             position = position_jitter(width = 0.1, seed = SEED),
             size = 3.5, shape = 21, stroke = 0.7,
             color = "white", alpha = 0.9) +
  geom_point(data = alfa_stats,
             aes(x = BeeSpecies, y = media, color = BeeSpecies),
             size = 5, shape = 18, show.legend = FALSE) +
  facet_wrap(~ Indice_label, scales = "free_y", nrow = 2, ncol = 2) +
  scale_fill_manual(values  = CORES_ESP, labels = LABELS_ESP,
                    name = "Esp\u00e9cie") +
  scale_color_manual(values = CORES_ESP) +
  scale_x_discrete(labels = LABELS_ESP) +
  TEMA_PUB +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1,
                                     face = "italic", size = 9),
        legend.position = "none") +
  labs(title    = "A \u2014 Diversidade Alfa",
       subtitle = subtitle_alpha,
       x = NULL, y = "Valor do \u00cdndice")


# G01b — Diversidade alfa por status nativo/introduzido
alfa_status_stats <- alfa_long |>
  dplyr::group_by(Indice_label, Nativo_Introduzido) |>
  dplyr::summarise(media = mean(Valor, na.rm = TRUE),
                   ep = sd(Valor, na.rm = TRUE) / sqrt(dplyr::n()),
                   .groups = "drop")

g01b_alpha_status <- ggplot() +
  geom_linerange(data = alfa_status_stats,
                 aes(x = Nativo_Introduzido, ymin = media - ep,
                     ymax = media + ep, color = Nativo_Introduzido),
                 linewidth = 1.2, show.legend = FALSE) +
  geom_point(data = alfa_long,
             aes(x = Nativo_Introduzido, y = Valor, fill = Nativo_Introduzido),
             position = position_jitter(width = 0.1, seed = SEED),
             size = 3.5, shape = 21, stroke = 0.7,
             color = "white", alpha = 0.9) +
  geom_point(data = alfa_status_stats,
             aes(x = Nativo_Introduzido, y = media, color = Nativo_Introduzido),
             size = 5, shape = 18, show.legend = FALSE) +
  facet_wrap(~ Indice_label, scales = "free_y", nrow = 2, ncol = 2) +
  scale_fill_manual(values = CORES_STATUS, name = "Origem") +
  scale_color_manual(values = CORES_STATUS) +
  TEMA_PUB +
  labs(title = "A2 — Diversidade Alfa por origem",
       subtitle = paste(
         "Descritivo/exploratório:",
         "origem é derivada de BeeSpecies e não representa efeito independente"
       ),
       x = NULL, y = "Valor do Índice")

###############################################################################
# G02 / G03 — PCoA (PC1×PC2 e PC1×PC3)
###############################################################################

cat("G02/G03 - PCoA...\n")

ord_pcoa <- ordinate(ps_c9, method = "PCoA", distance = dist_bray)
eig      <- ord_pcoa$values$Relative_eig * 100

pcoa_df           <- as.data.frame(ord_pcoa$vectors[, 1:3])
colnames(pcoa_df) <- c("PC1", "PC2", "PC3")
pcoa_df$BeeSpecies <- meta_c9[rownames(pcoa_df), "BeeSpecies"]
pcoa_df$Nativo_Introduzido <- meta_c9[rownames(pcoa_df), "Nativo_Introduzido"]
pcoa_df$SampleID   <- rownames(pcoa_df)
pcoa_df$SampleLabel <- meta_c9[rownames(pcoa_df), "SampleLabel"]

# Centróides e spiders para PC1×PC2 e PC1×PC3
cen12 <- pcoa_df |> dplyr::group_by(BeeSpecies) |>
  dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
cen13 <- pcoa_df |> dplyr::group_by(BeeSpecies) |>
  dplyr::summarise(PC1 = mean(PC1), PC3 = mean(PC3), .groups = "drop")

seg12 <- dplyr::left_join(pcoa_df, cen12, by = "BeeSpecies",
                           suffix = c("","_c")) |>
  dplyr::rename(x = PC1, y = PC2, xend = PC1_c, yend = PC2_c)
seg13 <- dplyr::left_join(pcoa_df, cen13, by = "BeeSpecies",
                           suffix = c("","_c")) |>
  dplyr::rename(x = PC1, y = PC3, xend = PC1_c, yend = PC3_c)

# Polígonos convexos exigem >=3 pontos não colineares. Grupos com n=2
# permanecem representados por pontos e spiders, sem polígono degenerado.
calcular_hull_seguro <- function(df, x, y) {
  df |>
    dplyr::group_by(BeeSpecies) |>
    dplyr::filter(
      dplyr::n() >= 3L,
      dplyr::n_distinct(.data[[x]]) >= 2L,
      dplyr::n_distinct(.data[[y]]) >= 2L
    ) |>
    dplyr::slice(grDevices::chull(.data[[x]], .data[[y]])) |>
    dplyr::ungroup()
}
hull12 <- calcular_hull_seguro(pcoa_df, "PC1", "PC2")
hull13 <- calcular_hull_seguro(pcoa_df, "PC1", "PC3")

# Função interna: monta um plot PCoA
fazer_pcoa <- function(df, seg, hull, cen, xv, yv, xcv, ycv,
                        xlab, ylab, titulo) {
  ggplot(df, aes(x = .data[[xv]], y = .data[[yv]],
                  color = BeeSpecies, shape = BeeSpecies)) +
    geom_polygon(data = hull,
                 aes(x = .data[[xv]], y = .data[[yv]],
                     fill = BeeSpecies, group = BeeSpecies),
                 alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_segment(data = seg,
                 aes(x = xend, y = yend, xend = x, yend = y,
                     color = BeeSpecies),
                 linewidth = 0.6, alpha = 0.55, show.legend = FALSE) +
    geom_point(size = 5, alpha = 0.95, stroke = 1.2) +
    geom_point(data = cen,
               aes(x = .data[[xcv]], y = .data[[ycv]], color = BeeSpecies),
               size = 3.5, shape = 3, stroke = 1.8, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(label = SampleLabel),
                              size = 3, fontface = "bold",
                              box.padding = 0.4, point.padding = 0.3,
                              segment.color = "grey60", segment.size = 0.4,
                              show.legend = FALSE) +
    scale_color_manual(values = CORES_ESP, labels = LABELS_ESP,
                       name = "Esp\u00e9cie") +
    scale_fill_manual( values = CORES_ESP, labels = LABELS_ESP,
                       name = "Esp\u00e9cie") +
    scale_shape_manual(values = c(16, 17, 15), labels = LABELS_ESP,
                       name = "Esp\u00e9cie") +
    TEMA_PUB +
    theme(legend.text = element_text(face = "italic")) +
    labs(title = titulo, subtitle = perm_label,
         x = xlab, y = ylab)
}

g02_pcoa12 <- fazer_pcoa(pcoa_df, seg12, hull12, cen12,
  "PC1","PC2","PC1","PC2",
  sprintf("PCoA1 [%.1f%%]", eig[1]),
  sprintf("PCoA2 [%.1f%%]", eig[2]),
  "B \u2014 PCoA Bray-Curtis (PC1 \u00d7 PC2)")

g03_pcoa13 <- fazer_pcoa(pcoa_df, seg13, hull13, cen13,
  "PC1","PC3","PC1","PC3",
  sprintf("PCoA1 [%.1f%%]", eig[1]),
  sprintf("PCoA3 [%.1f%%]", eig[3]),
  "B2 \u2014 PCoA Bray-Curtis (PC1 \u00d7 PC3)")


# G02b — PCoA por status nativo/introduzido
cen_status <- pcoa_df |>
  dplyr::group_by(Nativo_Introduzido) |>
  dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
seg_status <- dplyr::left_join(pcoa_df, cen_status, by = "Nativo_Introduzido",
                               suffix = c("", "_c")) |>
  dplyr::rename(x = PC1, y = PC2, xend = PC1_c, yend = PC2_c)

g02b_pcoa_status <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                                        color = Nativo_Introduzido,
                                        shape = Nativo_Introduzido)) +
  geom_segment(data = seg_status,
               aes(x = xend, y = yend, xend = x, yend = y,
                   color = Nativo_Introduzido),
               linewidth = 0.6, alpha = 0.55, show.legend = FALSE) +
  geom_point(size = 5, alpha = 0.95, stroke = 1.2) +
  geom_point(data = cen_status,
             aes(x = PC1, y = PC2, color = Nativo_Introduzido),
             size = 3.5, shape = 3, stroke = 1.8, show.legend = FALSE) +
  ggrepel::geom_text_repel(aes(label = SampleLabel), size = 3,
                            fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = CORES_STATUS, name = "Origem") +
  scale_shape_manual(values = c("Introduzida" = 17, "Nativa" = 16), name = "Origem") +
  TEMA_PUB +
  labs(title = "B3 — PCoA Bray-Curtis por origem",
       subtitle = paste(
         "Comparação exploratória; origem é derivada de BeeSpecies",
         "e não pode ser interpretada como efeito independente"
       ),
       x = sprintf("PCoA1 [%.1f%%]", eig[1]),
       y = sprintf("PCoA2 [%.1f%%]", eig[2]))

###############################################################################
# G04 — NMDS
###############################################################################

cat("G04 - NMDS...\n")

set.seed(SEED)
ord_nmds <- ordinate(ps_c9, method = "NMDS", distance = dist_bray,
                      trymax = 200, trace = FALSE)
cat(sprintf("  Stress: %.4f\n", ord_nmds$stress))

nmds_df           <- as.data.frame(ord_nmds$points)
colnames(nmds_df) <- c("NMDS1","NMDS2")
nmds_df$BeeSpecies <- meta_c9[rownames(nmds_df), "BeeSpecies"]
nmds_df$SampleLabel <- meta_c9[rownames(nmds_df), "SampleLabel"]

n_grupo_nmds <- table(nmds_df$BeeSpecies)
grupos_elipse <- names(n_grupo_nmds[n_grupo_nmds >= 3])

g04_nmds <- ggplot(nmds_df,
                   aes(x = NMDS1, y = NMDS2,
                       color = BeeSpecies, shape = BeeSpecies)) +
  geom_point(size = 5, alpha = 0.9, stroke = 1) +
  ggrepel::geom_text_repel(
    aes(label = SampleLabel), size = 3, fontface = "bold",
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "grey60", segment.size = 0.4,
    show.legend = FALSE
  ) +
  {if (length(grupos_elipse) > 0)
    stat_ellipse(data = nmds_df[nmds_df$BeeSpecies %in% grupos_elipse, ],
                 aes(group = BeeSpecies), level = 0.90,
                 linetype = 2, linewidth = 0.8, show.legend = FALSE)
  } +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           size = 3.5,
           label = paste0("Stress = ", round(ord_nmds$stress, 3),
                          "  |  ", perm_label)) +
  scale_color_manual(values = CORES_ESP, labels = LABELS_ESP,
                     name = "Esp\u00e9cie") +
  scale_shape_manual(values = c(16, 17, 15), labels = LABELS_ESP,
                     name = "Esp\u00e9cie") +
  TEMA_PUB +
  theme(legend.text = element_text(face = "italic")) +
  labs(title    = "C \u2014 NMDS Bray-Curtis",
       subtitle = "Elipses de confian\u00e7a 90% | PERMANOVA sob gr\u00e1fico",
       x = "NMDS1", y = "NMDS2")

###############################################################################
# G05 / G06 / G07 — COMPOSIÇÃO POR FILO E GÊNERO
###############################################################################

cat("G05/G06/G07 - Composição...\n")

# ── G05: Filo ──
ps_phy     <- tax_glom(ps_c9, taxrank = "Phylum", NArm = FALSE)
ps_phy_rel <- transform_sample_counts(ps_phy, function(x) x / sum(x) * 100)
filo_df    <- psmelt(ps_phy_rel) |>
  dplyr::mutate(
    Especie_label = LABELS_ESP[BeeSpecies],
    SampleLabel_x = rotulo_amostra(Sample),
    Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum)
  )

filos_uniq <- sort(unique(filo_df$Phylum))
cores_filo <- setNames(
  colorRampPalette(CORES_MULTI)(length(filos_uniq)), filos_uniq)

g05_filo <- ggplot(filo_df, aes(x = SampleLabel_x, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack",
           width = 0.85, color = "white", linewidth = 0.12) +
  facet_wrap(~ Especie_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cores_filo, name = "Filo") +
  scale_y_continuous(limits = c(0, 101), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  TEMA_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.text = element_text(size = 8, face = "italic"),
        panel.grid.major.x = element_blank()) +
  labs(title = "D \u2014 Composi\u00e7\u00e3o Bacteriana: Filo",
       x = NULL, y = "Abund\u00e2ncia Relativa (%)")

# ── G06 / G07: Gênero ──
ps_gen_noNA <- subset_taxa(ps_c9, !is.na(Genus) & Genus != "")
ps_gen      <- tax_glom(ps_gen_noNA, taxrank = "Genus", NArm = FALSE)
ps_gen_rel  <- transform_sample_counts(ps_gen, function(x) x / sum(x) * 100)

genus_df_all <- psmelt(ps_gen_rel) |>
  dplyr::filter(!is.na(Genus), Genus != "") |>
  dplyr::mutate(
    Especie_label = LABELS_ESP[BeeSpecies],
    SampleLabel_x = rotulo_amostra(Sample)
  )

ordem_gen <- genus_df_all |>
  dplyr::group_by(Genus) |>
  dplyr::summarise(media = mean(Abundance), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(media)) |>
  dplyr::pull(Genus)

genus_df_all$Genus <- factor(genus_df_all$Genus, levels = rev(ordem_gen))
n_gen <- length(levels(genus_df_all$Genus))
cores_gen <- setNames(rep_len(CORES_MULTI, n_gen), levels(genus_df_all$Genus))

g06_genus <- ggplot(genus_df_all,
                    aes(x = SampleLabel_x, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack",
           width = 0.85, color = "white", linewidth = 0.12) +
  facet_wrap(~ Especie_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cores_gen, name = "G\u00eanero") +
  scale_y_continuous(limits = c(0, 101), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  TEMA_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.key.size = grid::unit(0.38, "cm"),
        legend.text = element_text(size = 7, face = "italic"),
        panel.grid.major.x = element_blank()) +
  labs(title    = paste0("E \u2014 Composi\u00e7\u00e3o Bacteriana: G\u00eaneros (",
                          n_gen, " identificados)"),
       subtitle = "Abund\u00e2ncia relativa (%)",
       x = NULL, y = "Abund\u00e2ncia Relativa (%)")

# G07: Top 15 gêneros
top15_gen  <- ordem_gen[seq_len(min(15L, length(ordem_gen)))]
genus_df15 <- dplyr::filter(genus_df_all, Genus %in% top15_gen)
genus_df15$Genus <- factor(genus_df15$Genus, levels = rev(top15_gen))
cores_15 <- setNames(rep_len(CORES_MULTI, 15), rev(top15_gen))

g07_top15 <- ggplot(genus_df15,
                    aes(x = SampleLabel_x, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack",
           width = 0.9, color = "white", linewidth = 0.12) +
  facet_wrap(~ Especie_label, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = cores_15, name = "G\u00eanero") +
  scale_y_continuous(limits = c(0, 101), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  TEMA_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.key.size = grid::unit(0.42, "cm"),
        legend.text = element_text(size = 8, face = "italic"),
        panel.grid.major.x = element_blank()) +
  labs(title    = "E2 \u2014 Composi\u00e7\u00e3o Bacteriana: Top 15 G\u00eaneros",
       subtitle = "Os 15 g\u00eaneros de maior abundacia m\u00e9dia",
       x = NULL, y = "Abund\u00e2ncia Relativa (%)")

###############################################################################
# G08 — iNEXT
###############################################################################

cat("G08 - iNEXT...\n")

g08_inext <- ggiNEXT(rare_result, type = 1) +
  facet_wrap(~ Order.q, scales = "free_y",
             labeller = labeller(Order.q = c(
               "0" = "q = 0  Riqueza",
               "1" = "q = 1  Diversidade Shannon",
               "2" = "q = 2  Diversidade Simpson"))) +
  scale_color_manual(values = cores_inext, name = "Amostra (SampleLabel)") +
  scale_fill_manual( values = cores_inext, name = "Amostra (SampleLabel)") +
  TEMA_PUB +
  theme(legend.key.size = grid::unit(0.38, "cm"),
        legend.text = element_text(size = 8)) +
  labs(title    = "F \u2014 Curvas de Rarefa\u00e7\u00e3o e Extrapola\u00e7\u00e3o (iNEXT)",
       subtitle = "Linha s\u00f3lida = interpola\u00e7\u00e3o | Tracejada = extrapola\u00e7\u00e3o",
       x = "N\u00famero de reads", y = "Diversidade")

###############################################################################
# G09 — WHITTAKER (partição beta — recalcula inline)
###############################################################################

cat("G09 - Whittaker...\n")

beta_g <- betapart::beta.multi(asv_pa_c9, index.family = "sorensen")
especies_c9 <- sort(unique(as.character(meta_c9$BeeSpecies)))

whitt_sp <- dplyr::bind_rows(lapply(especies_c9, function(sp) {
  am <- rownames(meta_c9[meta_c9$BeeSpecies == sp, ])
  if (length(am) < 2) return(NULL)
  b <- betapart::beta.multi(asv_pa_c9[am, , drop = FALSE],
                             index.family = "sorensen")
  data.frame(Nivel = sp,
             Beta_total = round(b$beta.SOR, 3),
             Turnover   = round(b$beta.SIM, 3),
             Nestedness = round(b$beta.SNE, 3),
             stringsAsFactors = FALSE)
}))

whitt_df <- rbind(
  data.frame(Nivel = "Global",
             Beta_total = round(beta_g$beta.SOR, 3),
             Turnover   = round(beta_g$beta.SIM, 3),
             Nestedness = round(beta_g$beta.SNE, 3),
             stringsAsFactors = FALSE),
  whitt_sp[!is.na(whitt_sp$Beta_total), ]
)

whitt_df$Nivel <- factor(whitt_df$Nivel,
                          levels = c("Global", especies_c9))
whitt_long <- tidyr::pivot_longer(whitt_df,
                                   cols = c("Turnover","Nestedness"),
                                   names_to = "Componente", values_to = "Valor")
whitt_long$Nivel <- factor(whitt_long$Nivel, levels = levels(whitt_df$Nivel))

g09_whitt <- ggplot(whitt_long,
                    aes(x = Nivel, y = Valor, fill = Componente)) +
  geom_col(position = "stack", width = 0.55,
           color = "white", linewidth = 0.5) +
  geom_text(data = whitt_df,
            aes(x = Nivel, y = Beta_total + 0.03,
                label = round(Beta_total, 3)),
            size = 4.2, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(
    values = c("Turnover" = "#E31A1C", "Nestedness" = "#1F78B4"),
    name = "Componente") +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2),
                     expand = c(0, 0)) +
  scale_x_discrete(labels = c(
    "Global"               = "Global",
    "Melipona fasciculata" = "M. fasciculata",
    "Melipona scutellaris" = "M. scutellaris",
    "Melipona subnitida"   = "M. subnitida")) +
  TEMA_PUB +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")) +
  labs(title    = "G \u2014 Parti\u00e7\u00e3o da Diversidade Beta (S\u00f8rensen)",
       subtitle = sprintf("SOR = %.3f | SIM = %.3f | SNE = %.3f",
                           whitt_df$Beta_total[whitt_df$Nivel == "Global"],
                           whitt_df$Turnover[whitt_df$Nivel   == "Global"],
                           whitt_df$Nestedness[whitt_df$Nivel == "Global"]),
       x = NULL, y = "Diversidade Beta")

###############################################################################
# G10 — HEATMAP TOP 20 GÊNEROS
###############################################################################

cat("G10 - Heatmap...\n")

ps_g2   <- subset_taxa(ps_c9, !is.na(Genus) & Genus != "")
ps_gen2 <- tax_glom(ps_g2, taxrank = "Genus", NArm = FALSE)
ps_g2r  <- transform_sample_counts(ps_gen2, function(x) x / sum(x) * 100)

top20g    <- names(sort(taxa_sums(ps_g2r),
                        decreasing = TRUE)[seq_len(min(20L, ntaxa(ps_g2r)))])
ps_top20g <- prune_taxa(top20g, ps_g2r)

heat_df <- psmelt(ps_top20g) |>
  dplyr::mutate(
    Especie_label = LABELS_ESP[BeeSpecies],
    SampleLabel_x = rotulo_amostra(Sample),
    pct_lab = ifelse(Abundance >= 0.5, paste0(round(Abundance, 1), "%"), ""),
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)
  )

ord_heat <- heat_df |>
  dplyr::group_by(Genus) |>
  dplyr::summarise(media = mean(Abundance), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(media))
heat_df$Genus <- factor(heat_df$Genus, levels = rev(ord_heat$Genus))

g10_heat <- ggplot(heat_df,
                   aes(x = SampleLabel_x, y = Genus, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = pct_lab,
                color = ifelse(Abundance > 30, "white", "black")),
            size = 2.6, fontface = "bold", show.legend = FALSE) +
  facet_grid(~ Especie_label, scales = "free_x", space = "free") +
  scale_fill_gradientn(
    colours  = c("#003f5c","#2f4b7c","#665191","#a05195",
                 "#d45087","#f95d6a","#ff7c43","#ffa600"),
    values   = scales::rescale(c(0, 5, 15, 30, 50, 100)),
    name     = "Abund.\nRelativa (%)",
    na.value = "grey95") +
  scale_color_identity() +
  TEMA_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(face = "italic", size = 9),
        axis.title  = element_blank(),
        panel.grid  = element_blank()) +
  labs(title    = "H \u2014 Heatmap Top 20 G\u00eaneros",
       subtitle = "Percentual exibido nas c\u00e9lulas \u2265 0.5%")

###############################################################################
# G11 — DIAGRAMA DE VENN
###############################################################################

cat("G11 - Venn...\n")

asv_por_esp <- setNames(
  lapply(especies_c9, function(sp) {
    am <- rownames(meta_c9[meta_c9$BeeSpecies == sp, ])
    colnames(asv_c9)[colSums(asv_c9[am, , drop = FALSE]) > 0]
  }),
  paste0("M. ", sub("Melipona ", "", especies_c9))
)

g11_venn <- ggVennDiagram::ggVennDiagram(
  asv_por_esp, label_alpha = 0, edge_size = 1.2) +
  scale_fill_gradientn(
    colours = c("#2f4b7c","#a05195","#d45087","#f95d6a"),
    name    = "N de ASVs") +
  TEMA_PUB +
  theme(panel.border = element_blank(),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid   = element_blank()) +
  labs(title    = "I \u2014 Diagrama de Venn: ASVs compartilhadas",
       subtitle = "Compara\u00e7\u00e3o da diversidade bacteriana entre as tr\u00eas esp\u00e9cies")

###############################################################################
# G12 — VOLCANO DESeq2 (3 comparações)
###############################################################################

cat("G12 - Volcano DESeq2...\n")

comps_deseq <- c("scutellaris_vs_fasciculata",
                 "subnitida_vs_fasciculata",
                 "subnitida_vs_scutellaris")

volcano_list <- lapply(comps_deseq, function(nm) {
  arq <- file.path(deseq_path, paste0("core9_completo_", nm, ".csv"))
  if (!file.exists(arq)) {
    log_msg(paste("Arquivo DESeq2 ausente:", basename(arq)), "WARN")
    return(NULL)
  }
  res <- read.csv(arq, stringsAsFactors = FALSE)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
  res$padj_plot <- pmax(res$padj, .Machine$double.xmin)
  res$Dir <- factor(
    dplyr::case_when(
      res$padj < ALPHA & res$log2FoldChange >  LFC_THRESHOLD ~ "Up",
      res$padj < ALPHA & res$log2FoldChange < -LFC_THRESHOLD ~ "Down",
      TRUE ~ "NS"),
    levels = c("Up","Down","NS"))
  n_up   <- sum(res$Dir == "Up")
  n_down <- sum(res$Dir == "Down")
  titulo <- gsub("_"," vs ", nm)
  top5 <- res |> dplyr::filter(Dir != "NS") |>
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) |>
    dplyr::slice_head(n = 5) |>
    dplyr::mutate(rotulo = ifelse(!is.na(Genus) & Genus != "",
                                   Genus, ASV_ID))
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj_plot), color = Dir)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed",
               color = "grey55", linewidth = 0.7) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
               color = "grey55", linewidth = 0.7) +
    ggrepel::geom_text_repel(data = top5,
                              aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)),
                                  label = rotulo),
                              size = 2.8, fontface = "italic",
                              box.padding = 0.4, show.legend = FALSE,
                              inherit.aes = FALSE) +
    scale_color_manual(
      values = c("Up" = "#E31A1C", "Down" = "#1F78B4", "NS" = "grey72"),
      name = "Regula\u00e7\u00e3o",
      labels = c("Up"   = paste0("Up (n=",   n_up,   ")"),
                 "Down" = paste0("Down (n=", n_down, ")"),
                 "NS"   = "NS")) +
    TEMA_PUB +
    labs(title    = titulo,
         subtitle = sprintf("|log2FC| > %.1f e padj < %.2f | Up = %d | Down = %d",
                             LFC_THRESHOLD, ALPHA, n_up, n_down),
         x = "log2 Fold Change", y = "-log10(padj)")
})
names(volcano_list) <- comps_deseq
volcano_list <- volcano_list[!sapply(volcano_list, is.null)]

g12_volcano <- if (length(volcano_list) > 0) {
  patchwork::wrap_plots(volcano_list, ncol = 3) +
    patchwork::plot_annotation(
      title = "J \u2014 Volcano DESeq2 \u2014 Abund\u00e2ncia Diferencial",
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
} else {
  log_msg("Arquivos DESeq2 não encontrados; G12 omitido.", "WARN")
  NULL
}


# G12b — Volcano DESeq2 nativo/introduzido
arq_status_deseq <- file.path(deseq_path, "nativo_introduzido_prev2",
                              "core9_completo_nativa_vs_introduzida.csv")
g12b_volcano_status <- if (file.exists(arq_status_deseq)) {
  res <- read.csv(arq_status_deseq, stringsAsFactors = FALSE)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
  res$padj_plot <- pmax(res$padj, .Machine$double.xmin)
  res$Dir <- factor(dplyr::case_when(
    res$padj < ALPHA & res$log2FoldChange >  LFC_THRESHOLD ~ "Nativa ↑",
    res$padj < ALPHA & res$log2FoldChange < -LFC_THRESHOLD ~ "Introduzida ↑",
    TRUE ~ "NS"), levels = c("Nativa ↑", "Introduzida ↑", "NS"))
  top5 <- res |>
    dplyr::filter(Dir != "NS") |>
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) |>
    dplyr::slice_head(n = 5) |>
    dplyr::mutate(rotulo = ifelse(!is.na(Genus) & Genus != "", Genus, ASV_ID))
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj_plot), color = Dir)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "grey55") +
    ggrepel::geom_text_repel(data = top5,
                              aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)), label = rotulo),
                              size = 2.8, fontface = "italic", show.legend = FALSE,
                              inherit.aes = FALSE) +
    scale_color_manual(values = c("Nativa ↑" = "#2f4b7c", "Introduzida ↑" = "#FF6361", "NS" = "grey72"),
                       name = "Direção") +
    TEMA_PUB +
    labs(title = "J2 — DESeq2: Nativas vs Introduzida",
         subtitle = "Contraste composto e confundido com espécie; análise exploratória",
         x = "log2 Fold Change", y = "-log10(padj)")
} else {
  log_msg("Arquivo DESeq2 nativo/introduzido ausente; G12b omitido.", "WARN")
  NULL
}

###############################################################################
# G13 — SIMPER top 10 ASVs por comparação
###############################################################################

cat("G13 - SIMPER...\n")

arq_simper <- file.path(anal_path, "core9_simper_top10_por_par.csv")
g13_simper <- if (file.exists(arq_simper)) {
  simper_all <- read.csv(arq_simper, stringsAsFactors = FALSE)
  simper_all$Comparacao <- gsub("_", " vs ", simper_all$Comparacao)
  simper_all$Label <- ifelse(!is.na(simper_all$Genus) & simper_all$Genus != "",
                             simper_all$Genus, simper_all$ASV_ID)
  ggplot(simper_all,
         aes(x = reorder(Label, Contrib), y = Contrib * 100, fill = Comparacao)) +
    geom_col(position = "dodge", width = 0.7,
             color = "white", linewidth = 0.2) +
    coord_flip() +
    TEMA_PUB +
    labs(title = "K — SIMPER: Contribuição por Táxon",
         subtitle = "Top 10 ASVs por par de espécies; resultado lido do Script 08",
         x = "Táxon (Gênero)", y = "Contribuição média (%)")
} else {
  log_msg("core9_simper_top10_por_par.csv ausente; execute Script 08.", "WARN")
  NULL
}

###############################################################################
# G14 — IndVal ASVs indicadoras (FDR-BH < 0.05)
###############################################################################

cat("G14 - IndVal...\n")

arq_indval <- file.path(anal_path, "core9_indval_significativas.csv")
g14_indval <- if (file.exists(arq_indval)) {
  ind_sig <- read.csv(arq_indval, stringsAsFactors = FALSE)
  if (nrow(ind_sig) > 0) {
    if (!all(c("ASV_ID", "stat", "Grupo") %in% colnames(ind_sig))) {
      stop("core9_indval_significativas.csv sem colunas obrigatorias.", call. = FALSE)
    }
    col_p_ind <- if ("p_adj_BH" %in% colnames(ind_sig)) "p_adj_BH" else "p.value"
    if (!col_p_ind %in% colnames(ind_sig)) {
      stop("IndVal sem p_adj_BH ou p.value.", call. = FALSE)
    }
    ind_sig$Label <- ifelse(!is.na(ind_sig$Genus) & ind_sig$Genus != "",
                            ind_sig$Genus, ind_sig$ASV_ID)
    ind_sig$p_plot <- as.numeric(ind_sig[[col_p_ind]])
    ggplot(ind_sig, aes(x = reorder(Label, stat), y = stat, fill = Grupo)) +
      geom_col(width = 0.65, color = "white") +
      geom_text(aes(label = paste0(ifelse(col_p_ind == "p_adj_BH", "FDR=", "p="),
                                   signif(p_plot, 3))),
                hjust = -0.1, size = 3, color = "grey20") +
      coord_flip() +
      TEMA_PUB +
      labs(title = "L — IndVal: ASVs Indicadoras (FDR-BH < 0.05)",
           subtitle = "Resultado lido do Script 08",
           x = "Táxon (Gênero)", y = "Estatística IndVal")
  } else {
    log_msg("Nenhuma ASV indicadora significativa no Script 08; G14 omitido.", "INFO")
    NULL
  }
} else {
  log_msg("core9_indval_significativas.csv ausente; execute Script 08.", "WARN")
  NULL
}

###############################################################################
# SALVAR — PDF CONSOLIDADO + PNGs INDIVIDUAIS
###############################################################################

cat("\n=== SALVANDO FIGURAS ===\n\n")

graficos <- list(
  list(id = "G01_alpha",        g = g01_alpha,   w = 16, h = 10),
  list(id = "G01b_alpha_status", g = g01b_alpha_status, w = 12, h = 8),
  list(id = "G02_pcoa_PC1xPC2", g = g02_pcoa12,  w = 10, h =  7),
  list(id = "G02b_pcoa_status",  g = g02b_pcoa_status, w = 10, h = 7),
  list(id = "G03_pcoa_PC1xPC3", g = g03_pcoa13,  w = 10, h =  7),
  list(id = "G04_nmds",         g = g04_nmds,    w = 10, h =  7),
  list(id = "G05_filo",         g = g05_filo,    w = 14, h =  7),
  list(id = "G06_genus_todos",  g = g06_genus,   w = 16, h =  8),
  list(id = "G07_genus_top15",  g = g07_top15,   w = 12, h =  8),
  list(id = "G08_inext",        g = g08_inext,   w = 14, h =  6),
  list(id = "G09_whittaker",    g = g09_whitt,   w =  9, h =  7),
  list(id = "G10_heatmap",      g = g10_heat,    w = 12, h = 10),
  list(id = "G11_venn",         g = g11_venn,    w =  9, h =  8)
)

# Adicionar graficos condicionais. Objetos NULL nunca sao enviados ao ggsave,
# pois ggsave(plot = NULL) pode reutilizar silenciosamente o ultimo grafico.
if (!is.null(g13_simper))
  graficos <- c(graficos, list(list(id = "G13_simper", g = g13_simper, w = 14, h = 8)))
if (!is.null(g12_volcano))
  graficos <- c(graficos, list(list(id = "G12_volcano", g = g12_volcano, w = 18, h = 7)))
if (!is.null(g12b_volcano_status))
  graficos <- c(graficos, list(list(id = "G12b_volcano_status", g = g12b_volcano_status, w = 10, h = 7)))
if (!is.null(g14_indval))
  graficos <- c(graficos, list(list(id = "G14_indval",  g = g14_indval,  w = 10, h = 7)))

# PDF consolidado e PNGs individuais, com auditoria explícita de falhas.
erros_graficos <- list()
registrar_erro_grafico <- function(id, formato, e) {
  erros_graficos[[length(erros_graficos) + 1L]] <<- data.frame(
    Grafico = id,
    Formato = formato,
    Erro = conditionMessage(e),
    stringsAsFactors = FALSE
  )
  log_msg(paste(id, formato, "falhou:", conditionMessage(e)), "ERRO")
  invisible(NULL)
}

cat("  Gerando PDF consolidado...\n")
arq_pdf_completo <- file.path(fig_path, "figuras_completas.pdf")
grDevices::pdf(arq_pdf_completo, width = 16, height = 10, onefile = TRUE)
pdf_aberto <- TRUE
tryCatch(
  {
    for (item in graficos) {
      tryCatch(
        print(item$g),
        error = function(e) registrar_erro_grafico(item$id, "PDF_consolidado", e)
      )
    }
  },
  finally = {
    if (isTRUE(pdf_aberto) && grDevices::dev.cur() > 1L) grDevices::dev.off()
  }
)
if (!file.exists(arq_pdf_completo) || file.size(arq_pdf_completo) == 0L) {
  stop("Falha ao gerar figuras_completas.pdf.", call. = FALSE)
}
log_msg("figuras_completas.pdf salvo", "SAVE")

cat("  Gerando PNGs individuais...\n")
for (item in graficos) {
  ok <- tryCatch(
    {
      salvar_plot(item$g, paste0(item$id, ".png"), fig_path,
                  width = item$w, height = item$h)
      TRUE
    },
    error = function(e) {
      registrar_erro_grafico(item$id, "PNG", e)
      FALSE
    }
  )
  if (isTRUE(ok)) cat("  OK:", item$id, "\n")
}

erros_graficos_df <- if (length(erros_graficos)) {
  dplyr::bind_rows(erros_graficos)
} else {
  data.frame(Grafico = character(), Formato = character(), Erro = character())
}
write.csv(
  erros_graficos_df,
  file.path(fig_path, "auditoria_erros_graficos_script10.csv"),
  row.names = FALSE,
  quote = TRUE
)

###############################################################################
# PLUS10 — EXECUCAO GRAFICA PROTEGIDA
###############################################################################

graficos_plus10 <- list()
erros_p10_df <- data.frame(Grafico = character(), Formato = character(), Erro = character())
ps_p10 <- NULL
asv_p10 <- NULL
fig_path_p10 <- file.path(fig_path, "plus10_exploratorio")
plus10_figuras_ok <- tryCatch(
  {
    ###############################################################################
    # BLOCO PLUS10 - FIGURAS EXPLORATORIAS SEPARADAS
    #
    # O plus10 inclui uma unica amostra da segunda corrida. As figuras e os testes
    # associados sao salvos em diretorio proprio e nao substituem o core9.
    ###############################################################################

    cat("\n=== FIGURAS PLUS10 - SENSIBILIDADE EXPLORATORIA ===\n\n")

    fig_path_p10 <- file.path(fig_path, "plus10_exploratorio")
    dir.create(fig_path_p10, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(fig_path_p10)) {
      stop("Falha ao criar diretorio de figuras plus10: ", fig_path_p10, call. = FALSE)
    }

    arq_plus10 <- ctx$contracts[["phyloseq_plus10"]]
    if (!file.exists(arq_plus10) || is.na(file.size(arq_plus10)) || file.size(arq_plus10) == 0L) {
      stop("Objeto plus10 ausente ou vazio: ", arq_plus10, call. = FALSE)
    }
    ps_p10 <- tryCatch(
      readRDS(arq_plus10),
      error = function(e) stop("Falha ao ler plus10: ", conditionMessage(e), call. = FALSE)
    )
    if (!inherits(ps_p10, "phyloseq") || nsamples(ps_p10) != 10L) {
      stop("plus10 deve ser phyloseq com exatamente 10 amostras.", call. = FALSE)
    }
    meta_p10 <- validar_samplelabel_meta(as(sample_data(ps_p10), "data.frame"))
    meta_p10$SampleID <- rownames(meta_p10)
    meta_p10$BeeSpecies <- factor(
      meta_p10$BeeSpecies,
      levels = c("Melipona fasciculata", "Melipona scutellaris", "Melipona subnitida")
    )
    if (anyNA(meta_p10$BeeSpecies) ||
        !identical(as.integer(table(meta_p10$BeeSpecies)), c(3L, 3L, 4L))) {
      stop("plus10 deve apresentar distribuicao BeeSpecies 3/3/4.", call. = FALSE)
    }
    if (sum(as.character(meta_p10$Run) == "run_main") != 9L ||
        sum(as.character(meta_p10$Run) == "run_aux") != 1L) {
      stop("plus10 deve conter 9 amostras da corrida principal e 1 da corrida auxiliar.", call. = FALSE)
    }
    if (!"S10" %in% rownames(meta_p10)) {
      stop("Amostra S10 ausente no plus10.", call. = FALSE)
    }

    asv_p10 <- matriz_amostras(ps_p10)
    asv_p10 <- asv_p10[rownames(meta_p10), , drop = FALSE]
    if (anyNA(asv_p10) || any(!is.finite(asv_p10)) || any(asv_p10 < 0) ||
        any(rowSums(asv_p10) == 0)) {
      stop("Matriz plus10 contem valores invalidos ou amostras zeradas.", call. = FALSE)
    }
    asv_rel_p10 <- relativo(asv_p10)

    presenca_auxiliar <- asv_p10["S10", ] > 0
    presenca_demais <- colSums(
      asv_p10[setdiff(rownames(asv_p10), "S10"), , drop = FALSE] > 0
    ) > 0
    n_asvs_auxiliar <- sum(presenca_auxiliar)
    n_asvs_auxiliar_exclusivas <- sum(presenca_auxiliar & !presenca_demais)
    n_asvs_auxiliar_compartilhadas <- sum(presenca_auxiliar & presenca_demais)

    # Distancias oficiais produzidas pelo Script 06.
    arq_bray_p10 <- ctx$contracts[["plus10_dist_bray"]]
    arq_jac_p10 <- ctx$contracts[["plus10_dist_jaccard"]]
    distancias_ausentes <- c(arq_bray_p10, arq_jac_p10)[
      !file.exists(c(arq_bray_p10, arq_jac_p10))
    ]
    if (length(distancias_ausentes)) {
      stop(
        "Distancias plus10 ausentes: ",
        paste(distancias_ausentes, collapse = ", "),
        ". Execute o Script 06 antes do Script 10.",
        call. = FALSE
      )
    }
    dist_bray_p10 <- readRDS(arq_bray_p10)
    dist_jac_p10 <- readRDS(arq_jac_p10)
    dist_bray_p10 <- alinhar_distancia(dist_bray_p10, rownames(meta_p10), "plus10 Bray-Curtis")
    dist_jac_p10 <- alinhar_distancia(dist_jac_p10, rownames(meta_p10), "plus10 Jaccard")

    # G01P - diversidade alfa com p exato do Script 08.
    alfa_p10 <- phyloseq::estimate_richness(
      ps_p10, measures = c("Observed", "Shannon", "Simpson")
    )
    alfa_p10$Pielou <- ifelse(alfa_p10$Observed > 1,
                              alfa_p10$Shannon / log(alfa_p10$Observed), NA_real_)
    alfa_p10$SampleID <- rownames(alfa_p10)
    alfa_p10 <- cbind(
      alfa_p10,
      meta_p10[alfa_p10$SampleID,
               c("SampleLabel", "BeeSpecies", "Run"), drop = FALSE]
    )
    alfa_p10_long <- tidyr::pivot_longer(
      alfa_p10,
      cols = c("Observed", "Shannon", "Simpson", "Pielou"),
      names_to = "Metrica", values_to = "Valor"
    )
    kw_p10_path <- file.path(anal_path, "plus10_kruskal_wallis.csv")
    if (!file.exists(kw_p10_path)) {
      stop(
        "Resultado plus10 de diversidade alfa ausente: ", kw_p10_path,
        ". Execute o Script 08 antes do Script 10.",
        call. = FALSE
      )
    }
    kw_p10_plot <- read.csv(
      kw_p10_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    col_metrica_kw <- intersect(c("Medida", "Metrica"), colnames(kw_p10_plot))[1]
    col_p_kw <- intersect(c("p_adj_BH", "p_FDR_BH", "p_exato", "p_raw"), colnames(kw_p10_plot))[1]
    if (is.na(col_metrica_kw) || is.na(col_p_kw)) {
      stop("Tabela KW plus10 sem colunas de metrica/p-valor reconheciveis.", call. = FALSE)
    }
    kw_lab_p10 <- data.frame(
      Metrica = kw_p10_plot[[col_metrica_kw]],
      Rotulo_p = paste0("p exato/FDR = ", formatC(kw_p10_plot[[col_p_kw]], digits = 3, format = "fg")),
      stringsAsFactors = FALSE
    )
    alfa_p10_long <- dplyr::left_join(alfa_p10_long, kw_lab_p10, by = "Metrica")
    kw_annot_p10 <- dplyr::distinct(alfa_p10_long, Metrica, Rotulo_p)
    kw_annot_p10$BeeSpecies <- factor(
      levels(meta_p10$BeeSpecies)[1], levels = levels(meta_p10$BeeSpecies)
    )

    gp10_alpha <- ggplot(
      alfa_p10_long,
      aes(x = BeeSpecies, y = Valor, color = BeeSpecies)
    ) +
      geom_boxplot(outlier.shape = NA, alpha = 0.15, linewidth = 0.5) +
      geom_point(aes(shape = Run), size = 3,
                 position = position_jitter(width = 0.08, height = 0)) +
      ggrepel::geom_text_repel(
        aes(label = SampleLabel), size = 2.6, show.legend = FALSE,
        max.overlaps = Inf, box.padding = 0.25
      ) +
      facet_wrap(~ Metrica, scales = "free_y", ncol = 2) +
      scale_color_manual(values = CORES_ESP, labels = LABELS_ESP) +
      scale_x_discrete(labels = LABELS_ESP) +
      TEMA_PUB +
      theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
      labs(
        title = "Plus10 - diversidade alfa exploratoria",
        subtitle = paste(
          "Teste exato por enumeracao, grupos 3/3/4.",
          "Auxiliar pertence a segunda corrida; Run nao e estimavel separadamente."
        ),
        x = NULL, y = NULL, color = "Especie", shape = "Run"
      ) +
      geom_text(
        data = kw_annot_p10,
        aes(x = BeeSpecies, y = Inf, label = Rotulo_p),
        inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 3
      )

    # G02P - PCoA Bray-Curtis e Jaccard.
    criar_pcoa_p10 <- function(dist_obj, titulo) {
      ord <- phyloseq::ordinate(ps_p10, method = "PCoA", distance = dist_obj)
      df <- phyloseq::plot_ordination(ps_p10, ord, justDF = TRUE)
      df$SampleID <- rownames(df)
      df$SampleLabel <- meta_p10[df$SampleID, "SampleLabel"]
      df$BeeSpecies <- meta_p10[df$SampleID, "BeeSpecies"]
      df$Run <- meta_p10[df$SampleID, "Run"]
      eig <- ord$values$Relative_eig
      xlab <- if (length(eig) >= 1L) sprintf("PCoA1 (%.1f%%)", 100 * eig[1]) else "PCoA1"
      ylab <- if (length(eig) >= 2L) sprintf("PCoA2 (%.1f%%)", 100 * eig[2]) else "PCoA2"
      ggplot(df, aes(Axis.1, Axis.2, color = BeeSpecies, shape = Run)) +
        geom_point(size = 3.5) +
        ggrepel::geom_text_repel(aes(label = SampleLabel), size = 2.8,
                                  show.legend = FALSE, max.overlaps = Inf) +
        scale_color_manual(values = CORES_ESP, labels = LABELS_ESP) +
        TEMA_PUB +
        labs(
          title = titulo,
          subtitle = "Plus10 exploratorio; a segunda corrida possui uma unica amostra",
          x = xlab, y = ylab, color = "Especie", shape = "Run"
        )
    }
    gp10_pcoa_bray <- criar_pcoa_p10(dist_bray_p10, "Plus10 - PCoA Bray-Curtis")
    gp10_pcoa_jac <- criar_pcoa_p10(dist_jac_p10, "Plus10 - PCoA Jaccard")

    # G03P - composicao dos 15 generos mais abundantes.
    psg_p10 <- phyloseq::tax_glom(ps_p10, taxrank = "Genus", NArm = FALSE)
    psg_rel_p10 <- transform_sample_counts(psg_p10, function(x) x / sum(x))
    long_g_p10 <- phyloseq::psmelt(psg_rel_p10)
    long_g_p10$Genus <- as.character(long_g_p10$Genus)
    long_g_p10$Genus[is.na(long_g_p10$Genus) | trimws(long_g_p10$Genus) == ""] <- "Nao classificado"
    long_g_p10$SampleID <- as.character(long_g_p10$Sample)
    long_g_p10$SampleLabel <- meta_p10[long_g_p10$SampleID, "SampleLabel"]
    long_g_p10$BeeSpecies <- meta_p10[long_g_p10$SampleID, "BeeSpecies"]
    long_g_heat_p10 <- long_g_p10
    top15_p10 <- long_g_p10 |>
      dplyr::group_by(Genus) |>
      dplyr::summarise(Total = sum(Abundance), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(Total)) |>
      dplyr::slice_head(n = 15) |>
      dplyr::pull(Genus)
    long_g_p10$Genus_plot <- ifelse(long_g_p10$Genus %in% top15_p10,
                                    long_g_p10$Genus, "Outros")
    long_g_p10 <- long_g_p10 |>
      dplyr::group_by(SampleID, SampleLabel, BeeSpecies, Genus_plot) |>
      dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
    cores_g_p10 <- setNames(
      grDevices::colorRampPalette(CORES_MULTI)(length(unique(long_g_p10$Genus_plot))),
      sort(unique(long_g_p10$Genus_plot))
    )
    gp10_genus <- ggplot(long_g_p10,
                          aes(x = SampleLabel, y = Abundance * 100, fill = Genus_plot)) +
      geom_col(width = 0.8) +
      facet_grid(~ BeeSpecies, scales = "free_x", space = "free_x",
                 labeller = as_labeller(LABELS_ESP)) +
      scale_fill_manual(values = cores_g_p10) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
      TEMA_PUB +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.text = element_text(face = "italic", size = 8)) +
      labs(
        title = "Plus10 - composicao dos 15 generos mais abundantes",
        subtitle = sprintf(
          "Todas as %d ASVs presentes em S10 entram antes da agregacao taxonomica",
          n_asvs_auxiliar
        ),
        x = NULL, y = "Abundancia relativa (%)", fill = "Genero"
      )

    # G04P - heatmap dos 20 generos, com valores percentuais internos.
    top20_p10 <- long_g_heat_p10 |>
      dplyr::group_by(Genus) |>
      dplyr::summarise(Total = sum(Abundance), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(Total)) |>
      dplyr::slice_head(n = 20) |>
      dplyr::pull(Genus)
    heat_p10 <- long_g_heat_p10 |>
      dplyr::filter(Genus %in% top20_p10) |>
      dplyr::group_by(SampleLabel, BeeSpecies, Genus) |>
      dplyr::summarise(Percentual = sum(Abundance) * 100, .groups = "drop") |>
      dplyr::rename(Genus_plot = Genus)
    heat_p10$Rotulo <- ifelse(
      heat_p10$Percentual == 0, "",
      ifelse(heat_p10$Percentual < 0.01, "<0,01",
             formatC(heat_p10$Percentual, digits = 2, format = "f", decimal.mark = ","))
    )
    gp10_heat <- ggplot(heat_p10,
                         aes(x = SampleLabel, y = Genus_plot, fill = Percentual)) +
      geom_tile(color = "grey80", linewidth = 0.25) +
      geom_text(aes(label = Rotulo), size = 2.2, color = "black") +
      scale_fill_gradient(low = "white", high = "steelblue", name = "Abundancia\nrelativa (%)") +
      facet_grid(~ BeeSpecies, scales = "free_x", space = "free_x",
                 labeller = as_labeller(LABELS_ESP)) +
      TEMA_PUB +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(face = "italic", size = 8)) +
      labs(
        title = "Plus10 - heatmap dos generos dominantes",
        subtitle = "Valores internos em porcentagem; fundo branco para abundancias baixas",
        x = NULL, y = "Genero"
      )

    # G05P - volcanoes e comparacao de LFC core9 x plus10.
    criar_volcano_conjunto <- function(prefixo_arquivo, titulo_prefixo) {
      plots <- lapply(comps_deseq, function(nm) {
        arq <- file.path(deseq_path, paste0(prefixo_arquivo, nm, ".csv"))
        if (!file.exists(arq)) return(NULL)
        res <- read.csv(arq, stringsAsFactors = FALSE, check.names = FALSE)
        res <- res[is.finite(res$log2FoldChange) & !is.na(res$padj), , drop = FALSE]
        if (!nrow(res)) return(NULL)
        res$padj_plot <- pmax(res$padj, .Machine$double.xmin)
        res$Dir <- factor(dplyr::case_when(
          res$padj < ALPHA & res$log2FoldChange > LFC_THRESHOLD ~ "Up",
          res$padj < ALPHA & res$log2FoldChange < -LFC_THRESHOLD ~ "Down",
          TRUE ~ "NS"), levels = c("Up", "Down", "NS"))
        top5 <- res |>
          dplyr::filter(Dir != "NS") |>
          dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) |>
          dplyr::slice_head(n = 5) |>
          dplyr::mutate(rotulo = ifelse(!is.na(Genus) & Genus != "", Genus, ASV_ID))
        ggplot(res, aes(log2FoldChange, -log10(padj_plot), color = Dir)) +
          geom_point(alpha = 0.7, size = 1.8) +
          geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed") +
          geom_hline(yintercept = -log10(ALPHA), linetype = "dashed") +
          ggrepel::geom_text_repel(
            data = top5,
            aes(log2FoldChange, -log10(pmax(padj, .Machine$double.xmin)), label = rotulo),
            inherit.aes = FALSE, size = 2.5, show.legend = FALSE
          ) +
          scale_color_manual(values = c("Up" = "#E31A1C", "Down" = "#1F78B4", "NS" = "grey72")) +
          TEMA_PUB +
          labs(title = gsub("_", " vs ", nm), x = "log2 Fold Change", y = "-log10(padj)")
      })
      plots <- plots[!vapply(plots, is.null, logical(1))]
      if (!length(plots)) return(NULL)
      patchwork::wrap_plots(plots, ncol = 3) +
        patchwork::plot_annotation(
          title = titulo_prefixo,
          subtitle = paste(
            "Analise exploratoria plus10; design ~ BeeSpecies.",
            "Run nao foi ajustada porque a segunda corrida tem n=1."
          )
        )
    }
    gp10_volcano <- criar_volcano_conjunto(
      "plus10_completo_", "Plus10 - abundancia diferencial DESeq2"
    )

    lfc_plots_p10 <- lapply(comps_deseq, function(nm) {
      arq <- file.path(deseq_path, paste0("sensibilidade_core9_vs_plus10_", nm, ".csv"))
      if (!file.exists(arq)) return(NULL)
      z <- read.csv(arq, stringsAsFactors = FALSE, check.names = FALSE)
      z <- z[is.finite(z$log2FoldChange_core9) & is.finite(z$log2FoldChange_plus10), , drop = FALSE]
      if (nrow(z) < 3L) return(NULL)
      rho <- suppressWarnings(cor(z$log2FoldChange_core9, z$log2FoldChange_plus10,
                                  method = "spearman"))
      ggplot(z, aes(log2FoldChange_core9, log2FoldChange_plus10)) +
        geom_hline(yintercept = 0, color = "grey80") +
        geom_vline(xintercept = 0, color = "grey80") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_point(aes(color = Sig_plus10 | Sig_core9), alpha = 0.7) +
        scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B2182B"),
                           name = "Significativa\nem algum conjunto") +
        coord_equal() + TEMA_PUB +
        labs(title = gsub("_", " vs ", nm),
             subtitle = paste0("Spearman rho = ", signif(rho, 3)),
             x = "log2FC core9", y = "log2FC plus10")
    })
    lfc_plots_p10 <- lfc_plots_p10[!vapply(lfc_plots_p10, is.null, logical(1))]
    gp10_lfc_sens <- if (length(lfc_plots_p10)) {
      patchwork::wrap_plots(lfc_plots_p10, ncol = 3) +
        patchwork::plot_annotation(title = "Sensibilidade DESeq2: core9 versus plus10")
    } else NULL

    # G06P/G07P - SIMPER e IndVal do plus10.
    arq_simper_p10 <- file.path(anal_path, "plus10_simper_top10_por_par.csv")
    gp10_simper <- if (file.exists(arq_simper_p10)) {
      zz <- read.csv(arq_simper_p10, stringsAsFactors = FALSE, check.names = FALSE)
      zz$Comparacao <- gsub("_", " vs ", zz$Comparacao)
      zz$Label <- ifelse(!is.na(zz$Genus) & zz$Genus != "", zz$Genus, zz$ASV_ID)
      ggplot(zz, aes(reorder(Label, Contrib), Contrib * 100, fill = Comparacao)) +
        geom_col(position = "dodge") + coord_flip() + TEMA_PUB +
        labs(title = "Plus10 - SIMPER exploratorio", x = "Taxon", y = "Contribuicao media (%)")
    } else NULL

    arq_ind_p10 <- file.path(anal_path, "plus10_indval_significativas.csv")
    gp10_indval <- if (file.exists(arq_ind_p10)) {
      zz <- read.csv(arq_ind_p10, stringsAsFactors = FALSE, check.names = FALSE)
      if (!nrow(zz)) NULL else {
        zz$Label <- ifelse(!is.na(zz$Genus) & zz$Genus != "", zz$Genus, zz$ASV_ID)
        ggplot(zz, aes(reorder(Label, stat), stat, fill = Grupo)) +
          geom_col() + coord_flip() + TEMA_PUB +
          labs(title = "Plus10 - ASVs indicadoras", x = "Taxon", y = "IndVal")
      }
    } else NULL

    # G08P - Venn plus10.
    asv_por_esp_p10 <- setNames(
      lapply(levels(meta_p10$BeeSpecies), function(sp) {
        am <- rownames(meta_p10)[meta_p10$BeeSpecies == sp]
        colnames(asv_p10)[colSums(asv_p10[am, , drop = FALSE]) > 0]
      }),
      unname(LABELS_ESP[levels(meta_p10$BeeSpecies)])
    )
    gp10_venn <- ggVennDiagram::ggVennDiagram(
      asv_por_esp_p10, label_alpha = 0, edge_size = 1.1
    ) +
      scale_fill_gradientn(colours = c("#2f4b7c", "#a05195", "#d45087", "#f95d6a"),
                           name = "N de ASVs") +
      TEMA_PUB +
      theme(panel.border = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), panel.grid = element_blank()) +
      labs(title = "Plus10 - ASVs compartilhadas",
           subtitle = sprintf(
             "S10 inclui %d ASVs compartilhadas com core9 e %d exclusivas",
             n_asvs_auxiliar_compartilhadas,
             n_asvs_auxiliar_exclusivas
           ))

    graficos_plus10 <- list(
      list(id = "P10_G01_alpha", g = gp10_alpha, w = 15, h = 10),
      list(id = "P10_G02_pcoa_bray", g = gp10_pcoa_bray, w = 10, h = 7),
      list(id = "P10_G03_pcoa_jaccard", g = gp10_pcoa_jac, w = 10, h = 7),
      list(id = "P10_G04_genus_top15", g = gp10_genus, w = 15, h = 8),
      list(id = "P10_G05_heatmap_genus", g = gp10_heat, w = 14, h = 10),
      list(id = "P10_G08_venn", g = gp10_venn, w = 9, h = 8)
    )
    if (!is.null(gp10_volcano)) graficos_plus10 <- c(
      graficos_plus10, list(list(id = "P10_G06_volcano_deseq2", g = gp10_volcano, w = 18, h = 7))
    )
    if (!is.null(gp10_lfc_sens)) graficos_plus10 <- c(
      graficos_plus10, list(list(id = "P10_G07_lfc_core9_vs_plus10", g = gp10_lfc_sens, w = 18, h = 7))
    )
    if (!is.null(gp10_simper)) graficos_plus10 <- c(
      graficos_plus10, list(list(id = "P10_G09_simper", g = gp10_simper, w = 14, h = 8))
    )
    if (!is.null(gp10_indval)) graficos_plus10 <- c(
      graficos_plus10, list(list(id = "P10_G10_indval", g = gp10_indval, w = 10, h = 7))
    )

    # Saidas separadas; falha em um grafico e registrada sem reutilizar o ultimo plot.
    erros_graficos_p10 <- list()
    registrar_erro_p10 <- function(id, formato, e) {
      erros_graficos_p10[[length(erros_graficos_p10) + 1L]] <<- data.frame(
        Grafico = id, Formato = formato, Erro = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      log_msg(paste("plus10", id, formato, "falhou:", conditionMessage(e)), "ERRO")
    }

    pdf_p10 <- file.path(fig_path_p10, "figuras_plus10_exploratorias.pdf")
    grDevices::pdf(pdf_p10, width = 16, height = 10, onefile = TRUE)
    tryCatch(
      {
        for (item in graficos_plus10) {
          tryCatch(print(item$g), error = function(e) registrar_erro_p10(item$id, "PDF", e))
        }
      },
      finally = {
        if (grDevices::dev.cur() > 1L) grDevices::dev.off()
      }
    )
    for (item in graficos_plus10) {
      tryCatch(
        salvar_plot(item$g, paste0(item$id, ".png"), fig_path_p10,
                    width = item$w, height = item$h),
        error = function(e) registrar_erro_p10(item$id, "PNG", e)
      )
    }
    erros_p10_df <- if (length(erros_graficos_p10)) {
      dplyr::bind_rows(erros_graficos_p10)
    } else {
      data.frame(Grafico = character(), Formato = character(), Erro = character())
    }
    write.csv(erros_p10_df, file.path(fig_path_p10, "auditoria_erros_graficos_plus10.csv"),
              row.names = FALSE, quote = TRUE, na = "")
    write.csv(
      data.frame(
        Conjunto = "plus10",
        N_amostras = nsamples(ps_p10),
        N_ASVs_objeto = ntaxa(ps_p10),
        N_ASVs_auxiliar = n_asvs_auxiliar,
        N_ASVs_auxiliar_exclusivas = n_asvs_auxiliar_exclusivas,
        Inclui_segunda_run = TRUE,
        Run_ajustavel = FALSE,
        Interpretacao = paste(
          sprintf(
            "Sensibilidade exploratoria separada; todos os calculos usam as %d ASVs presentes em S10.",
            n_asvs_auxiliar
          ),
          "A unica amostra da segunda corrida impede separar Run de BeeSpecies."
        ),
        stringsAsFactors = FALSE
      ),
      file.path(fig_path_p10, "metadata_figuras_plus10.csv"),
      row.names = FALSE, quote = TRUE, na = ""
    )
    log_msg(paste("Figuras plus10 salvas em", fig_path_p10), "SAVE")

    TRUE
  },
  error = function(e) {
    dir.create(fig_path_p10, recursive = TRUE, showWarnings = FALSE)
    msg <- conditionMessage(e)
    log_msg(paste("Bloco grafico plus10 falhou:", msg), "ERRO")
    erros_p10_df <<- data.frame(
      Grafico = "BLOCO_PLUS10",
      Formato = "NA",
      Erro = msg,
      stringsAsFactors = FALSE
    )
    try(
      write.csv(
        erros_p10_df,
        file.path(fig_path_p10, "auditoria_erros_graficos_plus10.csv"),
        row.names = FALSE, quote = TRUE, na = ""
      ),
      silent = TRUE
    )
    FALSE
  }
)

###############################################################################
# METADADOS DE EXECUÇÃO
###############################################################################

perm_diag_core9_path <- file.path(
  anal_path, "core9_resolucao_permutacoes_simper_indval.csv"
)
perm_diag_core9 <- if (file.exists(perm_diag_core9_path)) {
  read.csv(perm_diag_core9_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}
indval_perm_core9 <- if (nrow(perm_diag_core9) > 0L &&
                           "Fator" %in% colnames(perm_diag_core9)) {
  perm_diag_core9[perm_diag_core9$Fator == "BeeSpecies", , drop = FALSE]
} else {
  data.frame()
}

run_meta <- data.frame(
  Script          = "10_graficos",
  Versao          = VERSAO,
  Data_execucao   = DATA_EXECUCAO,
  Core9_amostras  = nsamples(ps_c9),
  Core9_ASVs      = ntaxa(ps_c9),
  Plus10_amostras = if (!is.null(ps_p10)) nsamples(ps_p10) else NA_integer_,
  Plus10_ASVs     = if (!is.null(ps_p10)) ntaxa(ps_p10) else NA_integer_,
  Plus10_ASVs_Auxiliar = if (!is.null(asv_p10)) sum(asv_p10["S10", ] > 0) else NA_integer_,
  Plus10_graficos = length(graficos_plus10),
  Plus10_diretorio = fig_path_p10,
  Plus10_inferencia = "sensibilidade exploratoria separada; Run nao estimavel",
  Plus10_figuras_ok = plus10_figuras_ok,
  N_graficos      = length(graficos),
  PERMANOVA_fonte = "Script 08; testes exatos por enumeracao das rotulacoes unicas",
  IndVal_fonte = "Script 08; matriz customizada de alocacoes unicas dos rotulos",
  IndVal_alocacoes_totais_core9 = if (nrow(indval_perm_core9)) {
    indval_perm_core9$N_alocacoes_rotulos_totais[[1L]]
  } else NA_real_,
  IndVal_permutacoes_alternativas_core9 = if (nrow(indval_perm_core9)) {
    indval_perm_core9$N_permutacoes_alternativas[[1L]]
  } else NA_real_,
  IndVal_p_min_teorico_core9 = if (nrow(indval_perm_core9)) {
    indval_perm_core9$p_min_teorico[[1L]]
  } else NA_real_,
  LFC_threshold_volcano = LFC_THRESHOLD,
  IndVal_criterio = "p_adj_BH < 0.05 quando disponivel",
  NMDS_stress     = round(ord_nmds$stress, 4),
  PERMANOVA_R2    = perm_R2,
  PERMANOVA_p     = perm_p,
  Alpha_testes = paste(indices_kw, collapse = ";"),
  Chao1_ACE_uso = "suplementar_descritivo; nao testados apos filtro de ASVs raras",
  Alpha_profundidade_min = min(profundidade_core9),
  Alpha_profundidade_max = max(profundidade_core9),
  Alpha_razao_max_min = round(razao_profundidade_core9, 6),
  Alpha_razao_excede_10x = razao_profundidade_core9 > LIMIAR_RAZAO_PROFUNDIDADE,
  Alpha_decisao_profundidade = ifelse(
    razao_profundidade_core9 > LIMIAR_RAZAO_PROFUNDIDADE,
    "priorizar_iNEXT_cobertura_para_riqueza; sem_rarefacao_automatica",
    "estimate_richness_mantido; iNEXT_cobertura_suplementar"
  ),
  iNEXT_cobertura_comum = cobertura_comum,
  iNEXT_estimateD_gerado = !is.null(estimateD_cobertura),
  iNEXT_fonte = inext_fonte,
  iNEXT_ASVs = nrow(asv_inext),
  Falhas_graficos = nrow(erros_graficos_df),
  Falhas_graficos_plus10 = nrow(erros_p10_df),
  Origem_status = "derivada_de_BeeSpecies; nao_independente",
  stringsAsFactors = FALSE
)
write.csv(run_meta,
  file.path(fig_path, "metadata_execucao_script10.csv"),
  row.names = FALSE)

cat("\n=============================================================\n")
cat("SCRIPT 10 CONCLUÍDO\n")
cat(sprintf("Total de figuras previstas: %d | falhas registradas: %d\n",
            length(graficos), nrow(erros_graficos_df)))
if (nrow(erros_graficos_df) > 0L || nrow(erros_p10_df) > 0L) {
  warning(
    paste(
      "Uma ou mais figuras falharam. Consulte as auditorias do core9 e do",
      "plus10_exploratorio."
    ),
    call. = FALSE
  )
}
cat("Saídas em:", fig_path, "\n")
cat("=============================================================\n\n")

gc(verbose = FALSE)
log_msg("Finalizado", "FINAL")
})
