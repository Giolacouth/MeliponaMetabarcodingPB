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
run_pipeline_script("01a_dada2_multirun.R", "dada2", function(ctx) {
###############################################################################
# SCRIPT 01 — PIPELINE DADA2 (MULTI-RUN)
#
# Projeto: Microbioma de mel de Melipona spp. (abelhas sem ferrao)
# Regiao : 16S rRNA V3-V4 | Primers 341F/805R | Illumina MiSeq PE 2x301
#
# # PARAMETROS PRINCIPAIS:
#   truncLen       : definido por 1_1_parametro_trunclen.R (curva de retencao 99.9%)
#   min_prevalence : 2 (~22% das 9 amostras principais; logica E)
#   min_reads      : 10 (criterio operacional pre-especificado do projeto)
#   overlap ref    : 427 pb (inserto sem primers: 465 - 17 - 21)
#
# # DESIGN MULTI-RUN:
#   Corrida principal (n=9):
#     - learnErrors proprio;
#     - inferencia com pool = TRUE;
#     - filtro por abundancia e prevalencia.
#
#   Corrida auxiliar (n=1):
#     - learnErrors proprio;
#     - inferencia sem pooling entre amostras;
#     - filtro somente por abundancia minima;
#     - modelo da corrida principal aplicado apenas como sensibilidade.
#
#   Remocao de quimeras:
#     - executada sobre a tabela combinada das 10 amostras;
#     - cada cenario da sensibilidade deve repetir essa etapa.
#
#   Saidas:
#     - seqtab_nochim.rds: corrida principal, colunas = sequencias;
#     - seqtab_nochim_ids.rds: corrida principal, colunas = IDs;
#     - seqtab_auxiliar.rds: corrida auxiliar, colunas = sequencias;
#     - seqtab_nochim_prefiltro.rds: combinado, pos-quimera e pre-filtros.
#
# LIMITACOES:
#   Sem controle negativo.
#   O modelo da corrida auxiliar foi estimado com uma unica amostra;
#   sua estabilidade deve ser avaliada por analise de sensibilidade.
###############################################################################

options(encoding = "UTF-8")
tryCatch(
  Sys.setlocale("LC_ALL", "pt_BR.UTF-8"),
  warning = function(w) message("Locale pt_BR.UTF-8 indisponivel; seguindo com locale atual."),
  error   = function(e) message("Locale pt_BR.UTF-8 indisponivel; seguindo com locale atual.")
)

# Semente compartilhada para reprodutibilidade entre as etapas estocasticas.

SEED_GLOBAL <- 1234L
set.seed(SEED_GLOBAL)

###############################################################################
# 1. PACOTES
###############################################################################
pkgs <- c("dada2", "Biostrings", "ShortRead",
          "ggplot2", "dplyr", "tibble", "tidyr", "vegan")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Pacote nao instalado: ", p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("dada2 version:", as.character(packageVersion("dada2")), "\n")

###############################################################################
# 2. PARAMETROS EDITAVEIS
###############################################################################
pipeline_version <- ctx$version
base_path <- ctx$base_path
raw_path <- ctx$raw_path
cutadapt <- Sys.getenv("CUTADAPT_BIN", unset = "")
if (!nzchar(cutadapt)) {
  cutadapt <- Sys.which("cutadapt")
}

FWD <- "CCTACGGGNGGCWGCAG"
REV <- "GACTACHVGGGTATCTAATCC"

truncLen_by_run <- list("run_main" = c(283, 278), "run_aux" = c(283, 278))
default_truncLen <- c(283, 278)
maxEE_par  <- c(2, 2)
truncQ_par <- 2
minLen_par <- 50

nbases_par <- Inf
# Política explícita por corrida: pooling apenas na corrida principal.
pool_by_run <- c("run_main" = TRUE, "run_aux" = FALSE)
pool_fallback <- "pseudo"

minOverlap_par  <- 20
maxMismatch_par <- 0

expected_min_len <- 370
expected_max_len <- 480

min_total_reads <- 10   # limiar absoluto de leituras totais por ASV
min_prevalence  <- 2    # >= 2 de 9 amostras principais (~22%); logica E

# ---------------------------------------------------------------------------
# CHECKPOINT DE FILTRO (interrupcao obrigatoria antes de fixar o cenario)
# Fluxo operacional:
#   (1) Execute este script. Ele processa ate o filtro, grava a auditoria do
#       que seria removido e PARA com stop() no checkpoint.
#   (2) Inspecione ASVs_descartadas_principal_com_sequencia.tsv e, para
#       comparar cenarios, rode 01b_auditoria_quimeras_recuperacao_ASVs.R (le o
#       seqtab_nochim_prefiltro.rds, ja salvo antes do checkpoint).
#   (3) Defina CENARIO_FILTRO abaixo (edite esta linha) e EXECUTE O SCRIPT
#       NOVAMENTE. A decisao NAO e lida de nenhum arquivo: e voce quem a
#       insere aqui.
# Opcoes:
#   NA_character_ -> nao confirmado; o script para no checkpoint.
#   "aplicar"     -> aplica o filtro (reads < min_total_reads E prev < min_prevalence).
#   "sem_filtro"  -> nao aplica filtro de frequencia; mantem todas as ASVs nao-quimericas.
# Observacao: definir a escolha DEPOIS de inspecionar exige um criterio
# pre-especificado e independente dos desfechos a jusante (evita vies de
# selecao); trate cenarios alternativos como analise de sensibilidade.
if (!exists("CENARIO_FILTRO")) CENARIO_FILTRO <-  "aplicar"

###############################################################################
# 3. CAMINHOS E DIRETORIOS
###############################################################################
trimmed_path  <- file.path(ctx$stage$root, "fastq_trimmed")
filtered_path <- file.path(ctx$stage$root, "fastq_filtered")
quality_path  <- file.path(ctx$stage$root, "qualidade")
output_path   <- ctx$stage$root
log_path      <- ctx$stage$logs
plot_path     <- ctx$stage$figures
for (d in c(trimmed_path, filtered_path, quality_path, output_path, log_path, plot_path)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
}

parametros_execucao <- data.frame(
  Pipeline_version = pipeline_version,
  Run_principal = "run_main",
  Run_auxiliar = "run_aux",
  Min_total_reads_principal = min_total_reads,
  Min_prevalence_principal = min_prevalence,
  Min_total_reads_auxiliar = min_total_reads,
  Cenario_filtro = CENARIO_FILTRO,
  Min_overlap = minOverlap_par,
  Max_mismatch_merge = maxMismatch_par,
  Expected_min_len = expected_min_len,
  Expected_max_len = expected_max_len,
  Pool_principal = unname(pool_by_run["run_main"]),
  Pool_auxiliar = unname(pool_by_run["run_aux"]),
  stringsAsFactors = FALSE
)
write.csv(
  parametros_execucao,
  file.path(output_path, "pipeline_parametros_execucao.csv"),
  row.names = FALSE, quote = TRUE
)

if (!nzchar(cutadapt) || !file.exists(cutadapt)) {
  stop(
    "cutadapt nao localizado. Instale-o no PATH ou defina CUTADAPT_BIN.",
    call. = FALSE
  )
}

###############################################################################
# 4. METADADOS
# Apenas carregamento e validacao estrutural.
# Identificacao das corridas (run_principal/auxiliar) e feita na Secao 6,
# apos o cruzamento com os FASTQs disponiveis.
###############################################################################
cat("\n=== IMPORTANDO METADADOS ===\n")
metadata_path <- file.path(raw_path, "metadados.tsv")
if (!file.exists(metadata_path))
  stop("metadados.tsv nao encontrado em: ", raw_path)

metadata <- read.table(metadata_path, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("", "NA", "NaN"))

# Validacao das colunas obrigatorias ANTES de qualquer operacao sobre elas
# Aceita SampleLabels como alias legado, mas padroniza para SampleLabel.
if (!"SampleLabel" %in% colnames(metadata) && "SampleLabels" %in% colnames(metadata)) {
  colnames(metadata)[colnames(metadata) == "SampleLabels"] <- "SampleLabel"
}
if (!"SampleID"    %in% colnames(metadata)) stop("Coluna 'SampleID' ausente.")
if (!"SampleLabel" %in% colnames(metadata)) stop("Coluna 'SampleLabel' ausente.")
if (!"Run"         %in% colnames(metadata)) stop("Coluna 'Run' ausente.")

colnames(metadata) <- make.names(colnames(metadata), unique = TRUE)
char_cols <- sapply(metadata, is.character)
metadata[char_cols] <- lapply(metadata[char_cols], trimws)

metadata$SampleID <- trimws(as.character(metadata$SampleID))
metadata$Run <- trimws(as.character(metadata$Run))

if (anyNA(metadata$SampleID) || any(metadata$SampleID == "")) {
  stop("Existem SampleID ausentes ou vazios no metadado.")
}
if (anyDuplicated(metadata$SampleID) > 0L) {
  stop(
    "Existem SampleID duplicados: ",
    paste(
      unique(metadata$SampleID[duplicated(metadata$SampleID)]),
      collapse = ", "
    )
  )
}
if (anyNA(metadata$Run) || any(metadata$Run == "")) {
  stop("Existem valores Run ausentes ou vazios no metadado.")
}

if (anyNA(metadata$SampleLabel) || any(metadata$SampleLabel == "")) {
  stop("Existem SampleLabel ausentes ou vazios no metadado.")
}
if (anyDuplicated(metadata$SampleLabel) > 0) {
  stop("Existem SampleLabel duplicados: ",
       paste(unique(metadata$SampleLabel[duplicated(metadata$SampleLabel)]), collapse = ", "),
       ". Cada amostra deve possuir um rotulo unico para os graficos.")
}

cat("Amostras nos metadados:", nrow(metadata), "\n")
cat("Corridas:\n"); print(table(metadata$Run))

###############################################################################
# 5. LEITURA E PAREAMENTO DOS FASTQ
###############################################################################

cat("\n=== LEITURA DOS FASTQ ===\n")

fnFs <- sort(list.files(
  raw_path,
  pattern = "_R1_001.fastq.gz",
  full.names = TRUE
))

fnRs <- sort(list.files(
  raw_path,
  pattern = "_R2_001.fastq.gz",
  full.names = TRUE
))

# --------------------------------------------------------------------------
# Validacao e pareamento explicito dos arquivos R1 e R2
# --------------------------------------------------------------------------

if (length(fnFs) == 0L || length(fnRs) == 0L) {
  stop(
    "Nenhum FASTQ com padrao ",
    "_R1_001/_R2_001.fastq.gz em ",
    raw_path
  )
}

baseF <- sub(
  "_R1_001\\.fastq\\.gz$",
  "",
  basename(fnFs)
)

baseR <- sub(
  "_R2_001\\.fastq\\.gz$",
  "",
  basename(fnRs)
)

if (anyDuplicated(baseF)) {
  duplicados_R1 <- unique(
    baseF[duplicated(baseF)]
  )

  stop(
    "Mais de um arquivo R1 encontrado para: ", paste(duplicados_R1, collapse = ", ")
  )
}

if (anyDuplicated(baseR)) {duplicados_R2 <- unique( baseR[duplicated(baseR)]
  )

  stop( "Mais de um arquivo R2 encontrado para: ",paste(duplicados_R2, collapse = ", ")
  )
}

sem_R2 <- setdiff(baseF, baseR)
sem_R1 <- setdiff(baseR, baseF)

if (length(sem_R2) > 0L || length(sem_R1) > 0L) {
  stop(
    "Pares FASTQ incompletos.\n",
    "Sem R2: ",
    if (length(sem_R2) > 0L) {
      paste(sem_R2, collapse = ", ")
    } else {
      "nenhum"
    },
    "\nSem R1: ",
    if (length(sem_R1) > 0L) {
      paste(sem_R1, collapse = ", ")
    } else {
      "nenhum"
    }
  )
}

sample_names <- sort(
  intersect(baseF, baseR)
)

fnFs <- fnFs[
  match(sample_names, baseF)
]

fnRs <- fnRs[
  match(sample_names, baseR)
]

sample_names_R1 <- sub(
  "_R1_001\\.fastq\\.gz$",
  "",
  basename(fnFs)
)

sample_names_R2 <- sub(
  "_R2_001\\.fastq\\.gz$",
  "",
  basename(fnRs)
)

if (!identical(sample_names_R1, sample_names)) {
  stop(
    "Falha ao alinhar os arquivos R1 com os SampleID."
  )
}

if (!identical(sample_names_R2, sample_names)) {
  stop(
    "Falha ao alinhar os arquivos R2 com os SampleID."
  )
}

if (!identical(sample_names_R1, sample_names_R2)) {
  stop(
    "R1 e R2 permaneceram em ordens diferentes apos o pareamento."
  )
}

cat(
  "Amostras com par completo:",
  length(sample_names),
  "\n"
)

print(sample_names)
###############################################################################
# 5B. FUNCAO CONTAR_READS + CONTAGEM BRUTA (PRE-CUTADAPT)
###############################################################################
contar_reads <- function(f)
  tryCatch(as.integer(ShortRead::countFastq(f)$records),
           error = function(e) NA_integer_)

cat("\n=== CONTAGEM DE READS BRUTAS (PRE-CUTADAPT) ===\n")
reads_raw_R1 <- sapply(fnFs, contar_reads); names(reads_raw_R1) <- sample_names
reads_raw_R2 <- sapply(fnRs, contar_reads); names(reads_raw_R2) <- sample_names
cat("Reads brutas por amostra (R1):\n"); print(reads_raw_R1)
cat("Total de reads brutas (R1):",
    format(sum(reads_raw_R1, na.rm = TRUE), big.mark = ","), "\n")

if (any(is.na(reads_raw_R1)) || any(is.na(reads_raw_R2)))
  stop("Falha ao contar reads brutas em uma ou mais amostras. Verifique os FASTQ.")
if (any(reads_raw_R1 != reads_raw_R2))
  stop("Contagens brutas R1/R2 diferentes. Verifique pareamento dos FASTQ.")

###############################################################################
# 6. CRUZAMENTO FASTQ x METADADOS
# Identificacao das corridas feita aqui, sobre o metadata JA CRUZADO com os
# FASTQs disponiveis — evita resultado diferente se metadata tem amostras
# sem FASTQ correspondente.
###############################################################################
cat("\n=== CRUZAMENTO FASTQ x METADADOS ===\n")
missing_metadata <- setdiff(sample_names, metadata$SampleID)
missing_fastq    <- setdiff(metadata$SampleID, sample_names)
if (length(missing_metadata) > 0)
  stop("FASTQ sem metadados: ", paste(missing_metadata, collapse = ", "))
if (length(missing_fastq) > 0)
  warning("Metadados sem FASTQ: ", paste(missing_fastq, collapse = ", "))

metadata <- metadata[match(sample_names, metadata$SampleID), ]
rownames(metadata) <- metadata$SampleID
run_of_sample <- setNames(metadata$Run, metadata$SampleID)
cat("Mapa amostra -> corrida:\n"); print(run_of_sample)

write.table(metadata, file.path(output_path, "metadata_final.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Identificacao das corridas — feita UMA UNICA VEZ, aqui
# metadata ja esta cruzado com os FASTQs disponíveis
run_n_amostras <- table(metadata$Run)
if (length(run_n_amostras) != 2L) {
  stop("Este pipeline espera exatamente 2 corridas. Encontrado: ",
       paste(names(run_n_amostras), collapse = ", "))
}
run_principal <- "run_main"
run_auxiliar  <- "run_aux"
if (!run_principal %in% metadata$Run) stop("Corrida principal run_main ausente nos metadados.")
if (!run_auxiliar  %in% metadata$Run) stop("Corrida auxiliar run_aux ausente nos metadados.")
if (nrow(metadata) != 10L ||
    unname(run_n_amostras[run_principal]) != 9L ||
    unname(run_n_amostras[run_auxiliar]) != 1L) {
  stop(
    "Desenho experimental invalido: eram esperadas 10 amostras, com 9 em ",
    run_principal, " e 1 em ", run_auxiliar, ". Encontrado: ",
    paste(names(run_n_amostras), as.integer(run_n_amostras), sep = "=", collapse = "; ")
  )
}
runs_ordenados <- c(run_principal, run_auxiliar)

aux_ids <- metadata$SampleID[metadata$Run == run_auxiliar]
if (length(aux_ids) != 1L || !identical(unname(aux_ids), "S10")) {
  stop(
    "A corrida auxiliar deve conter exclusivamente S10; encontrado: ",
    paste(aux_ids, collapse = ", ")
  )
}

cat(sprintf("\nCorrida principal : %s (%d amostras) — learnErrors proprio\n",
            run_principal, run_n_amostras[run_principal]))
cat(sprintf("Corrida auxiliar  : %s (%d amostra)  — learnErrors proprio\n",
            run_auxiliar,  run_n_amostras[run_auxiliar]))
cat("Ordem de processamento:", paste(runs_ordenados, collapse = " -> "), "\n\n")

###############################################################################
# 7. ORIENTACOES E REV-COMP DOS PRIMERS
###############################################################################

validar_primer <- function(x, nome = "primer") {
  x <- toupper(trimws(as.character(x)))

  if (length(x) != 1 || is.na(x) || x == "") {
    stop(nome, " vazio ou invalido.")
  }

  if (!grepl("^[ACGTRYSWKMBDHVN]+$", x)) {
    chars <- unique(strsplit(x, "")[[1]])
    invalidos <- setdiff(chars, strsplit("ACGTRYSWKMBDHVN", "")[[1]])

    stop(
      nome, " contem caractere(s) invalido(s): ",
      paste(invalidos, collapse = ", "),
      "\nValor recebido: ", x,
      "\nVerifique se houve troca de letra O por numero 0, espaco oculto ou caractere copiado incorretamente."
    )
  }

  x
}

FWD <- validar_primer(FWD, "FWD")
REV <- validar_primer(REV, "REV")

allOrients <- function(primer) {
  primer <- validar_primer(primer, "primer")

  dna <- Biostrings::DNAString(primer)

  out <- c(
    Forward    = as.character(dna),
    Complement = as.character(Biostrings::complement(dna)),
    Reverse    = as.character(Biostrings::reverse(dna)),
    RevComp    = as.character(Biostrings::reverseComplement(dna))
  )

  out
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

FWD.RC <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(FWD)))
REV.RC <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(REV)))

primerHits <- function(primer, fn) {
  primer <- validar_primer(primer, "primer em primerHits")

  reads <- ShortRead::sread(ShortRead::readFastq(fn))

  sum(
    Biostrings::vcountPattern(
      pattern = Biostrings::DNAString(primer),
      subject = reads,
      fixed = FALSE
    ) > 0
  )
}

cat("\n=== PRIMERS — PRE-CUTADAPT (amostra 1) ===\n")

primer_diag <- tryCatch(
  {
    rbind(
      FWD.R1 = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
      REV.R1 = sapply(REV.orients, primerHits, fn = fnFs[[1]]),
      FWD.R2 = sapply(FWD.orients, primerHits, fn = fnRs[[1]]),
      REV.R2 = sapply(REV.orients, primerHits, fn = fnRs[[1]])
    )
  },
  error = function(e) {
    warning(
      "Diagnostico PRE-CUTADAPT de primers falhou, mas o pipeline seguira para o cutadapt.\n",
      "Erro original: ", conditionMessage(e)
    )
    NULL
  }
)

if (!is.null(primer_diag)) {
  print(primer_diag)
}
###############################################################################
# 8. CUTADAPT
###############################################################################
cat("\n=== CUTADAPT ===\n")
fnFs.cut <- file.path(trimmed_path, basename(fnFs))
fnRs.cut <- file.path(trimmed_path, basename(fnRs))

g_R1 <- paste0("^", FWD); G_R2 <- paste0("^", REV)
a_R1 <- REV.RC;            A_R2 <- FWD.RC

cut_flags <- c(
  "-g", g_R1, "-G", G_R2, "-a", a_R1, "-A", A_R2,
  "-e", "0.15",
  "--match-read-wildcards",
  "--max-n", "0",
  "--discard-untrimmed",
  "--pair-filter=any",
  "--minimum-length", as.character(minLen_par),
  "--overlap", "15",
  "--cores", "0"
)

for (i in seq_along(fnFs)) {
  cat("  Processando:", sample_names[i], "\n")
  log_i <- file.path(log_path, paste0("cutadapt_", sample_names[i], ".log"))
  status <- system2(cutadapt,
                    args = c(cut_flags, "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs[i], fnRs[i]),
                    stdout = log_i, stderr = log_i)
  if (!is.null(attr(status, "status")) || (is.numeric(status) && status != 0))
    stop("cutadapt retornou status != 0 para ", sample_names[i],
         ". Inspecione ", log_i)
  if (!file.exists(fnFs.cut[i]) || !file.exists(fnRs.cut[i]))
    stop("cutadapt nao gerou arquivos de saida para ", sample_names[i],
         ". Inspecione ", log_i)
}

###############################################################################
# 9. VERIFICACAO POS-CUTADAPT
###############################################################################
cat("\n=== VERIFICACAO POS-CUTADAPT ===\n")

reads_R1 <- sapply(fnFs.cut, contar_reads); names(reads_R1) <- sample_names
reads_R2 <- sapply(fnRs.cut, contar_reads); names(reads_R2) <- sample_names

if (any(is.na(reads_R1)) || any(is.na(reads_R2)))
  stop("Falha ao contar reads pos-cutadapt. Verifique logs em ", log_path)
if (any(reads_R1 == 0 | reads_R2 == 0))
  stop("Uma ou mais amostras ficaram sem reads apos cutadapt.")
if (any(reads_R1 != reads_R2))
  stop("Contagens R1/R2 diferentes apos cutadapt. Pares dessincronizados.")

pct_cutadapt <- round(100 * reads_R1 / reads_raw_R1[sample_names], 1)
cutadapt_summary <- data.frame(
  SampleID   = sample_names,
  Reads_raw  = reads_raw_R1[sample_names],
  Reads_cut  = reads_R1,
  Pct_retida = pct_cutadapt)
cat("Retencao apos cutadapt:\n"); print(cutadapt_summary)
write.csv(cutadapt_summary,
          file.path(output_path, "cutadapt_retention.csv"), row.names = FALSE)

cat(sprintf("\nTotal bruto        : %s reads\n",
            format(sum(reads_raw_R1[sample_names], na.rm = TRUE), big.mark = ",")))
cat(sprintf("Total pos-cutadapt : %s reads\n",
            format(sum(reads_R1, na.rm = TRUE), big.mark = ",")))
cat(sprintf("Retencao global    : %.1f%%\n",
            100 * sum(reads_R1, na.rm = TRUE) /
              sum(reads_raw_R1[sample_names], na.rm = TRUE)))
min_idx <- which.min(reads_R1)
cat(sprintf("Amostra com menos reads: %s (%s reads)\n",
            sample_names[min_idx], format(reads_R1[min_idx], big.mark = ",")))

vazias <- sample_names[is.na(reads_R1) | reads_R1 == 0 |
                         is.na(reads_R2) | reads_R2 == 0]
if (length(vazias) > 0)
  warning("Amostras vazias apos cutadapt: ", paste(vazias, collapse = ", "))

descasados <- sample_names[!is.na(reads_R1) & !is.na(reads_R2) &
                             reads_R1 != reads_R2]
if (length(descasados) > 0)
  warning("R1/R2 com contagens diferentes: ", paste(descasados, collapse = ", "))

# Primer residual
cat("\n=== PRIMER RESIDUAL ===\n")
check_primer_residual <- function(idx) {
  pos <- rbind(
    FWD.R1 = sapply(FWD.orients, primerHits, fn = fnFs.cut[[idx]]),
    REV.R2 = sapply(REV.orients, primerHits, fn = fnRs.cut[[idx]]))
  total   <- sum(pos); n_reads <- reads_R1[idx]
  data.frame(amostra = sample_names[idx], primer_hits = total,
             n_reads = n_reads,
             pct = ifelse(is.na(n_reads) || n_reads == 0,
                          NA_real_, round(100 * total / n_reads, 4)),
             row.names = NULL)
}
residual_tab <- do.call(rbind, lapply(seq_along(fnFs.cut), check_primer_residual))
print(residual_tab)
write.csv(residual_tab,
          file.path(output_path, "primer_residual_check.csv"), row.names = FALSE)
acima <- residual_tab$amostra[!is.na(residual_tab$pct) & residual_tab$pct > 0.1]
if (length(acima) > 0)
  warning("Primer residual >0.1% em: ", paste(acima, collapse = ", "))

# Comprimento pos-corte
cat("\n=== COMPRIMENTO POS-CORTE x truncLen ===\n")
runs <- unique(metadata$Run)
len_check <- do.call(rbind, lapply(runs, function(r) {
  idx <- which(run_of_sample[sample_names] == r)
  if (length(idx) == 0) return(NULL)
  f1 <- fnFs.cut[idx][reads_R1[idx] > 0][1]
  f2 <- fnRs.cut[idx][reads_R2[idx] > 0][1]
  if (is.na(f1) || is.na(f2)) return(NULL)
  medF <- median(width(sread(readFastq(f1))))
  medR <- median(width(sread(readFastq(f2))))
  tl <- truncLen_by_run[[r]]; if (is.null(tl)) tl <- default_truncLen
  data.frame(Run = r,
             mediana_R1 = medF, truncLen_R1 = tl[1], folga_R1 = medF - tl[1],
             mediana_R2 = medR, truncLen_R2 = tl[2], folga_R2 = medR - tl[2])
}))
print(len_check)
if (any(len_check$folga_R1 < 0 | len_check$folga_R2 < 0))
  warning("truncLen excede a mediana pos-corte em alguma corrida.")

# Perfis de qualidade — usa runs_ordenados (definido na Secao 6)
cat("\n=== PERFIS DE QUALIDADE ===\n")
plot_qualidade_seguro <- function(arquivos, titulo, pdf_file) {
  ok <- arquivos[sapply(arquivos, function(f) {
    n <- contar_reads(f); !is.na(n) && n > 0 })]
  if (length(ok) == 0) { warning("Nenhum arquivo: ", titulo); return(invisible()) }
  pdf(pdf_file, width = 10, height = 6)
  tryCatch(
    print(plotQualityProfile(ok, aggregate = TRUE) +
            ggtitle(paste(titulo, "- AGREGADO"))),
    error = function(e) warning("Falha no perfil agregado: ", conditionMessage(e)))
  for (f in ok)
    tryCatch(print(plotQualityProfile(f) + ggtitle(paste(titulo, "-", basename(f)))),
             error = function(e) warning("Falha: ", basename(f)))
  dev.off()
}
for (r in runs_ordenados) {
  idx <- which(run_of_sample[sample_names] == r)
  if (length(idx) == 0) next
  plot_qualidade_seguro(fnFs.cut[idx], paste("R1", r),
                        file.path(quality_path, paste0("qualidade_R1_", r, ".pdf")))
  plot_qualidade_seguro(fnRs.cut[idx], paste("R2", r),
                        file.path(quality_path, paste0("qualidade_R2_", r, ".pdf")))
}
cat("Perfis salvos em:", quality_path, "\n")

###############################################################################
# 10. FUNCAO DADA COM FALLBACK
###############################################################################
dada_with_fallback <- function(derep, err, pool_primary, pool_fb, label) {
  res <- tryCatch(
    {
      out <- dada(derep, err = err, pool = pool_primary, multithread = TRUE)
      attr(out, "pool_used") <- pool_primary
      out
    },
    error = function(e) {
      warning(label, ": pool=", pool_primary, " falhou. Usando '", pool_fb, "'.")
      out <- dada(derep, err = err, pool = pool_fb, multithread = TRUE)
      attr(out, "pool_used") <- pool_fb
      out
    })
  return(res)
}

###############################################################################
# 11. PROCESSAMENTO POR CORRIDA
# Ordem: corrida principal primeiro; corrida auxiliar aprende modelo proprio.
# O modelo principal e transferido para a auxiliar somente na analise de sensibilidade.
###############################################################################
seqtab_list <- list(); track_list  <- list()
filter_list <- list(); merge_list  <- list(); pool_used <- list()

# Modelos de erro da corrida principal — inicializados como NULL;
# a corrida auxiliar falha explicitamente se o principal nao foi processado
errF_principal <- NULL
errR_principal <- NULL

cat("\nMemoria disponivel:\n"); print(gc())

for (r in runs_ordenados) {

  cat("\n############################################################\n")
  cat("CORRIDA:", r, "\n")
  if (r == run_principal)
    cat("TIPO: PRINCIPAL — learnErrors proprio\n") else
      cat("TIPO: AUXILIAR  — learnErrors proprio (amostra unica)\n")
  cat("############################################################\n")

  idx  <- which(run_of_sample[sample_names] == r)
  sn_r <- sample_names[idx]
  if (length(idx) == 0L) stop("Nenhuma amostra encontrada para a corrida esperada: ", r)
  cat("Amostras (", length(sn_r), "):\n", sep = ""); print(sn_r)

  fnFs.cut.r  <- fnFs.cut[idx]; fnRs.cut.r <- fnRs.cut[idx]
  fnFs.filt.r <- file.path(filtered_path, paste0(sn_r, "_F_filt.fastq.gz"))
  fnRs.filt.r <- file.path(filtered_path, paste0(sn_r, "_R_filt.fastq.gz"))
  names(fnFs.filt.r) <- sn_r; names(fnRs.filt.r) <- sn_r

  tl <- truncLen_by_run[[r]]; if (is.null(tl)) tl <- default_truncLen

  overlap_esperado <- tl[1] + tl[2] - 427
  if (overlap_esperado < minOverlap_par)
    warning("Corrida ", r, ": sobreposicao esperada (~", overlap_esperado,
            " pb) < minOverlap (", minOverlap_par, ").")
  cat("truncLen:", tl[1], "/", tl[2],
      "| sobreposicao esperada ~", overlap_esperado, "pb\n")

  out_r <- filterAndTrim(
    fnFs.cut.r, fnFs.filt.r,
    fnRs.cut.r, fnRs.filt.r,
    truncLen = tl,
    maxN = 0,
    maxEE = maxEE_par,
    truncQ = truncQ_par,
    minLen = minLen_par,
    rm.phix = TRUE,
    rm.lowcomplex = c(8, 8),
    matchIDs = TRUE,
    compress = TRUE,
    multithread = TRUE,
    n = 5e4,
    verbose = TRUE
  )
  rownames(out_r) <- sn_r

  surv <- rownames(out_r)[out_r[, "reads.out"] > 0]
  if (length(surv) < length(sn_r)) {
    stop(
      "Corrida ", r, ": amostra(s) zerada(s) apos filterAndTrim; ",
      "o universo amostral nao sera reduzido silenciosamente: ",
      paste(setdiff(sn_r, surv), collapse = ", ")
    )
  }
  fnFs.filt.r <- fnFs.filt.r[surv]; fnRs.filt.r <- fnRs.filt.r[surv]
  sn_r_surv <- surv

  filter_list[[r]] <- data.frame(
    SampleID           = rownames(out_r), Run = r,
    Reads_Input        = out_r[, "reads.in"],
    Reads_Filtradas    = out_r[, "reads.out"],
    Retencao_Filtragem = ifelse(out_r[, "reads.in"] > 0,
                                round(out_r[, "reads.out"] / out_r[, "reads.in"] * 100, 1),
                                NA_real_))

  # --------------------------------------------------------------------------
  # learnErrors: proprio para corrida principal; substituido para auxiliar
  # --------------------------------------------------------------------------
  if (r == run_principal) {
    cat("Aprendendo erros R1 —", r, "(corrida principal, nbases=Inf)...\n")
    errF_r <- learnErrors(fnFs.filt.r, nbases = nbases_par, multithread = TRUE,
                          verbose = TRUE)
    cat("Aprendendo erros R2 —", r, "(corrida principal, nbases=Inf)...\n")
    errR_r <- learnErrors(fnRs.filt.r, nbases = nbases_par, multithread = TRUE,
                          verbose = TRUE)
    errF_principal <- errF_r
    errR_principal <- errR_r

    pdf(file.path(plot_path, paste0("modelo_erros_", r, ".pdf")), 10, 8)

    suppressWarnings({
      print(plotErrors(errF_r, nominalQ = TRUE) +
              ggtitle(paste("Erros R1 —", r, "(modelo proprio)")))
      print(plotErrors(errR_r, nominalQ = TRUE) +
              ggtitle(paste("Erros R2 —", r, "(modelo proprio)")))
    })
    dev.off()

  } else {

    # Corrida auxiliar: aprende o proprio modelo de erros (learnErrors).
    # Justificativa (Prodan 2020; Callahan 2016): o DADA2 estima o erro por
    # corrida, pois corridas distintas podem ter perfis de erro distintos. O
    # learnErrors utiliza as reads disponíveis, mas uma unica amostra reduz a
    # diversidade de perfis de qualidade representados. O ajuste e tecnicamente
    # executavel, nao automaticamente estavel; por isso deve ser inspecionado e
    # comparado ao modelo transferido da corrida principal.

    cat(sprintf(
      "Corrida %s (%d amostra): learnErrors PROPRIO (nbases=Inf).\n", r, length(sn_r)))
    cat("Aprendendo erros R1 —", r, "(corrida auxiliar, modelo proprio)...\n")
    errF_r <- learnErrors(fnFs.filt.r, nbases = nbases_par, multithread = TRUE,
                          verbose = TRUE)
    cat("Aprendendo erros R2 —", r, "(corrida auxiliar, modelo proprio)...\n")
    errR_r <- learnErrors(fnRs.filt.r, nbases = nbases_par, multithread = TRUE,
                          verbose = TRUE)

    # PDF do modelo proprio da auxiliar — INSPECIONAR convergencia dos bins
    pdf(file.path(plot_path, paste0("modelo_erros_", r, "_proprio.pdf")), 10, 8)
    suppressWarnings({
      print(plotErrors(errF_r, nominalQ = TRUE) +
              ggtitle(paste("Erros R1 —", r,
                            "(modelo proprio — amostra unica; inspecionar ajuste)")))
      print(plotErrors(errR_r, nominalQ = TRUE) +
              ggtitle(paste("Erros R2 —", r,
                            "(modelo proprio — amostra unica; inspecionar ajuste)")))
    })
    dev.off()
  }
  # --------------------------------------------------------------------------
  # Dereplicate, infer ASVs, merge
  # --------------------------------------------------------------------------
  derepFs_r <- derepFastq(fnFs.filt.r)
  derepRs_r <- derepFastq(fnRs.filt.r)
  if (inherits(derepFs_r, "derep")) {
    derepFs_r <- setNames(list(derepFs_r), sn_r_surv)
  } else {
    names(derepFs_r) <- sn_r_surv
  }
  if (inherits(derepRs_r, "derep")) {
    derepRs_r <- setNames(list(derepRs_r), sn_r_surv)
  } else {
    names(derepRs_r) <- sn_r_surv
  }

  n_uniq_F <- sum(sapply(derepFs_r, function(x) length(x$uniques)))
  n_uniq_R <- sum(sapply(derepRs_r, function(x) length(x$uniques)))
  cat("Sequencias unicas — R1:", n_uniq_F, "| R2:", n_uniq_R, "\n")

  pool_r <- unname(pool_by_run[r])
  if (length(pool_r) != 1L || is.na(pool_r)) {
    stop("Politica de pooling ausente para a corrida: ", r)
  }
  cat("Inferindo ASVs R1 —", r, "(pool=", pool_r, ")...\n")
  dadaFs_r <- dada_with_fallback(derepFs_r, errF_r, pool_r, pool_fallback,
                                 paste0("R1/", r))
  cat("Inferindo ASVs R2 —", r, "(pool=", pool_r, ")...\n")
  dadaRs_r <- dada_with_fallback(derepRs_r, errR_r, pool_r, pool_fallback,
                                 paste0("R2/", r))
  pool_F <- attr(dadaFs_r, "pool_used")
  pool_R <- attr(dadaRs_r, "pool_used")
  if (!identical(pool_F, pool_R)) {
    stop("Fallback de pooling divergente entre R1 e R2 na corrida ", r,
         ": R1=", pool_F, "; R2=", pool_R)
  }
  pool_used[[r]] <- pool_F

  if (length(sn_r_surv) == 1) {
    if (inherits(dadaFs_r, "dada"))
      dadaFs_r <- setNames(list(dadaFs_r), sn_r_surv)
    if (inherits(dadaRs_r, "dada"))
      dadaRs_r <- setNames(list(dadaRs_r), sn_r_surv)
  }

  mergers_r <- mergePairs(dadaFs_r, derepFs_r, dadaRs_r, derepRs_r,
                          justConcatenate = FALSE,
                          minOverlap = minOverlap_par,
                          maxMismatch = maxMismatch_par, verbose = TRUE)
  if (length(sn_r_surv) == 1 && is.data.frame(mergers_r))
    mergers_r <- setNames(list(mergers_r), sn_r_surv)

  merge_zero <- sn_r_surv[sapply(mergers_r, function(x)
    is.data.frame(x) && (nrow(x) == 0 || sum(x$abundance[x$accept]) == 0))]
  if (length(merge_zero) > 0L) {
    stop(
      "Corrida ", r, ": merge zerou em amostra(s); nenhuma sera descartada: ",
      paste(merge_zero, collapse = ", ")
    )
  }

  merge_list[[r]] <- data.frame(
    SampleID  = sn_r_surv, Run = r,
    Mergeadas = sapply(mergers_r, function(x) sum(x$abundance[x$accept])),
    Total     = sapply(dadaFs_r,  function(x) sum(getUniques(x))))
  merge_list[[r]]$Retencao_Merge <-
    round(merge_list[[r]]$Mergeadas / merge_list[[r]]$Total * 100, 1)

  seqtab_list[[r]] <- makeSequenceTable(mergers_r)

  getN <- function(x) sum(getUniques(x))
  out_surv <- out_r[sn_r_surv, , drop = FALSE]
  track_r  <- cbind(input       = out_surv[, "reads.in"],
                    filtradas   = out_surv[, "reads.out"],
                    denoised_R1 = sapply(dadaFs_r, getN),
                    denoised_R2 = sapply(dadaRs_r, getN),
                    merged      = sapply(mergers_r, getN))
  rownames(track_r) <- sn_r_surv
  track_list[[r]] <- track_r

  # Preserva dereps + modelo da auxiliar para a analise de sensibilidade (Secao 15B).
  # A corrida principal libera tudo normalmente.
  if (r == run_auxiliar) {
    derepFs_aux <- derepFs_r
    derepRs_aux <- derepRs_r
    sn_aux_surv <- sn_r_surv
  }
  rm(dadaFs_r, dadaRs_r, mergers_r, errF_r, errR_r)
  if (r == run_principal) rm(derepFs_r, derepRs_r)
  gc()
}

###############################################################################
# 12. UNIAO DAS TABELAS E COLAPSO DE ASVs (todas as corridas)
###############################################################################
cat("\n=== UNINDO TABELAS ===\n")
if (!setequal(names(seqtab_list), runs_ordenados) || length(seqtab_list) != 2L) {
  stop(
    "As duas corridas esperadas nao produziram tabelas de sequencias. Esperadas: ",
    paste(runs_ordenados, collapse = ", "), "; produzidas: ",
    paste(names(seqtab_list), collapse = ", ")
  )
}
seqtab_all <- dada2::mergeSequenceTables(
  tables = unname(seqtab_list[runs_ordenados])
)
cat("Tabela combinada:", nrow(seqtab_all), "amostras x", ncol(seqtab_all), "ASVs\n")

saveRDS(seqtab_all, file.path(output_path, "seqtab_pre_filtro_comprimento.rds"))
length_filter_audit <- data.frame(
  Sequence   = colnames(seqtab_all),
  Length     = nchar(colnames(seqtab_all)),
  TotalReads = as.numeric(colSums(seqtab_all)),
  Prevalence = as.integer(colSums(seqtab_all > 0)),
  stringsAsFactors = FALSE)
write.table(length_filter_audit,
            file.path(output_path, "ASVs_pre_filtro_comprimento.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


###############################################################################
# 13. DISTRIBUICAO DE COMPRIMENTOS E RASTREAMENTO POR LOTE
###############################################################################
cat("\n=== DISTRIBUICAO DE COMPRIMENTOS ===\n")
seq_lengths <- nchar(getSequences(seqtab_all))

# 1. Estatisticas Globais
length_distribution <- table(seq_lengths)
reads_por_comprimento <- tapply(colSums(seqtab_all), seq_lengths, sum)
pico     <- as.integer(names(which.max(reads_por_comprimento)))
pico_asv <- as.integer(names(which.max(length_distribution)))
cat(sprintf("Pico Global -> Por reads: %d pb | Por ASVs: %d pb\n", pico, pico_asv))

# 2. Aplicacao do Filtro de Janela Alvo
keep_len <- seq_lengths >= expected_min_len & seq_lengths <= expected_max_len
seqtab_all.filt <- seqtab_all[, keep_len, drop = FALSE]

if (ncol(seqtab_all.filt) == 0 || sum(seqtab_all.filt) == 0)
  stop("Filtro de comprimento removeu todas as ASVs/reads. Revise expected_min_len/expected_max_len.")

saveRDS(seqtab_all.filt, file.path(output_path, "seqtab_prechimera.rds"))
length_removed <- data.frame(
  Sequence   = colnames(seqtab_all)[!keep_len],
  Length     = nchar(colnames(seqtab_all)[!keep_len]),
  TotalReads = as.numeric(colSums(seqtab_all)[!keep_len]),
  Prevalence = as.integer(colSums(seqtab_all > 0)[!keep_len]),
  Reason     = "Fora_da_janela_de_comprimento",
  stringsAsFactors = FALSE)
write.table(length_removed,
            file.path(output_path, "ASVs_removidas_filtro_comprimento.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(length_removed) > 0L) {
  fasta_length <- Biostrings::DNAStringSet(length_removed$Sequence)
  names(fasta_length) <- paste0("LEN_REMOVED_",
                                sprintf("%06d", seq_len(nrow(length_removed))),
                                "|length=", length_removed$Length,
                                "|reads=", length_removed$TotalReads)
  Biostrings::writeXStringSet(
    fasta_length,
    filepath = file.path(output_path, "ASVs_removidas_filtro_comprimento.fa"),
    format = "fasta")
}

# collapseNoMismatch nao aplicado: ambas as corridas usam truncLen identico
# c(283, 278). Sem truncagem diferencial, nao ha sub-sequencias para colapsar.

# 3. VERIFICACAO DE RETENCAO LOCAL
# Avalia se a amostra isolada ou algum lote especifico sofreu descarte desproporcional
reads_antes_amostra  <- rowSums(seqtab_all)
reads_depois_amostra <- rowSums(seqtab_all.filt)
retencao_amostra     <- round((reads_depois_amostra / reads_antes_amostra) * 100, 1)

amostras_criticas <- names(retencao_amostra[retencao_amostra < 85])
if (length(amostras_criticas) > 0) {
  warning("ALERTA LOCAL: As seguintes amostras perderam >15% das reads no filtro de comprimento:\n  ",
          paste(amostras_criticas, collapse = ", "))
}

write.csv(data.frame(SampleID = names(retencao_amostra), Retencao_Pct = retencao_amostra, Run = run_of_sample[names(retencao_amostra)]),
          file.path(output_path, "retencao_comprimento_por_amostra.csv"), row.names = FALSE)

# 4. Visualizacao Aprimorada (Janela de Corte)
df_plot_len <- data.frame(
  Comprimento = as.numeric(names(reads_por_comprimento)),
  Reads       = as.numeric(reads_por_comprimento)
)

p_len <- ggplot(df_plot_len, aes(x = Comprimento, y = Reads)) +
  geom_col(fill = "#2C7FB8", width = 1) +
  geom_vline(
    xintercept = c(expected_min_len, expected_max_len),
    color = "#D95F02",
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = expected_min_len + 2,
    y = max(df_plot_len$Reads) * 0.92,
    label = paste0("Mínimo: ", expected_min_len, " pb"),
    color = "#D95F02",
    hjust = 0
  ) +
  annotate(
    "text",
    x = expected_max_len - 2,
    y = max(df_plot_len$Reads) * 0.92,
    label = paste0("Máximo: ", expected_max_len, " pb"),
    color = "#D95F02",
    hjust = 1
  ) +
  labs(
    title = "Distribuição de comprimento das ASVs",
    subtitle = sprintf(
      "Janela aceita: %d a %d pb",
      expected_min_len,
      expected_max_len
    ),
    x = "Comprimento (pb)",
    y = "Abundância total de leituras"
  ) +
  theme_classic()

ggsave(file.path(plot_path, "length_distribution.pdf"), p_len, width = 8, height = 5)

# 5. Resumo de Saida
pct_len_reads   <- round(sum(seqtab_all.filt) / sum(seqtab_all) * 100, 1)
pct_len_asv <- round(ncol(seqtab_all.filt) / ncol(seqtab_all) * 100, 1)

cat(sprintf("Filtro aplicado: %d a %d pb\n", expected_min_len, expected_max_len))
cat(sprintf("ASVs retidas globais   : %d de %d (%.1f%%)\n", ncol(seqtab_all.filt), ncol(seqtab_all), pct_len_asv))
cat(sprintf("Reads retidas globais  : %.1f%%\n", pct_len_reads))
if (pct_len_reads < 90) {
  cat("!! ALERTA GLOBAL: < 90% das reads totais foram retidas. Inspecione o PDF.\n")
}

################################################################################
# 14. REMOCAO E AUDITORIA DE QUIMERAS (tabela COMBINADA — 10 amostras)
###############################################################################
cat("\n=== REMOCAO E AUDITORIA DE QUIMERAS (abordagem combinada) ===\n")

# isBimeraDenovoTable: identificacao com vetor logico para auditoria
# Nao aceita "method" — esse parametro e de removeBimeraDenovo
bimeras_logico <- isBimeraDenovoTable(
  seqtab_all.filt,
  minSampleFraction          = 0.9,   # sequencia bimerica em >=90% das amostras
  ignoreNNegatives           = 1,
  minFoldParentOverAbundance = 1.5,   # padrao documentado
  minParentAbundance         = 2,
  allowOneOff                = FALSE,
  multithread                = FALSE,
  verbose                    = TRUE
)

names(bimeras_logico) <- colnames(seqtab_all.filt)
saveRDS(bimeras_logico, file.path(output_path, "chimera_flags_baseline.rds"))
saveRDS(seqtab_all.filt[, bimeras_logico, drop = FALSE],
        file.path(output_path, "seqtab_quimeras_removidas.rds"))
chimera_audit_complete <- data.frame(
  Sequence   = colnames(seqtab_all.filt),
  Length     = nchar(colnames(seqtab_all.filt)),
  TotalReads = as.numeric(colSums(seqtab_all.filt)),
  Prevalence = as.integer(colSums(seqtab_all.filt > 0)),
  IsChimera  = as.logical(bimeras_logico),
  stringsAsFactors = FALSE)
write.table(chimera_audit_complete,
            file.path(output_path, "auditoria_completa_quimeras.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

asvs_removed  <- sum(bimeras_logico)
seqs_quimeras <- colnames(seqtab_all.filt)[bimeras_logico]
seqtab.nochim <- seqtab_all.filt[, !bimeras_logico, drop = FALSE]

if (ncol(seqtab.nochim) == 0 || sum(seqtab.nochim) == 0)
  stop("Remocao de quimeras removeu todas as ASVs/reads. Revise parametros.")

asvs_before <- ncol(seqtab_all.filt)
asvs_after  <- ncol(seqtab.nochim)
pct_nochim  <- round(sum(seqtab.nochim) / sum(seqtab_all.filt) * 100, 1)
# Auditoria das quimeras removidas
if (asvs_removed > 0) {
  auditoria_quimeras <- data.frame(
    Sequence    = seqs_quimeras,
    Total_Reads = colSums(seqtab_all.filt)[seqs_quimeras],
    Prevalencia = colSums(seqtab_all.filt[, seqs_quimeras, drop = FALSE] > 0),
    stringsAsFactors = FALSE)
  write.csv(auditoria_quimeras,
            file.path(output_path, "auditoria_quimeras_removidas.csv"),
            row.names = FALSE)
  cat(sprintf("auditoria_quimeras_removidas.csv: %d sequencias.\n", asvs_removed))
} else {
  cat("Nenhuma quimera detectada.\n")
}

write.csv(data.frame(ASVs_before = asvs_before, ASVs_after = asvs_after,
                     ASVs_removed = asvs_removed,
                     Pct_removed = round(100 * asvs_removed / asvs_before, 1)),
          file.path(output_path, "chimera_stats.csv"), row.names = FALSE)

cat("Retencao pos-quimera:", pct_nochim, "% (esperado 75-95%)\n")

chimera_check <- data.frame(
  SampleID     = rownames(seqtab.nochim),
  Reads_before = rowSums(seqtab_all.filt)[rownames(seqtab.nochim)],
  Reads_after  = rowSums(seqtab.nochim))
chimera_check$Pct_retention <- round(100 * chimera_check$Reads_after / chimera_check$Reads_before, 1)
write.csv(chimera_check, file.path(output_path, "chimera_retention_per_sample.csv"), row.names = FALSE)


if (any(chimera_check$Pct_retention < 70))
  warning("Amostras com retencao pos-quimera <70%.")

sem_quimera <- chimera_check$SampleID[chimera_check$Pct_retention == 100]
if (length(sem_quimera) > 0) {
  cat(sprintf("ATENCAO: %d amostra(s) com 100%% retencao pos-quimera.\n",
              length(sem_quimera)))
  cat("Causa provavel: baixa diversidade unica ou amostra com poucos candidatos\n")
  cat("a sequencias parentais. Declarar como limitacao nos Metodos.\n")
}

# Salvar tabela pre-filtro (combinada, 10 amostras)
saveRDS(seqtab.nochim, file.path(output_path, "seqtab_nochim_prefiltro.rds"))
cat("seqtab_nochim_prefiltro.rds salvo (combinado, pre-filtro — para utilitarios)\n")


# --------------------------------------------------------------------------
# Separacao principal / auxiliar (pos-quimera, pre-filtro de frequencia)
# --------------------------------------------------------------------------
amostras_principal <- metadata$SampleID[metadata$Run == run_principal]
amostras_auxiliar  <- metadata$SampleID[metadata$Run == run_auxiliar]

amostras_main_pres <- intersect(amostras_principal, rownames(seqtab.nochim))
amostras_aux_pres  <- intersect(amostras_auxiliar,  rownames(seqtab.nochim))

seqtab_main_nochim <- seqtab.nochim[amostras_main_pres, , drop = FALSE]
seqtab_aux_nochim  <- seqtab.nochim[amostras_aux_pres,  , drop = FALSE]

seqtab_main_nochim <- seqtab_main_nochim[
  , colSums(seqtab_main_nochim) > 0, drop = FALSE]
seqtab_aux_nochim  <- seqtab_aux_nochim[
  , colSums(seqtab_aux_nochim)  > 0, drop = FALSE]

cat(sprintf("\nSeparacao pos-quimera:\n"))
cat(sprintf("  seqtab_main_nochim : %d amostras | %d ASVs\n",
            nrow(seqtab_main_nochim), ncol(seqtab_main_nochim)))
cat(sprintf("  seqtab_aux_nochim  : %d amostra  | %d ASVs\n",
            nrow(seqtab_aux_nochim),  ncol(seqtab_aux_nochim)))

if (nrow(seqtab_main_nochim) == 0 || ncol(seqtab_main_nochim) == 0)
  stop("Tabela principal vazia apos separacao pos-quimera.")
# A corrida auxiliar integra o universo canônico. Não é permitido seguir com
# tabela degenerada ou remover a amostra silenciosamente.
aux_vazia <- (nrow(seqtab_aux_nochim) != 1L || ncol(seqtab_aux_nochim) == 0L)
if (aux_vazia) {
  stop(
    "Tabela auxiliar inválida após a remoção de quimeras: amostras=",
    nrow(seqtab_aux_nochim), "; ASVs=", ncol(seqtab_aux_nochim),
    ". Esperado: 1 amostra e pelo menos 1 ASV.",
    call. = FALSE
  )
}

###############################################################################
# 14B. PROFUNDIDADE POR AMOSTRA (pos-quimera, pre-filtro)
#
# Usa seqtab.nochim:
#   - depois da remocao de quimeras;
#   - antes dos filtros de abundancia e prevalencia;
#   - inclui as 10 amostras das duas corridas.
###############################################################################

cat("\n=== PROFUNDIDADE POR AMOSTRA (pos-quimera, pre-filtro) ===\n")

# --------------------------------------------------------------------------
# 1. Calculo e validacao
# --------------------------------------------------------------------------

sample_depth <- rowSums(seqtab.nochim)

if (length(sample_depth) == 0L) {
  stop(
    "Nao foi possivel calcular a profundidade: ",
    "seqtab.nochim nao contem amostras."
  )
}

if (is.null(names(sample_depth))) {
  stop(
    "seqtab.nochim nao possui nomes de amostras nas linhas."
  )
}

if (anyNA(sample_depth)) {
  stop(
    "Foram encontrados valores NA na profundidade pos-quimera."
  )
}

if (any(!is.finite(sample_depth))) {
  stop(
    "Foram encontrados valores nao finitos na profundidade pos-quimera."
  )
}

if (any(sample_depth < 0)) {
  stop(
    "Foram encontrados valores negativos na profundidade pos-quimera."
  )
}

amostras_sem_reads <- names(sample_depth)[sample_depth == 0]

if (length(amostras_sem_reads) > 0L) {
  warning(
    "Amostras sem reads apos remocao de quimeras: ",
    paste(amostras_sem_reads, collapse = ", ")
  )
}

# --------------------------------------------------------------------------
# 2. Tabela de profundidade
# --------------------------------------------------------------------------

sample_depth_df <- data.frame(
  SampleID = names(sample_depth),
  SampleLabel = unname(setNames(metadata$SampleLabel, metadata$SampleID)[names(sample_depth)]),
  Run = unname(run_of_sample[names(sample_depth)]),
  Reads_pos_quimera = as.numeric(sample_depth),
  stringsAsFactors = FALSE
)

if (anyNA(sample_depth_df$SampleLabel) || any(sample_depth_df$SampleLabel == "")) {
  stop("SampleLabel nao encontrado para uma ou mais amostras da tabela de profundidade.")
}

if (anyNA(sample_depth_df$Run)) {
  amostras_sem_run <- sample_depth_df$SampleID[
    is.na(sample_depth_df$Run)
  ]

  stop(
    "Corrida nao encontrada para a(s) amostra(s): ",
    paste(amostras_sem_run, collapse = ", ")
  )
}

sample_depth_df <- sample_depth_df[
  order(sample_depth_df$Reads_pos_quimera),
  ,
  drop = FALSE
]

write.csv(
  sample_depth_df,
  file.path(output_path, "sample_depth_prefiltro.csv"),
  row.names = FALSE
)

cat("Profundidade pos-quimera por amostra:\n")
print(sample_depth_df, row.names = FALSE)

cat("\nResumo da profundidade pos-quimera:\n")
print(summary(sample_depth_df$Reads_pos_quimera))

cat(
  sprintf(
    paste0(
      "Minima: %d | Mediana: %.1f | Media: %.1f | Maxima: %d\n"
    ),
    min(sample_depth_df$Reads_pos_quimera),
    median(sample_depth_df$Reads_pos_quimera),
    mean(sample_depth_df$Reads_pos_quimera),
    max(sample_depth_df$Reads_pos_quimera)
  )
)

# --------------------------------------------------------------------------
# 3. Histograma geral
# --------------------------------------------------------------------------

hist_file <- file.path(
  plot_path,
  "sample_depth_hist.pdf"
)

grDevices::pdf(
  hist_file,
  width = 8,
  height = 6
)

tryCatch(
  {
    hist(
      sample_depth_df$Reads_pos_quimera,
      main = "Profundidade por amostra apos remocao de quimeras",
      xlab = "Reads pos-quimera",
      ylab = "Numero de amostras",
      breaks = "Sturges"
    )

    abline(
      v = median(sample_depth_df$Reads_pos_quimera),
      lty = 2,
      lwd = 2
    )

    legend(
      "topright",
      legend = paste0(
        "Mediana = ",
        format(
          median(sample_depth_df$Reads_pos_quimera),
          big.mark = ".",
          scientific = FALSE,
          trim = TRUE
        )
      ),
      lty = 2,
      lwd = 2,
      bty = "n"
    )
  },
  finally = {
    grDevices::dev.off()
  }
)

# --------------------------------------------------------------------------
# 4. Grafico por amostra
# --------------------------------------------------------------------------

barplot_file <- file.path(
  plot_path,
  "sample_depth_por_amostra.pdf"
)

grDevices::pdf(
  barplot_file,
  width = 13,
  height = 8
)

tryCatch(
  {
    # Aumenta a margem esquerda para comportar os nomes completos
    par(
      mar = c(5, 14, 4, 2) + 0.1,
      mgp = c(2.5, 0.8, 0)
    )

    barplot(
      height = sample_depth_df$Reads_pos_quimera,
      names.arg = sample_depth_df$SampleLabel,
      horiz = TRUE,
      las = 1,
      cex.names = 0.75,
      xlab = "Reads pós-quimera",
      ylab = "",
      main = paste0(
        "Profundidade por amostra ",
        "(pós-quimera, pré-filtro)"
      ),
      xlim = c(
        0,
        max(sample_depth_df$Reads_pos_quimera) * 1.08
      )
    )

    abline(
      v = median(sample_depth_df$Reads_pos_quimera),
      lty = 2,
      lwd = 2
    )

    legend(
      "bottomright",
      legend = paste0(
        "Mediana = ",
        format(
          median(sample_depth_df$Reads_pos_quimera),
          big.mark = ".",
          decimal.mark = ",",
          scientific = FALSE,
          trim = TRUE
        ),
        " reads"
      ),
      lty = 2,
      lwd = 2,
      bty = "n"
    )
  },
  finally = {
    grDevices::dev.off()
  }
)

cat("Arquivos gerados:\n")
cat(" - ", hist_file, "\n", sep = "")
cat(" - ", barplot_file, "\n", sep = "")
cat(
  " - ",
  file.path(output_path, "sample_depth_prefiltro.csv"),
  "\n",
  sep = ""
)
###############################################################################
# 14C. ANALISE DE SENSIBILIDADE DOS FILTROS
# Executada sobre seqtab_main_nochim (9 amostras) — os limiares sao
# selecionados para a analise primaria. min_prev avaliado dentro do n=9.
###############################################################################
cat("\n=== SENSIBILIDADE DOS FILTROS (corrida principal, n=9) ===\n")
evaluate_threshold <- function(seqtab, min_reads, min_prev) {
  total_reads <- colSums(seqtab)
  prevalence  <- colSums(seqtab > 0)   # n de amostras com reads > 0
  keep    <- !(total_reads < min_reads & prevalence < min_prev)
  seqtab2 <- seqtab[, keep, drop = FALSE]
  if (sum(seqtab) == 0) stop("Tabela vazia na analise de sensibilidade.")
  data.frame(
    min_reads     = min_reads,
    min_prev      = min_prev,
    ASVs          = ncol(seqtab2),
    Reads         = sum(seqtab2),
    Pct_reads     = round(sum(seqtab2) / sum(seqtab) * 100, 2),
    Mean_Richness = mean(rowSums(seqtab2 > 0)),
    Mean_Shannon  = mean(vegan::diversity(seqtab2, index = "shannon")))
}
grid <- expand.grid(min_reads = c(1, 5, 10, 20), min_prev = c(1, 2, 3))
sensitivity_results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  evaluate_threshold(seqtab_main_nochim, grid$min_reads[i], grid$min_prev[i])))
write.csv(sensitivity_results,
          file.path(output_path, "sensitivity_thresholds.csv"), row.names = FALSE)
cat("Sensibilidade sobre corrida principal (", nrow(seqtab_main_nochim),
    " amostras):\n")
print(sensitivity_results)

###############################################################################
# 15. FILTRO DE PREVALENCIA + ABUNDANCIA
#
# Corrida principal (seqtab_main_nochim, n=9):
#   logica E: descartar se (reads < min_total_reads) E (prevalencia < min_prevalence)
#   prevalencia = numero de amostras com reads > 0
#   min_prevalence = 2 aplicado sobre 9 amostras (~22%)
#
# Corrida auxiliar (seqtab_aux_nochim, n=1):
#   filtro por leitura minima apenas (prevalencia maxima = 1; min_prev invalido)
###############################################################################
cat("\n=== FILTRO DE PREVALENCIA + ABUNDANCIA ===\n")

# --- Corrida principal ---
cat(sprintf("Corrida principal (n=%d amostras):\n", nrow(seqtab_main_nochim)))
cat(sprintf("  Criterio: descartar se (reads < %d) E (prevalencia < %d)\n",
            min_total_reads, min_prevalence))

asv_reads_main <- colSums(seqtab_main_nochim)
asv_prev_main  <- colSums(seqtab_main_nochim > 0)   # > 0, nao > 2
descartar_main <- (asv_reads_main < min_total_reads) & (asv_prev_main < min_prevalence)
manter_main    <- !descartar_main

# Usa colnames() enquanto a tabela ainda contem todas as colunas (pre-filtro).
asvs_desc_principal_seq <- data.frame(
  Sequence    = colnames(seqtab_main_nochim)[descartar_main],
  Length      = nchar(colnames(seqtab_main_nochim)[descartar_main]),
  Reads_total = as.numeric(asv_reads_main[descartar_main]),
  Prevalencia = as.integer(asv_prev_main[descartar_main]),
  Reason      = rep(paste0("Reads_total < ", min_total_reads,
                           " E Prevalencia < ", min_prevalence),
                    sum(descartar_main)),
  stringsAsFactors = FALSE)
asvs_desc_principal_seq <- asvs_desc_principal_seq[
  order(asvs_desc_principal_seq$Reads_total, decreasing = TRUE), , drop = FALSE]
write.table(asvs_desc_principal_seq,
            file.path(output_path, "ASVs_descartadas_principal_com_sequencia.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- CHECKPOINT DE FILTRO: interrupcao obrigatoria -------------------------
if (is.na(CENARIO_FILTRO)) {
  cat("\n===================== CHECKPOINT DE FILTRO =====================\n")
  cat(sprintf("Limiares atuais: descartar se reads < %d E prevalencia < %d.\n",
              min_total_reads, min_prevalence))
  cat(sprintf("Seriam descartadas %d de %d ASVs na corrida principal.\n",
              sum(descartar_main), length(descartar_main)))
  cat("Auditoria: ASVs_descartadas_principal_com_sequencia.tsv\n")
  cat("Para comparar cenarios, rode 01b_auditoria_quimeras_recuperacao_ASVs.R\n")
  cat("(le seqtab_nochim_prefiltro.rds, ja salvo).\n")
  cat("Depois, edite CENARIO_FILTRO ('aplicar' ou 'sem_filtro') no topo\n")
  cat("deste script e execute-o novamente.\n")
  cat("================================================================\n")
  stop("Filtro nao confirmado. Defina CENARIO_FILTRO e re-execute.",
       call. = FALSE)
} else if (identical(CENARIO_FILTRO, "sem_filtro")) {
  manter_main    <- rep(TRUE, length(manter_main))
  descartar_main <- !manter_main
  cat("CENARIO_FILTRO = 'sem_filtro': nenhuma ASV removida por frequencia (principal).\n")
} else if (!identical(CENARIO_FILTRO, "aplicar")) {
  stop("CENARIO_FILTRO invalido: use 'aplicar' ou 'sem_filtro'.", call. = FALSE)
}
# --- fim CHECKPOINT ---

asvs_main_antes  <- rowSums(seqtab_main_nochim > 0)
reads_main_antes <- sum(seqtab_main_nochim)
seqtab_main_nochim <- seqtab_main_nochim[, manter_main, drop = FALSE]
pct_reads_main <- round(sum(seqtab_main_nochim) / reads_main_antes * 100, 1)

cat(sprintf("  ASVs antes: %d | descarte: %d | retidas: %d\n",
            length(manter_main), sum(descartar_main), sum(manter_main)))
cat(sprintf("  Reads retidas: %.1f%%\n", pct_reads_main))

amostras_vazias_main <- rownames(seqtab_main_nochim)[
  rowSums(seqtab_main_nochim) == 0]
if (length(amostras_vazias_main) > 0) {
  stop(
    "Amostra(s) da corrida principal ficaram com zero reads apos o filtro: ",
    paste(amostras_vazias_main, collapse = ", "),
    ". O pipeline nao removera amostras silenciosamente. Reveja os filtros."
  )
}

# --- Corrida auxiliar ---
cat(sprintf("\nCorrida auxiliar (n=%d amostra):\n", nrow(seqtab_aux_nochim)))
cat(sprintf("  Criterio: descartar se (reads < %d) [min_prevalence nao aplicavel]\n",
            min_total_reads))

if (nrow(seqtab_aux_nochim) > 0 && ncol(seqtab_aux_nochim) > 0 && sum(seqtab_aux_nochim) > 0) {
  asv_reads_aux <- colSums(seqtab_aux_nochim)
  descartar_aux <- asv_reads_aux < min_total_reads
  manter_aux    <- !descartar_aux
  if (identical(CENARIO_FILTRO, "sem_filtro")) {
    manter_aux    <- rep(TRUE, length(manter_aux))
    descartar_aux <- !manter_aux
    cat("CENARIO_FILTRO = 'sem_filtro': nenhuma ASV removida por frequencia (auxiliar).\n")
  }

  asvs_desc_aux_seq <- data.frame(
    Sequence    = colnames(seqtab_aux_nochim)[descartar_aux],
    Length      = nchar(colnames(seqtab_aux_nochim)[descartar_aux]),
    Reads_total = as.numeric(asv_reads_aux[descartar_aux]),
    Prevalencia = as.integer(colSums(seqtab_aux_nochim > 0)[descartar_aux]),
    Reason      = rep(paste0("Reads_total < ", min_total_reads,
                             " na corrida auxiliar"), sum(descartar_aux)),
    stringsAsFactors = FALSE)
  write.table(asvs_desc_aux_seq,
              file.path(output_path, "ASVs_descartadas_auxiliar_com_sequencia.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  reads_aux_antes <- sum(seqtab_aux_nochim)
  seqtab_aux_nochim <- seqtab_aux_nochim[, manter_aux, drop = FALSE]
  pct_reads_aux <- round(sum(seqtab_aux_nochim) / reads_aux_antes * 100, 1)
} else {
  stop(
    "Tabela auxiliar vazia antes do filtro de abundancia. ",
    "A amostra S10 nao pode ser removida silenciosamente.",
    call. = FALSE
  )
}

cat(sprintf("  ASVs antes: %d | descarte: %d | retidas: %d\n",
            length(manter_aux), sum(descartar_aux), sum(manter_aux)))
cat(sprintf("  Reads retidas: %.1f%%\n", pct_reads_aux))

amostras_vazias_aux <- rownames(seqtab_aux_nochim)[
  rowSums(seqtab_aux_nochim) == 0]
if (length(amostras_vazias_aux) > 0) {
  stop(
    "A amostra da corrida auxiliar ficou com zero reads apos o filtro: ",
    paste(amostras_vazias_aux, collapse = ", "),
    ". O pipeline nao removera essa amostra silenciosamente; reveja o limiar.",
    call. = FALSE
  )
}

# Auditoria das ASVs descartadas da corrida principal
asvs_desc_df <- data.frame(
  ASV_idx     = which(descartar_main),
  Reads_total = asv_reads_main[descartar_main],
  Prevalencia = asv_prev_main[descartar_main],
  stringsAsFactors = FALSE)
asvs_desc_df <- asvs_desc_df[order(asvs_desc_df$Reads_total, decreasing = TRUE), ]
limiar_alerta <- min_total_reads * 10
asvs_alerta   <- asvs_desc_df[asvs_desc_df$Reads_total >= limiar_alerta, ]

cat(sprintf("\n=== AUDITORIA DE ASVs DESCARTADAS (corrida principal) ===\n"))
cat(sprintf("Total descartadas: %d | Com reads >= %d: %d\n",
            nrow(asvs_desc_df), limiar_alerta, nrow(asvs_alerta)))
if (nrow(asvs_alerta) > 0) {
  cat(sprintf("ATENCAO: %d ASV(s) com >= %d reads — inspecionar manualmente:\n",
              nrow(asvs_alerta), limiar_alerta))
  print(head(asvs_alerta, 20))
}
write.csv(asvs_desc_df,
          file.path(output_path, "asvs_descartadas_auditoria.csv"),
          row.names = FALSE)

# Riqueza antes/depois por amostra (corrida principal)
asvs_main_depois <- rowSums(seqtab_main_nochim > 0)
amostras_comuns  <- intersect(names(asvs_main_antes), names(asvs_main_depois))
riqueza_check <- data.frame(
  SampleID      = amostras_comuns,
  ASVs_antes    = asvs_main_antes[amostras_comuns],
  ASVs_depois   = asvs_main_depois[amostras_comuns])
riqueza_check$ASVs_removidas <- riqueza_check$ASVs_antes - riqueza_check$ASVs_depois
riqueza_check$Pct_retido <- round(
  riqueza_check$ASVs_depois / riqueza_check$ASVs_antes * 100, 1)
print(riqueza_check, row.names = FALSE)
write.csv(riqueza_check,
          file.path(output_path, "riqueza_por_amostra_filtro.csv"),
          row.names = FALSE)
if (any(riqueza_check$Pct_retido < 50))
  warning("Amostra(s) perderam >50% das ASVs no filtro de frequencia.")

# Profundidade final — corrida principal (base para rarefacao no Script 8)
if (nrow(seqtab_main_nochim) == 0 || ncol(seqtab_main_nochim) == 0)
  stop("Tabela principal vazia apos filtro de prevalencia/abundancia.")

depth_final_main <- rowSums(seqtab_main_nochim)
cat(sprintf("\nProfundidade final — corrida principal (seqtab_nochim.rds):\n"))
print(sort(depth_final_main))
cat(sprintf("Minima: %d | Mediana: %d | Maxima: %d\n",
            min(depth_final_main), as.integer(median(depth_final_main)),
            max(depth_final_main)))
min_idx_fin <- which.min(depth_final_main)
cat(sprintf("Amostra limitante (rarefacao): %s (%d reads)\n",
            names(depth_final_main)[min_idx_fin], depth_final_main[min_idx_fin]))

depth_final_aux <- rowSums(seqtab_aux_nochim)
cat(sprintf("\nProfundidade final — corrida auxiliar (seqtab_auxiliar.rds):\n"))
print(depth_final_aux)


###############################################################################
# 15B. SENSIBILIDADE DO MODELO DE ERRO DA AUXILIAR
#
# Compara:
#   1. modelo proprio da corrida auxiliar;
#   2. modelo transferido da corrida principal.
#
# Os dois cenarios repetem:
#   dada -> merge -> comprimento -> quimeras -> filtro de abundancia.
#
# Nenhuma saida primaria e sobrescrita.
###############################################################################

cat(
  "\n=== SENSIBILIDADE COMPLETA: ",
  "MODELO PROPRIO vs. TRANSFERIDO ===\n",
  sep = ""
)

sens_ok <- (
  exists("derepFs_aux") &&
    exists("derepRs_aux") &&
    !is.null(errF_principal) &&
    !is.null(errR_principal) &&
    !isTRUE(aux_vazia)
)

processar_cenario_sensibilidade <- function(
    seqtab_aux_pre,
    nome_modelo) {

  seqtab_comb <- dada2::mergeSequenceTables(
    seqtab_list[[run_principal]],
    seqtab_aux_pre
  )

  comprimentos <- nchar(
    dada2::getSequences(seqtab_comb)
  )

  keep_len_sens <- (
    comprimentos >= expected_min_len &
      comprimentos <= expected_max_len
  )

  seqtab_len <- seqtab_comb[
    ,
    keep_len_sens,
    drop = FALSE
  ]

  if (ncol(seqtab_len) == 0L || sum(seqtab_len) == 0) {
    stop(
      "Cenario ",
      nome_modelo,
      ": tabela vazia apos filtro de comprimento."
    )
  }

  bimeras_sens <- dada2::isBimeraDenovoTable(
    seqtab_len,
    minSampleFraction = 0.9,
    ignoreNNegatives = 1,
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    multithread = FALSE,
    verbose = TRUE
  )

  seqtab_nochim_sens <- seqtab_len[
    ,
    !bimeras_sens,
    drop = FALSE
  ]

  main_sens <- seqtab_nochim_sens[
    intersect(
      amostras_principal,
      rownames(seqtab_nochim_sens)
    ),
    ,
    drop = FALSE
  ]

  aux_sens <- seqtab_nochim_sens[
    intersect(
      amostras_auxiliar,
      rownames(seqtab_nochim_sens)
    ),
    ,
    drop = FALSE
  ]

  main_sens <- main_sens[
    ,
    colSums(main_sens) > 0,
    drop = FALSE
  ]

  aux_sens <- aux_sens[
    ,
    colSums(aux_sens) > 0,
    drop = FALSE
  ]

  main_pre_filtro <- main_sens
  aux_pre_filtro  <- aux_sens

  if (ncol(main_sens) > 0L) {
    reads_main_sens <- colSums(main_sens)
    prev_main_sens  <- colSums(main_sens > 0)

    descartar_main_sens <- (
      reads_main_sens < min_total_reads &
        prev_main_sens < min_prevalence
    )

    main_sens <- main_sens[
      ,
      !descartar_main_sens,
      drop = FALSE
    ]
  }

  if (ncol(aux_sens) > 0L) {
    aux_sens <- aux_sens[
      ,
      colSums(aux_sens) >= min_total_reads,
      drop = FALSE
    ]
  }

  list(
    Modelo = nome_modelo,
    seqtab_combinada_prequimera = seqtab_comb,
    seqtab_main_final = main_sens,
    seqtab_aux_posquimera = aux_pre_filtro,
    seqtab_aux_final = aux_sens,
    N_ASVs_aux_prequimera = ncol(seqtab_aux_pre),
    Reads_aux_prequimera = sum(seqtab_aux_pre),
    N_ASVs_aux_posquimera = ncol(aux_pre_filtro),
    Reads_aux_posquimera = sum(aux_pre_filtro),
    N_ASVs_aux_final = ncol(aux_sens),
    Reads_aux_final = sum(aux_sens),
    N_ASVs_main_final = ncol(main_sens),
    Reads_main_final = sum(main_sens)
  )
}

if (sens_ok) {

  pool_aux <- unname(pool_by_run[run_auxiliar])

  dadaF_transf <- dada2::dada(
    derepFs_aux,
    err = errF_principal,
    pool = pool_aux,
    multithread = TRUE
  )

  dadaR_transf <- dada2::dada(
    derepRs_aux,
    err = errR_principal,
    pool = pool_aux,
    multithread = TRUE
  )

  if (inherits(dadaF_transf, "dada")) {
    dadaF_transf <- setNames(
      list(dadaF_transf),
      sn_aux_surv
    )
  }

  if (inherits(dadaR_transf, "dada")) {
    dadaR_transf <- setNames(
      list(dadaR_transf),
      sn_aux_surv
    )
  }

  mergers_transf <- dada2::mergePairs(
    dadaF_transf,
    derepFs_aux,
    dadaR_transf,
    derepRs_aux,
    minOverlap = minOverlap_par,
    maxMismatch = maxMismatch_par,
    justConcatenate = FALSE,
    verbose = TRUE
  )

  if (is.data.frame(mergers_transf)) {
    mergers_transf <- setNames(
      list(mergers_transf),
      sn_aux_surv
    )
  }

  seqtab_aux_proprio_pre <- seqtab_list[[run_auxiliar]]



  seqtab_aux_transf_pre <- dada2::makeSequenceTable(
    mergers_transf
  )

  cenario_proprio <- processar_cenario_sensibilidade(
    seqtab_aux_pre = seqtab_aux_proprio_pre,
    nome_modelo = "proprio"
  )

  cenario_transferido <- processar_cenario_sensibilidade(
    seqtab_aux_pre = seqtab_aux_transf_pre,
    nome_modelo = "transferido"
  )

  resumo_sens <- data.frame(
    Modelo = c(
      cenario_proprio$Modelo,
      cenario_transferido$Modelo
    ),
    ASVs_prequimera = c(
      cenario_proprio$N_ASVs_aux_prequimera,
      cenario_transferido$N_ASVs_aux_prequimera
    ),
    Reads_prequimera = c(
      cenario_proprio$Reads_aux_prequimera,
      cenario_transferido$Reads_aux_prequimera
    ),
    ASVs_posquimera = c(
      cenario_proprio$N_ASVs_aux_posquimera,
      cenario_transferido$N_ASVs_aux_posquimera
    ),
    Reads_posquimera = c(
      cenario_proprio$Reads_aux_posquimera,
      cenario_transferido$Reads_aux_posquimera
    ),
    ASVs_final = c(
      cenario_proprio$N_ASVs_aux_final,
      cenario_transferido$N_ASVs_aux_final
    ),
    Reads_final = c(
      cenario_proprio$Reads_aux_final,
      cenario_transferido$Reads_aux_final
    ),
    ASVs_main_final = c(
      cenario_proprio$N_ASVs_main_final,
      cenario_transferido$N_ASVs_main_final
    ),
    Reads_main_final = c(
      cenario_proprio$Reads_main_final,
      cenario_transferido$Reads_main_final
    ),
    stringsAsFactors = FALSE
  )

  seqs_proprio <- colnames(
    cenario_proprio$seqtab_aux_final
  )

  seqs_transferido <- colnames(
    cenario_transferido$seqtab_aux_final
  )

  seqs_uniao <- union(
    seqs_proprio,
    seqs_transferido
  )

  seqs_compartilhadas_sens <- intersect(
    seqs_proprio,
    seqs_transferido
  )

  jaccard_sens <- if (length(seqs_uniao) > 0L) {
    length(seqs_compartilhadas_sens) /
      length(seqs_uniao)
  } else {
    NA_real_
  }

  abund_proprio <- setNames(
    numeric(length(seqs_uniao)),
    seqs_uniao
  )

  abund_transferido <- setNames(
    numeric(length(seqs_uniao)),
    seqs_uniao
  )

  if (length(seqs_proprio) > 0L) {
    abund_proprio[seqs_proprio] <- colSums(
      cenario_proprio$seqtab_aux_final
    )
  }

  if (length(seqs_transferido) > 0L) {
    abund_transferido[seqs_transferido] <- colSums(
      cenario_transferido$seqtab_aux_final
    )
  }

  bray_curtis_sens <- if (
    sum(abund_proprio + abund_transferido) > 0
  ) {
    sum(abs(abund_proprio - abund_transferido)) /
      sum(abund_proprio + abund_transferido)
  } else {
    NA_real_
  }

  auditoria_sens <- data.frame(
    Sequence = seqs_uniao,
    Reads_modelo_proprio = as.numeric(
      abund_proprio[seqs_uniao]
    ),
    Reads_modelo_transferido = as.numeric(
      abund_transferido[seqs_uniao]
    ),
    Compartilhada = seqs_uniao %in%
      seqs_compartilhadas_sens,
    stringsAsFactors = FALSE
  )

  metricas_concordancia <- data.frame(
    ASVs_compartilhadas =
      length(seqs_compartilhadas_sens),
    ASVs_exclusivas_proprio =
      length(setdiff(
        seqs_proprio,
        seqs_transferido
      )),
    ASVs_exclusivas_transferido =
      length(setdiff(
        seqs_transferido,
        seqs_proprio
      )),
    Jaccard = jaccard_sens,
    Bray_Curtis = bray_curtis_sens,
    stringsAsFactors = FALSE
  )

  saveRDS(
    seqtab_aux_proprio_pre,
    file.path(
      output_path,
      "seqtab_aux_modelo_proprio_prequimera.rds"
    )
  )

  saveRDS(
    seqtab_aux_transf_pre,
    file.path(
      output_path,
      "seqtab_aux_modelo_transferido_prequimera.rds"
    )
  )

  saveRDS(
    cenario_proprio$seqtab_aux_final,
    file.path(
      output_path,
      "seqtab_aux_modelo_proprio_final.rds"
    )
  )

  saveRDS(
    cenario_transferido$seqtab_aux_final,
    file.path(
      output_path,
      "seqtab_aux_modelo_transferido_final.rds"
    )
  )

  write.csv(
    resumo_sens,
    file.path(
      output_path,
      "sensibilidade_modelo_auxiliar_resumo.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    metricas_concordancia,
    file.path(
      output_path,
      "sensibilidade_modelo_auxiliar_concordancia.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    auditoria_sens,
    file.path(
      output_path,
      "sensibilidade_modelo_auxiliar_ASVs.csv"
    ),
    row.names = FALSE
  )

  cat("\nResumo da sensibilidade:\n")
  print(resumo_sens, row.names = FALSE)

  cat("\nConcordancia final:\n")
  print(metricas_concordancia, row.names = FALSE)

  rm(
    dadaF_transf,
    dadaR_transf,
    mergers_transf,
    seqtab_aux_transf_pre
  )

} else {
  warning(
    "Sensibilidade nao executada: ",
    "dereps, modelos ou tabela auxiliar indisponiveis."
  )
}

rm(
  list = intersect(
    c("derepFs_aux", "derepRs_aux"),
    ls()
  )
)

gc()
###############################################################################
# 16. RASTREAMENTO GLOBAL DE READS
# "final" = reads na saida filtrada de cada amostra (main ou aux)
###############################################################################
cat("\n=== RASTREAMENTO DE READS ===\n")

final_per_sample <- c(rowSums(seqtab_main_nochim),
                      rowSums(seqtab_aux_nochim))

track_all <- do.call(rbind, track_list)
# Incluir todas as amostras que passaram pelo pipeline; 0 se filtradas inteiras
all_track_samples <- rownames(track_all)
final_col <- setNames(rep(0L, length(all_track_samples)), all_track_samples)
final_col[names(final_per_sample)] <- as.integer(final_per_sample)
track_all <- cbind(track_all, final = final_col)

raw_col   <- reads_raw_R1[rownames(track_all)]
track_all <- cbind(raw = raw_col, track_all)

track_pct <- round(sweep(track_all, 1, track_all[, "raw"], "/") * 100, 1)
cat("\nAbsoluto:\n"); print(track_all)
cat("\nPercentual (em relacao as reads brutas):\n"); print(track_pct)
global_retention <- round(sum(track_all[, "final"]) /
                            sum(track_all[, "raw"], na.rm = TRUE) * 100, 1)
cat("\nRetencao global (bruto->final, todas as amostras):", global_retention, "%\n")

filter_stats <- do.call(rbind, filter_list)
merge_stats  <- do.call(rbind, merge_list)

# Retencao sequencial: cada etapa sobre a IMEDIATAMENTE anterior.
# Usa track_all (denominadores encadeados), nao a subtracao de percentuais
# com bases distintas.
retencao_seq <- data.frame(
  SampleID            = rownames(track_all),
  Run                 = run_of_sample[rownames(track_all)],
  bruto_para_input    = round(100 * track_all[, "input"]     / track_all[, "raw"],       1),
  input_para_filtrado = round(100 * track_all[, "filtradas"] / track_all[, "input"],     1),
  filtrado_para_denoised = round(100 * track_all[, "denoised_R1"] / track_all[, "filtradas"], 1),
  denoised_para_merge = round(100 * track_all[, "merged"]    / track_all[, "denoised_R1"], 1),
  merge_para_final    = round(100 * track_all[, "final"]     / track_all[, "merged"],     1),
  bruto_para_final    = round(100 * track_all[, "final"]     / track_all[, "raw"],        1),
  row.names = NULL
)
write.csv(retencao_seq,
          file.path(output_path, "retencao_sequencial.csv"), row.names = FALSE)

# Alerta sobre a etapa de merge com denominador denoised -> merge.
quedas_merge <- retencao_seq$SampleID[
  !is.na(retencao_seq$denoised_para_merge) & retencao_seq$denoised_para_merge < 70]
if (length(quedas_merge) > 0)
  warning("Amostras com retencao denoised->merge < 70%: ",
          paste(quedas_merge, collapse = ", "))

###############################################################################
# 17. CONTAMINANTES CONHECIDOS
###############################################################################
contaminant_genera_ref <- c(
  "Ralstonia", "Pseudomonas", "Bradyrhizobium", "Cutibacterium",
  "Propionibacterium", "Acinetobacter", "Sphingomonas", "Stenotrophomonas",
  "Delftia", "Comamonas", "Methylobacterium", "Herbaspirillum",
  "Burkholderia", "Escherichia", "Shigella", "Staphylococcus")
writeLines(contaminant_genera_ref,
           file.path(output_path, "REF_generos_contaminacao_conhecida.txt"))
cat("\n[Sem controle negativo] Lista de contaminantes exportada.\n")

###############################################################################
# 18. DICIONARIO GLOBAL DE IDs, FASTA UNIFICADO E TABELAS
#
# DESIGN:
#   O universo de ASVs e particionado em tres subconjuntos:
#     (a) principal          : colnames(seqtab_main_nochim)  → IDs ASV_1…N
#     (b) compartilhadas     : intersecao entre principal e auxiliar
#     (c) exclusivas_aux     : colnames(seqtab_aux_nochim) \ principal → IDs ASV_AUX_1…M
#
#   O FASTA global (ASVs.fa) contem (a) + (c), garantindo que o Script 4
#   (rBLAST) processe a totalidade das sequencias de ambas as corridas.
#   O prefixo ASV_AUX_ distingue as exclusivas da auxiliar sem alterar os
#   IDs ASV_N da corrida principal. Os scripts a jusante preservam a unidade ASV.
#
#   Scripts 2a-2e e 4: operam sobre seqtab_global_nochim.rds (main + aux),
#   classificando TODAS as ASVs em uma unica passagem. O particionamento por
#   origem (principal vs. auxiliar_exclusiva) e preservado em ASV_sequences.tsv
#   e pode ser auditado por todos os scripts a jusante.
###############################################################################
cat("\n=== DICIONARIO GLOBAL DE IDs, FASTA E TABELAS ===\n")

asv_seqs_main    <- colnames(seqtab_main_nochim)
asv_seqs_aux_all <- colnames(seqtab_aux_nochim)

# Particionamento
asv_seqs_exclusivas_aux  <- setdiff(asv_seqs_aux_all, asv_seqs_main)
asv_seqs_compartilhadas  <- intersect(asv_seqs_aux_all, asv_seqs_main)

cat(sprintf("ASVs corrida principal             : %d\n", length(asv_seqs_main)))
cat(sprintf("ASVs corrida auxiliar (total)      : %d\n", length(asv_seqs_aux_all)))
cat(sprintf("  compartilhadas com principal     : %d\n", length(asv_seqs_compartilhadas)))
cat(sprintf("  exclusivas da auxiliar           : %d\n", length(asv_seqs_exclusivas_aux)))

# Esquema de IDs:
#   ASV_1 … ASV_N      : corrida principal (retrocompativel com Scripts 2-8)
#   ASV_AUX_1 … ASV_AUX_M : exclusivas da auxiliar
asv_ids_main         <- paste0("ASV_", seq_along(asv_seqs_main))
asv_ids_excl_aux     <- if (length(asv_seqs_exclusivas_aux) > 0)
  paste0("ASV_AUX_", seq_along(asv_seqs_exclusivas_aux)) else character(0)

# Dicionario global: sequencia -> ID
seq2id_global <- c(
  setNames(asv_ids_main,     asv_seqs_main),
  setNames(asv_ids_excl_aux, asv_seqs_exclusivas_aux))

# Renomear colunas dos seqtabs com IDs consistentes
seqtab_main_nochim_ids           <- seqtab_main_nochim
colnames(seqtab_main_nochim_ids) <- asv_ids_main

seqtab_aux_nochim_ids            <- seqtab_aux_nochim
colnames(seqtab_aux_nochim_ids)  <- seq2id_global[asv_seqs_aux_all]
# Blindagem: todo colname da auxiliar deve ter ID resolvido. NA aqui faria
# pivot_longer descartar a coluna silenciosamente, perdendo ASVs da tabela long.
if (ncol(seqtab_aux_nochim_ids) > 0 && anyNA(colnames(seqtab_aux_nochim_ids)))
  stop("IDs nao resolvidos para ", sum(is.na(colnames(seqtab_aux_nochim_ids))),
       " ASV(s) da auxiliar. Verifique seq2id_global.")

# --------------------------------------------------------------------------
# FASTA GLOBAL E LOCAL via writeFasta (dada2 p.63)
# --------------------------------------------------------------------------
all_seqs_global <- c(asv_seqs_main, asv_seqs_exclusivas_aux)
all_ids_global  <- c(asv_ids_main,  asv_ids_excl_aux)

# Fasta Global (Entrada Script 4 - rBLAST)
writeLines(c(rbind(paste0(">", all_ids_global), all_seqs_global)),
           file.path(output_path, "ASVs.fa"))

# FASTA restrito corrida principal
writeLines(c(rbind(paste0(">", asv_ids_main), asv_seqs_main)),
           file.path(output_path, "ASVs_principal.fa"))

# FASTA restrito corrida auxiliar (se houver exclusivas)
if (length(asv_seqs_exclusivas_aux) > 0)
  writeLines(c(rbind(paste0(">", asv_ids_excl_aux), asv_seqs_exclusivas_aux)),
             file.path(output_path, "ASVs_auxiliar_exclusivas.fa"))

# --------------------------------------------------------------------------
# Tabela de sequencias global com anotacao de origem
# --------------------------------------------------------------------------
asv_headers_global <- data.frame(
  ASV_ID   = all_ids_global,
  Sequence = all_seqs_global,
  Length   = nchar(all_seqs_global),
  Origem   = c(rep("principal",        length(asv_seqs_main)),
               rep("auxiliar_exclusiva", length(asv_seqs_exclusivas_aux))),
  stringsAsFactors = FALSE)
write.table(asv_headers_global,
            file.path(output_path, "ASV_sequences.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --------------------------------------------------------------------------
# Tabela de contagens de ASVs e tabela long — corrida principal
# --------------------------------------------------------------------------
write.table(seqtab_main_nochim_ids,
            file.path(output_path, "ASV_table.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

asv_long <- as.data.frame(seqtab_main_nochim_ids)
asv_long$SampleID <- rownames(asv_long)
asv_long <- tidyr::pivot_longer(asv_long, cols = starts_with("ASV_"),
                                names_to = "ASV_ID", values_to = "Reads")
asv_long <- asv_long[asv_long$Reads > 0, ]
write.table(asv_long, file.path(output_path, "ASV_table_long.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Tabela long da corrida auxiliar com IDs consistentes
if (nrow(seqtab_aux_nochim_ids) > 0) {
  asv_long_aux <- as.data.frame(seqtab_aux_nochim_ids)
  asv_long_aux$SampleID <- rownames(asv_long_aux)
  asv_long_aux <- tidyr::pivot_longer(
    asv_long_aux,
    cols = dplyr::starts_with("ASV_"),
    names_to = "ASV_ID", values_to = "Reads")
  asv_long_aux <- asv_long_aux[asv_long_aux$Reads > 0, ]
  write.table(asv_long_aux,
              file.path(output_path, "ASV_table_auxiliar_long.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("ASV_table_auxiliar_long.tsv: %d linhas (ASVs presentes na auxiliar)\n",
              nrow(asv_long_aux)))
}

# --------------------------------------------------------------------------
# seqtab_global_nochim: merge pos-filtro de main + aux
# Uso: entrada canonica dos Scripts 2a-2e e 4.
# Classifica TODAS as ASVs (main + exclusivas_aux) em uma unica passagem.
# O particionamento por origem e feito no Script 7 via ASV_sequences.tsv.
# --------------------------------------------------------------------------
# A junção global é obrigatória: 9 amostras principais + 1 auxiliar.
if (isTRUE(aux_vazia) || nrow(seqtab_aux_nochim) != 1L ||
    ncol(seqtab_aux_nochim) == 0L) {
  stop(
    "A tabela auxiliar final não satisfaz o contrato de 1 amostra com ASVs. ",
    "O seqtab global não será criado incompleto.",
    call. = FALSE
  )
}
seqtab_global_nochim <- mergeSequenceTables(
  seqtab_main_nochim,
  seqtab_aux_nochim
)
# A tabela global com IDs e derivada da tabela global por sequencia, evitando
# usar mergeSequenceTables() sobre rotulos que nao sao sequencias biologicas.
ids_globais <- unname(seq2id_global[colnames(seqtab_global_nochim)])
if (anyNA(ids_globais) || anyDuplicated(ids_globais) > 0L) {
  stop("Falha ao mapear as sequencias globais para ASV_ID unicos.", call. = FALSE)
}
seqtab_global_nochim_ids <- seqtab_global_nochim
colnames(seqtab_global_nochim_ids) <- ids_globais

amostras_esperadas <- metadata$SampleID
if (nrow(seqtab_global_nochim) != length(amostras_esperadas) ||
    !setequal(rownames(seqtab_global_nochim), amostras_esperadas)) {
  stop(
    "seqtab_global_nochim não preservou o universo amostral: esperadas=",
    length(amostras_esperadas), "; encontradas=", nrow(seqtab_global_nochim),
    call. = FALSE
  )
}
seqtab_global_nochim <- seqtab_global_nochim[
  amostras_esperadas, , drop = FALSE
]
seqtab_global_nochim_ids <- seqtab_global_nochim_ids[
  amostras_esperadas, , drop = FALSE
]
saveRDS(seqtab_global_nochim,     file.path(output_path, "seqtab_global_nochim.rds"))
saveRDS(seqtab_global_nochim_ids, file.path(output_path, "seqtab_global_nochim_ids.rds"))

cat(sprintf("seqtab_global_nochim.rds: %d amostras | %d ASVs (colunas = sequencias)\n",
            nrow(seqtab_global_nochim), ncol(seqtab_global_nochim)))
cat(sprintf("seqtab_global_nochim_ids.rds: %d amostras | %d ASVs (colunas = IDs)\n",
            nrow(seqtab_global_nochim_ids), ncol(seqtab_global_nochim_ids)))
# --------------------------------------------------------------------------
# Saidas primarias
# --------------------------------------------------------------------------
saveRDS(seqtab_main_nochim, file.path(output_path, "seqtab_nochim.rds"))
saveRDS(seqtab_aux_nochim,  file.path(output_path, "seqtab_auxiliar.rds"))

cat(sprintf("\nseqtab_nochim.rds   : %d amostras | %d ASVs (principal — auditoria e analise core9)\n",
            nrow(seqtab_main_nochim), ncol(seqtab_main_nochim)))
cat(sprintf("seqtab_auxiliar.rds : %d amostra  | %d ASVs (auxiliar  — analise suplementar)\n",
            nrow(seqtab_aux_nochim),  ncol(seqtab_aux_nochim)))
cat("AVISO: seqtab_auxiliar.rds usa modelo de erros PROPRIO da corrida auxiliar\n")
cat("       (amostra unica; ajuste inspecionado em modelo_erros_*_proprio.pdf).\n")
cat("       Para classificacao taxonomica completa da auxiliar (incluindo\n")
cat("       ASVs exclusivas), usar seqtab_global_nochim.rds nos Scripts 2a-2e e 4.\n")
cat("       Declarar ambas as limitacoes na secao de Material e Metodos.\n")

###############################################################################
# 19. ESTATISTICAS E RELATORIO FINAL
###############################################################################
write.csv(filter_stats, file.path(output_path, "filter_stats.csv"), row.names = FALSE)
write.csv(merge_stats,  file.path(output_path, "merge_stats.csv"),  row.names = FALSE)
write.csv(track_all,    file.path(output_path, "track_reads.csv"))
write.csv(track_pct,    file.path(output_path, "track_reads_percentual.csv"))
write.csv(as.data.frame(length_distribution),
          file.path(output_path, "length_distribution.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_path, "sessionInfo.txt"))

report_file <- file.path(output_path,
                         paste0("pipeline_report_", pipeline_version, ".txt"))
sink(report_file)
tryCatch({
  cat(
    "RELATORIO FINAL — PIPELINE DADA2 (multi-run) ",pipeline_version, "\n\n", sep = ""
  )
  cat("dada2:", as.character(packageVersion("dada2")), "\n\n")
  cat("Corridas processadas:\n"); print(table(metadata$Run))
  cat(sprintf("\nCorrida principal : %s (%d amostras) — learnErrors proprio\n",
              run_principal, run_n_amostras[run_principal]))
  cat(sprintf(  "Corrida auxiliar  : %s (%d amostra) — learnErrors proprio\n",
              run_auxiliar, run_n_amostras[run_auxiliar]))
  cat("\nPolitica de pooling por corrida:\n")
  print(pool_by_run)
  cat("Fallback de pooling:", pool_fallback, "\n")
  cat("Pool por corrida:\n")
  for (r in names(pool_used)) cat(" ", r, ":", pool_used[[r]], "\n")
  cat(sprintf("\nAmostras finais (principal): %d\n", nrow(seqtab_main_nochim)))
  cat(sprintf("ASVs finais     (principal): %d  [IDs ASV_1…ASV_%d]\n",
              ncol(seqtab_main_nochim), ncol(seqtab_main_nochim)))
  cat(sprintf("Amostras finais (auxiliar) : %d\n", nrow(seqtab_aux_nochim)))
  cat(sprintf("ASVs finais     (auxiliar) : %d  [compartilhadas: %d | exclusivas: %d]\n",
              ncol(seqtab_aux_nochim),
              length(asv_seqs_compartilhadas),
              length(asv_seqs_exclusivas_aux)))
  cat(sprintf("ASVs.fa (FASTA global)     : %d sequencias\n", length(all_ids_global)))
  cat("Reads brutas (pre-cutadapt):",
      format(sum(reads_raw_R1, na.rm = TRUE), big.mark = ","), "\n")
  cat("Retencao global:", global_retention, "%\n")
  cat("Retencao pos-quimera (combinado):", pct_nochim, "%\n")
  cat("Filtro comprimento:", expected_min_len, "-", expected_max_len, "pb\n")
  cat("Filtro principal: reads <", min_total_reads,
      "E prevalencia <", min_prevalence, "/ 9 amostras (logica E)\n")
  cat("Filtro auxiliar: reads <", min_total_reads, "apenas\n")
  cat(sprintf("Profundidade minima principal (rarefacao): %d reads — %s\n",
              min(depth_final_main), names(which.min(depth_final_main))))
  cat("Controle negativo: AUSENTE.\n")
  cat("\nLIMITACAO METODOLOGICA:\n")
  cat("  O modelo de erros da corrida auxiliar (", run_auxiliar, ") foi estimado\n")
  cat("  a partir de uma unica amostra. Embora a retencao de reads tenha sido\n")
  cat("  elevada, a menor diversidade de amostras pode reduzir a estabilidade\n")
  cat("  do ajuste em faixas de qualidade com poucas observacoes. A amostra foi\n")
  cat("  adicionalmente processada com o modelo da corrida principal como analise\n")
  cat("  de sensibilidade (ver sensibilidade_modelo_auxiliar_resumo.csv). A amostra\n")
  cat("  auxiliar nao foi incluida nos testes estatisticos principais entre grupos.\n")
  print(sessionInfo())
}, finally = {
  sink()
})

###############################################################################
# 20. FINALIZACAO
###############################################################################
cat("\n=== PIPELINE DADA2 (multi-run) FINALIZADO ===\n")
cat("Saidas em:", output_path, "\n")
cat("-- OBJETOS CANONICOS E ANALISE PRINCIPAL --\n")
cat(" - seqtab_nochim.rds            (principal, pos-filtro — colunas = sequencias)\n")
cat(" - ASVs.fa                      (principal + exclusivas_aux — FASTA para BLAST)\n")
cat(" - ASV_table.tsv / ASV_table_long.tsv (contagens por ASV; apenas corrida principal)\n")
cat(" - ASV_sequences.tsv            (todos os IDs com coluna Origem)\n")
cat("-- ANALISE SUPLEMENTAR (auxiliar) --\n")
cat(" - seqtab_auxiliar.rds          (auxiliar, pos-filtro; modelo de erros proprio; modelo principal apenas em sensibilidade)\n")
cat(" - seqtab_global_nochim.rds     (main + aux — colunas = sequencias; entrada Scripts 2-3)\n")
cat(" - seqtab_global_nochim_ids.rds (main + aux — colunas = IDs ASV/ASV_AUX)\n")
cat(" - ASVs_principal.fa            (FASTA so da corrida principal)\n")
cat(" - ASVs_auxiliar_exclusivas.fa  (FASTA das ASVs exclusivas da auxiliar)\n")
cat(" - ASV_table_auxiliar_long.tsv\n")
cat("-- UTILITARIOS --\n")
cat(" - seqtab_nochim_prefiltro.rds  (combinado, pre-filtro — auditoria/sensibilidade; nao e entrada do Script 9)\n")
cat(" - track_reads.csv              (raw->input->filtradas->denoised->merged->final)\n")
cat(" - sensitivity_thresholds.csv   (avaliado sobre corrida principal, n=9)\n")
cat(" - cutadapt_retention.csv / primer_residual_check.csv\n")
cat("Proximo passo: 01b_auditoria_quimeras_recuperacao_ASVs.R (validar limiares de filtro)\n")
})
