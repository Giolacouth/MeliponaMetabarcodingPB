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
run_pipeline_script("01b_auditoria_quimeras_recuperacao_ASVs.R", "auditoria_asvs", function(ctx) {
###############################################################################
# SCRIPT 01b — AUDITORIA E RECUPERACAO DE ASVs
#
# OBJETIVO
#   Auditar ASVs removidas após a remoção de quimeras e produzir cenários de
#   recuperação SEM sobrescrever os resultados primários do pipeline.
#
# O QUE ESTE SCRIPT CONSEGUE RECUPERAR
#   - ASVs presentes em seqtab_nochim_prefiltro.rds:
#       já passaram pelo filtro de comprimento;
#       foram consideradas NÃO QUIMÉRICAS pelo pipeline;
#       podem ter sido removidas apenas pelos filtros posteriores de
#       abundância/prevalência.
#
# O QUE ESTE SCRIPT NÃO CONSEGUE RECUPERAR SEM ARQUIVOS ADICIONAIS
#   - reads removidas por Cutadapt/filterAndTrim;
#   - pares rejeitados no merge;
#   - ASVs removidas pelo filtro de comprimento antes de serem salvas;
#   - ASVs classificadas como quimeras, caso seqtab_prechimera.rds não exista.
#
# SAÍDAS
#   - inventário completo das ASVs pós-quimera e pré-filtro;
#   - ASVs removidas pela regra atual, com sequência e contagens por amostra;
#   - cenários atual, relaxado e sem filtro de frequência;
#   - tabelas RDS de recuperação, sem modificar os objetos oficiais;
#   - FASTA das candidatas;
#   - auditoria opcional das quimeras, quando seqtab_prechimera.rds e
#     chimera_flags_baseline.rds estiverem disponiveis.
###############################################################################

options(
  encoding = "UTF-8",
  stringsAsFactors = FALSE,
  warn = 1
)

SEED_GLOBAL <- 1234L
set.seed(SEED_GLOBAL)

suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
})

###############################################################################
# 1. CONFIGURAÇÃO
###############################################################################

pipeline_version <- ctx$version
base_path <- ctx$base_path
raw_path <- ctx$raw_path
metadata_path <- file.path(raw_path, "metadados.tsv")
output_path <- ctx$output_root

audit_path <- ctx$stage$root
table_path <- ctx$stage$tables
fasta_path <- file.path(ctx$stage$root, "fasta")
plot_path <- ctx$stage$figures
rds_path <- ctx$stage$objects
log_path <- ctx$stage$logs
for (d in c(audit_path, table_path, fasta_path, plot_path, rds_path, log_path)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
}

# Entradas obrigatórias e opcionais
prefilter_file <- ctx$contracts[["seqtab_nochim_prefiltro"]]

main_final_file <- ctx$contracts[["seqtab_nochim"]]

aux_final_file <- ctx$contracts[["seqtab_auxiliar"]]

prelength_file <- ctx$contracts[["seqtab_pre_filtro_comprimento"]]

prechimera_file <- ctx$contracts[["seqtab_prechimera"]]

chimera_flags_file <- ctx$contracts[["chimera_flags_baseline"]]

###############################################################################
# 2. REGRAS ATUAIS E CENÁRIOS DE RECUPERAÇÃO
###############################################################################

# Regra aplicada no Script 01. Preferir o arquivo de parametros gravado pela
# mesma execucao para impedir divergencia silenciosa entre scripts.
param_file <- ctx$contracts[["pipeline_parametros"]]
if (file.exists(param_file)) {
  param_exec <- read.csv(param_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(param_exec) != 1L || !all(c(
    "Min_total_reads_principal", "Min_prevalence_principal",
    "Min_total_reads_auxiliar", "Cenario_filtro"
  ) %in% colnames(param_exec))) {
    stop("pipeline_parametros_execucao.csv invalido.", call. = FALSE)
  }
  current_min_total_reads_main <- as.integer(param_exec$Min_total_reads_principal[1L])
  current_min_prevalence_main <- as.integer(param_exec$Min_prevalence_principal[1L])
  current_min_total_reads_aux <- as.integer(param_exec$Min_total_reads_auxiliar[1L])
  cenario_filtro_executado <- as.character(param_exec$Cenario_filtro[1L])
} else {
  warning("pipeline_parametros_execucao.csv ausente; usando valores legados 10/2/10.")
  current_min_total_reads_main <- 10L
  current_min_prevalence_main <- 2L
  current_min_total_reads_aux <- 10L
  cenario_filtro_executado <- "desconhecido"
}

# Cenário relaxado, destinado à sensibilidade.
# Não substitui automaticamente a análise principal.
relaxed_min_total_reads_main <- 3L
relaxed_min_prevalence_main <- 2L
relaxed_min_total_reads_aux <- 3L

# Preservar no cenário relaxado uma sequência observada nas duas corridas,
# mesmo que cada corrida isoladamente tenha baixa contagem.
preserve_shared_between_runs <- TRUE

# Taxonomia opcional das ASVs descartadas.
# Deixe FALSE até confirmar o caminho do banco SILVA.
RUN_OPTIONAL_TAXONOMY <- FALSE
silva_train_set <- file.path(
  base_path,
  "bancodados",
  "Silva",
  "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)
taxonomy_min_boot <- 80L

# Auditoria opcional de quimeras.
# Só roda quando seqtab_prechimera.rds estiver disponível.
RUN_OPTIONAL_CHIMERA_SENSITIVITY <- TRUE

###############################################################################
# 3. FUNÇÕES
###############################################################################

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

pct <- function(n, d, digits = 4L) {
  if (length(d) == 0L || is.na(d) || d == 0L) {
    return(NA_real_)
  }
  round(100 * n / d, digits)
}

write_tsv <- function(x, filename) {
  destino <- file.path(table_path, filename)
  write.table(
    x,
    file = destino,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
  if (!file.exists(destino) ||
      is.na(file.size(destino)) ||
      file.size(destino) == 0L) {
    abort("Falha ao salvar TSV ou arquivo vazio: %s", destino)
  }
  invisible(destino)
}

is_dna_sequence <- function(x) {
  grepl("^[ACGT]+$", x)
}

validate_seqtab <- function(x, object_name) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    abort("%s não é uma matriz/data.frame.", object_name)
  }

  x <- as.matrix(x)

  if (nrow(x) == 0L) {
    abort("%s não possui amostras.", object_name)
  }

  if (ncol(x) == 0L) {
    abort("%s não possui ASVs.", object_name)
  }

  if (is.null(rownames(x)) || any(rownames(x) == "")) {
    abort("%s não possui SampleID válidos nas linhas.", object_name)
  }

  if (is.null(colnames(x)) || any(colnames(x) == "")) {
    abort("%s não possui sequências/IDs nas colunas.", object_name)
  }

  if (!is.numeric(x)) {
    abort("%s deve conter contagens numericas.", object_name)
  }

  if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    abort("%s possui contagens invalidas.", object_name)
  }

  if (max(abs(x - round(x))) > 1e-8) {
    abort("%s deve conter contagens inteiras.", object_name)
  }

  if (!all(is_dna_sequence(colnames(x)))) {
    abort(
      paste0(
        "%s não possui sequências de DNA nas colunas. ",
        "Este script requer o objeto com colunas = sequências."
      ),
      object_name
    )
  }

  # DADA2 (isBimeraDenovoTable) espera contagens inteiras; manter integer.
  # Ressalva: valores > 2^31 sofreriam overflow -> NA. Irrelevante para
  # contagens de ASV por amostra (muito abaixo desse limite), mas registrado.
  storage.mode(x) <- "integer"
  x
}

read_optional_seqtab <- function(path, object_name) {
  if (!file.exists(path)) {
    return(NULL)
  }

  obj <- tryCatch(
    readRDS(path),
    error = function(e) {
      abort(
        "Falha ao ler %s (%s): %s",
        object_name,
        path,
        conditionMessage(e)
      )
    }
  )

  validate_seqtab(obj, object_name)
}

matrix_by_samples <- function(seqtab, sample_ids) {
  present <- intersect(sample_ids, rownames(seqtab))

  if (length(present) == 0L) {
    return(
      matrix(
        numeric(),
        nrow = 0L,
        ncol = 0L
      )
    )
  }

  out <- seqtab[present, , drop = FALSE]
  out <- out[, colSums(out) > 0, drop = FALSE]
  out
}

vector_for_sequences <- function(values, sequences, default = 0) {
  out <- rep(default, length(sequences))
  names(out) <- sequences

  common <- intersect(names(values), sequences)
  out[common] <- values[common]
  out
}

logical_for_sequences <- function(sequence_names, sequences) {
  sequences %in% sequence_names
}

make_audit_ids <- function(sequences) {
  ordered <- sort(unique(sequences))
  ids <- sprintf("AUDIT_ASV_%06d", seq_along(ordered))
  names(ids) <- ordered
  ids[sequences]
}

write_fasta <- function(sequences, names_vector, filename) {
  if (length(sequences) == 0L) {
    writeLines(
      character(),
      file.path(fasta_path, filename)
    )
    return(invisible(NULL))
  }

  x <- DNAStringSet(sequences)
  names(x) <- names_vector

  writeXStringSet(
    x,
    filepath = file.path(fasta_path, filename),
    format = "fasta"
  )
}

safe_max <- function(x) {
  if (length(x) == 0L) 0 else max(x)
}

classification_priority <- function(
    total_global,
    prevalence_global,
    shared_between_runs) {

  if (shared_between_runs || prevalence_global >= 2L || total_global >= 5L) {
    return("PRIORIDADE_1_MAIOR_EVIDENCIA")
  }

  if (total_global >= 2L) {
    return("PRIORIDADE_2_BAIXA_EVIDENCIA")
  }

  "PRIORIDADE_3_READ_UNICA"
}

###############################################################################
# 4. CARREGAMENTO E VALIDAÇÃO
###############################################################################

if (!file.exists(metadata_path)) {
  abort("Metadado não encontrado: %s", metadata_path)
}

if (!file.exists(prefilter_file)) {
  abort(
    paste0(
      "Arquivo obrigatório não encontrado: %s\n",
      "Execute o DADA2 pelo menos até salvar seqtab_nochim_prefiltro.rds."
    ),
    prefilter_file
  )
}

for (arq in c(main_final_file, aux_final_file)) {
  if (!file.exists(arq)) {
    abort("Arquivo final obrigatorio ausente: %s. A auditoria nao inferira a membresia final.", arq)
  }
}

metadata <- read.delim(
  metadata_path,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_metadata <- c("SampleID", "Run")
missing_metadata_columns <- setdiff(
  required_metadata,
  names(metadata)
)

if (length(missing_metadata_columns) > 0L) {
  abort(
    "Colunas ausentes no metadado: %s",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

metadata$SampleID <- trimws(as.character(metadata$SampleID))
metadata$Run <- trimws(as.character(metadata$Run))

if (anyNA(metadata$SampleID) || any(metadata$SampleID == "")) {
  abort("Ha SampleID ausente ou vazio no metadado.")
}
if (anyDuplicated(metadata$SampleID)) {
  abort("Ha SampleID duplicado no metadado.")
}
if (anyNA(metadata$Run) || any(metadata$Run == "")) {
  abort("Ha valor Run ausente ou vazio no metadado.")
}

run_principal <- "run_main"
run_auxiliar <- "run_aux"
run_sizes <- table(metadata$Run)
if (!setequal(names(run_sizes), c(run_principal, run_auxiliar)) ||
    unname(run_sizes[run_principal]) != 9L ||
    unname(run_sizes[run_auxiliar]) != 1L) {
  abort(
    "Desenho de corridas divergente: esperado %s=9 e %s=1; observado: %s",
    run_principal, run_auxiliar,
    paste(names(run_sizes), as.integer(run_sizes), sep = "=", collapse = ";")
  )
}

main_samples <- metadata$SampleID[metadata$Run == run_principal]
aux_samples <- metadata$SampleID[metadata$Run == run_auxiliar]
if (!identical(aux_samples, "S10")) {
  abort("A corrida auxiliar deve conter apenas S10.")
}

seqtab_prefilter_obj <- tryCatch(
  readRDS(prefilter_file),
  error = function(e) {
    abort(
      "Falha ao ler seqtab_nochim_prefiltro.rds (%s): %s",
      prefilter_file,
      conditionMessage(e)
    )
  }
)
seqtab_prefilter <- validate_seqtab(
  seqtab_prefilter_obj,
  "seqtab_nochim_prefiltro.rds"
)
rm(seqtab_prefilter_obj)

missing_samples <- setdiff(metadata$SampleID, rownames(seqtab_prefilter))
extra_samples <- setdiff(rownames(seqtab_prefilter), metadata$SampleID)
if (length(missing_samples) > 0L || length(extra_samples) > 0L) {
  abort(
    "Universo amostral divergente no pre-filtro. Ausentes=%s; extras=%s",
    paste(missing_samples, collapse = ","),
    paste(extra_samples, collapse = ",")
  )
}

seqtab_main_pre <- matrix_by_samples(
  seqtab_prefilter,
  main_samples
)

seqtab_aux_pre <- matrix_by_samples(
  seqtab_prefilter,
  aux_samples
)

seqtab_main_final <- read_optional_seqtab(
  main_final_file,
  "seqtab_nochim.rds"
)

seqtab_aux_final <- read_optional_seqtab(
  aux_final_file,
  "seqtab_auxiliar.rds"
)

###############################################################################
# 5. RECONSTRUIR O EFEITO DOS FILTROS ATUAIS
###############################################################################

main_total <- colSums(seqtab_main_pre)
main_prevalence <- colSums(seqtab_main_pre > 0)
main_max_sample <- apply(seqtab_main_pre, 2L, max)

aux_total <- if (ncol(seqtab_aux_pre) > 0L) {
  colSums(seqtab_aux_pre)
} else {
  numeric()
}

aux_prevalence <- if (ncol(seqtab_aux_pre) > 0L) {
  colSums(seqtab_aux_pre > 0)
} else {
  numeric()
}

aux_max_sample <- if (ncol(seqtab_aux_pre) > 0L) {
  apply(seqtab_aux_pre, 2L, max)
} else {
  numeric()
}

all_sequences <- colnames(seqtab_prefilter)
audit_ids <- make_audit_ids(all_sequences)

main_total_all <- vector_for_sequences(
  main_total,
  all_sequences,
  default = 0
)

main_prev_all <- vector_for_sequences(
  main_prevalence,
  all_sequences,
  default = 0
)

main_max_all <- vector_for_sequences(
  main_max_sample,
  all_sequences,
  default = 0
)

aux_total_all <- vector_for_sequences(
  aux_total,
  all_sequences,
  default = 0
)

aux_prev_all <- vector_for_sequences(
  aux_prevalence,
  all_sequences,
  default = 0
)

aux_max_all <- vector_for_sequences(
  aux_max_sample,
  all_sequences,
  default = 0
)

global_total <- colSums(seqtab_prefilter)
global_prevalence <- colSums(seqtab_prefilter > 0)

present_main <- main_total_all > 0
present_aux <- aux_total_all > 0
shared_between_runs <- present_main & present_aux

# Regra atual reproduzida exatamente
discard_main_current_rule <- (
  present_main &
    main_total_all < current_min_total_reads_main &
    main_prev_all < current_min_prevalence_main
)

discard_aux_current_rule <- (
  present_aux &
    aux_total_all < current_min_total_reads_aux
)

# Conferência contra os arquivos finais, quando disponíveis
if (!is.null(seqtab_main_final)) {
  in_main_final <- logical_for_sequences(
    colnames(seqtab_main_final),
    all_sequences
  )
} else {
  abort("seqtab_nochim.rds ausente; membresia final nao sera inferida.")
}

if (!is.null(seqtab_aux_final)) {
  in_aux_final <- logical_for_sequences(
    colnames(seqtab_aux_final),
    all_sequences
  )
} else {
  abort("seqtab_auxiliar.rds ausente; membresia final nao sera inferida.")
}

discarded_main_observed <- present_main & !in_main_final
discarded_aux_observed <- present_aux & !in_aux_final
discarded_any_observed <- discarded_main_observed | discarded_aux_observed

rule_mismatch_main <- (
  present_main &
    (discarded_main_observed != discard_main_current_rule)
)

rule_mismatch_aux <- (
  present_aux &
    (discarded_aux_observed != discard_aux_current_rule)
)

###############################################################################
# 6. CENÁRIO RELAXADO
###############################################################################

keep_main_relaxed <- (
  present_main &
    (
      main_total_all >= relaxed_min_total_reads_main |
        main_prev_all >= relaxed_min_prevalence_main |
        (
          preserve_shared_between_runs &
            shared_between_runs
        )
    )
)

keep_aux_relaxed <- (
  present_aux &
    (
      aux_total_all >= relaxed_min_total_reads_aux |
        (
          preserve_shared_between_runs &
            shared_between_runs
        )
    )
)

recovered_main_relaxed <- (
  discarded_main_observed &
    keep_main_relaxed
)

recovered_aux_relaxed <- (
  discarded_aux_observed &
    keep_aux_relaxed
)

recovered_any_relaxed <- (
  recovered_main_relaxed |
    recovered_aux_relaxed
)

###############################################################################
# 7. INVENTÁRIO DETALHADO
###############################################################################

status_main <- rep("AUSENTE_NA_PRINCIPAL", length(all_sequences))
status_main[present_main & in_main_final] <- "RETIDA_FINAL_PRINCIPAL"
status_main[discarded_main_observed] <- "DESCARTADA_POS_QUIMERA_PRINCIPAL"

status_aux <- rep("AUSENTE_NA_AUXILIAR", length(all_sequences))
status_aux[present_aux & in_aux_final] <- "RETIDA_FINAL_AUXILIAR"
status_aux[discarded_aux_observed] <- "DESCARTADA_POS_QUIMERA_AUXILIAR"

priority <- vapply(
  seq_along(all_sequences),
  function(i) {
    classification_priority(
      total_global = global_total[i],
      prevalence_global = global_prevalence[i],
      shared_between_runs = shared_between_runs[i]
    )
  },
  character(1L)
)

inventory <- data.frame(
  AuditID = unname(audit_ids),
  Sequence = all_sequences,
  Length = nchar(all_sequences),
  NonChimericAccordingToPrimaryPipeline = TRUE,
  Main_TotalReads = as.numeric(main_total_all),
  Main_Prevalence = as.integer(main_prev_all),
  Main_MaxReadsInOneSample = as.numeric(main_max_all),
  Aux_TotalReads = as.numeric(aux_total_all),
  Aux_Prevalence = as.integer(aux_prev_all),
  Aux_MaxReadsInOneSample = as.numeric(aux_max_all),
  Global_TotalReads = as.numeric(global_total),
  Global_Prevalence = as.integer(global_prevalence),
  PresentMainBeforeFrequencyFilter = present_main,
  PresentAuxBeforeFrequencyFilter = present_aux,
  SharedBetweenRunsBeforeFrequencyFilter = shared_between_runs,
  InMainFinal = in_main_final,
  InAuxFinal = in_aux_final,
  StatusMain = status_main,
  StatusAux = status_aux,
  CurrentRuleWouldDiscardMain = discard_main_current_rule,
  CurrentRuleWouldDiscardAux = discard_aux_current_rule,
  ObservedDiscardedMain = discarded_main_observed,
  ObservedDiscardedAux = discarded_aux_observed,
  RuleMismatchMain = rule_mismatch_main,
  RuleMismatchAux = rule_mismatch_aux,
  KeepMainRelaxed = keep_main_relaxed,
  KeepAuxRelaxed = keep_aux_relaxed,
  RecoveredMainRelaxed = recovered_main_relaxed,
  RecoveredAuxRelaxed = recovered_aux_relaxed,
  RecoveryPriority = priority,
  stringsAsFactors = FALSE
)

inventory <- inventory[
  order(
    !discarded_any_observed,
    inventory$RecoveryPriority,
    -inventory$Global_TotalReads,
    inventory$AuditID
  ),
  ,
  drop = FALSE
]

discarded_inventory <- inventory[
  inventory$ObservedDiscardedMain |
    inventory$ObservedDiscardedAux,
  ,
  drop = FALSE
]

relaxed_inventory <- inventory[
  inventory$RecoveredMainRelaxed |
    inventory$RecoveredAuxRelaxed,
  ,
  drop = FALSE
]

###############################################################################
# 8. CONTAGENS POR AMOSTRA DAS ASVs DESCARTADAS
###############################################################################

discarded_sequences <- discarded_inventory$Sequence

long_counts <- data.frame()

if (length(discarded_sequences) > 0L) {
  discarded_matrix <- seqtab_prefilter[
    ,
    discarded_sequences,
    drop = FALSE
  ]

  long_counts <- as.data.frame(
    as.table(discarded_matrix),
    stringsAsFactors = FALSE
  )

  names(long_counts) <- c(
    "SampleID",
    "Sequence",
    "Reads"
  )

  long_counts <- long_counts[
    long_counts$Reads > 0,
    ,
    drop = FALSE
  ]

  id_map <- setNames(
    inventory$AuditID,
    inventory$Sequence
  )

  run_map <- setNames(
    metadata$Run,
    metadata$SampleID
  )

  long_counts$AuditID <- unname(
    id_map[long_counts$Sequence]
  )

  long_counts$Run <- unname(
    run_map[long_counts$SampleID]
  )

  long_counts <- long_counts[
    ,
    c(
      "AuditID",
      "Sequence",
      "SampleID",
      "Run",
      "Reads"
    )
  ]

  long_counts <- long_counts[
    order(
      long_counts$AuditID,
      -long_counts$Reads,
      long_counts$SampleID
    ),
    ,
    drop = FALSE
  ]
}

###############################################################################
# 9. CONSTRUIR MATRIZES DE SENSIBILIDADE
###############################################################################

main_relaxed_sequences <- all_sequences[
  keep_main_relaxed
]

aux_relaxed_sequences <- all_sequences[
  keep_aux_relaxed
]

global_relaxed_sequences <- union(
  main_relaxed_sequences,
  aux_relaxed_sequences
)

seqtab_main_relaxed <- seqtab_main_pre[
  ,
  intersect(
    main_relaxed_sequences,
    colnames(seqtab_main_pre)
  ),
  drop = FALSE
]

seqtab_aux_relaxed <- if (ncol(seqtab_aux_pre) > 0L) {
  seqtab_aux_pre[
    ,
    intersect(
      aux_relaxed_sequences,
      colnames(seqtab_aux_pre)
    ),
    drop = FALSE
  ]
} else {
  seqtab_aux_pre
}

seqtab_global_relaxed <- seqtab_prefilter[
  ,
  global_relaxed_sequences,
  drop = FALSE
]

# O objeto sem filtro de frequência é uma cópia explícita do prefiltro.
seqtab_all_nonchimeric <- seqtab_prefilter

saveRDS(
  seqtab_main_relaxed,
  file.path(
    rds_path,
    "seqtab_main_recuperacao_relaxada.rds"
  )
)

saveRDS(
  seqtab_aux_relaxed,
  file.path(
    rds_path,
    "seqtab_aux_recuperacao_relaxada.rds"
  )
)

saveRDS(
  seqtab_global_relaxed,
  file.path(
    rds_path,
    "seqtab_global_recuperacao_relaxada.rds"
  )
)

saveRDS(
  seqtab_all_nonchimeric,
  file.path(
    rds_path,
    "seqtab_todas_ASVs_nao_quimericas_sem_filtro_frequencia.rds"
  )
)

###############################################################################
# 10. FASTA
###############################################################################

if (nrow(discarded_inventory) > 0L) {
  fasta_names_discarded <- paste0(
    discarded_inventory$AuditID,
    "|global_reads=",
    discarded_inventory$Global_TotalReads,
    "|global_prev=",
    discarded_inventory$Global_Prevalence,
    "|priority=",
    discarded_inventory$RecoveryPriority
  )

  write_fasta(
    discarded_inventory$Sequence,
    fasta_names_discarded,
    "ASVs_descartadas_nao_quimericas.fa"
  )
} else {
  write_fasta(
    character(),
    character(),
    "ASVs_descartadas_nao_quimericas.fa"
  )
}

if (nrow(relaxed_inventory) > 0L) {
  fasta_names_relaxed <- paste0(
    relaxed_inventory$AuditID,
    "|global_reads=",
    relaxed_inventory$Global_TotalReads,
    "|global_prev=",
    relaxed_inventory$Global_Prevalence,
    "|priority=",
    relaxed_inventory$RecoveryPriority
  )

  write_fasta(
    relaxed_inventory$Sequence,
    fasta_names_relaxed,
    "ASVs_candidatas_recuperacao_relaxada.fa"
  )
} else {
  write_fasta(
    character(),
    character(),
    "ASVs_candidatas_recuperacao_relaxada.fa"
  )
}

###############################################################################
# 11. TAXONOMIA OPCIONAL DAS ASVs DESCARTADAS
###############################################################################

taxonomy_discarded <- data.frame()

if (RUN_OPTIONAL_TAXONOMY) {
  if (!requireNamespace("dada2", quietly = TRUE)) {
    abort("Pacote dada2 necessário para a taxonomia opcional.")
  }

  if (!file.exists(silva_train_set)) {
    abort(
      "Banco SILVA não encontrado: %s",
      silva_train_set
    )
  }

  if (nrow(discarded_inventory) > 0L) {
    tax <- dada2::assignTaxonomy(
      seqs = discarded_inventory$Sequence,
      refFasta = silva_train_set,
      minBoot = taxonomy_min_boot,
      multithread = TRUE,
      tryRC = TRUE,
      outputBootstraps = FALSE,
      verbose = TRUE
    )

    taxonomy_discarded <- data.frame(
      AuditID = discarded_inventory$AuditID,
      Sequence = discarded_inventory$Sequence,
      as.data.frame(
        tax,
        stringsAsFactors = FALSE
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

###############################################################################
# 12. AUDITORIA OPCIONAL DAS QUIMERAS
###############################################################################

chimera_audit <- data.frame()
chimera_uncertain_matrix <- NULL

if (
  RUN_OPTIONAL_CHIMERA_SENSITIVITY &&
    file.exists(prechimera_file)
) {
  if (!requireNamespace("dada2", quietly = TRUE)) {
    abort("Pacote dada2 necessário para auditar quimeras.")
  }

  seqtab_prechimera <- read_optional_seqtab(
    prechimera_file,
    "seqtab_prechimera.rds"
  )

  if (file.exists(chimera_flags_file)) {
    baseline_flags <- tryCatch(
      readRDS(chimera_flags_file),
      error = function(e) {
        abort(
          "Falha ao ler chimera_flags_baseline.rds: %s",
          conditionMessage(e)
        )
      }
    )
    if (!is.logical(baseline_flags) || anyNA(baseline_flags)) {
      abort("chimera_flags_baseline.rds deve ser vetor logico sem NA.")
    }

    if (is.null(names(baseline_flags))) {
      if (length(baseline_flags) != ncol(seqtab_prechimera)) {
        abort(
          "chimera_flags_baseline.rds possui comprimento incompatível."
        )
      }

      names(baseline_flags) <- colnames(seqtab_prechimera)
    }

    baseline_flags <- baseline_flags[
      colnames(seqtab_prechimera)
    ]
  } else {
    baseline_flags <- dada2::isBimeraDenovoTable(
      seqtab_prechimera,
      minSampleFraction = 0.9,
      ignoreNNegatives = 1,
      minFoldParentOverAbundance = 1.5,
      minParentAbundance = 2,
      allowOneOff = FALSE,
      multithread = FALSE,
      verbose = TRUE
    )

    names(baseline_flags) <- colnames(seqtab_prechimera)
  }

  # Análise de sensibilidade que exige evidência parental mais forte.
  stronger_parent_evidence_flags <- dada2::isBimeraDenovoTable(
    seqtab_prechimera,
    minSampleFraction = 0.9,
    ignoreNNegatives = 1,
    minFoldParentOverAbundance = 2,
    minParentAbundance = 4,
    allowOneOff = FALSE,
    multithread = FALSE,
    verbose = TRUE
  )

  names(stronger_parent_evidence_flags) <- colnames(
    seqtab_prechimera
  )

  chimera_sequences <- colnames(seqtab_prechimera)
  chimera_total <- colSums(seqtab_prechimera)
  chimera_prevalence <- colSums(seqtab_prechimera > 0)

  chimera_status <- ifelse(
    baseline_flags & stronger_parent_evidence_flags,
    "QUIMERA_CONCORDANTE_MAIOR_EVIDENCIA",
    ifelse(
      baseline_flags & !stronger_parent_evidence_flags,
      "QUIMERA_BASELINE_APENAS_REVISAR",
      ifelse(
        !baseline_flags & stronger_parent_evidence_flags,
        "SENSIBILIDADE_APENAS_REVISAR",
        "NAO_QUIMERA_EM_AMBOS_CENARIOS"
      )
    )
  )

  chimera_audit <- data.frame(
    Sequence = chimera_sequences,
    Length = nchar(chimera_sequences),
    TotalReads = as.numeric(chimera_total),
    Prevalence = as.integer(chimera_prevalence),
    BaselineChimera = as.logical(baseline_flags),
    StrongerParentEvidenceChimera = as.logical(
      stronger_parent_evidence_flags
    ),
    Status = chimera_status,
    stringsAsFactors = FALSE
  )

  chimera_audit <- chimera_audit[
    order(
      chimera_audit$Status,
      -chimera_audit$TotalReads
    ),
    ,
    drop = FALSE
  ]

  uncertain_sequences <- chimera_audit$Sequence[
    chimera_audit$Status ==
      "QUIMERA_BASELINE_APENAS_REVISAR"
  ]

  if (length(uncertain_sequences) > 0L) {
    chimera_uncertain_matrix <- seqtab_prechimera[
      ,
      uncertain_sequences,
      drop = FALSE
    ]

    saveRDS(
      chimera_uncertain_matrix,
      file.path(
        rds_path,
        "seqtab_quimeras_baseline_incertas_para_revisao.rds"
      )
    )

    uncertain_names <- paste0(
      "CHIMERA_REVIEW_",
      sprintf("%06d", seq_along(uncertain_sequences))
    )

    write_fasta(
      uncertain_sequences,
      uncertain_names,
      "ASVs_quimera_baseline_incertas_para_revisao.fa"
    )
  }
}

###############################################################################
# 13. AUDITORIA DO FILTRO DE COMPRIMENTO, QUANDO DISPONÍVEL
###############################################################################

length_filter_audit <- data.frame()

if (file.exists(prelength_file)) {
  seqtab_prelength <- read_optional_seqtab(
    prelength_file,
    "seqtab_pre_filtro_comprimento.rds"
  )

  before_sequences <- colnames(seqtab_prelength)
  after_sequences <- if (file.exists(prechimera_file)) {
    colnames(
      read_optional_seqtab(
        prechimera_file,
        "seqtab_prechimera.rds"
      )
    )
  } else {
    character()
  }

  length_filter_audit <- data.frame(
    Sequence = before_sequences,
    Length = nchar(before_sequences),
    TotalReads = as.numeric(colSums(seqtab_prelength)),
    Prevalence = as.integer(colSums(seqtab_prelength > 0)),
    RetainedByLengthFilter = before_sequences %in% after_sequences,
    stringsAsFactors = FALSE
  )

  length_filter_audit <- length_filter_audit[
    order(
      length_filter_audit$RetainedByLengthFilter,
      length_filter_audit$Length,
      -length_filter_audit$TotalReads
    ),
    ,
    drop = FALSE
  ]
}

###############################################################################
# 14. RESUMO DOS CENÁRIOS
###############################################################################

main_current_asvs <- sum(present_main & in_main_final)
aux_current_asvs <- sum(present_aux & in_aux_final)

scenario_summary <- data.frame(
  Scenario = c(
    "ATUAL_PRINCIPAL",
    "RELAXADO_PRINCIPAL",
    "SEM_FILTRO_FREQUENCIA_PRINCIPAL",
    "ATUAL_AUXILIAR",
    "RELAXADO_AUXILIAR",
    "SEM_FILTRO_FREQUENCIA_AUXILIAR",
    "GLOBAL_RELAXADO",
    "GLOBAL_TODAS_NAO_QUIMERICAS"
  ),
  ASVs = c(
    main_current_asvs,
    ncol(seqtab_main_relaxed),
    ncol(seqtab_main_pre),
    aux_current_asvs,
    ncol(seqtab_aux_relaxed),
    ncol(seqtab_aux_pre),
    ncol(seqtab_global_relaxed),
    ncol(seqtab_prefilter)
  ),
  Reads = c(
    if (!is.null(seqtab_main_final)) {
      sum(seqtab_main_final)
    } else {
      sum(
        seqtab_main_pre[
          ,
          !discard_main_current_rule[
            colnames(seqtab_main_pre)
          ],
          drop = FALSE
        ]
      )
    },
    sum(seqtab_main_relaxed),
    sum(seqtab_main_pre),
    if (!is.null(seqtab_aux_final)) {
      sum(seqtab_aux_final)
    } else {
      sum(
        seqtab_aux_pre[
          ,
          !discard_aux_current_rule[
            colnames(seqtab_aux_pre)
          ],
          drop = FALSE
        ]
      )
    },
    sum(seqtab_aux_relaxed),
    sum(seqtab_aux_pre),
    sum(seqtab_global_relaxed),
    sum(seqtab_prefilter)
  ),
  stringsAsFactors = FALSE
)

scenario_summary$PercentReadsRelativeToNonChimericPrefilter <- pct(
  scenario_summary$Reads,
  sum(seqtab_prefilter)
)

###############################################################################
# 15. GRÁFICOS
###############################################################################

plot_data <- inventory[
  inventory$PresentMainBeforeFrequencyFilter,
  ,
  drop = FALSE
]

plot_data$FilterStatus <- ifelse(
  plot_data$ObservedDiscardedMain,
  "Descartada pela regra atual",
  "Retida na tabela principal"
)

p_main <- ggplot(
  plot_data,
  aes(
    x = Main_TotalReads,
    y = Main_Prevalence,
    shape = FilterStatus
  )
) +
  geom_point(
    alpha = 0.75,
    size = 2.4
  ) +
  geom_vline(
    xintercept = current_min_total_reads_main,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = current_min_prevalence_main,
    linetype = "dashed"
  ) +
  scale_x_log10() +
  scale_y_continuous(
    breaks = seq_len(length(main_samples))
  ) +
  labs(
    title = "ASVs pós-quimera antes e depois do filtro atual",
    subtitle = paste0(
      "A regra remove quando reads < ",
      current_min_total_reads_main,
      " E prevalência < ",
      current_min_prevalence_main
    ),
    x = "Reads totais na corrida principal (log10)",
    y = "Prevalência na corrida principal",
    shape = "Situação"
  ) +
  theme_bw()

ggsave(
  file.path(
    plot_path,
    "01_ASVs_pos_quimera_filtro_atual.png"
  ),
  p_main,
  width = 9,
  height = 6,
  dpi = 300
)

p_scenarios <- ggplot(
  scenario_summary,
  aes(
    x = Scenario,
    y = ASVs
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Número de ASVs por cenário de recuperação",
    x = NULL,
    y = "ASVs"
  ) +
  theme_bw()

ggsave(
  file.path(
    plot_path,
    "02_ASVs_por_cenario.png"
  ),
  p_scenarios,
  width = 9,
  height = 6,
  dpi = 300
)

###############################################################################
# 16. SALVAR TABELAS
###############################################################################

write_tsv(
  inventory,
  "01_inventario_ASVs_pos_quimera_prefiltro.tsv"
)

write_tsv(
  discarded_inventory,
  "02_ASVs_descartadas_nao_quimericas.tsv"
)

write_tsv(
  long_counts,
  "03_contagens_por_amostra_ASVs_descartadas.tsv"
)

write_tsv(
  relaxed_inventory,
  "04_ASVs_recuperadas_cenario_relaxado.tsv"
)

write_tsv(
  scenario_summary,
  "05_resumo_cenarios_recuperacao.tsv"
)

write_tsv(
  taxonomy_discarded,
  "06_taxonomia_opcional_ASVs_descartadas.tsv"
)

write_tsv(
  chimera_audit,
  "07_auditoria_sensibilidade_quimeras.tsv"
)

write_tsv(
  length_filter_audit,
  "08_auditoria_filtro_comprimento.tsv"
)

###############################################################################
# 17. RELATÓRIO DE CHECKPOINTS
###############################################################################

checkpoint <- data.frame(
  Check = c(
    "seqtab_nochim_prefiltro disponível",
    "Tabela final principal disponível",
    "Tabela final auxiliar disponível",
    "Regra reconstruída coincide com tabela principal",
    "Regra reconstruída coincide com tabela auxiliar",
    "Tabela pré-filtro de comprimento disponível",
    "Tabela pré-quimera disponível",
    "Flags originais de quimera disponíveis"
  ),
  Status = c(
    file.exists(prefilter_file),
    file.exists(main_final_file),
    file.exists(aux_final_file),
    !any(rule_mismatch_main),
    !any(rule_mismatch_aux),
    file.exists(prelength_file),
    file.exists(prechimera_file),
    file.exists(chimera_flags_file)
  ),
  Details = c(
    prefilter_file,
    main_final_file,
    aux_final_file,
    paste(sum(rule_mismatch_main), "ASVs divergentes"),
    paste(sum(rule_mismatch_aux), "ASVs divergentes"),
    prelength_file,
    prechimera_file,
    chimera_flags_file
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  checkpoint,
  "09_checkpoints_recuperacao.tsv"
)

###############################################################################
# 18. RESUMO EM TEXTO
###############################################################################

summary_lines <- c(
  "============================================================",
  "AUDITORIA DE ASVs NÃO QUIMÉRICAS",
  "============================================================",
  "",
  paste0("Corrida principal: ", run_principal),
  paste0("Amostras principais: ", length(main_samples)),
  paste0("Corrida auxiliar: ", run_auxiliar),
  paste0("Amostras auxiliares: ", length(aux_samples)),
  "",
  paste0(
    "ASVs pós-quimera e pré-filtro: ",
    ncol(seqtab_prefilter)
  ),
  paste0(
    "ASVs presentes na principal antes do filtro: ",
    ncol(seqtab_main_pre)
  ),
  paste0(
    "ASVs presentes na auxiliar antes do filtro: ",
    ncol(seqtab_aux_pre)
  ),
  "",
  paste0(
    "ASVs descartadas observadas na principal: ",
    sum(discarded_main_observed)
  ),
  paste0(
    "ASVs descartadas observadas na auxiliar: ",
    sum(discarded_aux_observed)
  ),
  paste0(
    "ASVs candidatas recuperadas no cenário relaxado: ",
    sum(recovered_any_relaxed)
  ),
  "",
  "INTERPRETAÇÃO:",
  paste0(
    "- As ASVs de seqtab_nochim_prefiltro.rds já foram consideradas ",
    "não quiméricas pelo pipeline."
  ),
  paste0(
    "- Recuperá-las significa retirar ou relaxar apenas o filtro posterior ",
    "de frequência/prevalência."
  ),
  paste0(
    "- Isso não prova que sejam biologicamente verdadeiras; variantes com ",
    "uma ou poucas reads exigem taxonomia, controles e análise de sensibilidade."
  ),
  paste0(
    "- Os objetos de recuperação foram salvos separadamente e não ",
    "sobrescrevem seqtab_nochim.rds nem seqtab_auxiliar.rds."
  ),
  "",
  "LIMITAÇÃO:",
  if (!file.exists(prechimera_file)) {
    paste0(
      "- Não é possível reavaliar ASVs removidas como quimeras sem ",
      "seqtab_prechimera.rds. Aplique o patch ao DADA2 e reexecute."
    )
  } else {
    paste0(
      "- A auditoria de quimeras foi executada. ASVs discordantes são ",
      "candidatas a revisão, não recuperação automática."
    )
  }
)

writeLines(
  summary_lines,
  con = file.path(
    audit_path,
    "resumo_auditoria_ASVs.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_path,
    "sessionInfo.txt"
  )
)

###############################################################################
# 19. CONSOLE
###############################################################################

cat("\n============================================================\n")
cat("AUDITORIA CONCLUÍDA\n")
cat("============================================================\n")
cat("Saída:", audit_path, "\n")
cat(
  "ASVs pós-quimera/pré-filtro:",
  ncol(seqtab_prefilter),
  "\n"
)
cat(
  "Descartadas na principal:",
  sum(discarded_main_observed),
  "\n"
)
cat(
  "Descartadas na auxiliar:",
  sum(discarded_aux_observed),
  "\n"
)
cat(
  "Recuperadas no cenário relaxado:",
  sum(recovered_any_relaxed),
  "\n"
)
cat(
  "Tabela principal para inspeção:",
  file.path(
    table_path,
    "02_ASVs_descartadas_nao_quimericas.tsv"
  ),
  "\n"
)
cat(
  "Objeto sem filtro de frequência:",
  file.path(
    rds_path,
    "seqtab_todas_ASVs_nao_quimericas_sem_filtro_frequencia.rds"
  ),
  "\n"
)
cat("\nNenhum resultado oficial foi sobrescrito.\n")
})
