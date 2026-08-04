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
run_pipeline_script("02a_silva138.R", "silva", function(ctx) {
###############################################################################
#   SCRIPT 02a — CLASSIFICACAO TAXONOMICA: SILVA 138.2
#
#   CARACTERISTICAS:
#   Curadoria manual guiada por arvore filogenetica (Yilmaz et al., 2014)
#   Nomenclatura ICNP classica (Firmicutes, nao Bacillota)
#   Atribuicao probabilistica via algoritmo Naive Bayes (Callahan et al., 2016)
#
#   ENTRADA: seqtab_global_nochim.rds (Matriz filtrada, livre de quimeras,
#            10 amostras — corrida principal + auxiliar)
#
################################################################################

library(dada2)

VERSAO        <- "1"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
BANCO_NOME    <- "silva138"

cat("=============================================================\n")
cat("SILVA 138.2: Classificacao Taxonomica v", VERSAO, "\n", sep = "")
cat("Data/Hora:", DATA_EXECUCAO, "\n")
cat("=============================================================\n\n")

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$stage$root
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_path)) stop("Falha ao criar diretorio: ", output_path, call. = FALSE)

# ---- ENTRADA CANONICA: identica para todos os bancos kmer ----
arq_seqtab   <- ctx$contracts[["seqtab_global_nochim"]]
arq_silva    <- file.path(base_path,
                          "bancodados/Silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
arq_silva_sp <- file.path(base_path,
                          "bancodados/Silva/silva_v138.2_assignSpecies.fa.gz")

arq_asvmap   <- ctx$contracts[["asv_sequences"]]

MINBOOT         <- 80
MULTITHREAD     <- TRUE
CHECKPOINT      <- TRUE

# ---------------------------------------------------------------------------
# FUNCOES AUXILIARES
# ---------------------------------------------------------------------------
log_msg <- function(msg, tipo = "INFO") {
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))
}

validar_arquivo <- function(arq, desc, obrig = TRUE) {
  if (!file.exists(arq)) {
    log_msg(paste(desc, "nao encontrado:", arq), "ERRO")
    if (obrig) stop(paste(desc, "nao encontrado:", arq), call. = FALSE)
    return(FALSE)
  }
  tamanho <- file.size(arq)
  if (is.na(tamanho) || tamanho == 0L) {
    log_msg(paste(desc, "esta vazio:", arq), "ERRO")
    if (obrig) stop(paste(desc, "esta vazio:", arq), call. = FALSE)
    return(FALSE)
  }
  log_msg(paste0(desc, " OK (", round(file.size(arq)/1024/1024, 2), " MB)"),
          "OK")
  TRUE
}

limpeza_mem <- function() { gc(verbose = FALSE); invisible(NULL) }

chkpt <- function(obj, sufixo) {
  if (!CHECKPOINT) return(invisible(NULL))
  arq <- file.path(output_path,
                   paste0("taxa_", BANCO_NOME, "_", sufixo, "_checkpoint.rds"))
  saveRDS(obj, arq)
  log_msg(paste("Checkpoint:", basename(arq)), "CHKPT")
}

com_asv <- function(mat) {
  cbind(ASV = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE))
}


carregar_seq2id <- function(arq) {
  if (!file.exists(arq))
    stop("ASV_sequences.tsv ausente: ", arq,
         "\nExecute o Script 01 antes deste script.")
  m <- tryCatch(
    read.delim(
      arq,
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
  req <- c("ASV_ID", "Sequence", "Origem")
  if (!all(req %in% colnames(m))) {
    stop(
      "ASV_sequences.tsv deve conter: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }
  m$ASV_ID <- trimws(as.character(m$ASV_ID))
  m$Sequence <- toupper(trimws(as.character(m$Sequence)))
  m$Origem <- trimws(as.character(m$Origem))
  if (anyNA(m$ASV_ID) || any(m$ASV_ID == "")) {
    stop("ASV_sequences.tsv contem ASV_ID ausente ou vazio.", call. = FALSE)
  }
  if (anyNA(m$Sequence) || any(m$Sequence == "")) {
    stop("ASV_sequences.tsv contem Sequence ausente ou vazia.", call. = FALSE)
  }
  if (anyNA(m$Origem) || any(m$Origem == "")) {
    stop("ASV_sequences.tsv contem Origem ausente ou vazia.", call. = FALSE)
  }
  if (anyDuplicated(m$ASV_ID) > 0) {
    stop("ASV_sequences.tsv contem ASV_ID duplicados.", call. = FALSE)
  }
  if (anyDuplicated(m$Sequence) > 0) {
    stop("ASV_sequences.tsv contem sequencias duplicadas.", call. = FALSE)
  }
  if (!all(grepl("^[ACGTRYSWKMBDHVN]+$", m$Sequence))) {
    stop("ASV_sequences.tsv contem sequencia fora do alfabeto IUPAC.", call. = FALSE)
  }
  setNames(m$ASV_ID, m$Sequence)
}


com_asv_id <- function(mat, seq2id) {
  seqs <- rownames(mat)
  ids  <- unname(seq2id[seqs])
  if (any(is.na(ids)))
    log_msg(sprintf("%d sequencia(s) sem ASV_ID correspondente em ASV_sequences.tsv.",
                    sum(is.na(ids))), "WARN")
  cbind(ASV_ID = ids, ASV_seq = seqs,
        as.data.frame(mat, stringsAsFactors = FALSE))
}

normalizar_colunas <- function(mat) {
  cn <- paste0(toupper(substr(colnames(mat), 1, 1)),
               tolower(substr(colnames(mat), 2, 99)))
  cn[cn == "Domain"] <- "Kingdom"
  manter <- cn != "Root"
  mat <- mat[, manter, drop = FALSE]
  colnames(mat) <- cn[manter]
  mat
}

# ---------------------------------------------------------------------------
# PARAMETROS E VALIDACOES
# ---------------------------------------------------------------------------
cat("Parametros:\n")
cat(sprintf("  MINBOOT=%d  MULTITHREAD=%s\n\n", MINBOOT, MULTITHREAD))

cat("=== VALIDACOES ===\n\n")
validar_arquivo(arq_seqtab,   "seqtab_global_nochim.rds", TRUE)
validar_arquivo(arq_silva,    "SILVA trainset",            TRUE)
validar_arquivo(arq_asvmap,   "ASV_sequences.tsv",         TRUE)
tem_sp_db <- validar_arquivo(arq_silva_sp, "SILVA assignSpecies", FALSE)
cat("\n")

seq2id <- carregar_seq2id(arq_asvmap)
log_msg(sprintf("Mapa ASV_ID <-> sequencia carregado: %d entradas.", length(seq2id)), "OK")

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
if (nrow(seqtab) != 10L) {
  stop(
    "Esperadas 10 amostras no seqtab; encontradas: ",
    nrow(seqtab),
    call. = FALSE
  )
}
if (ncol(seqtab) == 0L) {
  stop("seqtab nao possui ASVs.", call. = FALSE)
}
if (!is.numeric(seqtab)) {
  stop("seqtab deve conter contagens numericas.", call. = FALSE)
}
if (is.null(rownames(seqtab)) ||
    anyNA(rownames(seqtab)) ||
    any(rownames(seqtab) == "") ||
    anyDuplicated(rownames(seqtab)) > 0L) {
  stop("seqtab deve possuir SampleID unicos e nao vazios nas linhas.", call. = FALSE)
}
if (is.null(colnames(seqtab)) ||
    anyNA(colnames(seqtab)) ||
    any(colnames(seqtab) == "") ||
    anyDuplicated(colnames(seqtab)) > 0L) {
  stop("seqtab deve possuir sequencias unicas e nao vazias nas colunas.", call. = FALSE)
}
if (anyNA(seqtab) ||
    any(!is.finite(seqtab)) ||
    any(seqtab < 0) ||
    any(rowSums(seqtab) == 0) ||
    any(colSums(seqtab) == 0)) {
  stop(
    "seqtab contem valor invalido, amostra sem reads ou ASV com soma zero.",
    call. = FALSE
  )
}
if (max(abs(seqtab - round(seqtab))) > 1e-8) {
  stop("seqtab deve conter contagens inteiras.", call. = FALSE)
}
if (!all(grepl("^[ACGTRYSWKMBDHVN]+$", colnames(seqtab), ignore.case = TRUE))) {
  stop("seqtab contem sequencia fora do alfabeto IUPAC.", call. = FALSE)
}
if (!setequal(toupper(colnames(seqtab)), names(seq2id))) {
  stop("Universo do seqtab diverge de ASV_sequences.tsv.", call. = FALSE)
}
colnames(seqtab) <- toupper(colnames(seqtab))
cat(sprintf("ASVs: %d  |  Amostras: %d\n\n", ncol(seqtab), nrow(seqtab)))

# ==============================================================================
# PASSO 1: assignTaxonomy (kmer Naive Bayes)
# ==============================================================================
cat("=== PASSO 1: assignTaxonomy (kmer Naive Bayes) ===\n\n")
log_msg("Iniciando assignTaxonomy...", "INICIO")
t0 <- proc.time()

taxa_boot_obj <- tryCatch(
  assignTaxonomy(
    seqtab, arq_silva,
    multithread      = MULTITHREAD,
    minBoot          = MINBOOT,
    tryRC            = TRUE,
    outputBootstraps = TRUE),
  error = function(e) stop("SILVA assignTaxonomy falhou: ", conditionMessage(e), call. = FALSE)
)
if (!is.list(taxa_boot_obj) || !all(c("tax", "boot") %in% names(taxa_boot_obj))) {
  stop("SILVA assignTaxonomy nao retornou componentes tax/boot.")
}

taxa_raw <- taxa_boot_obj$tax
boots    <- taxa_boot_obj$boot
tempo_assign <- (proc.time() - t0)["elapsed"]
log_msg(sprintf("assignTaxonomy em %.1f s (%.1f min)",
                tempo_assign, tempo_assign/60), "OK")

taxa_raw <- normalizar_colunas(taxa_raw)
boots    <- normalizar_colunas(boots)
if (nrow(taxa_raw) != ncol(seqtab) ||
    nrow(boots) != ncol(seqtab) ||
    !identical(rownames(taxa_raw), colnames(seqtab)) ||
    !identical(rownames(boots), colnames(seqtab)) ||
    !identical(dim(taxa_raw), dim(boots)) ||
    !identical(colnames(taxa_raw), colnames(boots))) {
  stop(
    "SILVA retornou taxonomia/bootstrap desalinhados do universo canonico.",
    call. = FALSE
  )
}
log_msg(paste("Colunas:", paste(colnames(taxa_raw), collapse = ", ")), "INFO")
chkpt(taxa_raw, "apos_assigntaxonomy")
chkpt(boots,    "boots_assigntaxonomy")
limpeza_mem()

cat("\nCobertura inicial (%):\n")
cob_antes <- apply(taxa_raw, 2, function(x)
  round(100 * sum(!is.na(x)) / nrow(taxa_raw), 1))
print(cob_antes)

cat("\nBootstrap medio por nivel:\n")
boot_medio <- apply(boots, 2, function(x)
  round(mean(x[!is.na(x)], na.rm = TRUE), 1))
print(boot_medio)

# ==============================================================================
# PASSO 2: addSpecies (correspondencia exata)
# ==============================================================================
species_pre <- if ("Species" %in% colnames(taxa_raw)) {
  taxa_raw[, "Species"]
} else {
  rep(NA_character_, nrow(taxa_raw))
}

addspecies_executado <- FALSE
addspecies_sucesso <- FALSE
addspecies_erro <- NA_character_
species_exact <- rep(NA_character_, nrow(taxa_raw))

if (tem_sp_db) {
  cat("\n=== PASSO 2: addSpecies (correspondencia exata) ===\n\n")
  taxa_base <- taxa_raw[, setdiff(colnames(taxa_raw), "Species"), drop = FALSE]
  addspecies_executado <- TRUE
  taxa <- tryCatch({
    taxa_sp <- addSpecies(
      taxa_base,
      arq_silva_sp,
      tryRC = TRUE,
      allowMultiple = FALSE
    )
    addspecies_sucesso <- TRUE
    log_msg("addSpecies concluido", "OK")
    chkpt(taxa_sp, "apos_addspecies")
    taxa_sp
  }, error = function(e) {
    addspecies_erro <<- conditionMessage(e)
    log_msg(
      paste(
        "addSpecies falhou; preservando o resultado de assignTaxonomy:",
        addspecies_erro
      ),
      "WARN"
    )
    taxa_raw
  })

  species_post <- if ("Species" %in% colnames(taxa)) {
    taxa[, "Species"]
  } else {
    rep(NA_character_, nrow(taxa))
  }

  if (addspecies_sucesso) {
    species_exact <- species_post
    species_fonte <- ifelse(
      !is.na(species_exact),
      "exact_match",
      ifelse(!is.na(species_pre), "kmer_sem_match_exato", NA_character_)
    )
  } else {
    species_fonte <- ifelse(
      !is.na(species_post),
      "kmer_assignTaxonomy_fallback",
      NA_character_
    )
  }
  limpeza_mem()
} else {
  taxa <- taxa_raw
  species_post <- species_pre
  species_fonte <- ifelse(
    !is.na(species_pre),
    "kmer_assignTaxonomy",
    NA_character_
  )
  log_msg("addSpecies nao disponivel. Usando apenas assignTaxonomy", "WARN")
}

taxa_limpa  <- taxa
boots_limpa <- boots
cob_depois  <- apply(taxa_limpa, 2, function(x)
  round(100 * sum(!is.na(x)) / nrow(taxa_limpa), 1))

# ==============================================================================
# SALVAMENTO
# ==============================================================================
cat("\n=== SALVANDO ===\n\n")
salvar <- function(obj, nome) {
  arq_r <- file.path(output_path, paste0(nome, ".rds"))
  arq_c <- file.path(output_path, paste0(nome, ".csv"))
  saveRDS(obj, arq_r)
  if (!file.exists(arq_r) || is.na(file.size(arq_r)) || file.size(arq_r) == 0L)
    stop("Falha ao salvar RDS ou arquivo vazio: ", arq_r, call. = FALSE)
  # CSV traz ASV_ID numerico (primeira coluna) + ASV_seq, em vez de usar a
  # sequencia de DNA como unico identificador (ponto de rastreabilidade
  # manual entre tabelas de bancos diferentes).
  write.csv(com_asv_id(obj, seq2id), arq_c, row.names = FALSE, quote = TRUE)
  if (!file.exists(arq_c) || is.na(file.size(arq_c)) || file.size(arq_c) == 0L)
    stop("Falha ao salvar CSV ou arquivo vazio: ", arq_c, call. = FALSE)
  log_msg(paste("OK:", basename(arq_r), "|", basename(arq_c)), "SAVE")
}

salvar(taxa,        paste0("taxa_",  BANCO_NOME))
salvar(taxa_limpa,  paste0("taxa_",  BANCO_NOME, "_limpa"))
salvar(boots,       paste0("boots_", BANCO_NOME))
salvar(boots_limpa, paste0("boots_", BANCO_NOME, "_limpa"))

species_df <- data.frame(
  ASV_ID        = unname(seq2id[rownames(taxa_limpa)]),
  ASV           = rownames(taxa_limpa),
  Genus         = if ("Genus" %in% colnames(taxa_limpa)) taxa_limpa[, "Genus"] else NA,
  Species_kmer  = species_pre,
  Species_exact = species_exact,
  Species       = if ("Species" %in% colnames(taxa_limpa)) taxa_limpa[, "Species"] else NA,
  Species_fonte = species_fonte,
  Boot_Genus    = if ("Genus" %in% colnames(boots_limpa))
    boots_limpa[, "Genus"] else NA,
  stringsAsFactors = FALSE)
# NOTA DE COMPATIBILIDADE: a coluna "ASV" (sequencia) e mantida porque o
# Script 5 le este CSV e a referencia pelo nome literal "ASV". A coluna
# ASV_ID e adicionada para leitura/cruzamento manual sem quebrar o contrato.
write.csv(species_df,
          file.path(output_path, paste0("species_fonte_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)
saveRDS(species_df,
        file.path(output_path, paste0("species_fonte_", BANCO_NOME, ".rds")))

cob_df <- data.frame(
  Nivel      = names(cob_antes),
  Antes_pct  = cob_antes,
  Depois_pct = cob_depois[names(cob_antes)],
  Boot_medio = boot_medio[names(cob_antes)])
write.csv(cob_df,
          file.path(output_path, paste0("cobertura_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)

run_metadata <- data.frame(
  Script         = paste0("2_", BANCO_NOME, "_V1"),
  Versao         = VERSAO,
  Data_execucao  = DATA_EXECUCAO,
  Banco          = "SILVA 138.2",
  Seqtab_input   = basename(arq_seqtab),
  Metodo_assign  = if (addspecies_sucesso) {
    "assignTaxonomy_kmer + addSpecies_exact_univoco"
  } else if (addspecies_executado) {
    "assignTaxonomy_kmer_only; addSpecies_falhou"
  } else {
    "assignTaxonomy_kmer_only"
  },
  MINBOOT        = MINBOOT,
  MULTITHREAD    = MULTITHREAD,
  ASVs_entrada   = ncol(seqtab),
  Amostras       = nrow(seqtab),
  ASVs_processadas = nrow(taxa),
  Nota_filtragem = "Filtragem centralizada no Script 5",
  AddSpecies_banco_disponivel = tem_sp_db,
  AddSpecies_executado = addspecies_executado,
  AddSpecies_sucesso = addspecies_sucesso,
  AddSpecies_erro = addspecies_erro,
  Com_addSpecies = addspecies_sucesso,
  Tempo_assign_s = round(tempo_assign, 0),
  stringsAsFactors = FALSE)
write.csv(run_metadata,
          file.path(output_path, paste0("metadata_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)
log_msg(paste("OK: metadata_", BANCO_NOME, ".csv"), "SAVE")

# ==============================================================================
# RESUMO
# ==============================================================================
cat("\n=============================================================\n")
cat("SILVA 138.2: CONCLUIDO\n")
cat(sprintf("ASVs retidas   : %d (%.1f%%)\n",
            nrow(taxa_limpa), 100 * nrow(taxa_limpa) / nrow(taxa)))
cat(sprintf("Tempo assign   : %.1f s\n", tempo_assign))
cat(sprintf("Boot medio Genus: %.1f%%  |  Cobertura Genus: %.1f%%\n",
            boot_medio["Genus"], cob_depois["Genus"]))
cat("PROXIMO PASSO: 02b_rdp19.R\n")
cat("=============================================================\n\n")

rm(taxa, taxa_limpa, taxa_raw, taxa_boot_obj, boots, boots_limpa, seqtab)
limpeza_mem()
log_msg("Finalizado", "FINAL")
})
