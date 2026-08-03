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
run_pipeline_script("02b_rdp19.R", "rdp", function(ctx) {
###############################################################################
#   SCRIPT 02b — CLASSIFICACAO TAXONOMICA — RDP 19
#   #
#   CARACTERISTICAS:
#   - Curadoria pelo RDP Consortium
#   - PASSO UNICO: assignTaxonomy (kmer) ate Species
#     AVISO METODOLOGICO: Species via kmer e menos confiavel que addSpecies
#     (SILVA). Comparacoes de nivel Species entre RDP e SILVA sao assimetricas
#     e nao devem ser tratadas como equivalentes. Restringir comparacoes ao
#     nivel de genero nas analises quantitativas.
#   - RDP pode conter "Root" como nivel extra — removido na normalizacao
#   - Nomenclatura ICNP classica; menos atualizado para LAB novos
#     ex: Apilactobacillus pode aparecer como Lactobacillus
#
#   ENTRADA: seqtab_global_nochim.rds
###############################################################################

library(dada2)

VERSAO        <- "1"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
BANCO_NOME    <- "rdp19"

cat("=============================================================\n")
cat("RDP 19 — Classificacao Taxonomica v", VERSAO, "\n", sep = "")
cat("Data/Hora:", DATA_EXECUCAO, "\n")
cat("=============================================================\n\n")

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$stage$root
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_path)) stop("Falha ao criar diretorio: ", output_path, call. = FALSE)

# ---- ENTRADA CANONICA: identica para todos os bancos kmer ----
arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
arq_rdp    <- file.path(base_path,
                        "bancodados/RDP/rdp_19_toSpecies_trainset.fa.gz")
# Mapa canonico ASV_ID <-> sequencia (Script 01); usado para identificar
# cada ASV por numero nos CSVs deste script (ver com_asv_id()).
arq_asvmap <- ctx$contracts[["asv_sequences"]]

MINBOOT     <- 80
MULTITHREAD <- TRUE
CHECKPOINT  <- TRUE

# ---------------------------------------------------------------------------
# FUNCOES AUXILIARES
# ---------------------------------------------------------------------------
log_msg <- function(msg, tipo = "INFO")
  cat(sprintf("[%s] <%s> %s\n", format(Sys.time(), "%H:%M:%S"), tipo, msg))

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

chkpt <- function(obj, suf) {
  if (!CHECKPOINT) return(invisible(NULL))
  arq <- file.path(output_path,
                   paste0("taxa_", BANCO_NOME, "_", suf, "_checkpoint.rds"))
  saveRDS(obj, arq)
  log_msg(paste("Checkpoint:", basename(arq)), "CHKPT")
}

com_asv <- function(mat)
  cbind(ASV = rownames(mat), as.data.frame(mat, stringsAsFactors = FALSE))

# Carrega o mapa canonico ASV_ID <-> sequencia (mesma logica do Script 02/04).
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

# CSV com ASV_ID (numerico) + ASV_seq, em vez de so a sequencia como ID.
com_asv_id <- function(mat, seq2id) {
  seqs <- rownames(mat)
  ids  <- unname(seq2id[seqs])
  if (any(is.na(ids)))
    log_msg(sprintf("%d sequencia(s) sem ASV_ID correspondente.", sum(is.na(ids))),
            "WARN")
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
cat("NOTA: RDP usa passo unico (kmer ate Species).\n")
cat("      Species-level menos confiavel que addSpecies —",
    "comparar apenas a nivel Genus.\n\n")

cat("=== VALIDACOES ===\n\n")
validar_arquivo(arq_seqtab, "seqtab_global_nochim.rds", TRUE)
validar_arquivo(arq_rdp,    "RDP 19 trainset",          TRUE)
validar_arquivo(arq_asvmap, "ASV_sequences.tsv",        TRUE)
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
# ATRIBUICAO TAXONOMICA
# ==============================================================================
cat("=== assignTaxonomy (kmer, ate Species) ===\n\n")
log_msg("Iniciando assignTaxonomy...", "INICIO")
t0 <- proc.time()

taxa_boot_obj <- tryCatch(
  assignTaxonomy(
    seqtab, arq_rdp,
    multithread      = MULTITHREAD,
    minBoot          = MINBOOT,
    tryRC            = TRUE,
    outputBootstraps = TRUE),
  error = function(e) stop("RDP assignTaxonomy falhou: ", conditionMessage(e), call. = FALSE)
)
if (!is.list(taxa_boot_obj) || !all(c("tax", "boot") %in% names(taxa_boot_obj))) {
  stop("RDP assignTaxonomy nao retornou componentes tax/boot.")
}

taxa  <- taxa_boot_obj$tax
boots <- taxa_boot_obj$boot
tempo_assign <- (proc.time() - t0)["elapsed"]
log_msg(sprintf("assignTaxonomy em %.1f s (%.1f min)",
                tempo_assign, tempo_assign/60), "OK")

taxa  <- normalizar_colunas(taxa)
boots <- normalizar_colunas(boots)
if (nrow(taxa) != ncol(seqtab) ||
    nrow(boots) != ncol(seqtab) ||
    !identical(rownames(taxa), colnames(seqtab)) ||
    !identical(rownames(boots), colnames(seqtab)) ||
    !identical(dim(taxa), dim(boots)) ||
    !identical(colnames(taxa), colnames(boots))) {
  stop(
    "RDP retornou taxonomia/bootstrap desalinhados do universo canonico.",
    call. = FALSE
  )
}
log_msg(paste("Colunas:", paste(colnames(taxa), collapse = ", ")), "INFO")
chkpt(taxa,  "apos_assigntaxonomy")
chkpt(boots, "boots_assigntaxonomy")
limpeza_mem()

cat("\nCobertura inicial (%):\n")
cob_antes <- apply(taxa, 2, function(x)
  round(100 * sum(!is.na(x)) / nrow(taxa), 1))
print(cob_antes)
cat("\nBootstrap medio por nivel:\n")
boot_medio <- apply(boots, 2, function(x)
  round(mean(x[!is.na(x)], na.rm = TRUE), 1))
print(boot_medio)

# ==============================================================================
# AVISO ESPECIFICO PARA MICROBIOMA DE ABELHAS — Apilactobacillus
# ==============================================================================
if ("Genus" %in% colnames(taxa)) {
  gen_values <- taxa[, "Genus"]
  n_lacto <- sum(!is.na(gen_values) & gen_values == "Lactobacillus")
  if (n_lacto > 0) {
    log_msg(paste0(
      "ATENCAO: ", n_lacto, " ASVs classificadas como Lactobacillus. ",
      "Em microbioma de abelhas, parte pode ser Apilactobacillus ",
      "(RDP 19 pode nao ter nomenclatura atualizada). ",
      "Validar com BLAST (Script 4) e comparar com SILVA."), "WARN")
  }
}

taxa_limpa  <- taxa
boots_limpa <- boots

# Cobertura taxonomica apos a limpeza.
cob_depois <- apply(taxa_limpa, 2, function(x)
  round(100 * sum(!is.na(x)) / nrow(taxa_limpa), 1))

species_rdp_df <- data.frame(
  ASV_ID        = unname(seq2id[rownames(taxa_limpa)]),
  ASV           = rownames(taxa_limpa),
  Genus         = if ("Genus" %in% colnames(taxa_limpa)) taxa_limpa[, "Genus"] else NA,
  Species       = if ("Species" %in% colnames(taxa_limpa)) taxa_limpa[, "Species"] else NA,
  Species_fonte = if ("Species" %in% colnames(taxa_limpa))
    ifelse(!is.na(taxa_limpa[, "Species"]), "kmer_only", NA_character_) else NA_character_,
  stringsAsFactors = FALSE)
write.csv(species_rdp_df,
          file.path(output_path, paste0("species_fonte_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)

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
  write.csv(com_asv_id(obj, seq2id), arq_c, row.names = FALSE, quote = TRUE)
  if (!file.exists(arq_c) || is.na(file.size(arq_c)) || file.size(arq_c) == 0L)
    stop("Falha ao salvar CSV ou arquivo vazio: ", arq_c, call. = FALSE)
  log_msg(paste("OK:", basename(arq_r), "|", basename(arq_c)), "SAVE")
}

salvar(taxa,        paste0("taxa_",  BANCO_NOME))
salvar(taxa_limpa,  paste0("taxa_",  BANCO_NOME, "_limpa"))
salvar(boots,       paste0("boots_", BANCO_NOME))
salvar(boots_limpa, paste0("boots_", BANCO_NOME, "_limpa"))

cob_df <- data.frame(
  Nivel      = names(cob_antes),
  Antes_pct  = cob_antes,
  Depois_pct = cob_depois[names(cob_antes)],
  Boot_medio = boot_medio[names(cob_antes)])
write.csv(cob_df,
          file.path(output_path, paste0("cobertura_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)


run_metadata <- data.frame(
  Script        = paste0("2b_", BANCO_NOME, "_V1"),
  Versao        = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Banco         = "RDP 19",
  Seqtab_input  = basename(arq_seqtab),
  Metodo_assign = "assignTaxonomy_kmer_toSpecies",
  MINBOOT       = MINBOOT,
  MULTITHREAD   = MULTITHREAD,
  ASVs_entrada  = ncol(seqtab),
  Amostras      = nrow(seqtab),
  ASVs_saida    = nrow(taxa_limpa),
  Pct_retidas   = round(100 * nrow(taxa_limpa) / nrow(taxa), 1),
  Com_addSpecies  = FALSE,
  Tempo_assign_s  = round(tempo_assign, 0),
  Nota_species    = "kmer_apenas_menos_confiavel_que_addSpecies",
  stringsAsFactors = FALSE)
write.csv(run_metadata,
          file.path(output_path, paste0("metadata_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)
log_msg(paste("OK: metadata_", BANCO_NOME, ".csv"), "SAVE")

cat("\n=============================================================\n")
cat("RDP 19 — CONCLUIDO\n")
cat(sprintf("ASVs retidas: %d (%.1f%%)  |  Tempo: %.1f s\n",
            nrow(taxa_limpa), 100 * nrow(taxa_limpa) / nrow(taxa),
            tempo_assign))
cat(sprintf("Boot Genus: %.1f%%  |  Cobertura Genus: %.1f%%\n",
            boot_medio["Genus"], cob_depois["Genus"]))
cat("AVISO: comparacoes de Species com outros bancos nao sao equivalentes.\n")
cat("PROXIMO PASSO: 02c_Greengenes2.R; 2d_BEExact.R; 02e_GSRDB_V3V4_QIIME2.R; depois 04_rblast.R\n")
cat("=============================================================\n\n")

rm(taxa, taxa_limpa, taxa_boot_obj, boots, boots_limpa, seqtab)
limpeza_mem()
log_msg("Finalizado", "FINAL")
})
