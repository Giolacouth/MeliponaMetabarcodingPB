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
run_pipeline_script("02e_GSRDB_V3V4_QIIME2.R", "gsr", function(ctx) {
###############################################################################
# SCRIPT 2e — CLASSIFICACAO TAXONOMICA GSR-DB


options(encoding = "UTF-8", stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Biostrings)
})

###############################################################################
# 0. PARAMETROS EDITAVEIS
###############################################################################

VERSAO        <- ""
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
BANCO_NOME    <- ""

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$output_root
db_path <- file.path(
  base_path,
  "bancodados",
  "GSR -Greengenes, Silva e RDP"
)

classifier_qza <- file.path(
  db_path,
  "classifier_GSR-DB_V3-V4.qza"
)

# Entradas canonicas geradas pelo Script 01.
arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
arq_asvmap <- ctx$contracts[["asv_sequences"]]

# Pasta exclusiva desta analise.
gsr_root <- ctx$stage$root
gsr_qiime <- file.path(gsr_root, "qiime2")
gsr_rds <- file.path(gsr_root, "rds")
gsr_tab <- file.path(gsr_root, "tabelas")
gsr_log <- file.path(gsr_root, "logs")
gsr_chk <- file.path(gsr_root, "checkpoints")

for (d in c(gsr_root, gsr_qiime, gsr_rds, gsr_tab, gsr_log, gsr_chk)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# O artigo avaliou confidence = disable. A classificacao com 0.7 e salva como
# sensibilidade, sem substituir o resultado principal.
CONFIDENCE_PRINCIPAL <- "disable"
EXECUTAR_SENSIBILIDADE_07 <- TRUE
CONFIDENCE_SENSIBILIDADE <- "0.7"

N_JOBS <- 1L
READ_ORIENTATION <- "auto"
SOBRESCREVER <- TRUE


QIIME2_ENV <- Sys.getenv("QIIME2_ENV", unset = "")

RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

cat("=============================================================\n")
cat("GSR-DB V3-V4 / QIIME 2 — v", VERSAO, "\n", sep = "")
cat("Data:", DATA_EXECUCAO, "\n")
cat("Saida exclusiva:", gsr_root, "\n")
cat("=============================================================\n\n")

###############################################################################
# 1. FUNCOES DE LOG, I/O E EXECUCAO EXTERNA
###############################################################################

log_file <- file.path(gsr_log, "GSR_DB_V3V4_execucao.log")
command_file <- file.path(gsr_log, "comandos_qiime2_executados.txt")

log_msg <- function(msg, tipo = "INFO") {
  linha <- sprintf("[%s] <%s> %s", format(Sys.time(), "%H:%M:%S"), tipo, msg)
  cat(linha, "\n")
  cat(linha, "\n", file = log_file, append = TRUE)
  invisible(linha)
}

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

validar_arquivo <- function(arq, desc = basename(arq)) {
  if (!file.exists(arq)) {
    abort("%s nao encontrado: %s", desc, arq)
  }
  if (is.na(file.size(arq)) || file.size(arq) == 0) {
    abort("%s existe, mas esta vazio: %s", desc, arq)
  }
  log_msg(sprintf("%s OK (%.2f MB)", desc, file.size(arq) / 1024^2), "OK")
  invisible(TRUE)
}

salvar_csv <- function(df, arq) {
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    df,
    arq,
    row.names = FALSE,
    quote = TRUE,
    fileEncoding = "UTF-8",
    na = "NA"
  )
  validar_arquivo(arq, basename(arq))
  invisible(arq)
}

salvar_rds <- function(obj, arq) {
  dir.create(dirname(arq), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, arq)
  validar_arquivo(arq, basename(arq))
  invisible(arq)
}

remover_se_existe <- function(path) {
  if (!file.exists(path) && !dir.exists(path)) return(invisible(FALSE))
  if (!SOBRESCREVER) {
    abort("Saida ja existe e SOBRESCREVER=FALSE: %s", path)
  }
  unlink(path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

resolver_qiime <- function() {
  qiime_path <- unname(Sys.which("qiime"))
  if (nzchar(qiime_path)) {
    return(list(command = qiime_path, prefix = character(), modo = "PATH"))
  }

  if (nzchar(QIIME2_ENV)) {
    conda_path <- unname(Sys.which("conda"))
    if (nzchar(conda_path)) {
      return(list(
        command = conda_path,
        prefix = c("run", "-n", shQuote(QIIME2_ENV), "qiime"),
        modo = paste0("conda:", QIIME2_ENV)
      ))
    }

    mamba_path <- unname(Sys.which("mamba"))
    if (nzchar(mamba_path)) {
      return(list(
        command = mamba_path,
        prefix = c("run", "-n", shQuote(QIIME2_ENV), "qiime"),
        modo = paste0("mamba:", QIIME2_ENV)
      ))
    }

    micromamba_path <- unname(Sys.which("micromamba"))
    if (nzchar(micromamba_path)) {
      return(list(
        command = micromamba_path,
        prefix = c("run", "-n", shQuote(QIIME2_ENV), "qiime"),
        modo = paste0("micromamba:", QIIME2_ENV)
      ))
    }
  }

  abort(
    paste0(
      "QIIME 2 nao encontrado. Ative o ambiente antes de executar o R, por exemplo:\n"
     
    )
  )
}

qiime_exec <- resolver_qiime()
log_msg(paste("QIIME 2 resolvido por", qiime_exec$modo), "OK")

run_qiime <- function(args, etapa) {
  args_full <- c(qiime_exec$prefix, args)
  comando_legivel <- paste(shQuote(qiime_exec$command), paste(args_full, collapse = " "))

  cat(sprintf("\n# [%s] %s\n%s\n", DATA_EXECUCAO, etapa, comando_legivel),
      file = command_file, append = TRUE)
  log_msg(paste("Executando:", etapa), "RUN")

  saida <- system2(
    command = qiime_exec$command,
    args = args_full,
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(saida, "status")
  if (is.null(status)) status <- 0L

  cat(paste(saida, collapse = "\n"), "\n",
      file = log_file, append = TRUE)

  if (!identical(as.integer(status), 0L)) {
    abort(
      "Falha na etapa QIIME 2 '%s' (status %d). Consulte: %s",
      etapa,
      as.integer(status),
      log_file
    )
  }

  log_msg(paste("Concluido:", etapa), "OK")
  saida
}

###############################################################################
# 2. VALIDACAO DAS ENTRADAS DO PIPELINE DADA2
###############################################################################

cat("=== VALIDACAO DAS ENTRADAS ===\n\n")

validar_arquivo(arq_seqtab, "seqtab_global_nochim.rds")
validar_arquivo(arq_asvmap, "ASV_sequences.tsv")
validar_arquivo(classifier_qza, "classifier_GSR-DB_V3-V4.qza")

seqtab <- tryCatch(
  as.matrix(readRDS(arq_seqtab)),
  error = function(e) abort("Falha ao ler seqtab_global_nochim.rds: %s", conditionMessage(e))
)

if (nrow(seqtab) == 0L || ncol(seqtab) == 0L) {
  abort("seqtab_global_nochim.rds esta vazio.")
}
if (!is.numeric(seqtab)) {
  abort("seqtab_global_nochim.rds deve conter contagens numericas.")
}
if (nrow(seqtab) != 10L) {
  abort(
    "seqtab_global_nochim.rds deve conter 10 amostras; encontrou %d.",
    nrow(seqtab)
  )
}
if (is.null(rownames(seqtab)) || any(rownames(seqtab) == "")) {
  abort("seqtab_global_nochim.rds sem SampleID validos nas linhas.")
}
if (is.null(colnames(seqtab)) || any(colnames(seqtab) == "")) {
  abort("seqtab_global_nochim.rds sem sequencias nas colunas.")
}
if (anyDuplicated(rownames(seqtab))) {
  abort("seqtab_global_nochim.rds possui SampleID duplicados.")
}
if (anyDuplicated(colnames(seqtab))) {
  abort("seqtab_global_nochim.rds possui sequencias duplicadas.")
}
if (anyNA(seqtab) || any(!is.finite(seqtab)) || any(seqtab < 0)) {
  abort("seqtab_global_nochim.rds possui contagens invalidas.")
}
if (max(abs(seqtab - round(seqtab))) > 1e-8) {
  abort("seqtab_global_nochim.rds deve conter contagens inteiras.")
}
if (any(colSums(seqtab) == 0)) {
  abort("seqtab_global_nochim.rds possui ASV(s) com soma zero.")
}
if (any(rowSums(seqtab) == 0)) {
  abort("seqtab_global_nochim.rds possui amostra(s) com zero reads.")
}
if (!all(grepl("^[ACGTRYSWKMBDHVN]+$", colnames(seqtab), ignore.case = TRUE))) {
  abort("Uma ou mais colunas do seqtab nao parecem sequencias de DNA IUPAC.")
}

asv_map <- read.delim(
  arq_asvmap,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

req_map <- c("ASV_ID", "Sequence", "Origem")
faltam_map <- setdiff(req_map, colnames(asv_map))
if (length(faltam_map) > 0L) {
  abort("ASV_sequences.tsv sem coluna(s): %s", paste(faltam_map, collapse = ", "))
}

asv_map$ASV_ID <- trimws(as.character(asv_map$ASV_ID))
asv_map$Sequence <- toupper(trimws(as.character(asv_map$Sequence)))
asv_map$Origem <- trimws(as.character(asv_map$Origem))

if (anyNA(asv_map$ASV_ID) || any(asv_map$ASV_ID == "")) {
  abort("ASV_sequences.tsv possui ASV_ID ausente ou vazio.")
}
if (anyNA(asv_map$Sequence) || any(asv_map$Sequence == "")) {
  abort("ASV_sequences.tsv possui Sequence ausente ou vazia.")
}
if (anyNA(asv_map$Origem) || any(asv_map$Origem == "")) {
  abort("ASV_sequences.tsv possui Origem ausente ou vazia.")
}
if (anyDuplicated(asv_map$ASV_ID)) {
  abort("ASV_sequences.tsv possui ASV_ID duplicados.")
}
if (anyDuplicated(asv_map$Sequence)) {
  abort("ASV_sequences.tsv possui sequencias duplicadas.")
}

seqs_canonicas <- toupper(colnames(seqtab))
idx_map <- match(seqs_canonicas, asv_map$Sequence)
if (anyNA(idx_map)) {
  abort(
    "%d sequencia(s) do seqtab estao ausentes em ASV_sequences.tsv.",
    sum(is.na(idx_map))
  )
}

map_ord <- asv_map[idx_map, , drop = FALSE]
if (anyDuplicated(map_ord$ASV_ID)) {
  abort("O alinhamento seqtab -> ASV_sequences.tsv gerou ASV_ID duplicados.")
}

cat(sprintf("Seqtab: %d amostras | %d ASVs\n", nrow(seqtab), ncol(seqtab)))
cat("Origem das ASVs:\n")
print(table(map_ord$Origem, useNA = "ifany"))
cat("\n")

salvar_rds(
  list(seqtab_dim = dim(seqtab), asv_map = map_ord),
  file.path(gsr_chk, "checkpoint_01_entradas_validadas.rds")
)

###############################################################################
# 3. FASTA EXCLUSIVO PARA O GSR
###############################################################################

cat("=== CONSTRUCAO DO FASTA GSR ===\n\n")

fasta_gsr <- file.path(gsr_qiime, "ASVs_GSR_V3V4.fasta")
repseqs_qza <- file.path(gsr_qiime, "ASVs_GSR_V3V4_rep_seqs.qza")

remover_se_existe(fasta_gsr)
remover_se_existe(repseqs_qza)

asv_dna <- Biostrings::DNAStringSet(seqs_canonicas)
names(asv_dna) <- map_ord$ASV_ID
Biostrings::writeXStringSet(asv_dna, fasta_gsr, format = "fasta", width = 80L)
validar_arquivo(fasta_gsr, "ASVs_GSR_V3V4.fasta")

# Confere que os IDs e sequencias gravados sao exatamente os canonicos.
fasta_check <- Biostrings::readDNAStringSet(fasta_gsr)
if (!identical(names(fasta_check), map_ord$ASV_ID)) {
  abort("Os headers do FASTA GSR nao coincidem com ASV_sequences.tsv.")
}
if (!identical(unname(as.character(fasta_check)), unname(seqs_canonicas))) {
  abort("As sequencias do FASTA GSR divergem do seqtab_global_nochim.rds.")
}

log_msg(sprintf("FASTA GSR criado com %d ASVs.", length(fasta_check)), "OK")

###############################################################################
# 4. VALIDACAO DO QIIME 2
###############################################################################

cat("=== VALIDACAO DO CLASSIFICADOR GSR ===\n\n")

qiime_info <- tryCatch(
  run_qiime(c("info"), "registrar versao do QIIME 2"),
  error = function(e) {
    log_msg(conditionMessage(e), "WARN")
    character()
  }
)
if (length(qiime_info) > 0L) {
  writeLines(qiime_info, file.path(gsr_log, "qiime_info.txt"))
}

peek_classifier <- run_qiime(
  c("tools", "peek", shQuote(classifier_qza)),
  "inspecionar classifier_GSR-DB_V3-V4.qza"
)
writeLines(peek_classifier, file.path(gsr_log, "classifier_peek.txt"))

linha_tipo <- grep("^Type:", trimws(peek_classifier), value = TRUE)
if (length(linha_tipo) != 1L || !grepl("TaxonomicClassifier", linha_tipo, fixed = TRUE)) {
  abort(
    paste0(
      "O arquivo .qza nao foi reconhecido como TaxonomicClassifier.\n",
      "Resultado de 'qiime tools peek':\n%s"
    ),
    paste(peek_classifier, collapse = "\n")
  )
}

run_qiime(
  c("tools", "validate", shQuote(classifier_qza)),
  "validar integridade do classificador GSR"
)

classifier_md5 <- unname(tools::md5sum(classifier_qza))
log_msg(paste("MD5 do classificador:", classifier_md5), "INFO")

###############################################################################
# 5. IMPORTAR ASVs COMO FeatureData[Sequence]
###############################################################################

run_qiime(
  c(
    "tools", "import",
    "--input-path", shQuote(fasta_gsr),
    "--output-path", shQuote(repseqs_qza),
    "--type", shQuote("FeatureData[Sequence]")
  ),
  "importar ASVs como FeatureData[Sequence]"
)

peek_repseqs <- run_qiime(
  c("tools", "peek", shQuote(repseqs_qza)),
  "inspecionar ASVs importadas"
)
if (!any(grepl("FeatureData[Sequence]", peek_repseqs, fixed = TRUE))) {
  abort("ASVs importadas nao foram reconhecidas como FeatureData[Sequence].")
}

###############################################################################
# 6. FUNCOES PARA CLASSIFICAR, EXPORTAR E IMPORTAR PARA O R
###############################################################################

PLACEHOLDER_REGEX_GSR <- paste0(
  "^(unassigned|unknown|unclassified|",
  "uncultured(?:\\s+.*)?|unidentified(?:\\s+.*)?|",
  "environmental(?:\\s+.*)?|metagenome(?:\\s+.*)?|",
  "na|n/a|none|null|root)$"
)

normalizar_nome_rank <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^[dkpcofgs]__", "", x, ignore.case = TRUE)
  x[x == ""] <- NA_character_
  idx_placeholder <- !is.na(x) & grepl(
    PLACEHOLDER_REGEX_GSR,
    x,
    ignore.case = TRUE,
    perl = TRUE
  )
  x[idx_placeholder] <- NA_character_
  x
}

# Parser detalhado de linhagens QIIME 2.
#

parse_lineage_detalhado <- function(lineage) {
  out <- setNames(rep(NA_character_, length(RANKS)), RANKS)

  auditoria_vazia <- function(status) {
    data.frame(
      Status_parse = status,
      N_partes = 0L,
      N_prefixadas = 0L,
      N_sem_prefixo = 0L,
      N_sem_prefixo_recuperadas = 0L,
      N_sem_prefixo_nao_recuperadas = 0L,
      N_prefixos_desconhecidos = 0L,
      N_conflitos = 0L,
      Posicoes_prefixadas_coerentes = NA,
      Root_extrarank_removido = FALSE,
      Tokens_sem_prefixo = NA_character_,
      Tokens_nao_recuperados = NA_character_,
      Prefixos_desconhecidos = NA_character_,
      Conflitos_parse = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (length(lineage) != 1L || is.na(lineage) || trimws(lineage) == "") {
    return(list(taxa = out, audit = auditoria_vazia("lineage_ausente")))
  }

  lineage_raw <- trimws(as.character(lineage))
  if (tolower(lineage_raw) %in% c("unassigned", "unknown", "na", "n/a")) {
    return(list(taxa = out, audit = auditoria_vazia("lineage_nao_atribuida")))
  }

  # Acrescentar um separador temporario preserva campos vazios finais.
  partes <- strsplit(paste0(lineage_raw, ";"), ";", fixed = TRUE)[[1L]]
  partes <- partes[-length(partes)]
  partes <- trimws(partes)

  # "Root" sem prefixo representa um nivel acima de Kingdom e nao deve
  # deslocar o fallback posicional.
  root_removido <- FALSE
  idx_primeiro <- which(nzchar(partes))[1L]
  if (!is.na(idx_primeiro) &&
      grepl("^Root$", partes[idx_primeiro], ignore.case = TRUE)) {
    partes <- partes[-idx_primeiro]
    root_removido <- TRUE
  }

  if (length(partes) == 0L || !any(nzchar(partes))) {
    aud <- auditoria_vazia("lineage_sem_token_util")
    aud$Root_extrarank_removido <- root_removido
    return(list(taxa = out, audit = aud))
  }

  prefix_map <- c(
    d = "Kingdom", k = "Kingdom", p = "Phylum", c = "Class",
    o = "Order", f = "Family", g = "Genus", s = "Species"
  )
  rank_pos <- setNames(seq_along(RANKS), RANKS)

  tipo <- rep("vazio", length(partes))
  rank_prefixo <- rep(NA_character_, length(partes))
  valor_token <- rep(NA_character_, length(partes))
  prefixos_desconhecidos <- character()
  tokens_sem_prefixo <- character()
  conflitos <- character()
  ranks_bloqueados <- character()

  for (i in seq_along(partes)) {
    token <- partes[i]
    if (!nzchar(token)) next

    m <- regexec(
      "^([dkpcofgs])__\\s*(.*)$",
      token,
      ignore.case = TRUE,
      perl = TRUE
    )
    g <- regmatches(token, m)[[1L]]

    if (length(g) == 3L) {
      tipo[i] <- "prefixado"
      rank_prefixo[i] <- unname(prefix_map[tolower(g[2L])])
      valor_token[i] <- normalizar_nome_rank(g[3L])
    } else if (grepl("^[A-Za-z]__", token, perl = TRUE)) {
      tipo[i] <- "prefixo_desconhecido"
      prefixos_desconhecidos <- c(
        prefixos_desconhecidos,
        sprintf("pos%d=%s", i, token)
      )
    } else {
      tipo[i] <- "sem_prefixo"
      valor_token[i] <- normalizar_nome_rank(token)
      tokens_sem_prefixo <- c(
        tokens_sem_prefixo,
        sprintf("pos%d=%s", i, token)
      )
    }
  }

  # Primeiro, os tokens prefixados: o prefixo e a fonte primaria do rank.
  idx_prefixados <- which(tipo == "prefixado")
  for (i in idx_prefixados) {
    rank <- rank_prefixo[i]
    valor <- valor_token[i]
    if (is.na(valor)) next

    if (rank %in% ranks_bloqueados) next

    if (is.na(out[rank])) {
      out[rank] <- valor
    } else if (!identical(out[rank], valor)) {
      conflitos <- c(
        conflitos,
        sprintf(
          "rank_%s_duplicado:valor1=%s;valor2=%s;pos=%d",
          rank, out[rank], valor, i
        )
      )
      out[rank] <- NA_character_
      ranks_bloqueados <- union(ranks_bloqueados, rank)
    }
  }

  posicoes_prefixadas_coerentes <- if (length(idx_prefixados) == 0L) {
    NA
  } else {
    all(vapply(
      idx_prefixados,
      function(i) unname(rank_pos[rank_prefixo[i]]) == i,
      logical(1)
    ))
  }

 
  idx_sem_prefixo <- which(tipo == "sem_prefixo")
  recuperadas <- 0L
  nao_recuperadas <- character()

  for (i in idx_sem_prefixo) {
    valor <- valor_token[i]
    if (is.na(valor)) next

    if (i > length(RANKS)) {
      nao_recuperadas <- c(
        nao_recuperadas,
        sprintf("pos%d=%s[motivo=fora_dos_7_ranks]", i, partes[i])
      )
      next
    }

    if (length(idx_prefixados) > 0L && !isTRUE(posicoes_prefixadas_coerentes)) {
      nao_recuperadas <- c(
        nao_recuperadas,
        sprintf("pos%d=%s[motivo=posicao_mista_insegura]", i, partes[i])
      )
      next
    }

    rank <- RANKS[i]
    if (rank %in% ranks_bloqueados) {
      nao_recuperadas <- c(
        nao_recuperadas,
        sprintf("pos%d=%s[motivo=rank_bloqueado_por_conflito]", i, partes[i])
      )
      next
    }

    if (is.na(out[rank])) {
      out[rank] <- valor
      recuperadas <- recuperadas + 1L
    } else if (!identical(out[rank], valor)) {
      conflitos <- c(
        conflitos,
        sprintf(
          "rank_%s_prefixo_vs_posicao:prefixado=%s;posicional=%s;pos=%d",
          rank, out[rank], valor, i
        )
      )
      out[rank] <- NA_character_
      ranks_bloqueados <- union(ranks_bloqueados, rank)
      nao_recuperadas <- c(
        nao_recuperadas,
        sprintf("pos%d=%s[motivo=conflito_com_prefixo]", i, partes[i])
      )
    }
  }

  n_prefixadas <- length(idx_prefixados)
  n_sem_prefixo <- length(idx_sem_prefixo)
  n_desconhecidos <- sum(tipo == "prefixo_desconhecido")
  n_conflitos <- length(conflitos)
  n_nao_recuperadas <- length(nao_recuperadas)

  base_status <- if (n_prefixadas > 0L && n_sem_prefixo > 0L) {
    if (isTRUE(posicoes_prefixadas_coerentes)) {
      "misto_recuperado_posicional"
    } else {
      "misto_posicao_insegura"
    }
  } else if (n_prefixadas > 0L) {
    "somente_prefixos"
  } else if (n_sem_prefixo > 0L) {
    "sem_prefixos_posicional"
  } else if (n_desconhecidos > 0L) {
    "somente_prefixos_desconhecidos"
  } else {
    "lineage_sem_token_util"
  }

  precisa_revisao <- n_conflitos > 0L ||
    n_desconhecidos > 0L ||
    n_nao_recuperadas > 0L

  status <- if (precisa_revisao) {
    paste0(base_status, "_REVISAR")
  } else {
    base_status
  }

  audit <- data.frame(
    Status_parse = status,
    N_partes = sum(nzchar(partes)),
    N_prefixadas = n_prefixadas,
    N_sem_prefixo = n_sem_prefixo,
    N_sem_prefixo_recuperadas = recuperadas,
    N_sem_prefixo_nao_recuperadas = n_nao_recuperadas,
    N_prefixos_desconhecidos = n_desconhecidos,
    N_conflitos = n_conflitos,
    Posicoes_prefixadas_coerentes = posicoes_prefixadas_coerentes,
    Root_extrarank_removido = root_removido,
    Tokens_sem_prefixo = if (length(tokens_sem_prefixo) > 0L) {
      paste(tokens_sem_prefixo, collapse = " | ")
    } else {
      NA_character_
    },
    Tokens_nao_recuperados = if (n_nao_recuperadas > 0L) {
      paste(nao_recuperadas, collapse = " | ")
    } else {
      NA_character_
    },
    Prefixos_desconhecidos = if (n_desconhecidos > 0L) {
      paste(prefixos_desconhecidos, collapse = " | ")
    } else {
      NA_character_
    },
    Conflitos_parse = if (n_conflitos > 0L) {
      paste(conflitos, collapse = " | ")
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )

  list(taxa = out, audit = audit)
}

# Wrapper mantido para compatibilidade com chamadas unitarias.
parse_lineage <- function(lineage) {
  parse_lineage_detalhado(lineage)$taxa
}

eh_placeholder <- function(x) {
  ifelse(
    is.na(x),
    TRUE,
    grepl(
      PLACEHOLDER_REGEX_GSR,
      trimws(x),
      ignore.case = TRUE,
      perl = TRUE
    )
  )
}


eh_composta <- function(x, rank) {
  if (is.na(x) || trimws(x) == "") return(FALSE)
  x <- trimws(x)

  if (grepl(":", x, fixed = TRUE)) return(TRUE)

  if (rank %in% c("Genus", "Species")) {
    # Dois nomes taxonomicos capitalizados unidos por hifen.
    if (grepl("[A-Z][A-Za-z]+-[A-Z][A-Za-z]+", x, perl = TRUE)) return(TRUE)
  }

  FALSE
}

limpar_taxonomia_gsr <- function(taxa_raw) {
  out <- taxa_raw

  for (rk in colnames(out)) {
    out[, rk] <- trimws(as.character(out[, rk]))
    out[eh_placeholder(out[, rk]), rk] <- NA_character_
  }

  # Nao harmonizar Phylum neste script. 

  comp_genus <- vapply(out[, "Genus"], eh_composta, logical(1), rank = "Genus")
  comp_species <- vapply(out[, "Species"], eh_composta, logical(1), rank = "Species")

  # Sem genero univoco, a especie tambem nao pode ser considerada univoca.
  out[comp_genus, "Genus"] <- NA_character_
  out[comp_genus | comp_species, "Species"] <- NA_character_

  species_sem_genus <- !is.na(out[, "Species"]) & is.na(out[, "Genus"])
  out[species_sem_genus, "Species"] <- NA_character_

  list(
    taxa = out,
    comp_genus = comp_genus,
    comp_species = comp_species,
    species_sem_genus = species_sem_genus
  )
}

cobertura_por_rank <- function(mat, modo) {
  data.frame(
    Modo = modo,
    Rank = colnames(mat),
    N_classificadas = colSums(!is.na(mat) & trimws(mat) != ""),
    N_total = nrow(mat),
    Cobertura_pct = round(
      100 * colSums(!is.na(mat) & trimws(mat) != "") / nrow(mat),
      2
    ),
    stringsAsFactors = FALSE
  )
}

classificar_gsr <- function(confidence_value, sufixo) {
  tax_qza <- file.path(gsr_qiime, paste0("taxonomy_GSR_V3V4_", sufixo, ".qza"))
  export_dir <- file.path(gsr_qiime, paste0("taxonomy_GSR_V3V4_", sufixo, "_export"))

  remover_se_existe(tax_qza)
  remover_se_existe(export_dir)

  run_qiime(
    c(
      "feature-classifier", "classify-sklearn",
      "--i-reads", shQuote(repseqs_qza),
      "--i-classifier", shQuote(classifier_qza),
      "--p-confidence", as.character(confidence_value),
      "--p-n-jobs", as.character(N_JOBS),
      "--p-read-orientation", READ_ORIENTATION,
      "--o-classification", shQuote(tax_qza)
    ),
    paste0("classificar GSR V3-V4 [", sufixo, "]")
  )

  run_qiime(
    c(
      "tools", "export",
      "--input-path", shQuote(tax_qza),
      "--output-path", shQuote(export_dir)
    ),
    paste0("exportar taxonomia GSR [", sufixo, "]")
  )

  taxonomy_candidates <- list.files(
    export_dir,
    pattern = "taxonomy\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(taxonomy_candidates) != 1L) {
    abort(
      "Esperado exatamente um taxonomy.tsv em %s; encontrados: %d",
      export_dir,
      length(taxonomy_candidates)
    )
  }
  taxonomy_tsv <- taxonomy_candidates[[1L]]

  tab <- read.delim(
    taxonomy_tsv,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  id_col <- intersect(c("Feature ID", "Feature.ID", "FeatureID"), colnames(tab))[1]
  tax_col <- intersect(c("Taxon", "Taxonomy"), colnames(tab))[1]
  conf_col <- intersect(c("Confidence", "confidence"), colnames(tab))[1]

  if (is.na(id_col) || is.na(tax_col)) {
    abort(
      "taxonomy.tsv sem colunas de ID/Taxon. Colunas encontradas: %s",
      paste(colnames(tab), collapse = ", ")
    )
  }

  tab[[id_col]] <- trimws(as.character(tab[[id_col]]))
  if (anyDuplicated(tab[[id_col]])) {
    abort("taxonomy.tsv possui Feature ID duplicados.")
  }

  faltam <- setdiff(map_ord$ASV_ID, tab[[id_col]])
  extras <- setdiff(tab[[id_col]], map_ord$ASV_ID)
  if (length(faltam) > 0L || length(extras) > 0L) {
    abort(
      paste0(
        "Universo do taxonomy.tsv diverge do mapa canonico. ",
        "Faltantes=%d; extras=%d."
      ),
      length(faltam),
      length(extras)
    )
  }

  tab <- tab[match(map_ord$ASV_ID, tab[[id_col]]), , drop = FALSE]
  stopifnot(identical(tab[[id_col]], map_ord$ASV_ID))

  parsed_lineages <- lapply(tab[[tax_col]], parse_lineage_detalhado)

  taxa_raw <- do.call(
    rbind,
    lapply(parsed_lineages, function(x) x$taxa)
  )
  taxa_raw <- as.matrix(taxa_raw)
  colnames(taxa_raw) <- RANKS
  rownames(taxa_raw) <- seqs_canonicas
  storage.mode(taxa_raw) <- "character"

  auditoria_parse_core <- do.call(
    rbind,
    lapply(parsed_lineages, function(x) x$audit)
  )
  rownames(auditoria_parse_core) <- NULL

  auditoria_parse <- data.frame(
    ASV_ID = map_ord$ASV_ID,
    ASV_seq = seqs_canonicas,
    Origem = map_ord$Origem,
    Taxon_QIIME2_raw = tab[[tax_col]],
    auditoria_parse_core,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  n_lineages_mistas <- sum(
    auditoria_parse$N_prefixadas > 0L &
      auditoria_parse$N_sem_prefixo > 0L
  )
  n_lineages_revisar <- sum(grepl("_REVISAR$", auditoria_parse$Status_parse))
  n_conflitos_parse <- sum(auditoria_parse$N_conflitos > 0L)
  n_tokens_recuperados <- sum(auditoria_parse$N_sem_prefixo_recuperadas)
  n_tokens_nao_recuperados <- sum(auditoria_parse$N_sem_prefixo_nao_recuperadas)

  if (n_lineages_mistas > 0L) {
    log_msg(
      sprintf(
        "%s: %d linhagem(ns) com prefixos mistos; ver auditoria de parse.",
        sufixo, n_lineages_mistas
      ),
      "WARN"
    )
  }
  if (n_lineages_revisar > 0L) {
    log_msg(
      sprintf(
        "%s: %d linhagem(ns) exigem revisao de parse; nenhum token foi descartado silenciosamente.",
        sufixo, n_lineages_revisar
      ),
      "WARN"
    )
  }

  limpeza <- limpar_taxonomia_gsr(taxa_raw)
  taxa_limpa <- limpeza$taxa
  rownames(taxa_limpa) <- seqs_canonicas

  conf_num <- rep(NA_real_, nrow(tab))
  if (!is.na(conf_col)) {
    conf_num <- suppressWarnings(as.numeric(tab[[conf_col]]))
  }
  names(conf_num) <- seqs_canonicas

  tabela_saida <- data.frame(
    ASV_ID = map_ord$ASV_ID,
    ASV_seq = seqs_canonicas,
    Origem = map_ord$Origem,
    Taxon_QIIME2_raw = tab[[tax_col]],
    Confidence_QIIME2 = conf_num,
    Parse_status = auditoria_parse$Status_parse,
    Parse_n_prefixadas = auditoria_parse$N_prefixadas,
    Parse_n_sem_prefixo = auditoria_parse$N_sem_prefixo,
    Parse_n_conflitos = auditoria_parse$N_conflitos,
    as.data.frame(taxa_raw, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tabela_limpa <- data.frame(
    ASV_ID = map_ord$ASV_ID,
    ASV_seq = seqs_canonicas,
    Origem = map_ord$Origem,
    Confidence_QIIME2 = conf_num,
    as.data.frame(taxa_limpa, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  idx_audit <- limpeza$comp_genus | limpeza$comp_species | limpeza$species_sem_genus
  n_audit <- sum(idx_audit)

 
  auditoria_composta <- data.frame(
    ASV_ID = map_ord$ASV_ID[idx_audit],
    ASV_seq = seqs_canonicas[idx_audit],
    Origem = map_ord$Origem[idx_audit],
    Genus_raw = taxa_raw[idx_audit, "Genus"],
    Species_raw = taxa_raw[idx_audit, "Species"],
    Genus_composto = limpeza$comp_genus[idx_audit],
    Species_composta = limpeza$comp_species[idx_audit],
    Species_sem_Genus_univoco = limpeza$species_sem_genus[idx_audit],
    Decisao = rep(
      "preservado_na_tabela_bruta_e_convertido_em_NA_na_matriz_limpa",
      n_audit
    ),
    stringsAsFactors = FALSE
  )

  salvar_rds(
    taxa_raw,
    file.path(gsr_rds, paste0("taxa_gsr_v3v4_", sufixo, ".rds"))
  )
  salvar_rds(
    taxa_limpa,
    file.path(gsr_rds, paste0("taxa_gsr_v3v4_", sufixo, "_limpa.rds"))
  )
  salvar_rds(
    conf_num,
    file.path(gsr_rds, paste0("confianca_gsr_v3v4_", sufixo, ".rds"))
  )

  salvar_csv(
    tabela_saida,
    file.path(gsr_tab, paste0("taxa_gsr_v3v4_", sufixo, ".csv"))
  )
  salvar_csv(
    tabela_limpa,
    file.path(gsr_tab, paste0("taxa_gsr_v3v4_", sufixo, "_limpa.csv"))
  )
  salvar_csv(
    auditoria_composta,
    file.path(gsr_tab, paste0("atribuicoes_compostas_gsr_v3v4_", sufixo, ".csv"))
  )
  salvar_csv(
    auditoria_parse,
    file.path(gsr_tab, paste0("auditoria_parse_lineage_gsr_v3v4_", sufixo, ".csv"))
  )

  cob_raw <- cobertura_por_rank(taxa_raw, paste0(sufixo, "_raw"))
  cob_limpa <- cobertura_por_rank(taxa_limpa, paste0(sufixo, "_limpa"))
  cobertura <- rbind(cob_raw, cob_limpa)
  salvar_csv(
    cobertura,
    file.path(gsr_tab, paste0("cobertura_gsr_v3v4_", sufixo, ".csv"))
  )

  salvar_rds(
    list(
      taxonomy_qza = tax_qza,
      taxonomy_tsv = taxonomy_tsv,
      taxa_raw = taxa_raw,
      taxa_limpa = taxa_limpa,
      confidence = conf_num,
      auditoria_composta = auditoria_composta,
      auditoria_parse = auditoria_parse,
      cobertura = cobertura
    ),
    file.path(gsr_chk, paste0("checkpoint_02_classificacao_", sufixo, ".rds"))
  )

  list(
    sufixo = sufixo,
    confidence_param = confidence_value,
    tax_qza = tax_qza,
    taxonomy_tsv = taxonomy_tsv,
    taxa_raw = taxa_raw,
    taxa_limpa = taxa_limpa,
    confidence = conf_num,
    cobertura = cobertura,
    auditoria_parse = auditoria_parse,
    n_compostas = nrow(auditoria_composta),
    n_lineages_mistas = n_lineages_mistas,
    n_lineages_revisar = n_lineages_revisar,
    n_conflitos_parse = n_conflitos_parse,
    n_tokens_sem_prefixo_recuperados = n_tokens_recuperados,
    n_tokens_sem_prefixo_nao_recuperados = n_tokens_nao_recuperados
  )
}

###############################################################################
# 7. CLASSIFICACAO PRINCIPAL E SENSIBILIDADE
###############################################################################

cat("=== CLASSIFICACAO PRINCIPAL GSR ===\n\n")

res_principal <- classificar_gsr(
  confidence_value = CONFIDENCE_PRINCIPAL,
  sufixo = "confidence_disable"
)

# Alias de contrato para facilitar comparacao futura no Script 5.
# Estes arquivos representam a classificacao principal do teste GSR.
salvar_rds(
  res_principal$taxa_raw,
  file.path(gsr_rds, "taxa_gsr_v3v4.rds")
)
salvar_rds(
  res_principal$taxa_limpa,
  file.path(gsr_rds, "taxa_gsr_v3v4_limpa.rds")
)
salvar_rds(
  res_principal$confidence,
  file.path(gsr_rds, "confianca_gsr_v3v4.rds")
)

res_sens <- NULL
if (isTRUE(EXECUTAR_SENSIBILIDADE_07)) {
  cat("\n=== SENSIBILIDADE GSR — CONFIDENCE 0.7 ===\n\n")
  res_sens <- classificar_gsr(
    confidence_value = CONFIDENCE_SENSIBILIDADE,
    sufixo = "confidence_0.7"
  )
}

###############################################################################
# 8. COMPARACAO ENTRE CONFIGURACOES DE CONFIANCA
###############################################################################

if (!is.null(res_sens)) {
  comparar_rank <- function(rank) {
    a <- res_principal$taxa_limpa[, rank]
    b <- res_sens$taxa_limpa[, rank]

    ambos <- !is.na(a) & !is.na(b)

    data.frame(
      Rank = rank,
      N_total = length(a),
      N_disable = sum(!is.na(a)),
      N_confidence_0.7 = sum(!is.na(b)),
      Apenas_disable = sum(!is.na(a) & is.na(b)),
      Apenas_confidence_0.7 = sum(is.na(a) & !is.na(b)),
      Concordantes_entre_informativos = sum(ambos & a == b),
      Divergentes_entre_informativos = sum(ambos & a != b),
      stringsAsFactors = FALSE
    )
  }

  comparacao_conf <- do.call(rbind, lapply(RANKS, comparar_rank))
  salvar_csv(
    comparacao_conf,
    file.path(gsr_tab, "comparacao_confidence_disable_vs_0.7_por_rank.csv")
  )

  divergencias <- do.call(
    rbind,
    lapply(RANKS, function(rank) {
      a <- res_principal$taxa_limpa[, rank]
      b <- res_sens$taxa_limpa[, rank]
      idx <- (!is.na(a) | !is.na(b)) & (is.na(a) != is.na(b) | (!is.na(a) & !is.na(b) & a != b))

      if (!any(idx)) return(NULL)

      data.frame(
        ASV_ID = map_ord$ASV_ID[idx],
        ASV_seq = seqs_canonicas[idx],
        Origem = map_ord$Origem[idx],
        Rank = rank,
        Confidence_disable = a[idx],
        Confidence_0.7 = b[idx],
        stringsAsFactors = FALSE
      )
    })
  )

  if (is.null(divergencias)) {
    divergencias <- data.frame(
      ASV_ID = character(), ASV_seq = character(), Origem = character(),
      Rank = character(), Confidence_disable = character(),
      Confidence_0.7 = character(), stringsAsFactors = FALSE
    )
  }

  salvar_csv(
    divergencias,
    file.path(gsr_tab, "divergencias_confidence_disable_vs_0.7_por_asv.csv")
  )
}

###############################################################################
# 9. METADADOS E CHECKPOINTS FINAIS
###############################################################################

classifier_peek_text <- paste(peek_classifier, collapse = " | ")
qiime_version_text <- if (length(qiime_info) > 0L) {
  paste(qiime_info, collapse = " | ")
} else {
  NA_character_
}

metadata_execucao <- data.frame(
  Script = "02e_GSRDB_V3V4_QIIME2.R",
  Versao = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Base_path = base_path,
  Output_GSR = gsr_root,
  Seqtab_entrada = arq_seqtab,
  ASV_map_entrada = arq_asvmap,
  Classifier_qza = classifier_qza,
  Classifier_MD5 = classifier_md5,
  Classifier_peek = classifier_peek_text,
  QIIME2_execucao = qiime_exec$modo,
  QIIME2_info = qiime_version_text,
  N_amostras = nrow(seqtab),
  N_ASVs = ncol(seqtab),
  Confidence_principal = CONFIDENCE_PRINCIPAL,
  Sensibilidade_0.7 = EXECUTAR_SENSIBILIDADE_07,
  N_jobs = N_JOBS,
  Read_orientation = READ_ORIENTATION,
  N_atribuicoes_compostas_principal = res_principal$n_compostas,
  N_lineages_mistas_principal = res_principal$n_lineages_mistas,
  N_lineages_revisar_parse_principal = res_principal$n_lineages_revisar,
  N_lineages_com_conflito_parse_principal = res_principal$n_conflitos_parse,
  N_tokens_sem_prefixo_recuperados_principal = res_principal$n_tokens_sem_prefixo_recuperados,
  N_tokens_sem_prefixo_nao_recuperados_principal = res_principal$n_tokens_sem_prefixo_nao_recuperados,
  Normalizacao_filo = "nao_aplicada_no_script_2e",
  Harmonizacao_filo = "delegada_ao_script_05_para_todos_os_bancos",
  Uso_recomendado = "analise_taxonomica_alternativa_de_sensibilidade",
  stringsAsFactors = FALSE
)

salvar_csv(
  metadata_execucao,
  file.path(gsr_tab, "metadata_execucao_gsr_v3v4.csv")
)

saveRDS(
  list(
    metadata = metadata_execucao,
    principal = res_principal,
    sensibilidade = res_sens
  ),
  file.path(gsr_chk, "checkpoint_03_GSR_concluido.rds")
)

cat("\n=============================================================\n")
cat("GSR-DB V3-V4 — CONCLUIDO\n")
cat(sprintf("Amostras: %d | ASVs: %d\n", nrow(seqtab), ncol(seqtab)))
cat("Classificacao principal: confidence=disable\n")
cat(sprintf(
  "Parse: %d linhagens mistas | %d para revisar | %d conflitos\n",
  res_principal$n_lineages_mistas,
  res_principal$n_lineages_revisar,
  res_principal$n_conflitos_parse
))
cat("Harmonizacao de Phylum: delegada ao Script 05\n")
cat("Taxonomia limpa principal:\n  ",
    file.path(gsr_rds, "taxa_gsr_v3v4_limpa.rds"), "\n", sep = "")
cat("Resultados completos:\n  ", gsr_root, "\n", sep = "")
cat("Este script NAO alterou taxa_consenso_final.rds nem o Script 5.\n")
cat("=============================================================\n\n")

log_msg("Classificacao GSR finalizada.", "FINAL")
})
