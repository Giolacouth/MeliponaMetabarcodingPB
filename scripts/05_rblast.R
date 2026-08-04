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
run_pipeline_script("04_rblast.R", "blast", function(ctx) {
###############################################################################
# SCRIPT 04 — BLASTN / NCBI 16S PARA ASVs
#
# Objetivo
#   Classificar as ASVs por alinhamento contra o banco local NCBI 16S e
#   produzir evidencias auditaveis para o Script 05.
#
# 
###############################################################################

options(encoding = "UTF-8", stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(Biostrings)
})

VERSAO        <- "2.1_reprocessamento_contrato_taxonomia"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

###############################################################################
# 0. CONFIGURACAO
###############################################################################

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$output_root

blast_root <- ctx$stage$root
blast_in   <- file.path(blast_root, "entrada")
blast_rds  <- file.path(blast_root, "rds")
blast_tab  <- file.path(blast_root, "tabelas")
blast_log  <- file.path(blast_root, "logs")
blast_chk  <- file.path(blast_root, "checkpoints")

for (d in c(output_path, blast_root, blast_in, blast_rds,
            blast_tab, blast_log, blast_chk)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
}

blast_db_path <- file.path(
  base_path,
  "bancodados",
  "ncbi_16S",
  "16S_ribosomal_RNA"
)

arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
arq_asvmap <- ctx$contracts[["asv_sequences"]]
arq_query  <- file.path(blast_in, "ASVs_BLAST_NCBI16S.fasta")

# Coleta ampla. As camadas finais usam criterios mais estritos abaixo.
PERC_IDENT_COLETA <- 80
QCOV_COLETA       <- 80

# Triagem de genero. A aceitacao final ocorre no Script 05.
PERC_IDENT_GENUS <- 97
QCOV_GENUS       <- 100

# Match integral para a camada de especie.
PERC_IDENT_SPECIES <- 100
QCOV_SPECIES       <- 100

# Auditoria de baixa identidade.
PERC_IDENT_BAIXA <- 92
QCOV_BAIXA       <- 100

# Escolha operacional para reduzir truncamento dos melhores hits.
# Nao e um limiar biologico.
MAX_TARGET_SEQS <- 500L
NUM_THREADS     <- 4L

# TRUE reutiliza o TSV bruto ja existente e refaz somente o pos-processamento.
# Isso evita repetir o BLAST quando o universo de ASVs e o banco nao mudaram.
# FALSE forca uma nova busca blastn.
REUTILIZAR_BLAST_BRUTO <- TRUE

###############################################################################
# 1. FUNCOES AUXILIARES
###############################################################################

log_file <- file.path(blast_log, "BLAST_NCBI16S_execucao.log")

log_msg <- function(msg, tipo = "INFO") {
  linha <- sprintf("[%s] <%s> %s", format(Sys.time(), "%H:%M:%S"), tipo, msg)
  cat(linha, "\n")
  cat(linha, "\n", file = log_file, append = TRUE)
  invisible(linha)
}

abort <- function(...) stop(sprintf(...), call. = FALSE)

validar_arquivo <- function(arq, desc = basename(arq), vazio_ok = FALSE) {
  if (!file.exists(arq)) abort("%s nao encontrado: %s", desc, arq)
  if (!vazio_ok && file.size(arq) == 0) abort("%s esta vazio: %s", desc, arq)
  log_msg(sprintf("%s OK (%.2f MB)", desc, file.size(arq) / 1024^2), "OK")
  invisible(TRUE)
}

salvar_rds <- function(obj, arq) {
  saveRDS(obj, arq)
  validar_arquivo(arq, basename(arq))
  invisible(arq)
}

salvar_csv <- function(obj, arq) {
  write.csv(
    obj,
    arq,
    row.names = FALSE,
    quote = TRUE,
    fileEncoding = "UTF-8",
    na = ""
  )
  validar_arquivo(arq, basename(arq), vazio_ok = TRUE)
  invisible(arq)
}

valor_texto <- function(x) {
  x <- trimws(as.character(x))

  invalido <- is.na(x) |
    x == "" |
    toupper(x) %in% c(
      "NA", "N/A", "N.A.", "NULL", "NONE",
      "NOT AVAILABLE", "-"
    )

  x[invalido] <- NA_character_
  x
}

extrair_genero_especie_unico <- function(nome) {
  if (length(nome) != 1L || is.na(nome) || trimws(nome) == "") {
    return(c(Genus = NA_character_, Species = NA_character_))
  }

  x <- trimws(as.character(nome))
  partes <- strsplit(x, "\\s+")[[1L]]

  termos_nao_especie <- c(
    "sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.",
    "bacterium", "archaeon", "uncultured", "unidentified",
    "metagenome", "environmental", "candidate"
  )

  if (length(partes) == 0L) {
    return(c(Genus = NA_character_, Species = NA_character_))
  }

  if (tolower(partes[1L]) == "candidatus") {
    if (length(partes) < 2L) {
      return(c(Genus = NA_character_, Species = NA_character_))
    }
    genero <- paste(partes[1:2], collapse = " ")
    especie <- if (
      length(partes) >= 3L &&
      !tolower(partes[3L]) %in% termos_nao_especie
    ) {
      paste(partes[1:3], collapse = " ")
    } else {
      NA_character_
    }
  } else {
    genero <- partes[1L]
    especie <- if (
      length(partes) >= 2L &&
      !tolower(partes[2L]) %in% termos_nao_especie
    ) {
      paste(partes[1:2], collapse = " ")
    } else {
      NA_character_
    }
  }

  genericos <- c(
    "bacterium", "archaeon", "uncultured", "unidentified",
    "environmental", "metagenome"
  )

  if (
    !grepl("^(Candidatus )?[A-Z][A-Za-z0-9.-]*$", genero) ||
    tolower(genero) %in% genericos
  ) {
    genero <- NA_character_
    especie <- NA_character_
  }

  c(Genus = genero, Species = especie)
}

# Um campo sscinames pode conter mais de um nome separado por ";". Todos
# entram na avaliacao de ambiguidade, em vez de selecionar apenas o primeiro.
extrair_taxa_campo <- function(x) {
  if (length(x) != 1L || is.na(x) || trimws(x) == "") {
    return(data.frame(
      Genus = NA_character_,
      Species = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  nomes <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  nomes <- nomes[nomes != ""]

  if (length(nomes) == 0L) {
    return(data.frame(
      Genus = NA_character_,
      Species = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  mat <- t(vapply(
    nomes,
    extrair_genero_especie_unico,
    FUN.VALUE = c(Genus = NA_character_, Species = NA_character_)
  ))

  as.data.frame(mat, stringsAsFactors = FALSE)
}

selecionar_empates_melhores <- function(x) {
  if (nrow(x) == 0L) return(x)

  x <- x[x$Bits == max(x$Bits, na.rm = TRUE), , drop = FALSE]
  x <- x[x$Perc.Ident == max(x$Perc.Ident, na.rm = TRUE), , drop = FALSE]
  x <- x[
    x$Perc.Query.Coverage == max(x$Perc.Query.Coverage, na.rm = TRUE),
    ,
    drop = FALSE
  ]
  x <- x[x$Mismatches == min(x$Mismatches, na.rm = TRUE), , drop = FALSE]
  x <- x[x$Gap.Openings == min(x$Gap.Openings, na.rm = TRUE), , drop = FALSE]
  x <- x[
    x$Alignment.Length == max(x$Alignment.Length, na.rm = TRUE),
    ,
    drop = FALSE
  ]

  x
}

resumir_hits <- function(df_filtrado, camada) {
  if (nrow(df_filtrado) == 0L) {
    return(data.frame())
  }

  grupos <- split(df_filtrado, df_filtrado$QueryID)

  linhas <- lapply(grupos, function(x) {
    best <- selecionar_empates_melhores(x)

    taxa_list <- lapply(best$Sci.Name, extrair_taxa_campo)
    taxa_df <- do.call(rbind, taxa_list)

    generos <- unique(valor_texto(taxa_df$Genus))
    generos <- generos[!is.na(generos)]

    especies <- unique(valor_texto(taxa_df$Species))
    especies <- especies[!is.na(especies)]

    primeira <- best[1L, , drop = FALSE]
    primeira$Camada <- camada
    primeira$Genus_blast <- if (length(generos) == 1L) generos else NA_character_
    primeira$Species_blast <- if (
      length(especies) == 1L && length(generos) == 1L &&
      startsWith(especies, paste0(generos, " "))
    ) especies else NA_character_
    primeira$N_hits_melhor_nota <- nrow(best)
    primeira$N_generos_melhores <- length(generos)
    primeira$N_especies_melhores <- length(especies)
    primeira$Ambiguo_Genus <- length(generos) > 1L
    primeira$Ambiguo_Species <- length(especies) > 1L
    primeira$Generos_melhores <- paste(generos, collapse = ";")
    primeira$Especies_melhores <- paste(especies, collapse = ";")
    primeira
  })

  out <- do.call(rbind, linhas)
  rownames(out) <- NULL
  out
}

criar_matriz_taxa <- function(hits, sequencias) {
  mat <- matrix(
    NA_character_,
    nrow = length(sequencias),
    ncol = 2L,
    dimnames = list(sequencias, c("Genus", "Species"))
  )

  if (nrow(hits) == 0L) return(mat)

  idx <- match(hits$ASV_seq, sequencias)
  ok <- !is.na(idx)

  if (any(ok)) {
    mat[idx[ok], "Genus"] <- hits$Genus_blast[ok]
    mat[idx[ok], "Species"] <- hits$Species_blast[ok]
  }

  mat
}

###############################################################################
# 2. VALIDACAO E FASTA CANONICO
###############################################################################

cat("=============================================================\n")
cat("BLASTN / NCBI 16S — v", VERSAO, "\n", sep = "")
cat("Data:", DATA_EXECUCAO, "\n")
cat("Saida:", blast_root, "\n")
cat("=============================================================\n\n")

validar_arquivo(arq_seqtab, "seqtab_global_nochim.rds")
validar_arquivo(arq_asvmap, "ASV_sequences.tsv")

seqtab <- tryCatch(
  as.matrix(readRDS(arq_seqtab)),
  error = function(e) {
    abort(
      "Falha ao ler seqtab_global_nochim.rds: %s",
      conditionMessage(e)
    )
  }
)
if (nrow(seqtab) != 10L) {
  abort("seqtab_global_nochim.rds deve conter 10 amostras; encontrou %d.", nrow(seqtab))
}
if (ncol(seqtab) == 0L) {
  abort("seqtab_global_nochim.rds nao possui ASVs.")
}
if (!is.numeric(seqtab)) {
  abort("seqtab_global_nochim.rds deve conter contagens numericas.")
}
if (is.null(rownames(seqtab)) ||
    anyNA(rownames(seqtab)) ||
    any(rownames(seqtab) == "") ||
    anyDuplicated(rownames(seqtab)) > 0L) {
  abort("seqtab_global_nochim.rds deve possuir SampleID unicos e nao vazios.")
}
if (anyNA(seqtab) || any(!is.finite(seqtab)) || any(seqtab < 0) ||
    any(rowSums(seqtab) == 0) || any(colSums(seqtab) == 0)) {
  abort(
    "seqtab_global_nochim.rds contem contagem invalida, amostra sem reads ou ASV zerada."
  )
}
if (max(abs(seqtab - round(seqtab))) > 1e-8) {
  abort("seqtab_global_nochim.rds deve conter contagens inteiras.")
}
asv_map <- tryCatch(
  read.delim(
  arq_asvmap,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
  ),
  error = function(e) {
    abort("Falha ao ler ASV_sequences.tsv: %s", conditionMessage(e))
  }
)

if (!all(c("ASV_ID", "Sequence", "Origem") %in% colnames(asv_map))) {
  abort("ASV_sequences.tsv deve conter ASV_ID, Sequence e Origem.")
}
asv_map$ASV_ID <- trimws(as.character(asv_map$ASV_ID))
asv_map$Sequence <- toupper(trimws(as.character(asv_map$Sequence)))
asv_map$Origem <- trimws(as.character(asv_map$Origem))
if (anyNA(asv_map$ASV_ID) || any(asv_map$ASV_ID == "")) {
  abort("ASV_sequences.tsv contem ASV_ID ausente ou vazio.")
}
if (anyNA(asv_map$Sequence) || any(asv_map$Sequence == "")) {
  abort("ASV_sequences.tsv contem Sequence ausente ou vazia.")
}
if (anyNA(asv_map$Origem) || any(asv_map$Origem == "")) {
  abort("ASV_sequences.tsv contem Origem ausente ou vazia.")
}
if (!all(grepl("^[ACGTRYSWKMBDHVN]+$", asv_map$Sequence))) {
  abort("ASV_sequences.tsv contem sequencia fora do alfabeto IUPAC.")
}

if (anyDuplicated(asv_map$ASV_ID) > 0L) abort("ASV_ID duplicado no mapa.")
if (anyDuplicated(asv_map$Sequence) > 0L) abort("Sequence duplicada no mapa.")
if (is.null(colnames(seqtab)) || any(colnames(seqtab) == "")) {
  abort("seqtab sem sequencias nas colunas.")
}
if (anyDuplicated(colnames(seqtab)) > 0L) abort("seqtab contem ASVs duplicadas.")

idx_map <- match(colnames(seqtab), asv_map$Sequence)
if (anyNA(idx_map)) {
  abort(
    "%d ASVs do seqtab estao ausentes em ASV_sequences.tsv.",
    sum(is.na(idx_map))
  )
}

map_ord <- asv_map[idx_map, c("ASV_ID", "Sequence", "Origem"), drop = FALSE]
if (!identical(map_ord$Sequence, colnames(seqtab))) {
  abort("Falha ao ordenar ASV_sequences.tsv segundo o seqtab.")
}

asv_dna <- Biostrings::DNAStringSet(map_ord$Sequence)
names(asv_dna) <- map_ord$ASV_ID

Biostrings::writeXStringSet(
  asv_dna,
  filepath = arq_query,
  format = "fasta",
  width = 80L
)
validar_arquivo(arq_query, "FASTA canonico do BLAST")

fasta_check <- Biostrings::readDNAStringSet(arq_query)
if (!identical(names(fasta_check), map_ord$ASV_ID)) {
  abort("Headers do FASTA nao preservaram ASV_ID.")
}
if (!identical(unname(as.character(fasta_check)), map_ord$Sequence)) {
  abort("Sequencias do FASTA divergem do seqtab.")
}

header2seq <- setNames(map_ord$Sequence, map_ord$ASV_ID)
seq2id <- setNames(map_ord$ASV_ID, map_ord$Sequence)
n_asvs <- nrow(map_ord)

log_msg(sprintf("FASTA criado: %d ASVs.", n_asvs), "OK")

###############################################################################
# 3. VALIDACAO DO BLAST E DO BANCO
###############################################################################

blastn_check <- tryCatch(
  system2("blastn", args = "-version", stdout = TRUE, stderr = TRUE),
  error = function(e) character()
)
if (length(blastn_check) == 0L) abort("blastn nao encontrado no PATH.")
log_msg(blastn_check[1L], "OK")

db_files <- list.files(
  dirname(blast_db_path),
  pattern = paste0("^", basename(blast_db_path), "\\."),
  full.names = TRUE
)
if (length(db_files) == 0L) {
  abort("Banco BLAST nao encontrado com prefixo: %s", blast_db_path)
}
log_msg(sprintf("Banco NCBI 16S: %d arquivo(s) do indice.", length(db_files)), "OK")

# Facilita a recuperacao de sscinames quando taxdb.btd/taxdb.bti estao na
# mesma pasta do banco. O fallback por Subject.Title continua obrigatorio.
blast_db_dir <- dirname(blast_db_path)
blastdb_anterior <- Sys.getenv("BLASTDB", unset = "")
blastdb_partes <- c(blast_db_dir, blastdb_anterior)
blastdb_partes <- blastdb_partes[nzchar(blastdb_partes)]
Sys.setenv(BLASTDB = paste(unique(blastdb_partes), collapse = .Platform$path.sep))

taxdb_ok <- all(file.exists(file.path(blast_db_dir, c("taxdb.btd", "taxdb.bti"))))
if (taxdb_ok) {
  log_msg("taxdb.btd/taxdb.bti encontrados e BLASTDB configurado.", "OK")
} else {
  log_msg(
    "taxdb.btd/taxdb.bti nao encontrados; nomes cientificos dependerao do Subject.Title.",
    "WARN"
  )
}

###############################################################################
# 4. EXECUCAO BLASTN OU REPROCESSAMENTO DO TSV EXISTENTE
###############################################################################

blast_out_tsv <- file.path(blast_tab, "blast_resultados_brutos.tsv")
blast_out_tmp <- paste0(blast_out_tsv, ".tmp")
blast_stdout  <- file.path(blast_log, "blastn_stdout.log")
blast_stderr  <- file.path(blast_log, "blastn_stderr.log")
blast_cmd_log <- file.path(blast_log, "blastn_comando.txt")

if (file.exists(blast_out_tmp)) unlink(blast_out_tmp)

outfmt_blast <- paste(
  "6",
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "qcovhsp",
  "qcovs",
  "sscinames",
  "stitle"
)

blast_args <- c(
  "-query", arq_query,
  "-db", blast_db_path,
  "-perc_identity", as.character(PERC_IDENT_COLETA),
  "-qcov_hsp_perc", as.character(QCOV_COLETA),
  "-max_target_seqs", as.character(MAX_TARGET_SEQS),
  "-max_hsps", "1",
  "-num_threads", as.character(NUM_THREADS),
  "-outfmt", shQuote(outfmt_blast),
  "-out", blast_out_tmp
)

usar_blast_existente <- isTRUE(REUTILIZAR_BLAST_BRUTO) &&
  file.exists(blast_out_tsv) &&
  file.size(blast_out_tsv) > 0L

if (usar_blast_existente) {

  tempo_blast <- 0
  writeLines(
    c(
      paste("Data:", DATA_EXECUCAO),
      "Modo: reprocessamento do TSV bruto existente",
      paste("Arquivo:", blast_out_tsv),
      paste("Comando original esperado: blastn", paste(blast_args, collapse = " "))
    ),
    blast_cmd_log
  )

  log_msg(
    paste0(
      "Reutilizando blast_resultados_brutos.tsv existente; ",
      "a busca blastn nao sera repetida."
    ),
    "REUSE"
  )

} else {

  writeLines(
    c(
      paste("Data:", DATA_EXECUCAO),
      "Modo: nova execucao blastn",
      paste("blastn", paste(blast_args, collapse = " "))
    ),
    blast_cmd_log
  )

  log_msg(
    sprintf(
      "Executando blastn: coleta >=%d%% identidade; >=%d%% cobertura; max_target_seqs=%d.",
      PERC_IDENT_COLETA,
      QCOV_COLETA,
      MAX_TARGET_SEQS
    ),
    "RUN"
  )

  t0 <- proc.time()

  status_blast <- system2(
    "blastn",
    args = blast_args,
    stdout = blast_stdout,
    stderr = blast_stderr
  )

  tempo_blast <- unname((proc.time() - t0)["elapsed"])

  if (length(status_blast) != 1L || is.na(status_blast) || status_blast != 0L) {
    if (file.exists(blast_out_tmp)) unlink(blast_out_tmp)
    abort(
      "blastn falhou com status %s. Consulte: %s",
      as.character(status_blast),
      blast_stderr
    )
  }

  if (!file.exists(blast_out_tmp)) {
    abort("blastn nao criou o arquivo temporario esperado.")
  }

  if (file.exists(blast_out_tsv)) unlink(blast_out_tsv)
  if (!file.rename(blast_out_tmp, blast_out_tsv)) {
    abort("Falha ao finalizar arquivo BLAST: %s", blast_out_tsv)
  }

  log_msg(sprintf("blastn concluido em %.1f min.", tempo_blast / 60), "OK")
}

###############################################################################
# 5. LEITURA E VALIDACAO DOS HITS
###############################################################################

blast_cols <- c(
  "QueryID", "SubjectID", "Perc.Ident", "Alignment.Length",
  "Mismatches", "Gap.Openings", "Q.Start", "Q.End",
  "S.Start", "S.End", "E.Value", "Bits", "Qcov.HSP",
  "Qcov.Query", "Subject.Scientific.Name", "Subject.Title"
)

if (file.size(blast_out_tsv) == 0L) {
  blast_raw <- as.data.frame(
    setNames(replicate(length(blast_cols), character(), simplify = FALSE), blast_cols),
    stringsAsFactors = FALSE
  )
  log_msg("BLAST nao retornou hits; matrizes serao salvas integralmente como NA.", "WARN")
} else {
  blast_raw <- read.delim(
    blast_out_tsv,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = FALSE,
    stringsAsFactors = FALSE,
    col.names = blast_cols,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A")
  )
}

if (ncol(blast_raw) != length(blast_cols)) {
  abort(
    "Resultado BLAST possui %d colunas; eram esperadas %d.",
    ncol(blast_raw),
    length(blast_cols)
  )
}

if (nrow(blast_raw) > 0L) {
  blast_raw$Sci.Name <- valor_texto(blast_raw$Subject.Scientific.Name)

  idx_vazio <- is.na(blast_raw$Sci.Name)
  if (any(idx_vazio)) {
    fallback <- trimws(blast_raw$Subject.Title[idx_vazio])
    fallback <- sub(
      "^[A-Z]{1,6}_[A-Z0-9]+(?:\\.[0-9]+)?\\s+",
      "",
      fallback,
      perl = TRUE
    )
    blast_raw$Sci.Name[idx_vazio] <- fallback
  }

  num_cols <- c(
    "Perc.Ident", "Alignment.Length", "Mismatches", "Gap.Openings",
    "Q.Start", "Q.End", "S.Start", "S.End", "E.Value", "Bits",
    "Qcov.HSP", "Qcov.Query"
  )
  blast_raw[num_cols] <- lapply(
    blast_raw[num_cols],
    function(x) suppressWarnings(as.numeric(x))
  )
  numericos_invalidos <- vapply(
    blast_raw[num_cols],
    function(x) anyNA(x) || any(!is.finite(x)),
    logical(1)
  )
  if (any(numericos_invalidos)) {
    abort(
      "Resultado BLAST possui valor numerico invalido em: %s",
      paste(names(numericos_invalidos)[numericos_invalidos], collapse = ", ")
    )
  }
  if (any(blast_raw$Perc.Ident < 0 | blast_raw$Perc.Ident > 100) ||
      any(blast_raw$Alignment.Length <= 0) ||
      any(blast_raw$Mismatches < 0) ||
      any(blast_raw$Gap.Openings < 0) ||
      any(blast_raw$Q.Start < 1) ||
      any(blast_raw$Q.End < blast_raw$Q.Start)) {
    abort("Resultado BLAST possui metrica ou coordenada fora do dominio esperado.")
  }

  blast_raw$Perc.Query.Coverage <- blast_raw$Qcov.Query
  blast_raw$Perc.HSP.Coverage <- blast_raw$Qcov.HSP
  blast_raw$ASV_seq <- unname(header2seq[blast_raw$QueryID])
  blast_raw$ASV_ID <- unname(seq2id[blast_raw$ASV_seq])
  blast_raw$Query.Length <- nchar(blast_raw$ASV_seq)

  if (anyNA(blast_raw$ASV_seq)) {
    abort("Ha QueryID do BLAST sem correspondencia no FASTA canonico.")
  }

  query_ids_observados <- unique(as.character(blast_raw$QueryID))
  query_ids_invalidos <- setdiff(query_ids_observados, map_ord$ASV_ID)
  if (length(query_ids_invalidos) > 0L) {
    abort(
      "TSV bruto incompativel com o FASTA atual; QueryID(s) desconhecido(s): %s",
      paste(head(query_ids_invalidos, 10L), collapse = ", ")
    )
  }
}

salvar_csv(blast_raw, file.path(blast_tab, "blast_resultados_brutos.csv"))

###############################################################################
# 6. CAMADAS DE EVIDENCIA
###############################################################################

if (nrow(blast_raw) > 0L) {
  filtro_exact <- (
    blast_raw$Perc.Ident == PERC_IDENT_SPECIES &
    blast_raw$Perc.Query.Coverage == QCOV_SPECIES &
    blast_raw$Mismatches == 0 &
    blast_raw$Gap.Openings == 0 &
    blast_raw$Alignment.Length == blast_raw$Query.Length &
    blast_raw$Q.Start == 1 &
    blast_raw$Q.End == blast_raw$Query.Length
  )

  filtro_genus <- (
    blast_raw$Perc.Ident >= PERC_IDENT_GENUS &
    blast_raw$Perc.Query.Coverage >= QCOV_GENUS &
    blast_raw$Alignment.Length == blast_raw$Query.Length &
    blast_raw$Q.Start == 1 &
    blast_raw$Q.End == blast_raw$Query.Length
  )

  hits_100 <- resumir_hits(
    blast_raw[!is.na(filtro_exact) & filtro_exact, , drop = FALSE],
    "species_exact_full"
  )
  hits_97 <- resumir_hits(
    blast_raw[!is.na(filtro_genus) & filtro_genus, , drop = FALSE],
    "genus_candidate"
  )
} else {
  hits_100 <- data.frame()
  hits_97 <- data.frame()
}

if (nrow(hits_100) > 0L) {
  hits_100$ASV_seq <- unname(header2seq[hits_100$QueryID])
  hits_100$ASV_ID <- unname(seq2id[hits_100$ASV_seq])
}
if (nrow(hits_97) > 0L) {
  hits_97$ASV_seq <- unname(header2seq[hits_97$QueryID])
  hits_97$ASV_ID <- unname(seq2id[hits_97$ASV_seq])
  # Species dessa camada nao entra no consenso principal.
  hits_97$Species_blast_candidata <- hits_97$Species_blast
  hits_97$Species_blast <- NA_character_
}

taxa_blast100 <- criar_matriz_taxa(hits_100, map_ord$Sequence)
taxa_blast97 <- criar_matriz_taxa(hits_97, map_ord$Sequence)

salvar_rds(taxa_blast100, file.path(blast_rds, "taxa_blast100.rds"))
salvar_rds(taxa_blast97, file.path(blast_rds, "taxa_blast97.rds"))

salvar_csv(
  cbind(
    ASV_ID = map_ord$ASV_ID,
    ASV_seq = map_ord$Sequence,
    as.data.frame(taxa_blast100, stringsAsFactors = FALSE)
  ),
  file.path(blast_tab, "taxa_blast100.csv")
)
salvar_csv(
  cbind(
    ASV_ID = map_ord$ASV_ID,
    ASV_seq = map_ord$Sequence,
    as.data.frame(taxa_blast97, stringsAsFactors = FALSE)
  ),
  file.path(blast_tab, "taxa_blast97.csv")
)

salvar_csv(hits_100, file.path(blast_tab, "blast_hits_100pct_integral.csv"))
salvar_csv(hits_97, file.path(blast_tab, "blast_hits_97pct_candidatos.csv"))

###############################################################################
# 7. EVIDENCIA POR ASV E SATURACAO DA BUSCA
###############################################################################

n_hits_raw <- if (nrow(blast_raw) > 0L) {
  table(blast_raw$QueryID)
} else {
  integer()
}

evidencias <- data.frame(
  ASV_ID = map_ord$ASV_ID,
  ASV_seq = map_ord$Sequence,
  Origem = map_ord$Origem,
  Query_length = nchar(map_ord$Sequence),
  N_hits_retornados = as.integer(n_hits_raw[map_ord$ASV_ID]),
  Atingiu_MAX_TARGET_SEQS = FALSE,
  Genus_blast_exato = NA_character_,
  Species_blast_exata = NA_character_,
  N_hits_exatos_melhor_nota = 0L,
  Ambiguo_Genus_exato = FALSE,
  Ambiguo_Species_exata = FALSE,
  Genus_blast_candidato = NA_character_,
  Melhor_identidade_genus = NA_real_,
  Melhor_cobertura_genus = NA_real_,
  N_hits_genus_melhor_nota = 0L,
  Ambiguo_Genus_candidato = FALSE,
  stringsAsFactors = FALSE
)

evidencias$N_hits_retornados[is.na(evidencias$N_hits_retornados)] <- 0L
evidencias$Atingiu_MAX_TARGET_SEQS <-
  evidencias$N_hits_retornados >= MAX_TARGET_SEQS

if (nrow(hits_100) > 0L) {
  idx <- match(hits_100$ASV_ID, evidencias$ASV_ID)
  ok <- !is.na(idx)
  evidencias$Genus_blast_exato[idx[ok]] <- hits_100$Genus_blast[ok]
  evidencias$Species_blast_exata[idx[ok]] <- hits_100$Species_blast[ok]
  evidencias$N_hits_exatos_melhor_nota[idx[ok]] <-
    hits_100$N_hits_melhor_nota[ok]
  evidencias$Ambiguo_Genus_exato[idx[ok]] <- hits_100$Ambiguo_Genus[ok]
  evidencias$Ambiguo_Species_exata[idx[ok]] <- hits_100$Ambiguo_Species[ok]
}

if (nrow(hits_97) > 0L) {
  idx <- match(hits_97$ASV_ID, evidencias$ASV_ID)
  ok <- !is.na(idx)
  evidencias$Genus_blast_candidato[idx[ok]] <- hits_97$Genus_blast[ok]
  evidencias$Melhor_identidade_genus[idx[ok]] <- hits_97$Perc.Ident[ok]
  evidencias$Melhor_cobertura_genus[idx[ok]] <-
    hits_97$Perc.Query.Coverage[ok]
  evidencias$N_hits_genus_melhor_nota[idx[ok]] <-
    hits_97$N_hits_melhor_nota[ok]
  evidencias$Ambiguo_Genus_candidato[idx[ok]] <-
    hits_97$Ambiguo_Genus[ok]
}

salvar_rds(
  evidencias,
  file.path(blast_rds, "blast_evidencias_por_asv.rds")
)
salvar_csv(
  evidencias,
  file.path(blast_tab, "blast_evidencias_por_asv.csv")
)

n_saturadas <- sum(evidencias$Atingiu_MAX_TARGET_SEQS)
if (n_saturadas > 0L) {
  salvar_csv(
    evidencias[evidencias$Atingiu_MAX_TARGET_SEQS, , drop = FALSE],
    file.path(blast_tab, "blast_consultas_atingiram_limite_hits.csv")
  )
  log_msg(
    sprintf(
      "%d ASV(s) atingiram max_target_seqs=%d; ambiguidade pode estar subestimada.",
      n_saturadas,
      MAX_TARGET_SEQS
    ),
    "WARN"
  )
}

###############################################################################
# 8. BAIXA IDENTIDADE — SOMENTE AUDITORIA
###############################################################################

if (nrow(blast_raw) > 0L) {
  baixo <- blast_raw[
    blast_raw$Perc.Query.Coverage >= QCOV_BAIXA &
    blast_raw$Q.Start == 1L &
    blast_raw$Q.End == blast_raw$Query.Length &
    blast_raw$Alignment.Length == blast_raw$Query.Length &
    blast_raw$Perc.Ident >= PERC_IDENT_COLETA &
    blast_raw$Perc.Ident < PERC_IDENT_BAIXA,
    ,
    drop = FALSE
  ]

  if (nrow(baixo) > 0L) {
    baixo <- do.call(
      rbind,
      lapply(split(baixo, baixo$QueryID), selecionar_empates_melhores)
    )
    baixo$Flag <- "cobertura_integral_identidade_abaixo_limiar_auditoria"
  }
} else {
  baixo <- data.frame()
}

salvar_csv(
  baixo,
  file.path(blast_tab, "taxa_blast_baixa_identidade.csv")
)

###############################################################################
# 9. METADADOS E CHECKPOINT FINAL
###############################################################################

metadata_blast <- data.frame(
  Script = "04_rblast.R",
  Versao = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Banco = "NCBI 16S_ribosomal_RNA",
  Metodo = "blastn_tabular",
  Pasta_saida = blast_root,
  Query_fasta = arq_query,
  Perc_ident_coleta = PERC_IDENT_COLETA,
  Qcov_coleta = QCOV_COLETA,
  Perc_ident_genus = PERC_IDENT_GENUS,
  Qcov_genus = QCOV_GENUS,
  Perc_ident_species = PERC_IDENT_SPECIES,
  Perc_ident_exact = PERC_IDENT_SPECIES,
  Qcov_species = QCOV_SPECIES,
  Species_mismatches = 0L,
  Species_gap_openings = 0L,
  Species_alignment_full_length = TRUE,
  Max_target_seqs = MAX_TARGET_SEQS,
  Num_threads = NUM_THREADS,
  Blast_bruto_reutilizado = usar_blast_existente,
  ASVs_entrada = n_asvs,
  ASVs_genus_candidato = sum(!is.na(taxa_blast97[, "Genus"])),
  ASVs_species_exata = sum(!is.na(taxa_blast100[, "Species"])),
  ASVs_busca_saturada = n_saturadas,
  Tempo_blast_s = round(tempo_blast, 1),
  Rownames_tipo = "sequencia_nucleotidica",
  Criterio_genus = "triagem_97pct; alinhamento qstart=1,qend=qlen e length=qlen; decisao_final_no_script05",
  Criterio_species = "100pct_identidade_e_cobertura; zero_mismatch_gap; full_length; nao_ambigua",
  stringsAsFactors = FALSE
)

salvar_csv(
  metadata_blast,
  file.path(blast_tab, "metadata_blast.csv")
)

salvar_rds(
  list(
    metadata = metadata_blast,
    taxa_blast97 = taxa_blast97,
    taxa_blast100 = taxa_blast100,
    evidencias = evidencias,
    hits_97 = hits_97,
    hits_100 = hits_100
  ),
  file.path(blast_chk, "checkpoint_04_blast_concluido.rds")
)

cat("\n=============================================================\n")
cat("BLASTN / NCBI 16S — CONCLUIDO\n")
cat(sprintf("ASVs: %d\n", n_asvs))
cat(sprintf(
  "Genus candidato (>=%.1f%%; cobertura %.1f%%): %d\n",
  PERC_IDENT_GENUS,
  QCOV_GENUS,
  sum(!is.na(taxa_blast97[, "Genus"]))
))
cat(sprintf(
  "Species exata integral: %d\n",
  sum(!is.na(taxa_blast100[, "Species"]))
))
cat("Saidas:", blast_root, "\n")
cat("Proximo passo: Script 05 revisado.\n")
cat("=============================================================\n\n")

log_msg("Script 04 finalizado.", "FINAL")
})
