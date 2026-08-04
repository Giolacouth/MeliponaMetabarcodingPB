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
run_pipeline_script("00_conferencia_primers.R", "primers", function(ctx) {
###############################################################################
# 00_conferencia_primers.R
#

###############################################################################

suppressPackageStartupMessages({
  library(ShortRead)
  library(Biostrings)
})

###############################################################################
# 1. CONFIGURAÇÃO
###############################################################################

base_path <- ctx$base_path
raw_path <- ctx$raw_path
metadata_path <- file.path(raw_path, "metadados.tsv")
pipeline_version <- ctx$version
trimmed_path <- Sys.getenv(
  "PRIMER_TRIMMED_PATH",
  unset = file.path(ctx$layout$stages$dada2$root, "fastq_trimmed")
)
output_path <- ctx$stage$root
# Primers V3–V4
FWD <- DNAString()       
REV <- DNAString()   

# Região inicial examinada para confirmar presença na extremidade 5'
window_5p <- 40L

# Busca tolerante na janela 5': exata, 1 mismatch e 2 mismatches
mismatch_5p <- c(0L, 1L, 2L)

# Em toda a read, a busca exata é suficiente para detectar ocorrências internas
mismatch_whole_read <- 0L

# Não permitir inserções/deleções durante a busca
allow_indels <- FALSE

# Por padrão, analisar todas as reads.
# Um valor inteiro em max_reads limita a quantidade analisada.
max_reads <- Inf

# Quando max_reads for finito:
#   "random" = seleção aleatória reprodutível dos mesmos pares R1/R2
#   "first"  = primeiras reads
sampling_method <- "random"
random_seed <- 1234L

# Limiares apenas para geração de alertas no relatório
minimum_expected_primer_percent_raw <- 90
maximum_residual_primer_percent_trimmed <- 1
maximum_N_percent_5p <- 1

# Geração de gráficos PNG
generate_plots <- TRUE

# Se TRUE, uma amostra sem correspondência no metadado é considerada erro.
strict_metadata <- TRUE

###############################################################################
# 2. PREPARAÇÃO
###############################################################################


# showWarnings = FALSE), esta função não oculta a causa da falha.
ensure_directory <- function(path, label = "diretório") {
  path <- path.expand(path)


  if (file.exists(path) && !dir.exists(path)) {
    stop(
      sprintf(
        "%s não pode ser criado porque já existe um arquivo nesse caminho: %s",
        label,
        path
      ),
      call. = FALSE
    )
  }

  if (!dir.exists(path)) {
    created <- dir.create(
      path,
      recursive = TRUE,
      showWarnings = TRUE
    )

    if (!isTRUE(created) && !dir.exists(path)) {
      parent <- dirname(path)

      stop(
        paste0(
          "Falha ao criar ", label, ": ", path, "\n",
          "Diretório pai: ", parent, "\n",
          "Pai existe: ", dir.exists(parent), "\n",
          "Permissão de escrita no pai: ",
          if (dir.exists(parent)) file.access(parent, 2) == 0 else FALSE
        ),
        call. = FALSE
      )
    }
  }

  if (!dir.exists(path)) {
    stop(
      sprintf("%s não existe após a tentativa de criação: %s", label, path),
      call. = FALSE
    )
  }

  if (file.access(path, 2) != 0) {
    stop(
      sprintf("Sem permissão de escrita em %s: %s", label, path),
      call. = FALSE
    )
  }

  normalizePath(path, mustWork = TRUE)
}

output_path <- ensure_directory(
  output_path,
  "diretório principal de saída"
)

table_path <- ensure_directory(
  file.path(output_path, "tabelas"),
  "diretório de tabelas"
)

plot_path <- ensure_directory(
  file.path(output_path, "graficos"),
  "diretório de gráficos"
)

log_path <- ensure_directory(
  file.path(output_path, "logs"),
  "diretório de logs"
)

cat("\nDiretórios de saída validados:\n")
cat("  output:", output_path, "\n")
cat("  tabelas:", table_path, "\n")
cat("  gráficos:", plot_path, "\n")
cat("  logs:", log_path, "\n")

execution_start <- Sys.time()
script_version <- "2.0.0"

primers <- list(
  `341F` = FWD,
  `341F_RC` = reverseComplement(FWD),
  `805R` = REV,
  `805R_RC` = reverseComplement(REV)
)

expected_primer <- c(
  R1 = "341F",
  R2 = "805R"
)

search_plan <- list(
  five_prime = mismatch_5p,
  whole_read = mismatch_whole_read
)

###############################################################################
# 3. FUNÇÕES UTILITÁRIAS
###############################################################################

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

pct <- function(numerator, denominator, digits = 3L) {
  if (is.na(denominator) || denominator == 0L) {
    return(NA_real_)
  }
  round(100 * numerator / denominator, digits)
}

safe_min <- function(x) {
  if (length(x) == 0L) NA_integer_ else min(x)
}

safe_max <- function(x) {
  if (length(x) == 0L) NA_integer_ else max(x)
}

safe_median <- function(x) {
  if (length(x) == 0L) NA_real_ else median(x)
}

mode_integer <- function(x) {
  if (length(x) == 0L) {
    return(NA_integer_)
  }
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])
}

normalize_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

normalize_read_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\s.*$", "", x)
  x <- sub("/[12]$", "", x)
  x
}

strip_fastq_extension <- function(x) {
  sub("\\.(fastq|fq)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}

parse_fastq_filename <- function(file) {
  stem <- strip_fastq_extension(file)

  patterns <- c(
    "^(.*)_R([12])_001(?:_.*)?$",
    "^(.*)_R([12])(?:_.*)?$",
    "^(.*)_([12])(?:_.*)?$",
    "^(.*)\\.R([12])(?:[._].*)?$"
  )

  for (pattern in patterns) {
    hit <- regexec(pattern, stem, perl = TRUE)
    groups <- regmatches(stem, hit)[[1L]]

    if (length(groups) == 3L) {
      return(data.frame(
        RawFastqID = groups[2L],
        Direction = paste0("R", groups[3L]),
        File = normalizePath(file, mustWork = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  abort(
    "Não foi possível identificar R1/R2 no arquivo: %s. ",
    basename(file)
  )
}

list_fastq_files <- function(path) {
  if (!dir.exists(path)) {
    return(character())
  }

  list.files(
    path,
    pattern = "\\.(fastq|fq)(\\.gz)?$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
}

discover_fastq_pairs <- function(stage, path) {
  files <- list_fastq_files(path)

  if (length(files) == 0L) {
    return(data.frame())
  }

  parsed <- do.call(
    rbind,
    lapply(files, parse_fastq_filename)
  )

  key <- paste(parsed$RawFastqID, parsed$Direction, sep = "::")
  if (anyDuplicated(key)) {
    duplicated_files <- parsed[key %in% key[duplicated(key)], , drop = FALSE]
    abort(
      "Existem arquivos duplicados para uma mesma amostra/direção em %s:\n%s",
      path,
      paste(duplicated_files$File, collapse = "\n")
    )
  }

  r1 <- parsed[parsed$Direction == "R1", c("RawFastqID", "File"), drop = FALSE]
  r2 <- parsed[parsed$Direction == "R2", c("RawFastqID", "File"), drop = FALSE]
  names(r1)[2L] <- "R1_File"
  names(r2)[2L] <- "R2_File"

  pairs <- merge(r1, r2, by = "RawFastqID", all = TRUE, sort = TRUE)

  missing_r1 <- is.na(pairs$R1_File)
  missing_r2 <- is.na(pairs$R2_File)

  if (any(missing_r1 | missing_r2)) {
    incomplete <- pairs[missing_r1 | missing_r2, , drop = FALSE]
    abort(
      "Foram encontrados pares FASTQ incompletos em %s:\n%s",
      path,
      paste(capture.output(print(incomplete, row.names = FALSE)), collapse = "\n")
    )
  }

  pairs$Stage <- stage
  pairs$InputPath <- normalizePath(path, mustWork = TRUE)

  pairs[, c("Stage", "InputPath", "RawFastqID", "R1_File", "R2_File")]
}

read_metadata <- function(path) {
  if (!file.exists(path)) {
    abort("Metadado não encontrado: %s", path)
  }

  metadata <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )

  # Aceita o nome alternativo SampleLabels, mas padroniza internamente como
  # SampleLabel para manter um unico contrato em todo o pipeline.
  if (!"SampleLabel" %in% names(metadata) && "SampleLabels" %in% names(metadata)) {
    names(metadata)[names(metadata) == "SampleLabels"] <- "SampleLabel"
  }

  required <- c("SampleID", "SampleLabel", "Run")
  missing_columns <- setdiff(required, names(metadata))

  if (length(missing_columns) > 0L) {
    abort(
      "O metadado não possui as colunas obrigatórias: %s",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (anyNA(metadata$SampleID) || any(metadata$SampleID == "")) {
    abort("Existem SampleID ausentes no metadado.")
  }

  if (anyDuplicated(metadata$SampleID)) {
    abort("Existem SampleID duplicados no metadado.")
  }

  metadata$SampleLabel <- trimws(as.character(metadata$SampleLabel))
  if (anyNA(metadata$SampleLabel) || any(metadata$SampleLabel == "")) {
    abort("Existem SampleLabel ausentes ou vazios no metadado.")
  }
  if (anyDuplicated(metadata$SampleLabel)) {
    abort("Existem SampleLabel duplicados no metadado; os rotulos dos graficos devem ser unicos.")
  }

  metadata
}

resolve_metadata_row <- function(raw_fastq_id, metadata) {
  target <- normalize_key(raw_fastq_id)

  candidate_columns <- c("SampleID")
  if ("RawSampleID" %in% names(metadata)) {
    candidate_columns <- c(candidate_columns, "RawSampleID")
  }

  hits <- integer()

  for (column in candidate_columns) {
    values <- normalize_key(metadata[[column]])
    hits <- union(hits, which(values == target))
  }

  # Permite sufixos técnicos, por exemplo "_trimmed" ou "_cutadapt",
  # desde que apenas uma amostra seja identificada.
  if (length(hits) == 0L) {
    for (column in candidate_columns) {
      values <- normalize_key(metadata[[column]])
      prefix_hit <- vapply(
        values,
        function(value) {
          startsWith(target, paste0(value, "_")) ||
            startsWith(value, paste0(target, "_"))
        },
        logical(1L)
      )
      hits <- union(hits, which(prefix_hit))
    }
  }

  if (length(hits) == 1L) {
    return(metadata[hits, , drop = FALSE])
  }

  if (length(hits) > 1L) {
    abort(
      "O identificador FASTQ '%s' corresponde a mais de uma linha do metadado.",
      raw_fastq_id
    )
  }

  if (strict_metadata) {
    abort(
      paste0(
        "Não foi possível associar o FASTQ '%s' ao metadado. ",
        "Confirme se o nome-base do FASTQ corresponde a SampleID ou RawSampleID."
      ),
      raw_fastq_id
    )
  }

  data.frame(
    SampleID = raw_fastq_id,
    SampleLabel = raw_fastq_id,
    Run = NA_character_,
    stringsAsFactors = FALSE
  )
}

select_pair_indices <- function(n_reads, raw_fastq_id, stage) {
  indices <- seq_len(n_reads)

  if (!is.finite(max_reads) || n_reads <= max_reads) {
    return(indices)
  }

  n_select <- as.integer(max_reads)

  if (sampling_method == "first") {
    return(indices[seq_len(n_select)])
  }

  if (sampling_method != "random") {
    abort("sampling_method deve ser 'random' ou 'first'.")
  }

  stage_offset <- sum(utf8ToInt(stage))
  sample_offset <- sum(utf8ToInt(raw_fastq_id))
  set.seed(random_seed + stage_offset + sample_offset)

  sort(sample(indices, size = n_select, replace = FALSE))
}

extract_search_region <- function(seqs, region_name) {
  if (region_name == "whole_read") {
    return(seqs)
  }

  if (region_name == "five_prime") {
    region_width <- pmin(width(seqs), window_5p)
    return(subseq(seqs, start = 1L, width = region_width))
  }

  abort("Região de busca desconhecida: %s", region_name)
}

file_md5 <- function(path) {
  unname(tools::md5sum(path))
}

file_size <- function(path) {
  unname(file.info(path)$size)
}

as_character_or_na <- function(x) {
  if (length(x) == 0L || is.null(x)) NA_character_ else as.character(x)
}

rbind_or_empty <- function(items) {
  items <- Filter(
    function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L,
    items
  )

  if (length(items) == 0L) {
    return(data.frame())
  }

  do.call(rbind, items)
}

write_tsv <- function(x, filename) {
  # Revalida a pasta imediatamente antes de cada gravação. Isso também
  # protege contra exclusão ou desmontagem acidental durante uma execução longa.
  current_table_path <- ensure_directory(
    table_path,
    "diretório de tabelas"
  )

  destination <- file.path(
    current_table_path,
    filename
  )

  write.table(
    x,
    file = destination,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA",
    fileEncoding = "UTF-8"
  )

  if (!file.exists(destination)) {
    stop(
      sprintf("A tabela não foi criada: %s", destination),
      call. = FALSE
    )
  }

  invisible(destination)
}

add_issue <- function(
    issues,
    severity,
    stage,
    sample_id,
    code,
    message
) {
  issues[[length(issues) + 1L]] <- data.frame(
    Severity = severity,
    Stage = stage,
    SampleID = sample_id,
    Code = code,
    Message = message,
    stringsAsFactors = FALSE
  )
  issues
}

###############################################################################
# 4. ANÁLISE DE UMA AMOSTRA PAREADA
###############################################################################

analyze_fastq_pair <- function(pair_row, metadata) {
  stage <- pair_row$Stage
  raw_fastq_id <- pair_row$RawFastqID
  r1_file <- pair_row$R1_File
  r2_file <- pair_row$R2_File

  meta <- resolve_metadata_row(raw_fastq_id, metadata)
  sample_id <- meta$SampleID[1L]
  sample_label <- meta$SampleLabel[1L]
  run_id <- meta$Run[1L]

  fq1 <- readFastq(r1_file)
  fq2 <- readFastq(r2_file)

  total_r1 <- length(fq1)
  total_r2 <- length(fq2)

  if (total_r1 != total_r2) {
    abort(
      "%s (%s): R1 possui %d reads e R2 possui %d reads.",
      sample_id, stage, total_r1, total_r2
    )
  }

  if (total_r1 == 0L) {
    abort("%s (%s): os arquivos FASTQ estão vazios.", sample_id, stage)
  }

  ids_r1 <- normalize_read_id(id(fq1))
  ids_r2 <- normalize_read_id(id(fq2))
  ids_aligned <- identical(ids_r1, ids_r2)

  if (!ids_aligned) {
    first_mismatch <- which(ids_r1 != ids_r2)[1L]
    abort(
      paste0(
        "%s (%s): identificadores R1/R2 não estão alinhados. ",
        "Primeira divergência no par %d."
      ),
      sample_id, stage, first_mismatch
    )
  }

  selected <- select_pair_indices(total_r1, raw_fastq_id, stage)
  fq1_selected <- fq1[selected]
  fq2_selected <- fq2[selected]

  seqs_by_direction <- list(
    R1 = sread(fq1_selected),
    R2 = sread(fq2_selected)
  )

  inventory <- data.frame(
    Stage = stage,
    SampleID = sample_id,
    SampleLabel = sample_label,
    Run = run_id,
    RawFastqID = raw_fastq_id,
    R1_File = basename(r1_file),
    R2_File = basename(r2_file),
    R1_Path = r1_file,
    R2_Path = r2_file,
    R1_FileSizeBytes = file_size(r1_file),
    R2_FileSizeBytes = file_size(r2_file),
    R1_MD5 = file_md5(r1_file),
    R2_MD5 = file_md5(r2_file),
    R1_TotalReads = total_r1,
    R2_TotalReads = total_r2,
    TotalReadPairs = total_r1,
    ReadIDsAligned = ids_aligned,
    ReadPairsAnalyzed = length(selected),
    SamplingMethod = if (length(selected) == total_r1) "all" else sampling_method,
    R1_MinLength = min(width(seqs_by_direction$R1)),
    R1_MedianLength = median(width(seqs_by_direction$R1)),
    R1_MaxLength = max(width(seqs_by_direction$R1)),
    R2_MinLength = min(width(seqs_by_direction$R2)),
    R2_MedianLength = median(width(seqs_by_direction$R2)),
    R2_MaxLength = max(width(seqs_by_direction$R2)),
    stringsAsFactors = FALSE
  )

  primer_rows <- list()
  position_rows <- list()
  presence_cache <- list()
  N_cache <- list()

  for (direction in names(seqs_by_direction)) {
    seqs <- seqs_by_direction[[direction]]

    for (region_name in names(search_plan)) {
      region_seqs <- extract_search_region(seqs, region_name)

      N_count_per_read <- vcountPattern(
        DNAString("N"),
        region_seqs,
        fixed = TRUE
      )
      N_cache[[paste(direction, region_name, sep = "|")]] <- N_count_per_read > 0L

      for (max_mismatch in search_plan[[region_name]]) {
        for (primer_name in names(primers)) {
          primer <- primers[[primer_name]]

          matches <- vmatchPattern(
            pattern = primer,
            subject = region_seqs,
            max.mismatch = max_mismatch,
            min.mismatch = 0L,
            with.indels = allow_indels,
            fixed = FALSE
          )

          counts <- elementNROWS(matches)
          present <- counts > 0L
          multiple <- counts > 1L
          starts <- unlist(start(matches), use.names = FALSE)

          cache_key <- paste(
            direction,
            region_name,
            max_mismatch,
            primer_name,
            sep = "|"
          )
          presence_cache[[cache_key]] <- present

          starts_at_1 <- if (length(starts) == 0L) {
            0L
          } else {
            sum(starts == 1L)
          }

          N_present <- N_cache[[paste(direction, region_name, sep = "|")]]

          primer_rows[[length(primer_rows) + 1L]] <- data.frame(
            Stage = stage,
            SampleID = sample_id,
            SampleLabel = sample_label,
            Run = run_id,
            RawFastqID = raw_fastq_id,
            Direction = direction,
            File = if (direction == "R1") basename(r1_file) else basename(r2_file),
            SearchRegion = region_name,
            SearchWindowNt = if (region_name == "five_prime") window_5p else NA_integer_,
            PrimerName = primer_name,
            PrimerSequence = as.character(primer),
            PrimerLength = length(primer),
            ExpectedForDirection = primer_name == expected_primer[[direction]],
            MaxMismatch = max_mismatch,
            AllowIndels = allow_indels,
            TotalReadsInFile = if (direction == "R1") total_r1 else total_r2,
            ReadsAnalyzed = length(region_seqs),
            ReadsWithPrimer = sum(present),
            PercentReadsWithPrimer = pct(sum(present), length(region_seqs)),
            TotalPrimerOccurrences = sum(counts),
            ReadsWithMultipleOccurrences = sum(multiple),
            PercentReadsWithMultipleOccurrences = pct(sum(multiple), length(region_seqs)),
            OccurrencesStartingAtPosition1 = starts_at_1,
            PercentOccurrencesStartingAtPosition1 = pct(starts_at_1, length(starts)),
            MostFrequentStartPosition = mode_integer(starts),
            MedianStartPosition = safe_median(starts),
            MinimumStartPosition = safe_min(starts),
            MaximumStartPosition = safe_max(starts),
            ReadsWithNInRegion = sum(N_present),
            PercentReadsWithNInRegion = pct(sum(N_present), length(region_seqs)),
            stringsAsFactors = FALSE
          )

          if (length(starts) > 0L) {
            position_table <- as.data.frame(table(starts), stringsAsFactors = FALSE)
            names(position_table) <- c("StartPosition", "Occurrences")
            position_table$StartPosition <- as.integer(as.character(position_table$StartPosition))

            position_rows[[length(position_rows) + 1L]] <- data.frame(
              Stage = stage,
              SampleID = sample_id,
              SampleLabel = sample_label,
              Run = run_id,
              Direction = direction,
              SearchRegion = region_name,
              PrimerName = primer_name,
              MaxMismatch = max_mismatch,
              StartPosition = position_table$StartPosition,
              Occurrences = position_table$Occurrences,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  get_presence <- function(direction, mismatch, primer_name) {
    key <- paste(
      direction,
      "five_prime",
      mismatch,
      primer_name,
      sep = "|"
    )
    presence_cache[[key]]
  }

  pair_rows <- list()

  for (max_mismatch in mismatch_5p) {
    r1_expected <- get_presence("R1", max_mismatch, "341F")
    r2_expected <- get_presence("R2", max_mismatch, "805R")

    r1_inverted <- get_presence("R1", max_mismatch, "805R")
    r2_inverted <- get_presence("R2", max_mismatch, "341F")

    r1_any_rc <- (
      get_presence("R1", max_mismatch, "341F_RC") |
        get_presence("R1", max_mismatch, "805R_RC")
    )
    r2_any_rc <- (
      get_presence("R2", max_mismatch, "341F_RC") |
        get_presence("R2", max_mismatch, "805R_RC")
    )

    both_expected <- r1_expected & r2_expected
    only_r1_expected <- r1_expected & !r2_expected
    only_r2_expected <- !r1_expected & r2_expected
    neither_expected <- !r1_expected & !r2_expected
    both_inverted_direct <- r1_inverted & r2_inverted
    any_unexpected <- (
      r1_inverted |
        r2_inverted |
        r1_any_rc |
        r2_any_rc
    )

    n_pairs <- length(r1_expected)

    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      Stage = stage,
      SampleID = sample_id,
      SampleLabel = sample_label,
      Run = run_id,
      RawFastqID = raw_fastq_id,
      SearchRegion = "five_prime",
      SearchWindowNt = window_5p,
      MaxMismatch = max_mismatch,
      ReadPairsAnalyzed = n_pairs,
      PairsWithBothExpected = sum(both_expected),
      PercentPairsWithBothExpected = pct(sum(both_expected), n_pairs),
      PairsWithOnlyR1Expected = sum(only_r1_expected),
      PercentPairsWithOnlyR1Expected = pct(sum(only_r1_expected), n_pairs),
      PairsWithOnlyR2Expected = sum(only_r2_expected),
      PercentPairsWithOnlyR2Expected = pct(sum(only_r2_expected), n_pairs),
      PairsWithNeitherExpected = sum(neither_expected),
      PercentPairsWithNeitherExpected = pct(sum(neither_expected), n_pairs),
      PairsWithDirectInvertedOrientation = sum(both_inverted_direct),
      PercentPairsWithDirectInvertedOrientation = pct(sum(both_inverted_direct), n_pairs),
      PairsWithAnyReverseComplement = sum(r1_any_rc | r2_any_rc),
      PercentPairsWithAnyReverseComplement = pct(sum(r1_any_rc | r2_any_rc), n_pairs),
      PairsWithAnyUnexpectedPrimer = sum(any_unexpected),
      PercentPairsWithAnyUnexpectedPrimer = pct(sum(any_unexpected), n_pairs),
      stringsAsFactors = FALSE
    )
  }

  rm(fq1, fq2, fq1_selected, fq2_selected, seqs_by_direction)
  gc(verbose = FALSE)

  list(
    inventory = inventory,
    primer = rbind_or_empty(primer_rows),
    positions = rbind_or_empty(position_rows),
    pairs = rbind_or_empty(pair_rows)
  )
}

###############################################################################
# 5. DESCOBERTA DAS ETAPAS
###############################################################################

metadata <- read_metadata(metadata_path)

stage_pairs <- list(
  discover_fastq_pairs("antes_remocao", raw_path)
)

if (dir.exists(trimmed_path) && length(list_fastq_files(trimmed_path)) > 0L) {
  stage_pairs[[length(stage_pairs) + 1L]] <- discover_fastq_pairs(
    "apos_remocao",
    trimmed_path
  )
}

all_pairs <- rbind_or_empty(stage_pairs)

if (nrow(all_pairs) == 0L) {
  abort(
    "Nenhum par FASTQ foi encontrado. Verifique raw_path e trimmed_path."
  )
}

###############################################################################
# 6. EXECUÇÃO COM REGISTRO DE ERROS
###############################################################################

inventory_list <- list()
primer_list <- list()
position_list <- list()
pair_summary_list <- list()
issues <- list()

for (i in seq_len(nrow(all_pairs))) {
  pair_row <- all_pairs[i, , drop = FALSE]

  cat(
    sprintf(
      "[%d/%d] %s | %s\n",
      i,
      nrow(all_pairs),
      pair_row$Stage,
      pair_row$RawFastqID
    )
  )

  result <- tryCatch(
    analyze_fastq_pair(pair_row, metadata),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    issues <- add_issue(
      issues = issues,
      severity = "ERROR",
      stage = pair_row$Stage,
      sample_id = pair_row$RawFastqID,
      code = "ANALYSIS_FAILED",
      message = conditionMessage(result)
    )
    next
  }

  inventory_list[[length(inventory_list) + 1L]] <- result$inventory
  primer_list[[length(primer_list) + 1L]] <- result$primer
  position_list[[length(position_list) + 1L]] <- result$positions
  pair_summary_list[[length(pair_summary_list) + 1L]] <- result$pairs
}

inventory <- rbind_or_empty(inventory_list)
primer_report <- rbind_or_empty(primer_list)
position_report <- rbind_or_empty(position_list)
pair_report <- rbind_or_empty(pair_summary_list)

###############################################################################
# 7. ALERTAS DE QUALIDADE
###############################################################################

if (nrow(primer_report) > 0L) {
  raw_expected <- primer_report[
    primer_report$Stage == "antes_remocao" &
      primer_report$SearchRegion == "five_prime" &
      primer_report$MaxMismatch == 0L &
      primer_report$ExpectedForDirection,
    ,
    drop = FALSE
  ]

  for (i in seq_len(nrow(raw_expected))) {
    row <- raw_expected[i, , drop = FALSE]

    if (
      !is.na(row$PercentReadsWithPrimer) &&
        row$PercentReadsWithPrimer < minimum_expected_primer_percent_raw
    ) {
      issues <- add_issue(
        issues,
        "WARNING",
        row$Stage,
        row$SampleID,
        "LOW_EXPECTED_PRIMER_RAW",
        sprintf(
          "%s em %s foi detectado em %.3f%% das reads; limiar: %.1f%%.",
          row$PrimerName,
          row$Direction,
          row$PercentReadsWithPrimer,
          minimum_expected_primer_percent_raw
        )
      )
    }

    if (
      !is.na(row$PercentReadsWithNInRegion) &&
        row$PercentReadsWithNInRegion > maximum_N_percent_5p
    ) {
      issues <- add_issue(
        issues,
        "WARNING",
        row$Stage,
        row$SampleID,
        "HIGH_N_5P",
        sprintf(
          "%.3f%% das reads de %s possuem N nos primeiros %d nt.",
          row$PercentReadsWithNInRegion,
          row$Direction,
          window_5p
        )
      )
    }
  }

  trimmed_expected <- primer_report[
    primer_report$Stage == "apos_remocao" &
      primer_report$SearchRegion == "five_prime" &
      primer_report$MaxMismatch == 0L &
      primer_report$ExpectedForDirection,
    ,
    drop = FALSE
  ]

  for (i in seq_len(nrow(trimmed_expected))) {
    row <- trimmed_expected[i, , drop = FALSE]

    if (
      !is.na(row$PercentReadsWithPrimer) &&
        row$PercentReadsWithPrimer > maximum_residual_primer_percent_trimmed
    ) {
      issues <- add_issue(
        issues,
        "WARNING",
        row$Stage,
        row$SampleID,
        "RESIDUAL_PRIMER_AFTER_TRIMMING",
        sprintf(
          "%s residual em %s: %.3f%%; limiar: %.1f%%.",
          row$PrimerName,
          row$Direction,
          row$PercentReadsWithPrimer,
          maximum_residual_primer_percent_trimmed
        )
      )
    }
  }
}

issues_report <- rbind_or_empty(issues)

if (nrow(issues_report) == 0L) {
  issues_report <- data.frame(
    Severity = character(),
    Stage = character(),
    SampleID = character(),
    Code = character(),
    Message = character(),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 8. COMPARAÇÃO ANTES × DEPOIS
###############################################################################

before_after <- data.frame()

if (
  nrow(inventory) > 0L &&
    all(c("antes_remocao", "apos_remocao") %in% unique(inventory$Stage))
) {
  sample_ids <- intersect(
    inventory$SampleID[inventory$Stage == "antes_remocao"],
    inventory$SampleID[inventory$Stage == "apos_remocao"]
  )

  comparison_rows <- list()

  for (sample_id in sample_ids) {
    inv_before <- inventory[
      inventory$SampleID == sample_id &
        inventory$Stage == "antes_remocao",
      ,
      drop = FALSE
    ][1L, , drop = FALSE]

    inv_after <- inventory[
      inventory$SampleID == sample_id &
        inventory$Stage == "apos_remocao",
      ,
      drop = FALSE
    ][1L, , drop = FALSE]

    for (max_mismatch in mismatch_5p) {
      pair_before <- pair_report[
        pair_report$SampleID == sample_id &
          pair_report$Stage == "antes_remocao" &
          pair_report$MaxMismatch == max_mismatch,
        ,
        drop = FALSE
      ]

      pair_after <- pair_report[
        pair_report$SampleID == sample_id &
          pair_report$Stage == "apos_remocao" &
          pair_report$MaxMismatch == max_mismatch,
        ,
        drop = FALSE
      ]

      get_expected_percent <- function(stage, direction, primer_name) {
        row <- primer_report[
          primer_report$SampleID == sample_id &
            primer_report$Stage == stage &
            primer_report$Direction == direction &
            primer_report$SearchRegion == "five_prime" &
            primer_report$PrimerName == primer_name &
            primer_report$MaxMismatch == max_mismatch,
          ,
          drop = FALSE
        ]

        if (nrow(row) == 0L) NA_real_ else row$PercentReadsWithPrimer[1L]
      }

      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        SampleID = sample_id,
        SampleLabel = inv_before$SampleLabel,
        Run = inv_before$Run,
        MaxMismatch = max_mismatch,
        ReadPairsBefore = inv_before$TotalReadPairs,
        ReadPairsAfter = inv_after$TotalReadPairs,
        ReadPairsRemovedBetweenStages =
          inv_before$TotalReadPairs - inv_after$TotalReadPairs,
        PercentReadPairsRetained = pct(
          inv_after$TotalReadPairs,
          inv_before$TotalReadPairs
        ),
        PercentBothExpectedBefore =
          if (nrow(pair_before) == 0L) NA_real_
          else pair_before$PercentPairsWithBothExpected[1L],
        PercentBothExpectedAfter =
          if (nrow(pair_after) == 0L) NA_real_
          else pair_after$PercentPairsWithBothExpected[1L],
        Percent341F_R1_Before =
          get_expected_percent("antes_remocao", "R1", "341F"),
        Percent341F_R1_After =
          get_expected_percent("apos_remocao", "R1", "341F"),
        Percent805R_R2_Before =
          get_expected_percent("antes_remocao", "R2", "805R"),
        Percent805R_R2_After =
          get_expected_percent("apos_remocao", "R2", "805R"),
        stringsAsFactors = FALSE
      )
    }
  }

  before_after <- rbind_or_empty(comparison_rows)
}

###############################################################################
# 9. PARÂMETROS E REPRODUTIBILIDADE
###############################################################################

parameter_report <- data.frame(
  Parameter = c(
    "ScriptVersion",
    "ExecutionStart",
    "BasePath",
    "MetadataPath",
    "RawPath",
    "TrimmedPath",
    "OutputPath",
    "Primer341F",
    "Primer341F_RC",
    "Primer805R",
    "Primer805R_RC",
    "FivePrimeWindowNt",
    "MismatchFivePrime",
    "MismatchWholeRead",
    "AllowIndels",
    "MaxReads",
    "SamplingMethod",
    "RandomSeed",
    "MinimumExpectedPrimerPercentRaw",
    "MaximumResidualPrimerPercentTrimmed",
    "MaximumNPercentFivePrime",
    "RVersion",
    "ShortReadVersion",
    "BiostringsVersion"
  ),
  Value = c(
    script_version,
    format(execution_start, "%Y-%m-%d %H:%M:%S %z"),
    base_path,
    metadata_path,
    raw_path,
    trimmed_path,
    output_path,
    as.character(primers[["341F"]]),
    as.character(primers[["341F_RC"]]),
    as.character(primers[["805R"]]),
    as.character(primers[["805R_RC"]]),
    as.character(window_5p),
    paste(mismatch_5p, collapse = ","),
    paste(mismatch_whole_read, collapse = ","),
    as.character(allow_indels),
    if (is.finite(max_reads)) as.character(max_reads) else "all",
    sampling_method,
    as.character(random_seed),
    as.character(minimum_expected_primer_percent_raw),
    as.character(maximum_residual_primer_percent_trimmed),
    as.character(maximum_N_percent_5p),
    R.version.string,
    as.character(packageVersion("ShortRead")),
    as.character(packageVersion("Biostrings"))
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# 10. SALVAR TABELAS
###############################################################################

write_tsv(inventory, "01_inventario_fastq.tsv")
write_tsv(primer_report, "02_presenca_primers_detalhada.tsv")
write_tsv(pair_report, "03_resumo_primers_por_par.tsv")
write_tsv(position_report, "04_distribuicao_posicoes_primers.tsv")
write_tsv(before_after, "05_comparacao_antes_depois.tsv")
write_tsv(issues_report, "06_alertas_qc.tsv")
write_tsv(parameter_report, "07_parametros_execucao.tsv")

capture.output(
  sessionInfo(),
  file = file.path(log_path, "sessionInfo.txt")
)

###############################################################################
# 11. GRÁFICOS
###############################################################################

plot_expected_primers <- function(stage, filename, title_text) {
  data <- primer_report[
    primer_report$Stage == stage &
      primer_report$SearchRegion == "five_prime" &
      primer_report$MaxMismatch == 0L &
      primer_report$ExpectedForDirection,
    ,
    drop = FALSE
  ]

  if (nrow(data) == 0L) {
    return(invisible(NULL))
  }

  samples <- unique(data$SampleLabel)
  matrix_values <- matrix(
    NA_real_,
    nrow = 2L,
    ncol = length(samples),
    dimnames = list(c("341F em R1", "805R em R2"), samples)
  )

  for (sample_label in samples) {
    r1 <- data[
      data$SampleLabel == sample_label &
        data$Direction == "R1" &
        data$PrimerName == "341F",
      "PercentReadsWithPrimer"
    ]
    r2 <- data[
      data$SampleLabel == sample_label &
        data$Direction == "R2" &
        data$PrimerName == "805R",
      "PercentReadsWithPrimer"
    ]

    if (length(r1) > 0L) matrix_values["341F em R1", sample_label] <- r1[1L]
    if (length(r2) > 0L) matrix_values["805R em R2", sample_label] <- r2[1L]
  }

  png(
    filename = file.path(plot_path, filename),
    width = 1800,
    height = 1000,
    res = 160
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(10, 5, 4, 2) + 0.1)
  barplot(
    matrix_values,
    beside = TRUE,
    ylim = c(0, 100),
    las = 2,
    ylab = "Reads com primer (%)",
    main = title_text,
    legend.text = rownames(matrix_values),
    args.legend = list(x = "bottom", inset = c(0, -0.35), horiz = TRUE, bty = "n")
  )
  abline(h = c(maximum_residual_primer_percent_trimmed,
               minimum_expected_primer_percent_raw), lty = 3)

  invisible(NULL)
}

plot_pairs_expected <- function(stage, filename, title_text) {
  data <- pair_report[
    pair_report$Stage == stage &
      pair_report$MaxMismatch == 0L,
    ,
    drop = FALSE
  ]

  if (nrow(data) == 0L) {
    return(invisible(NULL))
  }

  png(
    filename = file.path(plot_path, filename),
    width = 1600,
    height = 900,
    res = 160
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(10, 5, 4, 2) + 0.1)
  barplot(
    height = data$PercentPairsWithBothExpected,
    names.arg = data$SampleLabel,
    ylim = c(0, 100),
    las = 2,
    ylab = "Pares com 341F em R1 e 805R em R2 (%)",
    main = title_text
  )
  abline(h = 90, lty = 3)

  invisible(NULL)
}

plot_retention <- function() {
  data <- before_after[
    before_after$MaxMismatch == 0L,
    ,
    drop = FALSE
  ]

  if (nrow(data) == 0L) {
    return(invisible(NULL))
  }

  png(
    filename = file.path(plot_path, "04_retencao_pares_antes_depois.png"),
    width = 1600,
    height = 900,
    res = 160
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mar = c(10, 5, 4, 2) + 0.1)
  barplot(
    height = data$PercentReadPairsRetained,
    names.arg = data$SampleLabel,
    ylim = c(0, 100),
    las = 2,
    ylab = "Pares retidos (%)",
    main = "Retenção de pares entre as etapas"
  )

  invisible(NULL)
}

if (generate_plots && nrow(primer_report) > 0L) {
  try(
    plot_expected_primers(
      "antes_remocao",
      "01_primers_esperados_antes_remocao.png",
      "Primers esperados antes da remoção — correspondência exata"
    ),
    silent = TRUE
  )

  try(
    plot_pairs_expected(
      "antes_remocao",
      "02_pares_com_ambos_primers_antes.png",
      "Pares com ambos os primers esperados antes da remoção"
    ),
    silent = TRUE
  )

  if ("apos_remocao" %in% primer_report$Stage) {
    try(
      plot_expected_primers(
        "apos_remocao",
        "03_primers_residuais_apos_remocao.png",
        "Presença residual dos primers após a remoção"
      ),
      silent = TRUE
    )
    try(plot_retention(), silent = TRUE)
  }
}

###############################################################################
# 12. RELATÓRIO HTML
###############################################################################

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

format_html_value <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      is.numeric(x),
      format(round(x, 3), trim = TRUE, scientific = FALSE),
      html_escape(x)
    )
  )
}

dataframe_to_html <- function(df, max_rows = Inf) {
  if (nrow(df) == 0L) {
    return("<p><em>Nenhum registro.</em></p>")
  }

  if (is.finite(max_rows) && nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }

  header <- paste0(
    "<tr>",
    paste0("<th>", html_escape(names(df)), "</th>", collapse = ""),
    "</tr>"
  )

  body <- apply(
    df,
    1L,
    function(row) {
      paste0(
        "<tr>",
        paste0("<td>", format_html_value(row), "</td>", collapse = ""),
        "</tr>"
      )
    }
  )

  paste0(
    "<table>",
    "<thead>", header, "</thead>",
    "<tbody>", paste(body, collapse = "\n"), "</tbody>",
    "</table>"
  )
}

build_overview_table <- function() {
  if (nrow(pair_report) == 0L || nrow(primer_report) == 0L) {
    return(data.frame())
  }

  base <- pair_report[
    pair_report$MaxMismatch == 0L,
    c(
      "Stage",
      "SampleID",
      "SampleLabel",
      "Run",
      "ReadPairsAnalyzed",
      "PercentPairsWithBothExpected",
      "PercentPairsWithAnyUnexpectedPrimer"
    ),
    drop = FALSE
  ]

  base$Percent341F_R1 <- NA_real_
  base$Percent805R_R2 <- NA_real_
  base$PercentN_R1_5p <- NA_real_
  base$PercentN_R2_5p <- NA_real_

  for (i in seq_len(nrow(base))) {
    sample_id <- base$SampleID[i]
    stage <- base$Stage[i]

    r1 <- primer_report[
      primer_report$SampleID == sample_id &
        primer_report$Stage == stage &
        primer_report$Direction == "R1" &
        primer_report$SearchRegion == "five_prime" &
        primer_report$PrimerName == "341F" &
        primer_report$MaxMismatch == 0L,
      ,
      drop = FALSE
    ]

    r2 <- primer_report[
      primer_report$SampleID == sample_id &
        primer_report$Stage == stage &
        primer_report$Direction == "R2" &
        primer_report$SearchRegion == "five_prime" &
        primer_report$PrimerName == "805R" &
        primer_report$MaxMismatch == 0L,
      ,
      drop = FALSE
    ]

    if (nrow(r1) > 0L) {
      base$Percent341F_R1[i] <- r1$PercentReadsWithPrimer[1L]
      base$PercentN_R1_5p[i] <- r1$PercentReadsWithNInRegion[1L]
    }

    if (nrow(r2) > 0L) {
      base$Percent805R_R2[i] <- r2$PercentReadsWithPrimer[1L]
      base$PercentN_R2_5p[i] <- r2$PercentReadsWithNInRegion[1L]
    }
  }

  base
}

overview <- build_overview_table()

inventory_html <- if (nrow(inventory) == 0L) {
  data.frame()
} else {
  inventory[
    ,
    c(
      "Stage", "SampleID", "SampleLabel", "Run", "RawFastqID",
      "TotalReadPairs", "ReadIDsAligned", "ReadPairsAnalyzed",
      "SamplingMethod", "R1_MD5", "R2_MD5"
    ),
    drop = FALSE
  ]
}

execution_end <- Sys.time()
elapsed_seconds <- round(
  as.numeric(difftime(execution_end, execution_start, units = "secs")),
  1
)

stages_text <- paste(unique(all_pairs$Stage), collapse = ", ")
n_errors <- sum(issues_report$Severity == "ERROR")
n_warnings <- sum(issues_report$Severity == "WARNING")

html_lines <- c(
  "<!doctype html>",
  "<html lang='pt-BR'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>Relatório detalhado de primers 341F/805R</title>",
  "<style>",
  "body{font-family:Arial,Helvetica,sans-serif;margin:32px;line-height:1.45;color:#222;}",
  "h1,h2,h3{color:#17365d;}",
  "table{border-collapse:collapse;width:100%;margin:12px 0 24px 0;font-size:13px;}",
  "th,td{border:1px solid #bbb;padding:6px 8px;text-align:left;vertical-align:top;}",
  "th{background:#d9eaf7;position:sticky;top:0;}",
  ".ok{padding:10px;background:#e8f5e9;border-left:5px solid #2e7d32;}",
  ".warn{padding:10px;background:#fff8e1;border-left:5px solid #f9a825;}",
  ".error{padding:10px;background:#ffebee;border-left:5px solid #c62828;}",
  "code{background:#f2f2f2;padding:2px 4px;}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>Relatório detalhado de verificação dos primers 341F/805R</h1>",
  sprintf(
    "<p>Gerado em %s. Duração: %.1f segundos.</p>",
    format(execution_end, "%d/%m/%Y %H:%M:%S"),
    elapsed_seconds
  ),
  "<h2>1. Escopo</h2>",
  sprintf(
    paste0(
      "<p>Etapas analisadas: <strong>%s</strong>. ",
      "A busca na extremidade 5' utilizou uma janela de %d nt, ",
      "com 0, 1 e 2 mismatches; a busca em toda a read utilizou ",
      "correspondência exata. Inserções e deleções: %s.</p>"
    ),
    html_escape(stages_text),
    window_5p,
    if (allow_indels) "permitidas" else "não permitidas"
  ),
  sprintf(
    "<p>Reads analisadas por arquivo: <strong>%s</strong>.</p>",
    if (is.finite(max_reads)) as.character(max_reads) else "todas"
  ),
  "<h2>2. Checkpoints</h2>",
  sprintf(
    "<p class='%s'>Amostras/etapas concluídas: %d. Erros: %d. Alertas: %d.</p>",
    if (n_errors > 0L) "error" else if (n_warnings > 0L) "warn" else "ok",
    nrow(inventory),
    n_errors,
    n_warnings
  ),
  "<h2>3. Visão geral por amostra</h2>",
  dataframe_to_html(overview),
  "<h2>4. Inventário e integridade dos FASTQ</h2>",
  dataframe_to_html(inventory_html),
  "<h2>5. Presença conjunta dos primers por par</h2>",
  dataframe_to_html(pair_report),
  "<h2>6. Comparação antes e depois da remoção</h2>",
  dataframe_to_html(before_after),
  "<h2>7. Alertas e erros</h2>",
  dataframe_to_html(issues_report),
  "<h2>8. Interpretação dos campos</h2>",
  paste0(
    "<p><code>PairsWithBothExpected</code> representa pares nos quais ",
    "341F foi encontrado em R1 e 805R em R2 dentro da janela 5'. ",
    "<code>PairsWithNeitherExpected</code> representa pares sem os dois ",
    "primers esperados sob a tolerância indicada.</p>"
  ),
  paste0(
    "<p><code>PercentReadsWithPrimer</code> é calculado por arquivo e direção. ",
    "O valor por par não pode ser obtido pela média dos percentuais de R1 e R2; ",
    "por isso foi calculado diretamente com os pares correspondentes.</p>"
  ),
  paste0(
    "<p><code>ReadPairsRemovedBetweenStages</code> é a diferença de contagem ",
    "entre as pastas antes e depois. Essa diferença não identifica, isoladamente, ",
    "o motivo do descarte quando outras etapas de filtragem foram aplicadas na ",
    "mesma pasta de saída.</p>"
  ),
  "<h2>9. Arquivos detalhados</h2>",
  "<ul>",
  "<li>01_inventario_fastq.tsv</li>",
  "<li>02_presenca_primers_detalhada.tsv</li>",
  "<li>03_resumo_primers_por_par.tsv</li>",
  "<li>04_distribuicao_posicoes_primers.tsv</li>",
  "<li>05_comparacao_antes_depois.tsv</li>",
  "<li>06_alertas_qc.tsv</li>",
  "<li>07_parametros_execucao.tsv</li>",
  "<li>logs/sessionInfo.txt</li>",
  "</ul>",
  "</body>",
  "</html>"
)

writeLines(
  html_lines,
  con = file.path(output_path, "relatorio_detalhado_primers.html"),
  useBytes = TRUE
)

###############################################################################
# 13. RESUMO NO CONSOLE E CHECKPOINT FINAL
###############################################################################

cat("\n============================================================\n")
cat("VERIFICAÇÃO DE PRIMERS CONCLUÍDA\n")
cat("============================================================\n")
cat("Saída:", output_path, "\n")
cat("Etapas:", stages_text, "\n")
cat("Amostras/etapas concluídas:", nrow(inventory), "\n")
cat("Erros:", n_errors, "\n")
cat("Alertas:", n_warnings, "\n")
cat("Relatório HTML:",
    file.path(output_path, "relatorio_detalhado_primers.html"),
    "\n")

if (!"apos_remocao" %in% unique(all_pairs$Stage)) {
  cat(
    "\n[AVISO] A pasta pós-remoção não foi analisada.\n",
    "Use PRIMER_TRIMMED_PATH para indicar outra pasta de FASTQ sem primers.\n",
    sep = ""
  )
}

if (n_errors > 0L) {
  stop(
    sprintf(
      "A execução terminou com %d erro(s). Consulte 06_alertas_qc.tsv.",
      n_errors
    ),
    call. = FALSE
  )
}

cat("\nCheckpoint final: OK\n")
})
