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
run_pipeline_script("05_comparacao_banco_de_dados.R", "taxonomia", function(ctx) {
###############################################################################
# SCRIPT 05 — CLASSIFICACAO TAXONOMICA INTEGRADA POR PRIORIDADE HIERARQUICA
#
# ORDEM DE PRIORIDADE
#   Kingdom -> Genus:
#     GSR-DB (confidence=disable) > BEExact > SILVA 138.2 > RDP 19 >
#     Greengenes2.
#   Species:
#     BLAST NCBI 16S exato e nao ambiguo > GSR-DB > BEExact > SILVA >
#     RDP > Greengenes2.
#
# REGRA DE DECISAO (Kingdom -> Species)
#   1. Apenas valores informativos e compativeis com a linhagem final ja
#      definida nos ranks anteriores podem participar da decisao.
#   2. Se nenhum banco fornecer valor valido, o rank permanece NA.
#   3. Se apenas um banco fornecer valor valido, sua classificacao e mantida.
#   4. Se dois ou mais bancos fornecerem valores e todos concordarem, o taxon
#      comum e mantido; a fonte registrada e o banco de maior prioridade.
#   5. Se dois ou mais bancos discordarem, vence o banco de maior prioridade
#      entre os bancos validos naquele rank.
#   6. A maioria numerica nao supera um banco de maior prioridade.
#
# SPECIES
#   - BLAST exato e nao ambiguo tem prioridade maxima quando satisfaz:
#       identidade=100%; cobertura=100%; mismatches=0; gaps=0;
#       alinhamento integral da ASV; melhores hits informativos concordantes.
#   - Quando o Genus do BLAST exato difere do Genus integrado, o par
#     Genus+Species do BLAST e promovido para preservar coerencia hierarquica.
#   - Sem Species exata por BLAST, aplica-se a hierarquia dos cinco
#     classificadores.
#   - Fonte unica e aceita, desde que o nome seja binomial, nao composto e
#     compativel com o Genus final e com a linhagem superior.
#   - O addSpecies exato do SILVA substitui o Species k-mer do proprio SILVA,
#     mas permanece um unico resultado SILVA.
#
###############################################################################

gc()
options(encoding = "UTF-8", stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(Biostrings)
})

VERSAO        <- "8.1_BLAST_species_exata_prioridade_maxima"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

###############################################################################
# 0. CAMINHOS E PARAMETROS
###############################################################################

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$stage$root

dir_comp <- file.path(output_path, "comparacao_bancos")
dir_aud  <- file.path(dir_comp, "auditorias")
dir_fig  <- file.path(dir_comp, "figuras")
dir_fa   <- file.path(dir_comp, "nao_identificadas_fasta")
dir_gsr_sens <- file.path(dir_comp, "sensibilidade_gsr")

for (d in c(output_path, dir_comp, dir_aud, dir_fig, dir_fa, dir_gsr_sens)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Falha ao criar diretorio: ", d, call. = FALSE)
  if (file.access(d, 2L) != 0L) {
    stop("Sem permissao de escrita no diretorio: ", d, call. = FALSE)
  }
}

NIVEIS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
NIVEIS_SUPERIORES <- c("Kingdom", "Phylum", "Class", "Order", "Family")
RANKS_CLASSIFICACAO <- NIVEIS
NIVEL_MINIMO_PRINCIPAL <- "Order"

CONT_KINGDOM <- c("Eukaryota", "Euk")
CONT_ORDER   <- c("Chloroplast")
CONT_FAMILY  <- c("Mitochondria")

TERMOS_VAGOS_REGEX <- paste(
  c(
    "uncultured", "unidentified", "unknown", "metagenome",
    "environmental", "uncharacterized", "Incertae Sedis",
    "incertae sedis", "gut metagenome"
  ),
  collapse = "|"
)

arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
arq_asvmap <- ctx$contracts[["asv_sequences"]]

arq_silva <- ctx$contracts[["taxa_silva"]]
arq_rdp   <- ctx$contracts[["taxa_rdp"]]
arq_gg2   <- ctx$contracts[["taxa_gg2"]]
arq_bee   <- ctx$contracts[["taxa_beexact"]]

# GSR entra como UM unico banco votante usando a classificacao principal
# confidence=disable, coerente com o Script 02e. O alias legado e aceito somente
# quando o arquivo explicito nao estiver disponivel.
arq_gsr_disable_principal <- ctx$contracts[["taxa_gsr_disable"]]
arq_gsr_disable_alias <- ctx$contracts[["taxa_gsr_disable_alias"]]
arq_gsr_disable <- if (file.exists(arq_gsr_disable_principal)) {
  arq_gsr_disable_principal
} else {
  arq_gsr_disable_alias
}
arq_gsr07 <- ctx$contracts[["taxa_gsr07"]]

arq_silva_species <- ctx$contracts[["species_silva"]]

blast_root <- ctx$layout$stages$blast$root
arq_blast97   <- file.path(blast_root, "rds", "taxa_blast97.rds")
arq_blast100  <- file.path(blast_root, "rds", "taxa_blast100.rds")
arq_blast_evid <- file.path(blast_root, "rds", "blast_evidencias_por_asv.rds")
arq_meta_blast <- file.path(blast_root, "tabelas", "metadata_blast.csv")

###############################################################################
# 1. LOG, I/O E VALIDACOES
###############################################################################

log_file <- file.path(dir_comp, "classificacao_hierarquica_log.txt")
if (file.exists(log_file)) unlink(log_file, force = TRUE)

log_msg <- function(msg, tipo = "INFO") {
  linha <- sprintf("[%s] <%s> %s", format(Sys.time(), "%H:%M:%S"), tipo, msg)
  cat(linha, "\n")
  cat(linha, "\n", file = log_file, append = TRUE)
  invisible(linha)
}

abort <- function(...) stop(sprintf(...), call. = FALSE)

valor_valido <- function(x) {
  length(x) == 1L &&
    !is.na(x) &&
    trimws(as.character(x)) != ""
}

salvar_csv <- function(x, arq) {
  write.csv(
    x, arq,
    row.names = FALSE, quote = TRUE,
    fileEncoding = "UTF-8", na = ""
  )
  if (!file.exists(arq) || is.na(file.size(arq)) || file.size(arq) == 0L) {
    abort("Falha ao salvar CSV ou arquivo vazio: %s", arq)
  }
  invisible(arq)
}

salvar_rds <- function(x, arq) {
  saveRDS(x, arq)
  if (!file.exists(arq) || file.size(arq) == 0L) {
    abort("Falha ao salvar RDS: %s", arq)
  }
  invisible(arq)
}

ler_matriz_obrigatoria <- function(arq, nome) {
  if (!file.exists(arq)) abort("%s ausente: %s", nome, arq)
  x <- tryCatch(
    as.matrix(readRDS(arq)),
    error = function(e) abort("Falha ao ler %s: %s", nome, conditionMessage(e))
  )
  if (nrow(x) == 0L || ncol(x) == 0L) {
    abort("%s esta vazio.", nome)
  }
  if (is.null(rownames(x)) || any(is.na(rownames(x))) || any(rownames(x) == "")) {
    abort("%s sem rownames validos.", nome)
  }
  if (anyDuplicated(rownames(x)) > 0L) {
    abort("%s possui rownames duplicados.", nome)
  }
  if (is.null(colnames(x)) || any(is.na(colnames(x))) || any(colnames(x) == "")) {
    abort("%s sem colnames validos.", nome)
  }
  x
}

garantir_ranks <- function(mat, ranks = NIVEIS) {
  x <- as.matrix(mat)
  faltam <- setdiff(ranks, colnames(x))
  if (length(faltam) > 0L) {
    add <- matrix(
      NA_character_,
      nrow = nrow(x), ncol = length(faltam),
      dimnames = list(rownames(x), faltam)
    )
    x <- cbind(x, add)
  }
  x[, ranks, drop = FALSE]
}

validar_universo <- function(mat, nome, universo) {
  faltam <- setdiff(universo, rownames(mat))
  extras <- setdiff(rownames(mat), universo)
  if (length(faltam) > 0L || length(extras) > 0L) {
    abort(
      "%s nao corresponde ao universo canonico: faltantes=%d; extras=%d.",
      nome, length(faltam), length(extras)
    )
  }
  invisible(TRUE)
}

com_asv_id <- function(mat, seq2id) {
  seqs <- rownames(mat)
  data.frame(
    ASV_ID = unname(seq2id[seqs]),
    ASV_seq = seqs,
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

###############################################################################
# 2. UNIVERSO CANONICO
###############################################################################

cat("=============================================================\n")
cat("CLASSIFICACAO TAXONOMICA POR PRIORIDADE — v", VERSAO, "\n", sep = "")
cat("Data:", DATA_EXECUCAO, "\n")
cat("=============================================================\n\n")

if (!file.exists(arq_seqtab)) abort("seqtab ausente: %s", arq_seqtab)
if (!file.exists(arq_asvmap)) abort("ASV_sequences.tsv ausente: %s", arq_asvmap)

seqtab <- tryCatch(
  as.matrix(readRDS(arq_seqtab)),
  error = function(e) abort("Falha ao ler seqtab_global_nochim.rds: %s", conditionMessage(e))
)
if (nrow(seqtab) == 0L || ncol(seqtab) == 0L) {
  abort("seqtab vazio.")
}
if (
  is.null(rownames(seqtab)) || any(is.na(rownames(seqtab))) ||
  any(rownames(seqtab) == "") || anyDuplicated(rownames(seqtab)) > 0L
) {
  abort("seqtab deve possuir SampleIDs unicos e nao vazios nas linhas.")
}
if (
  is.null(colnames(seqtab)) || any(is.na(colnames(seqtab))) ||
  any(colnames(seqtab) == "") || anyDuplicated(colnames(seqtab)) > 0L
) {
  abort("seqtab deve possuir sequencias ASV unicas e nao vazias nas colunas.")
}
if (!is.numeric(seqtab)) abort("seqtab deve conter contagens numericas.")
if (anyNA(seqtab) || any(!is.finite(seqtab)) || any(seqtab < 0)) {
  abort("seqtab contem NA, valor nao finito ou contagem negativa.")
}
if (max(abs(seqtab - round(seqtab))) > 1e-8) {
  abort("seqtab contem contagens nao inteiras.")
}
if (any(rowSums(seqtab) == 0)) abort("seqtab contem amostra com zero reads.")
if (any(colSums(seqtab) == 0)) abort("seqtab contem ASV com soma zero.")
if (!all(grepl("^[ACGTURYKMSWBDHVNacgturykmswbdhvn.-]+$", colnames(seqtab)))) {
  abort("seqtab possui colunas que nao sao sequencias 16S/IUPAC validas.")
}

asv_map <- read.delim(
  arq_asvmap,
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)
if (!all(c("ASV_ID", "Sequence", "Origem") %in% colnames(asv_map))) {
  abort("ASV_sequences.tsv deve conter ASV_ID, Sequence e Origem.")
}
if (anyNA(asv_map$ASV_ID) || any(trimws(as.character(asv_map$ASV_ID)) == "")) {
  abort("ASV_sequences.tsv contem ASV_ID ausente ou vazio.")
}
if (anyNA(asv_map$Sequence) || any(trimws(as.character(asv_map$Sequence)) == "")) {
  abort("ASV_sequences.tsv contem Sequence ausente ou vazia.")
}
if (anyDuplicated(asv_map$ASV_ID) > 0L) abort("ASV_ID duplicado.")
if (anyDuplicated(asv_map$Sequence) > 0L) abort("Sequence duplicada.")

all_asvs <- colnames(seqtab)
if (!setequal(all_asvs, asv_map$Sequence)) {
  abort(
    "ASV_sequences.tsv e seqtab possuem universos diferentes: faltantes=%d; extras=%d.",
    length(setdiff(all_asvs, asv_map$Sequence)),
    length(setdiff(asv_map$Sequence, all_asvs))
  )
}
idx_map <- match(all_asvs, asv_map$Sequence)
if (anyNA(idx_map)) abort("ASVs do seqtab ausentes em ASV_sequences.tsv.")

map_ord <- asv_map[idx_map, c("ASV_ID", "Sequence", "Origem"), drop = FALSE]
seq2id <- setNames(map_ord$ASV_ID, map_ord$Sequence)
id2seq <- setNames(map_ord$Sequence, map_ord$ASV_ID)
n_asvs_total <- length(all_asvs)

reads_por_seq <- colSums(seqtab)
prev_por_seq <- colSums(seqtab > 0)
amostras_por_seq <- apply(
  seqtab > 0, 2,
  function(x) paste(rownames(seqtab)[x], collapse = ";")
)

log_msg(
  sprintf("Universo canonico: %d ASVs | %d amostras.", n_asvs_total, nrow(seqtab)),
  "OK"
)

###############################################################################
# 3. CARREGAMENTO DOS BANCOS
#    5 fontes hierarquizadas + GSR 0.7 sensibilidade + BLAST conferencia
###############################################################################

silva_raw <- garantir_ranks(ler_matriz_obrigatoria(arq_silva, "SILVA 138.2"))
rdp_raw   <- garantir_ranks(ler_matriz_obrigatoria(arq_rdp,   "RDP 19"))
gg2_raw   <- garantir_ranks(ler_matriz_obrigatoria(arq_gg2,   "Greengenes2"))
bee_raw   <- garantir_ranks(ler_matriz_obrigatoria(arq_bee,   "BEExact"))

for (nm in c("silva_raw", "rdp_raw", "gg2_raw", "bee_raw")) {
  validar_universo(get(nm), nm, all_asvs)
}

matriz_taxa_na <- function(universo) {
  matrix(
    NA_character_,
    nrow = length(universo), ncol = length(NIVEIS),
    dimnames = list(universo, NIVEIS)
  )
}

gsr_disable_disponivel <- file.exists(arq_gsr_disable)
gsr07_disponivel <- file.exists(arq_gsr07)

if (!gsr_disable_disponivel) {
  abort(
    paste0(
      "GSR principal confidence=disable ausente. Esperado: ",
      arq_gsr_disable_principal,
      " (ou alias legado: ", arq_gsr_disable_alias, ")."
    )
  )
}
if (!identical(arq_gsr_disable, arq_gsr_disable_principal)) {
  log_msg(
    paste0(
      "Arquivo explicito GSR confidence=disable ausente; usando alias legado: ",
      arq_gsr_disable
    ),
    "WARN"
  )
}

gsr_disable_raw <- garantir_ranks(
  ler_matriz_obrigatoria(arq_gsr_disable, "GSR V3-V4 confidence=disable")
)
validar_universo(gsr_disable_raw, "gsr_disable_raw", all_asvs)

if (gsr07_disponivel) {
  gsr07_raw <- garantir_ranks(
    ler_matriz_obrigatoria(arq_gsr07, "GSR V3-V4 confidence=0.7")
  )
  validar_universo(gsr07_raw, "gsr07_raw", all_asvs)
} else {
  gsr07_raw <- matriz_taxa_na(all_asvs)
  log_msg(
    "GSR confidence=0.7 ausente; a sensibilidade taxonomica do Script 06 sera omitida.",
    "WARN"
  )
}

# Fixar a mesma ordem canonica em todas as matrizes antes da normalizacao.
# Isso evita desalinhamento posicional de atributos auxiliares (por exemplo,
# flags de Species composta) quando um banco entrega as mesmas ASVs em outra ordem.
reordenar_universo <- function(mat, universo) {
  mat[universo, , drop = FALSE]
}
silva_raw <- reordenar_universo(silva_raw, all_asvs)
rdp_raw <- reordenar_universo(rdp_raw, all_asvs)
gg2_raw <- reordenar_universo(gg2_raw, all_asvs)
bee_raw <- reordenar_universo(bee_raw, all_asvs)
gsr_disable_raw <- reordenar_universo(gsr_disable_raw, all_asvs)
gsr07_raw <- reordenar_universo(gsr07_raw, all_asvs)

blast97_raw  <- garantir_ranks(ler_matriz_obrigatoria(arq_blast97,  "BLAST Genus candidato"))
blast100_raw <- garantir_ranks(ler_matriz_obrigatoria(arq_blast100, "BLAST exato integral"))
validar_universo(blast97_raw,  "blast97_raw",  all_asvs)
validar_universo(blast100_raw, "blast100_raw", all_asvs)

if (!file.exists(arq_blast_evid)) {
  abort("blast_evidencias_por_asv.rds ausente: %s", arq_blast_evid)
}
blast_evid <- readRDS(arq_blast_evid)
if (!is.data.frame(blast_evid)) abort("blast_evidencias_por_asv.rds nao e data.frame.")

req_evid <- c(
  "ASV_ID", "ASV_seq", "Genus_blast_exato", "Species_blast_exata",
  "Genus_blast_candidato", "Melhor_identidade_genus",
  "Melhor_cobertura_genus", "Ambiguo_Genus_exato",
  "Ambiguo_Species_exata", "Ambiguo_Genus_candidato",
  "Atingiu_MAX_TARGET_SEQS"
)
faltam_evid <- setdiff(req_evid, colnames(blast_evid))
if (length(faltam_evid) > 0L) {
  abort(
    "blast_evidencias_por_asv.rds sem coluna(s): %s",
    paste(faltam_evid, collapse = ", ")
  )
}

if (!file.exists(arq_meta_blast)) abort("metadata_blast.csv ausente.")
meta_blast <- read.csv(
  arq_meta_blast,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

# Remover BOM/espacos acidentais dos nomes das colunas.
bom_utf8 <- intToUtf8(0xFEFF)
colnames(meta_blast) <- trimws(
  sub(paste0("^", bom_utf8), "", colnames(meta_blast))
)

if (nrow(meta_blast) != 1L) {
  abort(
    "metadata_blast.csv deve conter exatamente uma linha; encontradas %d.",
    nrow(meta_blast)
  )
}

# Compatibilidade controlada:
# versões anteriores do Script 04 usavam Perc_ident_exact para a mesma camada.
# O alias só corrige o nome do campo; os critérios biológicos continuam sendo
# validados abaixo.
if (
  !"Perc_ident_species" %in% colnames(meta_blast) &&
  "Perc_ident_exact" %in% colnames(meta_blast)
) {
  meta_blast$Perc_ident_species <- meta_blast$Perc_ident_exact
  log_msg(
    paste0(
      "metadata_blast.csv usa o nome legado Perc_ident_exact; ",
      "interpretado como Perc_ident_species."
    ),
    "WARN"
  )
}

req_meta_blast <- c(
  "Perc_ident_species", "Perc_ident_genus",
  "Qcov_species", "Qcov_genus", "Max_target_seqs",
  "ASVs_entrada"
)

faltam_meta_blast <- setdiff(req_meta_blast, colnames(meta_blast))
if (length(faltam_meta_blast) > 0L) {
  abort(
    paste0(
      "metadata_blast.csv incompativel; faltam: %s. ",
      "Substitua o Script 04 pela versao revisada e reprocese o TSV bruto."
    ),
    paste(faltam_meta_blast, collapse = ", ")
  )
}

num_meta_blast <- lapply(
  meta_blast[1L, req_meta_blast, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

invalidos_meta <- vapply(
  num_meta_blast,
  function(x) length(x) != 1L || !is.finite(x),
  logical(1)
)

if (any(invalidos_meta)) {
  abort(
    "metadata_blast.csv contem valor(es) nao numerico(s): %s.",
    paste(names(num_meta_blast)[invalidos_meta], collapse = ", ")
  )
}

if (
  num_meta_blast$Perc_ident_species != 100 ||
  num_meta_blast$Qcov_species != 100
) {
  abort("Camada BLAST de Species deve usar identidade=100 e cobertura=100.")
}

req_contrato_species <- c(
  "Species_mismatches",
  "Species_gap_openings",
  "Species_alignment_full_length",
  "Criterio_species"
)
faltam_contrato_species <- setdiff(
  req_contrato_species,
  colnames(meta_blast)
)
if (length(faltam_contrato_species) > 0L) {
  abort(
    paste0(
      "metadata_blast.csv sem o contrato completo de Species: %s. ",
      "Reprocese o TSV bruto com o Script 04 revisado."
    ),
    paste(faltam_contrato_species, collapse = ", ")
  )
}

species_mismatches_meta <- suppressWarnings(as.numeric(
  as.character(meta_blast$Species_mismatches[1L])
))
species_gaps_meta <- suppressWarnings(as.numeric(
  as.character(meta_blast$Species_gap_openings[1L])
))
species_full_meta <- toupper(trimws(as.character(
  meta_blast$Species_alignment_full_length[1L]
))) %in% c("TRUE", "T", "1")
criterio_species_meta <- tolower(as.character(
  meta_blast$Criterio_species[1L]
))

if (
  !is.finite(species_mismatches_meta) ||
  !is.finite(species_gaps_meta) ||
  species_mismatches_meta != 0 ||
  species_gaps_meta != 0 ||
  !isTRUE(species_full_meta)
) {
  abort(
    paste0(
      "Camada BLAST de Species nao satisfaz o contrato: ",
      "mismatches=0, gaps=0 e alinhamento integral."
    )
  )
}

tokens_criterio_species <- c(
  "100pct_identidade_e_cobertura",
  "zero_mismatch_gap",
  "full_length",
  "nao_ambigua"
)
if (!all(vapply(
  tokens_criterio_species,
  function(token) grepl(token, criterio_species_meta, fixed = TRUE),
  logical(1)
))) {
  abort(
    paste0(
      "Criterio_species de metadata_blast.csv nao comprova a camada ",
      "exata, integral e nao ambigua exigida para prioridade maxima."
    )
  )
}

if (
  num_meta_blast$Perc_ident_genus < 97 ||
  num_meta_blast$Qcov_genus != 100
) {
  abort("Camada BLAST candidata de Genus deve usar identidade>=97 e cobertura=100.")
}

if (num_meta_blast$Max_target_seqs < 1) {
  abort("Max_target_seqs invalido em metadata_blast.csv.")
}

if (num_meta_blast$ASVs_entrada != n_asvs_total) {
  abort(
    "metadata_blast.csv registra %d ASVs, mas o universo canonico possui %d.",
    as.integer(num_meta_blast$ASVs_entrada),
    n_asvs_total
  )
}

# A cobertura de 100% por qcovs nao substitui a checagem explícita de
# qstart=1, qend=qlen e alignment.length=qlen usada na versão revisada.
if (!"Criterio_genus" %in% colnames(meta_blast)) {
  abort(
    paste0(
      "metadata_blast.csv sem Criterio_genus. ",
      "Reprocese o TSV bruto com o Script 04 revisado."
    )
  )
}

criterio_genus_meta <- tolower(as.character(meta_blast$Criterio_genus[1L]))
tokens_criterio_genus <- c("qstart=1", "qend=qlen", "length=qlen")
criterio_genus_ok <- all(vapply(
  tokens_criterio_genus,
  function(token) grepl(token, criterio_genus_meta, fixed = TRUE),
  logical(1)
))

if (!criterio_genus_ok) {
  abort(
    paste0(
      "Os resultados BLAST foram gerados por uma versao anterior do Script 04, ",
      "sem a validacao explicita de alinhamento integral para Genus. ",
      "Reprocese o TSV bruto com REUTILIZAR_BLAST_BRUTO=TRUE; ",
      "nao e necessario repetir a busca blastn."
    )
  )
}

log_msg(
  sprintf(
    "Contrato BLAST validado: Species=%g/%g; Genus>=%g/%g; %d ASVs.",
    num_meta_blast$Perc_ident_species,
    num_meta_blast$Qcov_species,
    num_meta_blast$Perc_ident_genus,
    num_meta_blast$Qcov_genus,
    as.integer(num_meta_blast$ASVs_entrada)
  ),
  "OK"
)

if (
  anyDuplicated(blast_evid$ASV_seq) > 0L ||
  !setequal(blast_evid$ASV_seq, all_asvs)
) {
  abort("blast_evidencias_por_asv.rds nao corresponde ao universo canonico.")
}
blast_evid <- blast_evid[match(all_asvs, blast_evid$ASV_seq), , drop = FALSE]

###############################################################################
# 4. SPECIES EXATA DO SILVA (addSpecies)
###############################################################################

silva_species_exact <- setNames(rep(NA_character_, n_asvs_total), all_asvs)

if (file.exists(arq_silva_species)) {
  silva_species_audit <- read.csv(
    arq_silva_species,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  col_chave <- intersect(
    c("ASV_seq", "ASV", "Sequence"),
    colnames(silva_species_audit)
  )[1L]
  req <- c("Species", "Species_fonte")
  if (is.na(col_chave) || !all(req %in% colnames(silva_species_audit))) {
    abort(
      "species_fonte_silva138.csv precisa de chave ASV/ASV_seq/Sequence, Species e Species_fonte."
    )
  }
  idx_exact <- silva_species_audit$Species_fonte == "exact_match"
  chaves <- as.character(silva_species_audit[[col_chave]][idx_exact])
  vals   <- as.character(silva_species_audit$Species[idx_exact])
  chaves_seq <- chaves
  idx_ids <- chaves %in% names(id2seq)
  chaves_seq[idx_ids] <- unname(id2seq[chaves[idx_ids]])
  manter <- chaves_seq %in% all_asvs
  silva_species_exact[chaves_seq[manter]] <- vals[manter]
} else {
  log_msg("species_fonte_silva138.csv ausente; addSpecies exata sera ignorada.", "WARN")
}

###############################################################################
# 5. NORMALIZACAO TAXONOMICA (sinônimos + underscore + binômio)
###############################################################################

MAPA_FILO_NCBI <- c(
  "Bacteroidetes" = "Bacteroidota",
  "Planctomycetes" = "Planctomycetota",
  "Acidobacteria" = "Acidobacteriota",
  "Actinobacteria" = "Actinomycetota",
  "Actinobacteriota" = "Actinomycetota",
  "Aquificae" = "Aquificota",
  "Armatimonadetes" = "Armatimonadota",
  "Balneolaeota" = "Balneolota",
  "Caldiserica" = "Caldisericota",
  "Calditrichaeota" = "Calditrichota",
  "Chlamydiae" = "Chlamydiota",
  "Chlorobi" = "Chlorobiota",
  "Chloroflexi" = "Chloroflexota",
  "Chrysiogenetes" = "Chrysiogenota",
  "Crenarchaeota" = "Thermoproteota",
  "Deferribacteres" = "Deferribacterota",
  "Deinococcus-Thermus" = "Deinococcota",
  "Dictyoglomi" = "Dictyoglomota",
  "Elusimicrobia" = "Elusimicrobiota",
  "Fibrobacteres" = "Fibrobacterota",
  "Firmicutes" = "Bacillota",
  "Fusobacteria" = "Fusobacteriota",
  "Gemmatimonadetes" = "Gemmatimonadota",
  "Ignavibacteriae" = "Ignavibacteriota",
  "Kiritimatiellaeota" = "Kiritimatiellota",
  "Lentisphaerae" = "Lentisphaerota",
  "Nitrospinae" = "Nitrospinota",
  "Nitrospirae" = "Nitrospirota",
  "Proteobacteria" = "Pseudomonadota",
  "Rhodothermaeota" = "Rhodothermota",
  "Spirochaetes" = "Spirochaetota",
  "Synergistetes" = "Synergistota",
  "Tenericutes" = "Mycoplasmatota",
  "Thaumarchaeota" = "Nitrososphaerota",
  "Thermodesulfobacteria" = "Thermodesulfobacteriota",
  "Thermotogae" = "Thermotogota",
  "Verrucomicrobia" = "Verrucomicrobiota",
  "Epsilonproteobacteria" = "Campylobacterota",
  "Thermomicrobia" = "Thermomicrobiota",
  "Oligoflexia" = "Bdellovibrionota",
  "Deltaproteobacteria" = "Desulfobacterota"
)

# Normalizacao basica: remove prefixos "d__/k__/p__/...", underscores,
# espacos multiplos e placeholders. A remocao de sufixos auxiliares GTDB e
# aplicada SOMENTE nas fontes que usam essa nomenclatura (Greengenes2/GSR),
# impedindo alteracao indevida de nomes validos de SILVA, RDP ou BEExact.
normalizar_taxa_basica <- function(
  x,
  nivel = NULL,
  harmonizar_sufixos_gtdb = FALSE
) {
  y <- as.character(x)
  ok <- !is.na(y)
  y[ok] <- trimws(y[ok])
  y[ok] <- gsub("^[a-z]__+", "", y[ok], ignore.case = TRUE)
  y[ok] <- gsub("_", " ", y[ok], fixed = TRUE)
  y[ok] <- trimws(gsub("\\s{2,}", " ", y[ok]))
  y[!is.na(y) & y == ""] <- NA_character_

  if (
    isTRUE(harmonizar_sufixos_gtdb) &&
    !is.null(nivel) &&
    nivel %in% c("Phylum", "Class", "Order", "Family", "Genus")
  ) {
    idx <- !is.na(y)
    # Enterobacteriaceae A 725029 -> Enterobacteriaceae
    # Enterobacter B 681665       -> Enterobacter
    y[idx] <- sub("^(.+?)\\s+[A-Z]\\s+[0-9]+$", "\\1", y[idx], perl = TRUE)
    # Enterobacterales 737866     -> Enterobacterales
    # Restrito a nomes de uma palavra; preserva "Clostridium sensu stricto 1".
    y[idx] <- sub("^(\\S+)\\s+[0-9]+$", "\\1", y[idx], perl = TRUE)
    # Bacillota I / Bacilli A     -> Bacillota / Bacilli
    y[idx] <- sub("^(.+?)\\s+[A-Z]$", "\\1", y[idx], perl = TRUE)
    y[idx] <- trimws(y[idx])
  }

  if (any(!is.na(y))) {
    idx_vago <- !is.na(y) & grepl(TERMOS_VAGOS_REGEX, y, ignore.case = TRUE)
    y[idx_vago] <- NA_character_
  }

  if (!is.null(nivel) && nivel == "Genus") {
    y[!is.na(y) & grepl("^[0-9]+$", y)] <- NA_character_
  }

  y
}

normalizar_filo <- function(x, harmonizar_sufixos_gtdb = FALSE) {
  y <- normalizar_taxa_basica(
    x,
    nivel = "Phylum",
    harmonizar_sufixos_gtdb = harmonizar_sufixos_gtdb
  )
  idx <- !is.na(y) & y %in% names(MAPA_FILO_NCBI)
  y[idx] <- unname(MAPA_FILO_NCBI[y[idx]])
  y
}

# Reconstroi Species como binomio canonico. Quando o banco fornece uma
# linhagem com informacao adicional de cepa, apenas o binomio e mantido.
# Nomes "Candidatus Genus species" sao preservados como trinômio nomenclatural
# (Candidatus + genero + epiteto), e nao tratados como classificacao vaga.
normalizar_species_binomial <- function(
  mat,
  harmonizar_sufixos_gtdb = FALSE
) {
  species <- normalizar_taxa_basica(
    mat[, "Species"],
    nivel = "Species",
    harmonizar_sufixos_gtdb = FALSE
  )
  genus_l <- normalizar_taxa_basica(
    mat[, "Genus"],
    nivel = "Genus",
    harmonizar_sufixos_gtdb = harmonizar_sufixos_gtdb
  )

  separadores_composto <- c(":", ";", "/", "[", "]")
  tem_separador <- !is.na(species) & Reduce("|", lapply(
    separadores_composto,
    function(p) grepl(p, species, fixed = TRUE)
  ))
  tem_hifen_com_maisc <- !is.na(species) &
    grepl("-", species, fixed = TRUE) &
    grepl("[A-Z].*[A-Z]", species)
  eh_composta <- tem_separador | tem_hifen_com_maisc

  # Marcadores que nao representam especie resolvida.
  nao_especifica <- !is.na(species) & grepl(
    "(^|\\s)(sp\\.?|spp\\.?|cf\\.?|aff\\.?|bacterium|archaeon)($|\\s)",
    species,
    ignore.case = TRUE,
    perl = TRUE
  )
  species[nao_especifica] <- NA_character_

  out <- rep(NA_character_, length(species))

  candidatus_bin <- !is.na(species) & grepl(
    "^Candidatus\\s+[A-Z][A-Za-z.-]+\\s+[a-z0-9][A-Za-z0-9._-]*",
    species,
    perl = TRUE
  )
  if (any(candidatus_bin)) {
    out[candidatus_bin] <- vapply(
      strsplit(species[candidatus_bin], "\\s+"),
      function(z) paste(z[seq_len(min(3L, length(z)))], collapse = " "),
      character(1)
    )
  }

  binom_natural <- !is.na(species) & !candidatus_bin & grepl(
    "^[A-Z][A-Za-z.-]+\\s+[a-z0-9][A-Za-z0-9._-]*",
    species,
    perl = TRUE
  )
  if (any(binom_natural)) {
    out[binom_natural] <- vapply(
      strsplit(species[binom_natural], "\\s+"),
      function(z) paste(z[seq_len(min(2L, length(z)))], collapse = " "),
      character(1)
    )
  }

  epiteto_puro <- !is.na(species) & !candidatus_bin &
    !binom_natural & !grepl("\\s", species)
  ok_recon <- epiteto_puro & !is.na(genus_l) & nzchar(genus_l)
  out[ok_recon] <- paste(genus_l[ok_recon], species[ok_recon])

  # Preservar o texto composto apenas para auditoria; ele sera invalidado antes
  # da votacao de Species.
  out[eh_composta & !is.na(species)] <- species[eh_composta & !is.na(species)]

  attr(out, "composta") <- eh_composta
  out
}

normalizar_mat <- function(
  mat,
  origem,
  harmonizar_sufixos_gtdb = FALSE
) {
  x <- garantir_ranks(mat)
  for (rk in colnames(x)) {
    x[, rk] <- if (rk == "Phylum") {
      normalizar_filo(
        x[, rk],
        harmonizar_sufixos_gtdb = harmonizar_sufixos_gtdb
      )
    } else {
      normalizar_taxa_basica(
        x[, rk],
        nivel = rk,
        harmonizar_sufixos_gtdb = harmonizar_sufixos_gtdb
      )
    }
  }
  storage.mode(x) <- "character"

  sp_bin <- normalizar_species_binomial(
    x,
    harmonizar_sufixos_gtdb = harmonizar_sufixos_gtdb
  )
  x <- cbind(x, Species_binomial = sp_bin)

  attr(x, "origem") <- origem
  comp <- attr(sp_bin, "composta")
  names(comp) <- rownames(x)
  attr(x, "species_composta") <- comp
  x
}

###############################################################################
# 6. FUNCOES DE COMPARACAO TAXONOMICA
###############################################################################

# Os sufixos auxiliares ja foram tratados de forma especifica por fonte na
# normalizacao. Aqui apenas removemos espacos residuais.
raiz_taxon <- function(x) {
  if (!valor_valido(x)) return(NA_character_)
  trimws(as.character(x))
}

# expandir_hifenizados: separa "Escherichia-Shigella" em duas partes.
# Nao separa por espaço (para preservar binomios).
expandir_hifenizados <- function(x) {
  if (!valor_valido(x)) return(character())
  y <- raiz_taxon(x)
  partes <- strsplit(y, "-", fixed = TRUE)[[1L]]
  partes <- trimws(partes)
  partes[nzchar(partes)]
}

# mesmo_taxon: comparacao tolerante a sufixos GTDB, hifens e case.
mesmo_taxon <- function(a, b) {
  if (!valor_valido(a) || !valor_valido(b)) return(FALSE)
  ra <- tolower(expandir_hifenizados(a))
  rb <- tolower(expandir_hifenizados(b))
  length(intersect(ra, rb)) > 0L
}

# Species exige o mesmo binomio canonico; compatibilidade apenas de Genus nao
# e suficiente para criar concordancia no nivel de especie.
mesma_species_binomial <- function(a, b) {
  if (!valor_valido(a) || !valor_valido(b)) return(FALSE)
  chave <- function(x) {
    tolower(trimws(gsub("\\s{2,}", " ", as.character(x))))
  }
  identical(chave(a), chave(b))
}

# Extrai o Genus de um binomio canonico. Para "Candidatus Genus species",
# o Genus nomenclatural e representado por "Candidatus Genus".
extrair_genus_do_binomio <- function(sp_bin) {
  if (!valor_valido(sp_bin)) return(NA_character_)
  partes <- strsplit(trimws(as.character(sp_bin)), "\\s+")[[1L]]
  if (
    length(partes) >= 2L && identical(tolower(partes[1L]), "candidatus")
  ) {
    return(paste(partes[1:2], collapse = " "))
  }
  if (length(partes) >= 1L) partes[1L] else NA_character_
}

# linhagem_compativel: uma fonte so pode decidir o rank atual se sua linhagem
# nos ranks anteriores nao contradizer a classificacao integrada ja fixada.
linhagem_compativel <- function(mat, seq_asv, taxa_final, ranks_anteriores) {
  if (length(ranks_anteriores) == 0L) return(TRUE)
  for (anc in ranks_anteriores) {
    c_final <- taxa_final[seq_asv, anc]
    c_banco <- if (anc %in% colnames(mat)) mat[seq_asv, anc] else NA_character_
    if (
      valor_valido(c_final) && valor_valido(c_banco) &&
      !mesmo_taxon(c_final, c_banco)
    ) return(FALSE)
  }
  TRUE
}

###############################################################################
# 7. APLICACAO DA NORMALIZACAO E DEFINICAO DOS VOTANTES
###############################################################################

silva_norm <- normalizar_mat(
  silva_raw, "SILVA 138.2", harmonizar_sufixos_gtdb = FALSE
)
rdp_norm <- normalizar_mat(
  rdp_raw, "RDP 19", harmonizar_sufixos_gtdb = FALSE
)
gg2_norm <- normalizar_mat(
  gg2_raw, "Greengenes2", harmonizar_sufixos_gtdb = TRUE
)
bee_norm <- normalizar_mat(
  bee_raw, "BEExact", harmonizar_sufixos_gtdb = FALSE
)
gsr_disable_norm <- normalizar_mat(
  gsr_disable_raw,
  "GSR V3-V4 confidence=disable",
  harmonizar_sufixos_gtdb = TRUE
)
gsr07_norm <- normalizar_mat(
  gsr07_raw,
  "GSR V3-V4 confidence=0.7",
  harmonizar_sufixos_gtdb = TRUE
)

# O voto GSR corresponde integralmente a uma unica execucao do classificador:
# confidence=disable, definida como principal no Script 02e. Nao e permitido
# combinar 0.7 e disable em ranks diferentes da mesma ASV, pois isso criaria uma
# linhagem hibrida que nao foi produzida por nenhum classificador.
gsr_voto <- gsr_disable_norm
gsr_fonte_voto <- matrix(
  NA_character_,
  nrow = nrow(gsr_voto), ncol = ncol(gsr_voto),
  dimnames = dimnames(gsr_voto)
)
gsr_fonte_voto[!is.na(gsr_voto) & trimws(gsr_voto) != ""] <-
  "confidence_disable_principal"
attr(gsr_voto, "origem") <- "GSR V3-V4 confidence=disable (voto principal)"

blast97_raw <- reordenar_universo(blast97_raw, all_asvs)
blast100_raw <- reordenar_universo(blast100_raw, all_asvs)
blast97_norm <- normalizar_mat(
  blast97_raw, "BLAST candidato Genus", harmonizar_sufixos_gtdb = FALSE
)
blast100_norm <- normalizar_mat(
  blast100_raw, "BLAST exato integral", harmonizar_sufixos_gtdb = FALSE
)

# Species exata do SILVA (addSpecies): reconstruir binomial a partir do valor
# original + Genus normalizado do SILVA.
silva_species_exact_bin <- rep(NA_character_, n_asvs_total)
names(silva_species_exact_bin) <- all_asvs

canonizar_species_exata <- function(sp, genus) {
  if (!valor_valido(sp)) return(NA_character_)
  sp_clean <- normalizar_taxa_basica(sp, nivel = "Species")
  if (!valor_valido(sp_clean)) return(NA_character_)
  if (grepl(
    "(^|\\s)(sp\\.?|spp\\.?|cf\\.?|aff\\.?|bacterium|archaeon)($|\\s)",
    sp_clean, ignore.case = TRUE, perl = TRUE
  )) return(NA_character_)

  partes <- strsplit(sp_clean, "\\s+")[[1L]]
  if (
    length(partes) >= 3L && identical(tolower(partes[1L]), "candidatus")
  ) {
    return(paste(partes[1:3], collapse = " "))
  }
  if (length(partes) >= 2L && grepl("^[A-Z]", partes[1L])) {
    return(paste(partes[1:2], collapse = " "))
  }
  if (length(partes) == 1L && valor_valido(genus)) {
    return(paste(genus, partes[1L]))
  }
  NA_character_
}

for (s in all_asvs) {
  silva_species_exact_bin[s] <- canonizar_species_exata(
    silva_species_exact[s], silva_norm[s, "Genus"]
  )
}

# O addSpecies exato substitui apenas a coluna Species_binomial usada pelo voto
# SILVA. Nao cria um sexto voto e nao altera os demais ranks do SILVA.
silva_voto <- silva_norm
idx_silva_exact <- !is.na(silva_species_exact_bin) &
  trimws(silva_species_exact_bin) != ""
silva_voto[idx_silva_exact, "Species_binomial"] <-
  silva_species_exact_bin[idx_silva_exact]
silva_composta <- attr(silva_voto, "species_composta")
if (!is.null(silva_composta)) {
  silva_composta[idx_silva_exact] <- FALSE
  attr(silva_voto, "species_composta") <- silva_composta
}

# Cinco bancos integram a classificacao. A ordem do objeto e a ordem de
# prioridade decisoria e deve permanecer explicita.
ORDEM_PRIORIDADE <- c(
  "gsr", "beexact", "silva138", "rdp19", "greengenes2"
)
BANCOS <- list(
  gsr         = gsr_voto,
  beexact     = bee_norm,
  silva138    = silva_voto,
  rdp19       = rdp_norm,
  greengenes2 = gg2_norm
)
ORDEM_BANCOS <- ORDEM_PRIORIDADE

# O BLAST exato entra somente na decisao de Species. Nao e contado como um
# sexto classificador nos ranks superiores.
ORDEM_SPECIES <- c("blast100_exact", ORDEM_BANCOS)

# Camadas GSR preservadas para auditoria; somente "gsr" participa da hierarquia.
BANCOS_AUDITORIA <- c(
  BANCOS,
  list(
    gsr_07 = gsr07_norm,
    gsr_disable = gsr_disable_norm
  )
)

log_msg(
  sprintf("Ordem de prioridade: %s.", paste(ORDEM_BANCOS, collapse = " > ")),
  "OK"
)
log_msg(
  "GSR confidence=disable definido como primeira prioridade; GSR 0.7 permanece sensibilidade.",
  "OK"
)
log_msg(
  sprintf("Prioridade de Species: %s.", paste(ORDEM_SPECIES, collapse = " > ")),
  "OK"
)

if (
  length(ORDEM_BANCOS) != 5L ||
  !identical(ORDEM_BANCOS, c("gsr", "beexact", "silva138", "rdp19", "greengenes2")) ||
  !identical(ORDEM_SPECIES, c(
    "blast100_exact", "gsr", "beexact", "silva138", "rdp19", "greengenes2"
  )) ||
  !setequal(names(BANCOS), ORDEM_BANCOS)
) {
  abort(
    paste0(
      "Contrato de prioridade invalido. Kingdom-Genus: ",
      "gsr > beexact > silva138 > rdp19 > greengenes2; Species: ",
      "blast100_exact > gsr > beexact > silva138 > rdp19 > greengenes2."
    )
  )
}

# Auditoria da fonte efetivamente usada pelo voto GSR em cada rank.
auditoria_fonte_gsr <- do.call(rbind, lapply(NIVEIS, function(rk) {
  fonte <- gsr_fonte_voto[, rk]
  tab <- as.data.frame(table(
    Fonte = ifelse(is.na(fonte), "sem_valor", fonte),
    useNA = "no"
  ), stringsAsFactors = FALSE)
  data.frame(
    Rank = rk,
    Fonte = tab$Fonte,
    N_ASVs = tab$Freq,
    stringsAsFactors = FALSE
  )
}))
salvar_csv(
  auditoria_fonte_gsr,
  file.path(dir_aud, "auditoria_fonte_voto_gsr.csv")
)

# Auditoria das alteracoes de grafia/harmonizacao aplicadas ao Greengenes2.
auditoria_norm_gg2 <- do.call(rbind, lapply(NIVEIS, function(rk) {
  bruto <- as.character(gg2_raw[all_asvs, rk])
  norm <- as.character(gg2_norm[all_asvs, rk])
  mudou <- (
    (is.na(bruto) & !is.na(norm)) |
    (!is.na(bruto) & is.na(norm)) |
    (!is.na(bruto) & !is.na(norm) & trimws(bruto) != trimws(norm))
  )
  if (!any(mudou)) return(NULL)
  data.frame(
    ASV_ID = unname(seq2id[all_asvs[mudou]]),
    ASV_seq = all_asvs[mudou],
    Rank = rk,
    Valor_original = bruto[mudou],
    Valor_harmonizado = norm[mudou],
    stringsAsFactors = FALSE
  )
}))
if (!is.null(auditoria_norm_gg2) && nrow(auditoria_norm_gg2) > 0L) {
  salvar_csv(
    auditoria_norm_gg2,
    file.path(dir_aud, "auditoria_harmonizacao_greengenes2.csv")
  )
}

###############################################################################
# 8. CLASSIFICACAO KINGDOM -> GENUS POR PRIORIDADE HIERARQUICA
###############################################################################

taxa_integrada_todas <- matrix(
  NA_character_,
  nrow = n_asvs_total, ncol = length(NIVEIS),
  dimnames = list(all_asvs, NIVEIS)
)
fonte_todas <- taxa_integrada_todas
n_concordam_todas <- matrix(
  0L,
  nrow = n_asvs_total, ncol = length(NIVEIS),
  dimnames = list(all_asvs, NIVEIS)
)

# Numero maximo de bancos que apresentam um taxon compativel com o valor de
# qualquer banco. Campo apenas descritivo; nao participa da decisao.
calcular_apoio_maximo <- function(votos, bancos_validos, comparador) {
  bancos_validos <- as.character(bancos_validos)
  if (length(bancos_validos) == 0L) return(0L)
  max(vapply(
    bancos_validos,
    function(b) sum(vapply(
      bancos_validos,
      function(outro) comparador(votos[b], votos[outro]),
      logical(1)
    )),
    integer(1)
  ))
}

# Agrupamento descritivo deterministico. A decisao nao depende deste numero.
contar_grupos_taxonomicos <- function(votos, bancos_validos, comparador) {
  bancos_validos <- as.character(bancos_validos)
  if (length(bancos_validos) == 0L) return(0L)
  representantes <- character()
  for (b in bancos_validos) {
    val <- votos[b]
    pertence <- length(representantes) > 0L && any(vapply(
      representantes,
      function(rep_val) comparador(val, rep_val),
      logical(1)
    ))
    if (!pertence) representantes <- c(representantes, val)
  }
  length(representantes)
}

resolver_por_prioridade <- function(
  votos,
  bancos_validos,
  comparador = mesmo_taxon,
  n_bancos_com_valor = length(bancos_validos),
  ordem_prioridade = ORDEM_BANCOS
) {
  ordem_prioridade <- as.character(ordem_prioridade)
  bancos_validos <- ordem_prioridade[
    ordem_prioridade %in% as.character(bancos_validos)
  ]
  n_info <- length(bancos_validos)

  base <- list(
    valor = NA_character_,
    fonte = NA_character_,
    banco_decisor = NA_character_,
    concordantes = character(),
    discordantes = character(),
    n_com_valor = as.integer(n_bancos_com_valor),
    n_informativos = n_info,
    n_grupos = 0L,
    apoio_maximo = 0L,
    proporcao_concordancia = NA_real_,
    houve_conflito = FALSE,
    empate_real = FALSE,
    regra = "nenhum_valor_valido",
    status = "sem_classificacao"
  )

  if (n_info == 0L) {
    if (n_bancos_com_valor > 0L) {
      base$regra <- "todos_os_valores_foram_incompativeis_com_a_linhagem"
      base$status <- "sem_classificacao_compativel"
    }
    return(base)
  }

  # A ordem recebida e normativa; a primeira fonte valida decide.
  banco_decisor <- bancos_validos[1L]
  valor_final <- unname(votos[banco_decisor])
  concordantes <- bancos_validos[vapply(
    bancos_validos,
    function(b) comparador(valor_final, votos[b]),
    logical(1)
  )]
  discordantes <- setdiff(bancos_validos, concordantes)

  base$valor <- valor_final
  base$fonte <- banco_decisor
  base$banco_decisor <- banco_decisor
  base$concordantes <- concordantes
  base$discordantes <- discordantes
  base$n_grupos <- contar_grupos_taxonomicos(votos, bancos_validos, comparador)
  base$apoio_maximo <- calcular_apoio_maximo(votos, bancos_validos, comparador)
  base$proporcao_concordancia <- length(concordantes) / n_info
  base$houve_conflito <- length(discordantes) > 0L

  if (n_info == 1L) {
    base$regra <- "fonte_unica_mantida"
    base$status <- if (n_bancos_com_valor == 1L) {
      "fonte_unica_aceita"
    } else {
      "fonte_unica_aceita_apos_filtro_hierarquico"
    }
  } else if (!base$houve_conflito) {
    base$regra <- "concordancia_entre_bancos"
    base$status <- "concordancia_entre_bancos"
  } else {
    base$regra <- paste0("prioridade_hierarquica_", banco_decisor)
    base$status <- "prioridade_hierarquica_por_conflito"
  }

  base
}

# Detalhamento por ASV x rank.
detalhes_por_rank <- vector(
  "list", length(setdiff(RANKS_CLASSIFICACAO, "Species"))
)
names(detalhes_por_rank) <- setdiff(RANKS_CLASSIFICACAO, "Species")

resolver_rank_prioridade <- function(seq_asv, rk, ranks_anteriores) {
  votos <- setNames(rep(NA_character_, length(ORDEM_BANCOS)), ORDEM_BANCOS)
  compativel <- setNames(rep(FALSE, length(ORDEM_BANCOS)), ORDEM_BANCOS)

  for (banco in ORDEM_BANCOS) {
    mat <- BANCOS[[banco]]
    val <- mat[seq_asv, rk]
    if (!valor_valido(val)) next
    votos[banco] <- val
    compativel[banco] <- linhagem_compativel(
      mat, seq_asv, taxa_integrada_todas, ranks_anteriores
    )
  }

  tem_valor <- !is.na(votos)
  bancos_validos <- ORDEM_BANCOS[tem_valor & compativel]
  res <- resolver_por_prioridade(
    votos = votos,
    bancos_validos = bancos_validos,
    comparador = mesmo_taxon,
    n_bancos_com_valor = sum(tem_valor)
  )
  res$ausentes <- ORDEM_BANCOS[!tem_valor]
  res$incompativeis <- ORDEM_BANCOS[tem_valor & !compativel]
  res$votos_texto <- votos
  res$compatibilidade <- compativel
  res
}

for (rk in setdiff(RANKS_CLASSIFICACAO, "Species")) {
  pos_rk <- match(rk, NIVEIS)
  ranks_anteriores <- if (pos_rk > 1L) {
    NIVEIS[seq_len(pos_rk - 1L)]
  } else {
    character()
  }
  detalhes_rk <- vector("list", n_asvs_total)

  for (i in seq_along(all_asvs)) {
    seq_asv <- all_asvs[i]
    r <- resolver_rank_prioridade(seq_asv, rk, ranks_anteriores)

    taxa_integrada_todas[seq_asv, rk] <- r$valor
    fonte_todas[seq_asv, rk] <- r$fonte
    n_concordam_todas[seq_asv, rk] <- length(r$concordantes)

    detalhes_rk[[i]] <- data.frame(
      ASV_ID = unname(seq2id[seq_asv]),
      ASV_seq = seq_asv,
      Rank = rk,
      Classificacao_final = r$valor,
      Consenso = r$valor,
      Banco_decisor = r$banco_decisor,
      Fonte_grafia = r$fonte,
      N_bancos_com_valor = r$n_com_valor,
      N_bancos_informativos = r$n_informativos,
      N_bancos_concordam_com_resultado = length(r$concordantes),
      N_bancos_concordam = length(r$concordantes),
      Apoio_maximo = r$apoio_maximo,
      Proporcao_concordancia_resultado = r$proporcao_concordancia,
      Proporcao_apoio = r$proporcao_concordancia,
      N_grupos_taxonomicos = r$n_grupos,
      Houve_conflito = r$houve_conflito,
      Empate_real = FALSE,
      Bancos_concordam_com_resultado = paste(r$concordantes, collapse = ";"),
      Bancos_concordam = paste(r$concordantes, collapse = ";"),
      Bancos_discordam = paste(r$discordantes, collapse = ";"),
      Bancos_incompativeis = paste(r$incompativeis, collapse = ";"),
      Bancos_sem_valor = paste(r$ausentes, collapse = ";"),
      Regra_aplicada = r$regra,
      Status = r$status,
      Valor_gsr = unname(r$votos_texto["gsr"]),
      Valor_beexact = unname(r$votos_texto["beexact"]),
      Valor_silva138 = unname(r$votos_texto["silva138"]),
      Valor_rdp19 = unname(r$votos_texto["rdp19"]),
      Valor_greengenes2 = unname(r$votos_texto["greengenes2"]),
      Fonte_gsr_voto = unname(gsr_fonte_voto[seq_asv, rk]),
      Valor_gsr07_auditoria = unname(gsr07_norm[seq_asv, rk]),
      stringsAsFactors = FALSE
    )
  }

  detalhes_por_rank[[rk]] <- do.call(rbind, detalhes_rk)
  log_msg(sprintf("Rank %s resolvido por prioridade hierarquica.", rk), "OK")
}

detalhes_kingdom_genus <- do.call(rbind, detalhes_por_rank)
salvar_csv(
  detalhes_kingdom_genus,
  file.path(dir_aud, "detalhes_classificacao_kingdom_a_genus.csv")
)
# Nome legado preservado.
salvar_csv(
  detalhes_kingdom_genus,
  file.path(dir_aud, "detalhes_consenso_kingdom_a_genus.csv")
)

###############################################################################
# 8B. VALIDACAO DA REGRA DE PRIORIDADE
###############################################################################

validar_decisoes_rank <- function() {
  violacoes <- list()
  k <- 0L

  for (rk in setdiff(RANKS_CLASSIFICACAO, "Species")) {
    det <- detalhes_por_rank[[rk]]
    det <- det[match(all_asvs, det$ASV_seq), , drop = FALSE]

    for (i in seq_along(all_asvs)) {
      seq_asv <- all_asvs[i]
      pos_rk <- match(rk, NIVEIS)
      ranks_anteriores <- if (pos_rk > 1L) {
        NIVEIS[seq_len(pos_rk - 1L)]
      } else {
        character()
      }

      votos <- setNames(rep(NA_character_, length(ORDEM_BANCOS)), ORDEM_BANCOS)
      comp <- setNames(rep(FALSE, length(ORDEM_BANCOS)), ORDEM_BANCOS)
      for (banco in ORDEM_BANCOS) {
        val <- BANCOS[[banco]][seq_asv, rk]
        if (!valor_valido(val)) next
        votos[banco] <- val
        comp[banco] <- linhagem_compativel(
          BANCOS[[banco]], seq_asv, taxa_integrada_todas, ranks_anteriores
        )
      }
      validos <- ORDEM_BANCOS[!is.na(votos) & comp]
      esperado_banco <- if (length(validos)) validos[1L] else NA_character_
      esperado_valor <- if (length(validos)) unname(votos[esperado_banco]) else NA_character_

      banco_ok <- (
        (is.na(esperado_banco) && is.na(det$Banco_decisor[i])) ||
        identical(as.character(esperado_banco), as.character(det$Banco_decisor[i]))
      )
      valor_ok <- (
        (is.na(esperado_valor) && is.na(det$Classificacao_final[i])) ||
        identical(as.character(esperado_valor), as.character(det$Classificacao_final[i]))
      )

      if (!banco_ok || !valor_ok) {
        k <- k + 1L
        violacoes[[k]] <- data.frame(
          ASV_ID = unname(seq2id[seq_asv]),
          ASV_seq = seq_asv,
          Rank = rk,
          Banco_esperado = esperado_banco,
          Banco_obtido = det$Banco_decisor[i],
          Valor_esperado = esperado_valor,
          Valor_obtido = det$Classificacao_final[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(violacoes)) {
    out <- do.call(rbind, violacoes)
    arq <- file.path(dir_aud, "ERRO_validacao_prioridade_hierarquica.csv")
    salvar_csv(out, arq)
    abort(
      "Validacao da prioridade falhou em %d decisao(oes). Consulte %s.",
      nrow(out), arq
    )
  }

  log_msg("Validacao da prioridade Kingdom-Genus concluida sem violacoes.", "OK")
  invisible(TRUE)
}

validar_decisoes_rank()

###############################################################################
# 9. CLASSIFICACAO DE SPECIES — BLAST EXATO COMO PRIORIDADE MAXIMA
###############################################################################

# A matriz taxa_blast100.rds ja contem somente a camada exata produzida pelo
# Script 04. Aqui o contrato e reforcado por ASV: Species e Genus devem ser
# informativos, coerentes entre si e marcados como nao ambiguos.
blast_species_exata <- setNames(
  blast100_norm[all_asvs, "Species_binomial"],
  all_asvs
)
blast_genus_exato <- setNames(
  blast100_norm[all_asvs, "Genus"],
  all_asvs
)

blast_species_valida <- vapply(seq_along(all_asvs), function(i) {
  seq_asv <- all_asvs[i]
  sp <- blast_species_exata[seq_asv]
  gen <- blast_genus_exato[seq_asv]
  gen_sp <- extrair_genus_do_binomio(sp)

  valor_valido(sp) &&
    valor_valido(gen) &&
    valor_valido(gen_sp) &&
    mesmo_taxon(gen, gen_sp) &&
    identical(blast_evid$Ambiguo_Species_exata[i], FALSE) &&
    identical(blast_evid$Ambiguo_Genus_exato[i], FALSE)
}, logical(1))
names(blast_species_valida) <- all_asvs

# Registrar a decisao original de Genus antes de qualquer ajuste induzido por
# uma Species exata do BLAST.
genus_antes_blast <- setNames(
  taxa_integrada_todas[all_asvs, "Genus"],
  all_asvs
)
fonte_genus_antes_blast <- setNames(
  fonte_todas[all_asvs, "Genus"],
  all_asvs
)

genus_ajustado_por_blast <- blast_species_valida & vapply(
  all_asvs,
  function(seq_asv) {
    gen_atual <- genus_antes_blast[seq_asv]
    gen_blast <- blast_genus_exato[seq_asv]
    !valor_valido(gen_atual) || !mesmo_taxon(gen_atual, gen_blast)
  },
  logical(1)
)

if (any(genus_ajustado_por_blast)) {
  seqs_ajuste <- all_asvs[genus_ajustado_por_blast]

  taxa_integrada_todas[seqs_ajuste, "Genus"] <-
    blast_genus_exato[seqs_ajuste]
  fonte_todas[seqs_ajuste, "Genus"] <- "blast100_species_exact"

  det_gen <- detalhes_por_rank[["Genus"]]
  idx_det <- match(seqs_ajuste, det_gen$ASV_seq)
  if (anyNA(idx_det)) {
    abort("Falha ao localizar ASV(s) no detalhamento de Genus para ajuste BLAST.")
  }

  n_conc_gen <- vapply(seqs_ajuste, function(seq_asv) {
    gen_blast <- blast_genus_exato[seq_asv]
    sum(vapply(
      ORDEM_BANCOS,
      function(banco) {
        val <- BANCOS[[banco]][seq_asv, "Genus"]
        valor_valido(val) && mesmo_taxon(val, gen_blast)
      },
      logical(1)
    ))
  }, integer(1))

  n_info_gen <- vapply(seqs_ajuste, function(seq_asv) {
    sum(vapply(
      ORDEM_BANCOS,
      function(banco) valor_valido(BANCOS[[banco]][seq_asv, "Genus"]),
      logical(1)
    ))
  }, integer(1))

  bancos_conc_gen <- vapply(seqs_ajuste, function(seq_asv) {
    gen_blast <- blast_genus_exato[seq_asv]
    ok <- ORDEM_BANCOS[vapply(
      ORDEM_BANCOS,
      function(banco) {
        val <- BANCOS[[banco]][seq_asv, "Genus"]
        valor_valido(val) && mesmo_taxon(val, gen_blast)
      },
      logical(1)
    )]
    paste(ok, collapse = ";")
  }, character(1))

  bancos_disc_gen <- vapply(seqs_ajuste, function(seq_asv) {
    gen_blast <- blast_genus_exato[seq_asv]
    bad <- ORDEM_BANCOS[vapply(
      ORDEM_BANCOS,
      function(banco) {
        val <- BANCOS[[banco]][seq_asv, "Genus"]
        valor_valido(val) && !mesmo_taxon(val, gen_blast)
      },
      logical(1)
    )]
    paste(bad, collapse = ";")
  }, character(1))

  det_gen$Classificacao_final[idx_det] <- blast_genus_exato[seqs_ajuste]
  det_gen$Consenso[idx_det] <- blast_genus_exato[seqs_ajuste]
  det_gen$Banco_decisor[idx_det] <- "blast100_species_exact"
  det_gen$Fonte_grafia[idx_det] <- "blast100_species_exact"
  det_gen$N_bancos_concordam_com_resultado[idx_det] <- n_conc_gen
  det_gen$N_bancos_concordam[idx_det] <- n_conc_gen
  det_gen$Proporcao_concordancia_resultado[idx_det] <-
    ifelse(n_info_gen > 0L, n_conc_gen / n_info_gen, NA_real_)
  det_gen$Proporcao_apoio[idx_det] <-
    ifelse(n_info_gen > 0L, n_conc_gen / n_info_gen, NA_real_)
  det_gen$Houve_conflito[idx_det] <- TRUE
  det_gen$Bancos_concordam_com_resultado[idx_det] <- bancos_conc_gen
  det_gen$Bancos_concordam[idx_det] <- bancos_conc_gen
  det_gen$Bancos_discordam[idx_det] <- bancos_disc_gen
  det_gen$Regra_aplicada[idx_det] <-
    "genus_promovido_por_species_blast_exata"
  det_gen$Status[idx_det] <- "blast_species_exata_ajuste_genus"

  detalhes_por_rank[["Genus"]] <- det_gen
  n_concordam_todas[seqs_ajuste, "Genus"] <- n_conc_gen

  auditoria_ajuste_genus_blast <- data.frame(
    ASV_ID = unname(seq2id[seqs_ajuste]),
    ASV_seq = seqs_ajuste,
    Genus_antes = unname(genus_antes_blast[seqs_ajuste]),
    Fonte_antes = unname(fonte_genus_antes_blast[seqs_ajuste]),
    Genus_BLAST = unname(blast_genus_exato[seqs_ajuste]),
    Species_BLAST = unname(blast_species_exata[seqs_ajuste]),
    Busca_saturada = blast_evid$Atingiu_MAX_TARGET_SEQS[
      match(seqs_ajuste, blast_evid$ASV_seq)
    ],
    Regra = "Species BLAST exata tem prioridade maxima; Genus ajustado para coerencia",
    stringsAsFactors = FALSE
  )
} else {
  auditoria_ajuste_genus_blast <- data.frame(
    ASV_ID = character(), ASV_seq = character(),
    Genus_antes = character(), Fonte_antes = character(),
    Genus_BLAST = character(), Species_BLAST = character(),
    Busca_saturada = logical(), Regra = character(),
    stringsAsFactors = FALSE
  )
}

salvar_csv(
  auditoria_ajuste_genus_blast,
  file.path(dir_aud, "genus_ajustado_por_species_blast_exata.csv")
)

# Regravar o detalhamento Kingdom-Genus apos os ajustes induzidos por BLAST.
detalhes_kingdom_genus <- do.call(rbind, detalhes_por_rank)
salvar_csv(
  detalhes_kingdom_genus,
  file.path(dir_aud, "detalhes_classificacao_kingdom_a_genus.csv")
)
salvar_csv(
  detalhes_kingdom_genus,
  file.path(dir_aud, "detalhes_consenso_kingdom_a_genus.csv")
)

detalhes_species <- vector("list", n_asvs_total)

for (i in seq_along(all_asvs)) {
  seq_asv <- all_asvs[i]
  genus_final <- taxa_integrada_todas[seq_asv, "Genus"]

  votos_classificadores <- setNames(
    rep(NA_character_, length(ORDEM_BANCOS)),
    ORDEM_BANCOS
  )
  ok_classificadores <- setNames(
    rep(FALSE, length(ORDEM_BANCOS)),
    ORDEM_BANCOS
  )
  motivo_invalido_classificadores <- setNames(
    rep(NA_character_, length(ORDEM_BANCOS)),
    ORDEM_BANCOS
  )

  for (banco in ORDEM_BANCOS) {
    mat <- BANCOS[[banco]]
    sp_bin <- mat[seq_asv, "Species_binomial"]
    if (!valor_valido(sp_bin)) next
    votos_classificadores[banco] <- sp_bin

    composta_flag <- attr(mat, "species_composta")[[seq_asv]]
    if (isTRUE(composta_flag)) {
      motivo_invalido_classificadores[banco] <- "species_composta"
      next
    }

    genus_sp <- extrair_genus_do_binomio(sp_bin)
    if (!valor_valido(genus_final)) {
      motivo_invalido_classificadores[banco] <- "genus_final_ausente"
      next
    }
    if (!valor_valido(genus_sp) || !mesmo_taxon(genus_sp, genus_final)) {
      motivo_invalido_classificadores[banco] <- "genus_species_incompativel"
      next
    }

    if (!linhagem_compativel(
      mat, seq_asv, taxa_integrada_todas,
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    )) {
      motivo_invalido_classificadores[banco] <-
        "linhagem_superior_incompativel"
      next
    }

    ok_classificadores[banco] <- TRUE
  }

  bancos_ok_classificadores <- ORDEM_BANCOS[ok_classificadores]
  r_classificadores <- resolver_por_prioridade(
    votos = votos_classificadores,
    bancos_validos = bancos_ok_classificadores,
    comparador = mesma_species_binomial,
    n_bancos_com_valor = sum(!is.na(votos_classificadores)),
    ordem_prioridade = ORDEM_BANCOS
  )

  votos_sp <- setNames(
    rep(NA_character_, length(ORDEM_SPECIES)),
    ORDEM_SPECIES
  )
  votos_sp[ORDEM_BANCOS] <- votos_classificadores
  if (blast_species_valida[seq_asv]) {
    votos_sp["blast100_exact"] <- blast_species_exata[seq_asv]
  }

  ok_sp <- setNames(rep(FALSE, length(ORDEM_SPECIES)), ORDEM_SPECIES)
  ok_sp[ORDEM_BANCOS] <- ok_classificadores
  ok_sp["blast100_exact"] <- blast_species_valida[seq_asv]

  motivo_invalido <- setNames(
    rep(NA_character_, length(ORDEM_SPECIES)),
    ORDEM_SPECIES
  )
  motivo_invalido[ORDEM_BANCOS] <- motivo_invalido_classificadores

  if (!blast_species_valida[seq_asv] &&
      valor_valido(blast100_norm[seq_asv, "Species_binomial"])) {
    motivo_invalido["blast100_exact"] <-
      "blast_species_ausente_ambiguo_ou_incoerente"
  }

  fontes_ok <- ORDEM_SPECIES[ok_sp]
  r_sp <- resolver_por_prioridade(
    votos = votos_sp,
    bancos_validos = fontes_ok,
    comparador = mesma_species_binomial,
    n_bancos_com_valor = sum(!is.na(votos_sp)),
    ordem_prioridade = ORDEM_SPECIES
  )

  blast_usado <- identical(r_sp$banco_decisor, "blast100_exact")
  if (blast_usado) {
    r_sp$regra <- "blast_exato_100pct_prioridade_maxima"
    r_sp$status <- "blast_exato_prioridade_maxima"
  }

  taxa_integrada_todas[seq_asv, "Species"] <- r_sp$valor
  fonte_todas[seq_asv, "Species"] <- r_sp$fonte
  n_concordam_todas[seq_asv, "Species"] <- length(r_sp$concordantes)

  detalhes_species[[i]] <- data.frame(
    ASV_ID = unname(seq2id[seq_asv]),
    ASV_seq = seq_asv,
    Rank = "Species",
    Classificacao_final = r_sp$valor,
    Consenso = r_sp$valor,
    Banco_decisor = r_sp$banco_decisor,
    Fonte_grafia = r_sp$fonte,
    N_bancos_com_valor = r_sp$n_com_valor,
    N_bancos_informativos = r_sp$n_informativos,
    N_bancos_concordam_com_resultado = length(r_sp$concordantes),
    N_bancos_concordam = length(r_sp$concordantes),
    Apoio_maximo = r_sp$apoio_maximo,
    Proporcao_concordancia_resultado = r_sp$proporcao_concordancia,
    Proporcao_apoio = r_sp$proporcao_concordancia,
    N_grupos_taxonomicos = r_sp$n_grupos,
    Houve_conflito = r_sp$houve_conflito,
    Empate_real = FALSE,
    Bancos_concordam_com_resultado = paste(r_sp$concordantes, collapse = ";"),
    Bancos_concordam = paste(r_sp$concordantes, collapse = ";"),
    Bancos_discordam = paste(r_sp$discordantes, collapse = ";"),
    Bancos_validos = paste(fontes_ok, collapse = ";"),
    Bancos_invalidos = paste(
      ORDEM_SPECIES[!is.na(motivo_invalido)],
      collapse = ";"
    ),
    Motivos_invalidos = paste(
      paste0(
        ORDEM_SPECIES[!is.na(motivo_invalido)], "=",
        motivo_invalido[!is.na(motivo_invalido)]
      ),
      collapse = ";"
    ),
    Bancos_sem_valor = paste(ORDEM_SPECIES[is.na(votos_sp)], collapse = ";"),
    Regra_aplicada = r_sp$regra,
    Status = r_sp$status,
    BLAST_exato_prioritario = blast_usado,
    BLAST_busca_saturada = blast_evid$Atingiu_MAX_TARGET_SEQS[i],
    Genus_ajustado_por_BLAST = genus_ajustado_por_blast[seq_asv],
    Genus_antes_BLAST = genus_antes_blast[seq_asv],
    Fonte_Genus_antes_BLAST = fonte_genus_antes_blast[seq_asv],
    Species_classificadores_sem_BLAST = r_classificadores$valor,
    Fonte_classificadores_sem_BLAST = r_classificadores$fonte,
    Valor_blast100_exact = unname(votos_sp["blast100_exact"]),
    Valor_gsr = unname(votos_sp["gsr"]),
    Valor_beexact = unname(votos_sp["beexact"]),
    Valor_silva138 = unname(votos_sp["silva138"]),
    Valor_rdp19 = unname(votos_sp["rdp19"]),
    Valor_greengenes2 = unname(votos_sp["greengenes2"]),
    Fonte_gsr_voto = unname(
      gsr_fonte_voto[seq_asv, "Species_binomial"]
    ),
    Valor_gsr07_auditoria = unname(
      gsr07_norm[seq_asv, "Species_binomial"]
    ),
    stringsAsFactors = FALSE
  )
}

detalhes_species_df <- do.call(rbind, detalhes_species)
salvar_csv(
  detalhes_species_df,
  file.path(dir_aud, "detalhes_classificacao_species.csv")
)
salvar_csv(
  detalhes_species_df,
  file.path(dir_aud, "detalhes_consenso_species.csv")
)

# Auditoria especifica das Species promovidas pelo BLAST exato.
salvar_csv(
  detalhes_species_df[
    detalhes_species_df$BLAST_exato_prioritario,
    ,
    drop = FALSE
  ],
  file.path(dir_aud, "species_blast_exata_prioridade_maxima.csv")
)

# Validacao independente da prioridade em Species.
for (i in seq_len(nrow(detalhes_species_df))) {
  fontes_validas <- strsplit(
    detalhes_species_df$Bancos_validos[i], ";", fixed = TRUE
  )[[1L]]
  fontes_validas <- fontes_validas[nzchar(fontes_validas)]
  esperado <- if (length(fontes_validas)) fontes_validas[1L] else NA_character_
  obtido <- detalhes_species_df$Banco_decisor[i]
  ok <- (is.na(esperado) && is.na(obtido)) || identical(esperado, obtido)
  if (!ok) {
    arq <- file.path(dir_aud, "ERRO_validacao_prioridade_species.csv")
    salvar_csv(detalhes_species_df[i, , drop = FALSE], arq)
    abort("Validacao da prioridade em Species falhou. Consulte %s.", arq)
  }
}

# Validacao de coerencia Genus-Species final.
idx_sp_final <- which(!is.na(taxa_integrada_todas[, "Species"]))
incoerentes_sp <- idx_sp_final[!vapply(
  idx_sp_final,
  function(i) {
    seq_asv <- rownames(taxa_integrada_todas)[i]
    mesmo_taxon(
      taxa_integrada_todas[seq_asv, "Genus"],
      extrair_genus_do_binomio(taxa_integrada_todas[seq_asv, "Species"])
    )
  },
  logical(1)
)]
if (length(incoerentes_sp) > 0L) {
  seqs_bad <- rownames(taxa_integrada_todas)[incoerentes_sp]
  arq <- file.path(dir_aud, "ERRO_incoerencia_genus_species_final.csv")
  salvar_csv(
    data.frame(
      ASV_ID = unname(seq2id[seqs_bad]),
      ASV_seq = seqs_bad,
      Genus = taxa_integrada_todas[seqs_bad, "Genus"],
      Species = taxa_integrada_todas[seqs_bad, "Species"],
      stringsAsFactors = FALSE
    ),
    arq
  )
  abort("Genus e Species finais incoerentes. Consulte %s.", arq)
}

log_msg(
  sprintf(
    paste0(
      "Species concluida: %d ASVs por BLAST exato prioritario; ",
      "%d ajustes de Genus para coerencia."
    ),
    sum(detalhes_species_df$BLAST_exato_prioritario),
    sum(genus_ajustado_por_blast)
  ),
  "OK"
)

###############################################################################
# 10. CONTAMINANTES
###############################################################################

eh_contaminante <- function(mat) {
  out <- rep(FALSE, nrow(mat))
  names(out) <- rownames(mat)
  out <- out |
    (!is.na(mat[, "Kingdom"]) & mat[, "Kingdom"] %in% CONT_KINGDOM) |
    (!is.na(mat[, "Order"])   & mat[, "Order"]   %in% CONT_ORDER)   |
    (!is.na(mat[, "Family"])  & mat[, "Family"]  %in% CONT_FAMILY)
  out
}

# Remocao contratual baseada na classificacao integrada final.
cont_flag <- eh_contaminante(taxa_integrada_todas)
seqs_contaminantes <- names(cont_flag)[cont_flag]
seqs_analise <- setdiff(all_asvs, seqs_contaminantes)

motivo_cont <- vapply(seqs_contaminantes, function(s) {
  motivos <- character()
  if (valor_valido(taxa_integrada_todas[s, "Kingdom"]) &&
      taxa_integrada_todas[s, "Kingdom"] %in% CONT_KINGDOM) {
    motivos <- c(motivos, paste0("Kingdom=", taxa_integrada_todas[s, "Kingdom"]))
  }
  if (valor_valido(taxa_integrada_todas[s, "Order"]) &&
      taxa_integrada_todas[s, "Order"] %in% CONT_ORDER) {
    motivos <- c(motivos, paste0("Order=", taxa_integrada_todas[s, "Order"]))
  }
  if (valor_valido(taxa_integrada_todas[s, "Family"]) &&
      taxa_integrada_todas[s, "Family"] %in% CONT_FAMILY) {
    motivos <- c(motivos, paste0("Family=", taxa_integrada_todas[s, "Family"]))
  }
  paste(motivos, collapse = ";")
}, character(1))

asvs_contaminantes_excluir <- data.frame(
  ASV_ID = unname(seq2id[seqs_contaminantes]),
  ASV_seq = seqs_contaminantes,
  Motivo = if (length(seqs_contaminantes) == 0L) character(0) else motivo_cont,
  stringsAsFactors = FALSE
)

# Auditoria conservadora: registra qualquer fonte que tenha sinalizado
# Eukaryota, Chloroplast ou Mitochondria, mesmo quando a prioridade final escolhe
# outra classificacao. Esta tabela nao altera automaticamente a lista de remocao.
aud_cont <- do.call(rbind, lapply(ORDEM_BANCOS, function(banco) {
  mat <- BANCOS[[banco]][, NIVEIS, drop = FALSE]
  flag <- eh_contaminante(mat)
  if (!any(flag)) return(NULL)
  data.frame(
    ASV_ID = unname(seq2id[rownames(mat)[flag]]),
    ASV_seq = rownames(mat)[flag],
    Banco = banco,
    Kingdom = mat[flag, "Kingdom"],
    Order = mat[flag, "Order"],
    Family = mat[flag, "Family"],
    Removida_pela_classificacao_final = rownames(mat)[flag] %in% seqs_contaminantes,
    stringsAsFactors = FALSE
  )
}))
if (is.null(aud_cont)) {
  aud_cont <- data.frame(
    ASV_ID = character(), ASV_seq = character(), Banco = character(),
    Kingdom = character(), Order = character(), Family = character(),
    Removida_pela_classificacao_final = logical(),
    stringsAsFactors = FALSE
  )
}
salvar_csv(
  aud_cont,
  file.path(dir_aud, "auditoria_contaminantes_por_banco.csv")
)

taxa_final <- taxa_integrada_todas[seqs_analise, NIVEIS, drop = FALSE]
fonte_final <- fonte_todas[seqs_analise, NIVEIS, drop = FALSE]

# Contrato estrutural antes da escrita.
if (!identical(colnames(taxa_final), NIVEIS)) {
  abort("taxa_final nao possui os sete ranks canonicos na ordem esperada.")
}
if (anyDuplicated(rownames(taxa_final)) > 0L) {
  abort("taxa_final possui sequencias duplicadas nos rownames.")
}
if (!identical(rownames(taxa_final), seqs_analise)) {
  abort("taxa_final perdeu a ordem canonica apos remover contaminantes.")
}

# GSR 0.7 alinhado ao universo analitico para o Script 06.
if (gsr07_disponivel) {
  taxa_gsr07_analise <- gsr07_norm[seqs_analise, NIVEIS, drop = FALSE]
  taxa_gsr07_analise[, "Species"] <-
    gsr07_norm[seqs_analise, "Species_binomial"]

  if (!identical(rownames(taxa_gsr07_analise), rownames(taxa_final))) {
    abort("GSR 0.7 e taxa_final possuem universos ou ordens diferentes.")
  }

  salvar_rds(
    taxa_gsr07_analise,
    file.path(dir_gsr_sens, "taxa_gsr07_alinhada_analise.rds")
  )
  salvar_csv(
    com_asv_id(taxa_gsr07_analise, seq2id),
    file.path(dir_gsr_sens, "taxa_gsr07_alinhada_analise.csv")
  )

  sens_gsr <- do.call(rbind, lapply(NIVEIS, function(rk) {
    principal <- taxa_integrada_todas[all_asvs, rk]
    gsr07 <- if (rk == "Species") {
      gsr07_norm[all_asvs, "Species_binomial"]
    } else {
      gsr07_norm[all_asvs, rk]
    }
    data.frame(
      ASV_ID = unname(seq2id[all_asvs]),
      ASV_seq = all_asvs,
      Rank = rk,
      Classificacao_principal = principal,
      Consenso_principal = principal,
      GSR_07 = gsr07,
      Ambos_informativos = !is.na(principal) & !is.na(gsr07),
      Concordam = mapply(
        function(a, b) {
          if (!valor_valido(a) || !valor_valido(b)) return(NA)
          if (rk == "Species") {
            mesma_species_binomial(a, b)
          } else {
            mesmo_taxon(a, b)
          }
        },
        principal, gsr07,
        USE.NAMES = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }))
  salvar_csv(
    sens_gsr,
    file.path(dir_gsr_sens, "comparacao_classificacao_principal_vs_gsr07.csv")
  )
  # Nome legado preservado.
  salvar_csv(
    sens_gsr,
    file.path(dir_gsr_sens, "comparacao_consenso_principal_vs_gsr07.csv")
  )
}

# Arquivos contratuais — nomes legados preservados para os Scripts 06-12.
salvar_rds(taxa_final, file.path(output_path, "taxa_consenso_final.rds"))
salvar_csv(
  com_asv_id(taxa_final, seq2id),
  file.path(output_path, "taxa_consenso_final.csv")
)
salvar_csv(
  asvs_contaminantes_excluir,
  file.path(output_path, "asvs_contaminantes_excluir.csv")
)
salvar_rds(fonte_final, file.path(dir_comp, "fonte_consenso_por_nivel.rds"))
salvar_csv(
  com_asv_id(fonte_final, seq2id),
  file.path(dir_comp, "fonte_consenso_por_nivel.csv")
)
# Nomes metodologicamente explicitos, adicionais aos contratos legados.
salvar_rds(taxa_final, file.path(dir_comp, "taxa_integrada_prioridade_final.rds"))
salvar_csv(
  com_asv_id(taxa_final, seq2id),
  file.path(dir_comp, "taxa_integrada_prioridade_final.csv")
)

###############################################################################
# 11. TABELA UNICA DE CLASSIFICACAO E CONFERENCIA BLAST
###############################################################################

compor_coluna_banco <- function(banco_nome, rk) {
  BANCOS_AUDITORIA[[banco_nome]][all_asvs, rk]
}

compor_coluna_species_bin <- function(banco_nome) {
  BANCOS_AUDITORIA[[banco_nome]][all_asvs, "Species_binomial"]
}

comparar_opcional <- function(a, b, comparador = mesmo_taxon) {
  mapply(
    function(x, y) {
      if (!valor_valido(x) || !valor_valido(y)) return(NA)
      comparador(x, y)
    },
    a, b,
    USE.NAMES = FALSE
  )
}

tabela_final <- data.frame(
  ASV_ID = map_ord$ASV_ID,
  ASV_seq = all_asvs,
  Origem = map_ord$Origem,
  Reads_totais = unname(reads_por_seq[all_asvs]),
  N_amostras = unname(prev_por_seq[all_asvs]),
  Amostras = unname(amostras_por_seq[all_asvs]),
  Contaminante = all_asvs %in% seqs_contaminantes,
  stringsAsFactors = FALSE
)

for (rk in setdiff(RANKS_CLASSIFICACAO, "Species")) {
  tabela_final[[paste0(rk, "_silva138")]] <- compor_coluna_banco("silva138", rk)
  tabela_final[[paste0(rk, "_rdp19")]] <- compor_coluna_banco("rdp19", rk)
  tabela_final[[paste0(rk, "_greengenes2")]] <- compor_coluna_banco("greengenes2", rk)
  tabela_final[[paste0(rk, "_beexact")]] <- compor_coluna_banco("beexact", rk)
  tabela_final[[paste0(rk, "_gsr_disable")]] <- compor_coluna_banco("gsr", rk)
  tabela_final[[paste0(rk, "_gsr07")]] <- compor_coluna_banco("gsr_07", rk)
  tabela_final[[paste0(rk, "_gsr_voto")]] <- compor_coluna_banco("gsr", rk)
  tabela_final[[paste0(rk, "_gsr_fonte_voto")]] <- gsr_fonte_voto[all_asvs, rk]

  # Campos explicitos.
  tabela_final[[paste0(rk, "_final")]] <- taxa_integrada_todas[all_asvs, rk]
  tabela_final[[paste0(rk, "_banco_decisor")]] <- fonte_todas[all_asvs, rk]
  # Aliases legados.
  tabela_final[[paste0(rk, "_consenso")]] <- taxa_integrada_todas[all_asvs, rk]
  tabela_final[[paste0(rk, "_fonte")]] <- fonte_todas[all_asvs, rk]
  tabela_final[[paste0(rk, "_n_concordam")]] <- n_concordam_todas[all_asvs, rk]

  det <- detalhes_por_rank[[rk]]
  det <- det[match(all_asvs, det$ASV_seq), , drop = FALSE]
  tabela_final[[paste0(rk, "_n_bancos_com_valor")]] <- det$N_bancos_com_valor
  tabela_final[[paste0(rk, "_n_informativos")]] <- det$N_bancos_informativos
  tabela_final[[paste0(rk, "_apoio_maximo")]] <- det$Apoio_maximo
  tabela_final[[paste0(rk, "_proporcao_apoio")]] <- det$Proporcao_apoio
  tabela_final[[paste0(rk, "_houve_conflito")]] <- det$Houve_conflito
  tabela_final[[paste0(rk, "_empate_real")]] <- FALSE
  tabela_final[[paste0(rk, "_bancos_concordam")]] <- det$Bancos_concordam
  tabela_final[[paste0(rk, "_bancos_discordam")]] <- det$Bancos_discordam
  tabela_final[[paste0(rk, "_bancos_incompativeis")]] <- det$Bancos_incompativeis
  tabela_final[[paste0(rk, "_bancos_sem_valor")]] <- det$Bancos_sem_valor
  tabela_final[[paste0(rk, "_regra_aplicada")]] <- det$Regra_aplicada
  tabela_final[[paste0(rk, "_status")]] <- det$Status
}

# Species em formato binomial.
tabela_final[["Species_silva138"]] <- compor_coluna_species_bin("silva138")
tabela_final[["Species_rdp19"]] <- compor_coluna_species_bin("rdp19")
tabela_final[["Species_greengenes2"]] <- compor_coluna_species_bin("greengenes2")
tabela_final[["Species_beexact"]] <- compor_coluna_species_bin("beexact")
tabela_final[["Species_gsr_disable"]] <- compor_coluna_species_bin("gsr")
tabela_final[["Species_gsr07"]] <- compor_coluna_species_bin("gsr_07")
tabela_final[["Species_gsr_voto"]] <- compor_coluna_species_bin("gsr")
tabela_final[["Species_gsr_fonte_voto"]] <- gsr_fonte_voto[all_asvs, "Species_binomial"]
tabela_final[["Species_silva_addSpecies"]] <- silva_species_exact_bin[all_asvs]

sp_det <- detalhes_species_df[
  match(all_asvs, detalhes_species_df$ASV_seq),
  ,
  drop = FALSE
]

tabela_final[["Species_blast100_exact"]] <- blast_species_exata[all_asvs]
tabela_final[["Species_blast100_valida"]] <- blast_species_valida[all_asvs]
tabela_final[["Species_blast100_prioritaria"]] <-
  sp_det$BLAST_exato_prioritario
tabela_final[["Species_classificadores_sem_BLAST"]] <-
  sp_det$Species_classificadores_sem_BLAST
tabela_final[["Species_fonte_classificadores_sem_BLAST"]] <-
  sp_det$Fonte_classificadores_sem_BLAST
tabela_final[["Genus_antes_BLAST_species"]] <- genus_antes_blast[all_asvs]
tabela_final[["Genus_fonte_antes_BLAST_species"]] <-
  fonte_genus_antes_blast[all_asvs]
tabela_final[["Genus_ajustado_por_BLAST_species"]] <-
  genus_ajustado_por_blast[all_asvs]
tabela_final[["Species_final"]] <- taxa_integrada_todas[all_asvs, "Species"]
tabela_final[["Species_banco_decisor"]] <- fonte_todas[all_asvs, "Species"]
tabela_final[["Species_consenso"]] <- taxa_integrada_todas[all_asvs, "Species"]
tabela_final[["Species_fonte"]] <- fonte_todas[all_asvs, "Species"]
tabela_final[["Species_n_concordam"]] <- n_concordam_todas[all_asvs, "Species"]

tabela_final[["Species_n_bancos_com_valor"]] <- sp_det$N_bancos_com_valor
tabela_final[["Species_n_informativos"]] <- sp_det$N_bancos_informativos
tabela_final[["Species_apoio_maximo"]] <- sp_det$Apoio_maximo
tabela_final[["Species_proporcao_apoio"]] <- sp_det$Proporcao_apoio
tabela_final[["Species_houve_conflito"]] <- sp_det$Houve_conflito
tabela_final[["Species_empate_real"]] <- FALSE
tabela_final[["Species_bancos_concordam"]] <- sp_det$Bancos_concordam
tabela_final[["Species_bancos_discordam"]] <- sp_det$Bancos_discordam
tabela_final[["Species_bancos_invalidos"]] <- sp_det$Bancos_invalidos
tabela_final[["Species_bancos_sem_valor"]] <- sp_det$Bancos_sem_valor
tabela_final[["Species_regra_aplicada"]] <- sp_det$Regra_aplicada
tabela_final[["Species_status"]] <- sp_det$Status

# BLAST — Species exata participa da decisao; Genus candidato 97% permanece auditoria.
tabela_final[["BLAST_Genus_exato"]] <- blast100_norm[all_asvs, "Genus"]
tabela_final[["BLAST_Species_exata_bin"]] <- blast100_norm[all_asvs, "Species_binomial"]
tabela_final[["BLAST_Genus_candidato_97"]] <- blast97_norm[all_asvs, "Genus"]
tabela_final[["BLAST_identidade_%"]] <- blast_evid$Melhor_identidade_genus
tabela_final[["BLAST_cobertura_%"]] <- blast_evid$Melhor_cobertura_genus
tabela_final[["BLAST_ambiguo_Genus"]] <- blast_evid$Ambiguo_Genus_exato
tabela_final[["BLAST_ambiguo_Species"]] <- blast_evid$Ambiguo_Species_exata
tabela_final[["BLAST_busca_saturada"]] <- blast_evid$Atingiu_MAX_TARGET_SEQS

tabela_final[["BLAST_confirma_Genus_consenso"]] <- comparar_opcional(
  taxa_integrada_todas[all_asvs, "Genus"],
  blast100_norm[all_asvs, "Genus"],
  mesmo_taxon
)
tabela_final[["BLAST_confirma_Species_consenso"]] <- comparar_opcional(
  taxa_integrada_todas[all_asvs, "Species"],
  blast100_norm[all_asvs, "Species_binomial"],
  mesma_species_binomial
)

###############################################################################
# 12. AUDITORIAS COMPLEMENTARES
###############################################################################

comparar_par <- function(mat_a, mat_b, nome_a, nome_b) {
  do.call(rbind, lapply(NIVEIS, function(rk) {
    coluna <- if (rk == "Species") "Species_binomial" else rk
    a <- mat_a[all_asvs, coluna]
    b <- mat_b[all_asvs, coluna]
    ambos <- !is.na(a) & !is.na(b) & trimws(a) != "" & trimws(b) != ""
    comparador <- if (rk == "Species") mesma_species_binomial else mesmo_taxon
    conc <- ambos & mapply(comparador, a, b, USE.NAMES = FALSE)
    data.frame(
      Banco_A = nome_a,
      Banco_B = nome_b,
      Rank = rk,
      N_ambos = sum(ambos),
      N_concordam = sum(conc),
      Pct_concordancia = if (sum(ambos) > 0L) {
        round(100 * sum(conc) / sum(ambos), 2)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
}

pares_bancos <- combn(ORDEM_BANCOS, 2L, simplify = FALSE)
conc_bancos <- do.call(rbind, lapply(pares_bancos, function(p) {
  comparar_par(BANCOS[[p[1L]]], BANCOS[[p[2L]]], p[1L], p[2L])
}))
salvar_csv(
  conc_bancos,
  file.path(dir_comp, "concordancia_par_a_par_bancos_votantes.csv")
)
salvar_csv(
  conc_bancos,
  file.path(dir_comp, "concordancia_par_a_par_bancos.csv")
)

if (gsr07_disponivel) {
  conc_gsr07 <- do.call(rbind, lapply(ORDEM_BANCOS, function(b) {
    comparar_par(BANCOS[[b]], gsr07_norm, b, "gsr_07_sensibilidade")
  }))
  salvar_csv(
    conc_gsr07,
    file.path(dir_gsr_sens, "concordancia_gsr07_vs_bancos.csv")
  )
}

# ASVs decididas por cada fonte.
for (banco in ORDEM_BANCOS) {
  idx <- rowSums(fonte_final == banco, na.rm = TRUE) > 0L
  if (!any(idx)) next
  linhas <- data.frame(
    ASV_ID = unname(seq2id[rownames(fonte_final)[idx]]),
    ASV_seq = rownames(fonte_final)[idx],
    Kingdom = taxa_final[idx, "Kingdom"],
    Phylum = taxa_final[idx, "Phylum"],
    Class = taxa_final[idx, "Class"],
    Order = taxa_final[idx, "Order"],
    Family = taxa_final[idx, "Family"],
    Genus = taxa_final[idx, "Genus"],
    Species = taxa_final[idx, "Species"],
    Ranks_decididos_por_este_banco = apply(
      fonte_final[idx, , drop = FALSE],
      1L,
      function(r) paste(NIVEIS[!is.na(r) & r == banco], collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  salvar_csv(
    linhas,
    file.path(dir_comp, paste0("asvs_resolvidas_por_", banco, ".csv"))
  )
}

# Auditorias por status. As tabelas de Species possuem colunas adicionais de
# validade binomial; o preenchimento por nome evita rbind posicional ou erro por
# conjuntos de colunas diferentes.
rbind_fill <- function(...) {
  dfs <- list(...)
  todos_nomes <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    faltam <- setdiff(todos_nomes, names(df))
    for (nm in faltam) df[[nm]] <- NA
    df[, todos_nomes, drop = FALSE]
  })
  do.call(rbind, dfs)
}

detalhes_todos <- rbind_fill(detalhes_kingdom_genus, detalhes_species_df)
status_fonte_unica <- detalhes_todos$Status %in% c(
  "fonte_unica_aceita",
  "fonte_unica_aceita_apos_filtro_hierarquico"
)
status_conflito <- detalhes_todos$Status == "prioridade_hierarquica_por_conflito"

salvar_csv(
  detalhes_todos[status_fonte_unica, , drop = FALSE],
  file.path(dir_aud, "classificacoes_fonte_unica_ACEITAS.csv")
)
salvar_csv(
  detalhes_species_df[detalhes_species_df$Status %in% c(
    "fonte_unica_aceita", "fonte_unica_aceita_apos_filtro_hierarquico"
  ), , drop = FALSE],
  file.path(dir_aud, "species_fonte_unica_ACEITA.csv")
)
salvar_csv(
  detalhes_todos[status_conflito, , drop = FALSE],
  file.path(dir_aud, "conflitos_resolvidos_por_prioridade.csv")
)

# Remover arquivos legados que contradizem a regra atual, evitando interpretacao
# de resultados antigos como se fossem produzidos por esta execucao.
legados_obsoletos <- c(
  file.path(dir_aud, "species_fonte_unica_REJEITADA.csv"),
  file.path(dir_aud, "species_empate_real_sem_consenso.csv"),
  file.path(dir_aud, "ERRO_pares_concordantes_nao_reconhecidos.csv"),
  file.path(dir_aud, "ERRO_pares_species_nao_reconhecidos.csv")
)
unlink(legados_obsoletos[file.exists(legados_obsoletos)], force = TRUE)

nivel_mais_profundo <- function(row_vals) {
  valido <- !is.na(row_vals) & trimws(row_vals) != ""
  if (!any(valido)) return("Nao_identificada")
  tail(names(row_vals)[valido], 1L)
}

tabela_final[["Nivel_atribuicao"]] <- apply(
  taxa_integrada_todas[all_asvs, NIVEIS, drop = FALSE],
  1L,
  nivel_mais_profundo
)

# FASTA das ASVs sem classificacao ate Order.
idx_min <- match(NIVEL_MINIMO_PRINCIPAL, NIVEIS)
idx_nivel <- match(tabela_final$Nivel_atribuicao, NIVEIS)
nao_id <- (is.na(idx_nivel) | idx_nivel < idx_min) & !tabela_final$Contaminante

if (any(nao_id)) {
  dna_nao <- Biostrings::DNAStringSet(tabela_final$ASV_seq[nao_id])
  names(dna_nao) <- tabela_final$ASV_ID[nao_id]
  Biostrings::writeXStringSet(
    dna_nao,
    filepath = file.path(dir_fa, "nao_identificadas.fasta"),
    format = "fasta",
    width = 80L
  )
  salvar_csv(
    tabela_final[nao_id, , drop = FALSE],
    file.path(output_path, "tabela_nao_identificadas.csv")
  )
} else {
  arq_fasta_antigo <- file.path(dir_fa, "nao_identificadas.fasta")
  if (file.exists(arq_fasta_antigo)) unlink(arq_fasta_antigo, force = TRUE)
  salvar_csv(
    tabela_final[FALSE, , drop = FALSE],
    file.path(output_path, "tabela_nao_identificadas.csv")
  )
}

salvar_csv(
  tabela_final,
  file.path(dir_comp, "tabela_concordancia_final.csv")
)
log_msg("tabela_concordancia_final.csv salva com campos de prioridade.", "OK")

###############################################################################
# 13. FIGURAS E COBERTURA
###############################################################################

fonte_genus <- table(fonte_final[, "Genus"], useNA = "no")
if (length(fonte_genus) > 0L) {
  df_fonte <- data.frame(
    Banco_decisor = names(fonte_genus),
    N_ASVs = as.integer(fonte_genus),
    stringsAsFactors = FALSE
  )
  p_fonte <- ggplot(
    df_fonte,
    aes(x = reorder(Banco_decisor, N_ASVs), y = N_ASVs)
  ) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(
      title = "Banco decisor do Genus na classificacao integrada",
      x = NULL,
      y = "Numero de ASVs"
    )
  ggsave(
    file.path(dir_fig, "fig_fonte_genus.png"),
    p_fonte,
    width = 9,
    height = 5,
    dpi = 200
  )
}

status_por_rank <- c(detalhes_por_rank, list(Species = detalhes_species_df))
df_cobertura <- do.call(rbind, lapply(NIVEIS, function(rk) {
  det <- status_por_rank[[rk]]
  data.frame(
    Rank = rk,
    N_classificadas = sum(!is.na(taxa_integrada_todas[, rk])),
    N_fonte_unica = sum(det$Status %in% c(
      "fonte_unica_aceita", "fonte_unica_aceita_apos_filtro_hierarquico"
    )),
    N_concordancia_entre_bancos = sum(det$Status == "concordancia_entre_bancos"),
    N_conflito_resolvido_prioridade = sum(
      det$Status == "prioridade_hierarquica_por_conflito"
    ),
    N_blast_exato_prioritario = sum(
      det$Status %in% c(
        "blast_exato_prioridade_maxima",
        "blast_species_exata_ajuste_genus"
      )
    ),
    N_sem_classificacao = sum(grepl("^sem_classificacao", det$Status)),
    stringsAsFactors = FALSE
  )
}))

salvar_csv(
  df_cobertura,
  file.path(dir_comp, "cobertura_classificacao_por_rank.csv")
)
# Nome legado preservado.
salvar_csv(
  df_cobertura,
  file.path(dir_comp, "cobertura_consenso_por_rank.csv")
)

p_cob <- ggplot(
  df_cobertura,
  aes(x = factor(Rank, levels = NIVEIS), y = N_classificadas)
) +
  geom_col() +
  theme_classic(base_size = 11) +
  labs(
    title = "ASVs classificadas por rank",
    x = NULL,
    y = "Numero de ASVs"
  )
ggsave(
  file.path(dir_fig, "fig_cobertura_por_rank.png"),
  p_cob,
  width = 9,
  height = 5,
  dpi = 200
)

###############################################################################
# 14. METADADOS, CHECKPOINT E VALIDACAO DOS CONTRATOS
###############################################################################

descricao_regra <- paste(
  "Classificacao integrada rank a rank por prioridade hierarquica:",
  "GSR confidence=disable > BEExact > SILVA 138.2 > RDP 19 > Greengenes2.",
  "Fonte unica e aceita.",
  "Quando dois ou mais bancos discordam, vence o banco valido de maior prioridade.",
  "Quando todos os bancos validos concordam, o taxon comum e mantido.",
  "Compatibilidade com a linhagem superior e obrigatoria.",
  "Species BLAST exata, integral e nao ambigua tem prioridade maxima.",
  "Quando necessario, o Genus e ajustado para manter coerencia com a Species BLAST.",
  "Sem BLAST exato, Species segue a hierarquia dos cinco classificadores.",
  "GSR 0.7 e sensibilidade; BLAST Genus >=97% permanece conferencia."
)

metadata_consenso <- data.frame(
  Script = "05_comparacao_banco_de_dados.R",
  Versao = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Descricao_regra = descricao_regra,
  Descricao_hierarquia = descricao_regra,
  Unidade_analitica = "ASV",
  Nome_metodologico = "classificacao_taxonomica_integrada_por_prioridade",
  Arquivo_contratual_legado = "taxa_consenso_final.rds",
  Bancos_votantes = paste(ORDEM_BANCOS, collapse = ";"),
  Ordem_prioridade = paste(ORDEM_BANCOS, collapse = " > "),
  Ordem_prioridade_Species = paste(ORDEM_SPECIES, collapse = " > "),
  Numero_bancos_votantes = length(ORDEM_BANCOS),
  Numero_fontes_Species = length(ORDEM_SPECIES),
  Minimo_concordancia = 1L,
  Fonte_unica_aceita = TRUE,
  Regra_conflito = "prioridade_hierarquica",
  Regra_empate_real = "nao_aplicavel; conflito resolvido por prioridade",
  GSR_papel = "primeira_prioridade_confidence_disable",
  GSR_dependencia = paste(
    "GSR integra referencias relacionadas a SILVA, RDP e Greengenes;",
    "a hierarquia e uma regra decisoria pre-especificada, nao evidencia independente."
  ),
  GSR_disable_arquivo = arq_gsr_disable,
  GSR07_arquivo = arq_gsr07,
  GSR07_saida_script6 = file.path(
    dir_gsr_sens, "taxa_gsr07_alinhada_analise.rds"
  ),
  BLAST_papel = "Species_exata_prioridade_maxima; Genus_97pct_auditoria",
  BLAST_pasta = blast_root,
  ASVs_totais = n_asvs_total,
  ASVs_contaminantes = length(seqs_contaminantes),
  ASVs_classificacao_final = nrow(taxa_final),
  ASVs_consenso_final = nrow(taxa_final),
  Genus_classificado = sum(!is.na(taxa_final[, "Genus"])),
  Species_classificada = sum(!is.na(taxa_final[, "Species"])),
  Species_BLAST_exata_prioritaria = sum(
    detalhes_species_df$BLAST_exato_prioritario
  ),
  Genus_ajustado_por_BLAST_species = sum(genus_ajustado_por_blast),
  Classificacoes_fonte_unica_aceitas = sum(status_fonte_unica),
  Species_fonte_unica_aceita = sum(
    detalhes_species_df$Status %in% c(
      "fonte_unica_aceita", "fonte_unica_aceita_apos_filtro_hierarquico"
    )
  ),
  Conflitos_resolvidos_por_prioridade = sum(status_conflito),
  Validacao_declarada = paste(
    "Integracao taxonomica auditavel; nao substitui benchmark externo",
    "com sequencias V3-V4 de taxonomia conhecida."
  ),
  stringsAsFactors = FALSE
)

salvar_csv(
  metadata_consenso,
  file.path(output_path, "metadata_consenso_taxonomico.csv")
)

salvar_rds(
  list(
    metadata = metadata_consenso,
    taxa_consenso_final = taxa_final,
    taxa_integrada_prioridade = taxa_final,
    fonte = fonte_final,
    n_apoio = n_concordam_todas,
    tabela_final = tabela_final,
    concordancia_par_a_par = conc_bancos
  ),
  file.path(dir_comp, "checkpoint_05_consenso_concluido.rds")
)

# Validacao final dos tres contratos consumidos pelo Script 06.
arq_taxa_saida <- file.path(output_path, "taxa_consenso_final.rds")
arq_cont_saida <- file.path(output_path, "asvs_contaminantes_excluir.csv")
arq_meta_saida <- file.path(output_path, "metadata_consenso_taxonomico.csv")
for (arq in c(arq_taxa_saida, arq_cont_saida, arq_meta_saida)) {
  if (!file.exists(arq) || is.na(file.size(arq)) || file.size(arq) == 0L) {
    abort("Contrato do Script 06 ausente ou vazio: %s", arq)
  }
}

taxa_check <- as.matrix(readRDS(arq_taxa_saida))
if (
  !identical(colnames(taxa_check), NIVEIS) ||
  !identical(rownames(taxa_check), seqs_analise)
) {
  abort("Contrato taxa_consenso_final.rds invalido apos releitura.")
}

cont_check <- read.csv(arq_cont_saida, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("ASV_ID", "ASV_seq", "Motivo") %in% colnames(cont_check))) {
  abort("Contrato asvs_contaminantes_excluir.csv invalido apos releitura.")
}

meta_check <- read.csv(arq_meta_saida, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(meta_check) != 1L || !"Descricao_hierarquia" %in% colnames(meta_check)) {
  abort("Contrato metadata_consenso_taxonomico.csv invalido apos releitura.")
}

cat("\n=============================================================\n")
cat("CLASSIFICACAO TAXONOMICA POR PRIORIDADE — CONCLUIDA\n")
cat("Prioridade: ", paste(ORDEM_BANCOS, collapse = " > "), "\n", sep = "")
cat(sprintf("ASVs totais: %d\n", n_asvs_total))
cat(sprintf("Contaminantes: %d\n", length(seqs_contaminantes)))
cat(sprintf("ASVs na classificacao final: %d\n", nrow(taxa_final)))
cat(sprintf(
  "Genus: %d | Species: %d\n",
  sum(!is.na(taxa_final[, "Genus"])),
  sum(!is.na(taxa_final[, "Species"]))
))
cat(sprintf(
  "Classificacoes por fonte unica aceitas: %d\n",
  sum(status_fonte_unica)
))
cat(sprintf(
  "Conflitos resolvidos por prioridade entre classificadores: %d\n",
  sum(status_conflito)
))
cat("Tabela de auditoria: ", file.path(dir_comp, "tabela_concordancia_final.csv"), "\n", sep = "")
cat("Contratos para o Script 06:\n")
cat("  ", arq_taxa_saida, "\n", sep = "")
cat("  ", arq_cont_saida, "\n", sep = "")
cat("  ", arq_meta_saida, "\n", sep = "")
if (gsr07_disponivel) {
  cat(
    "  ",
    file.path(dir_gsr_sens, "taxa_gsr07_alinhada_analise.rds"),
    "\n",
    sep = ""
  )
}
cat("=============================================================\n\n")

log_msg("Script 05 finalizado.", "FINAL")
})
