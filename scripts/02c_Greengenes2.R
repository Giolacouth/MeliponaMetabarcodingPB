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
run_pipeline_script("02c_Greengenes2.R", "greengenes2", function(ctx) {
###############################################################################
#   SCRIPT 02c — CLASSIFICACAO TAXONOMICA — GREENGENES2 2024.09
#   Executar apos o Script 02b (RDP) e antes do Script 04 (BLAST)
#
#   CARACTERISTICAS:
#   - Banco mais recente (2024.09), integra genomas (Web of Life + GTDB + LTP)
#   - Trainset formatado para DADA2 (assignTaxonomy / Naive Bayes), mantido
#     oficialmente pela equipe do DADA2 (Callahan). Comparacao METODOLOGICAMENTE
#     JUSTA com SILVA e RDP: mesmo classificador, mesma regiao (full-length
#     recortada pelo kmer).
#
#   AVISO METODOLOGICO CRITICO — NOMENCLATURA GTDB:
#   - GG2 herda a taxonomia do GTDB, que renomeia filos e generos em relacao
#     ao SILVA/RDP (ICNP). Exemplos:
#       Firmicutes      -> Bacillota
#       Proteobacteria  -> Pseudomonadota
#       Actinobacteria  -> Actinomycetota
#       Bacteroidetes   -> Bacteroidota
#     Generos tambem podem ser divididos/fundidos (ex.: Lactobacillus s.l.).
#   - No Script 5, GG2 integra a evidência classificatória independente de GSR.
#   - O valor de GG2 é harmonizado antes da comparação por rank.
#   - Diferenças de nomenclatura GTDB permanecem registradas na auditoria.
#
#   NIVEL DE ESPECIE:
#   - O arquivo toGenus (recomendado) classifica ate Genus. O toSpecies vai
#     ate Species, mas a documentacao alerta que nao se deve esperar atribuicao
#     confiavel a especie. Este script usa toGenus por padrao; se apontado para
#     toSpecies, a coluna Species e marcada como kmer_gtdb (baixa confianca).
#
#   ENTRADA: seqtab_global_nochim.rds  (mesmo universo dos demais bancos)
#   ARQUIVOS GG2 (em bancodados/greengenes):
#     gg2_2024_09_toGenus_trainset.fa.gz    (PADRAO — recomendado)
#     gg2_2024_09_toSpecies_trainset.fa.gz  (opcional)
###############################################################################

library(dada2)

VERSAO        <- "1"
DATA_EXECUCAO <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
BANCO_NOME    <- "gg2"

cat("=============================================================\n")
cat("GREENGENES2 2024.09 — Classificacao Taxonomica v", VERSAO, "\n", sep = "")
cat("Data/Hora:", DATA_EXECUCAO, "\n")
cat("=============================================================\n\n")

base_path <- ctx$base_path
pipeline_version <- ctx$version
output_path <- ctx$stage$root
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_path)) stop("Falha ao criar diretorio: ", output_path, call. = FALSE)

# ---- ENTRADA CANONICA: identica para todos os bancos kmer ----
arq_seqtab <- ctx$contracts[["seqtab_global_nochim"]]
# Mapa canonico ASV_ID <-> sequencia (Script 01); usado para identificar
# cada ASV por numero nos CSVs deste script (ver com_asv_id()).
arq_asvmap <- ctx$contracts[["asv_sequences"]]

# ---- ESCOLHA DO TRAINSET GG2 ----
# Padrao: toGenus (recomendado para microbioma de abelhas).
# Para usar toSpecies, trocar USAR_TOSPECIES para TRUE.
USAR_TOSPECIES <- FALSE

dir_gg2     <- file.path(base_path, "bancodados/greengenes")
arq_gg2_gen <- file.path(dir_gg2, "gg2_2024_09_toGenus_trainset.fa.gz")
arq_gg2_spc <- file.path(dir_gg2, "gg2_2024_09_toSpecies_trainset.fa.gz")
arq_gg2     <- if (USAR_TOSPECIES) arq_gg2_spc else arq_gg2_gen

MINBOOT     <- 80
MULTITHREAD <- TRUE
CHECKPOINT  <- TRUE

# ---------------------------------------------------------------------------
# FUNCOES AUXILIARES (identicas as dos Scripts 2 e 3 para consistencia)
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
  log_msg(paste0(desc, " OK (", round(file.size(arq)/1024/1024, 2), " MB)"), "OK")
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

# Carrega o mapa canonico ASV_ID <-> sequencia (mesma logica dos Scripts 02/03/04).
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

# Normalizacao de colunas — IDENTICA a dos outros bancos, com uma diferenca:
# o trainset GG2 do DADA2 ja vem SEM prefixos (k__, p__...) porque a funcao
# makeTaxonomyFasta_GG2 os remove. Mantemos a capitalizacao e o tratamento de
# Domain->Kingdom e remocao de Root por seguranca.
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
cat(sprintf("  MINBOOT=%d  MULTITHREAD=%s  TRAINSET=%s\n\n",
            MINBOOT, MULTITHREAD,
            if (USAR_TOSPECIES) "toSpecies" else "toGenus"))
cat("NOTA METODOLOGICA: GG2 usa taxonomia GTDB.\n")
cat("      No Script 5 ocupa a terceira prioridade e preenche apenas lacunas\n")
cat("      compativeis apos SILVA e RDP.\n\n")

cat("=== VALIDACOES ===\n\n")
validar_arquivo(arq_seqtab, "seqtab_global_nochim.rds", TRUE)
validar_arquivo(arq_gg2,
                paste0("GG2 trainset (",
                       if (USAR_TOSPECIES) "toSpecies" else "toGenus", ")"),
                TRUE)
validar_arquivo(arq_asvmap, "ASV_sequences.tsv", TRUE)
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
# ATRIBUICAO TAXONOMICA (kmer Naive Bayes)
# ==============================================================================
cat("=== assignTaxonomy (kmer Naive Bayes, GG2/GTDB) ===\n\n")
log_msg("Iniciando assignTaxonomy...", "INICIO")
t0 <- proc.time()

# GG2 trainset segue a ordem padrao de 6 ou 7 niveis do DADA2:
# Kingdom, Phylum, Class, Order, Family, Genus [, Species].
# Nao e necessario taxLevels customizado (a formatacao makeTaxonomyFasta_GG2
# usa os mesmos rotulos do default do DADA2).
taxa_boot_obj <- tryCatch(
  assignTaxonomy(
    seqtab, arq_gg2,
    multithread      = MULTITHREAD,
    minBoot          = MINBOOT,
    tryRC            = TRUE,
    outputBootstraps = TRUE),
  error = function(e)
    stop("GG2 assignTaxonomy falhou: ", conditionMessage(e), call. = FALSE)
)
if (!is.list(taxa_boot_obj) || !all(c("tax", "boot") %in% names(taxa_boot_obj)))
  stop("GG2 assignTaxonomy nao retornou componentes tax/boot.")

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
    "GG2 retornou taxonomia/bootstrap desalinhados do universo canonico.",
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
# NOTA SOBRE NOMENCLATURA GTDB (diagnostico informativo)
# ==============================================================================
# Detecta filos com nomenclatura GTDB para alertar o usuario de que a
# comparacao com SILVA/RDP exige harmonizacao de nomes (feita no Script 5).
if ("Phylum" %in% colnames(taxa)) {
  filos_gtdb <- c("Bacillota", "Pseudomonadota", "Actinomycetota",
                  "Bacteroidota", "Campylobacterota", "Fusobacteriota")
  phy_vals <- taxa[, "Phylum"]
  n_gtdb <- sum(!is.na(phy_vals) & phy_vals %in% filos_gtdb)
  if (n_gtdb > 0)
    log_msg(paste0(
      n_gtdb, " ASVs com filo em nomenclatura GTDB (ex.: Bacillota=Firmicutes). ",
      "Harmonizacao de nomes ocorre no Script 5 para comparacao com SILVA/RDP."),
      "WARN")
}

# Aviso especifico para microbioma de abelhas (Lactobacillus s.l. no GTDB)
if ("Genus" %in% colnames(taxa)) {
  gen_values <- taxa[, "Genus"]
  n_lacto <- sum(!is.na(gen_values) &
                   grepl("Lactobacillus|Apilactobacillus|Bombilactobacillus",
                         gen_values))
  if (n_lacto > 0)
    log_msg(paste0(
      n_lacto, " ASVs em generos de bacterias laticas (Lactobacillaceae). ",
      "GTDB pode usar nomenclatura revisada (ex.: Apilactobacillus, ",
      "Bombilactobacillus) — comparar com SILVA/RDP/BLAST no Script 5."), "WARN")
}

taxa_limpa  <- taxa
boots_limpa <- boots
cob_depois  <- apply(taxa_limpa, 2, function(x)
  round(100 * sum(!is.na(x)) / nrow(taxa_limpa), 1))

# ==============================================================================
# AUDITORIA DE ESPECIE (apenas se toSpecies)
# ==============================================================================
tem_species <- "Species" %in% colnames(taxa_limpa)
species_gg2_df <- data.frame(
  ASV_ID        = unname(seq2id[rownames(taxa_limpa)]),
  ASV           = rownames(taxa_limpa),
  Genus         = if ("Genus" %in% colnames(taxa_limpa)) taxa_limpa[, "Genus"] else NA,
  Species       = if (tem_species) taxa_limpa[, "Species"] else NA,
  Species_fonte = if (tem_species)
    ifelse(!is.na(taxa_limpa[, "Species"]), "kmer_gtdb_baixa_confianca", NA_character_)
  else NA_character_,
  stringsAsFactors = FALSE)
write.csv(species_gg2_df,
          file.path(output_path, paste0("species_fonte_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)

# ==============================================================================
# SALVAMENTO (mesmo padrao de nomes dos Scripts 2 e 3)
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
  Script        = paste0("2c_", BANCO_NOME, "_V1"),
  Versao        = VERSAO,
  Data_execucao = DATA_EXECUCAO,
  Banco         = "Greengenes2 2024.09 (GTDB)",
  Trainset      = if (USAR_TOSPECIES) "toSpecies" else "toGenus",
  Seqtab_input  = basename(arq_seqtab),
  Metodo_assign = "assignTaxonomy_kmer",
  MINBOOT       = MINBOOT,
  MULTITHREAD   = MULTITHREAD,
  ASVs_entrada  = ncol(seqtab),
  Amostras      = nrow(seqtab),
  ASVs_saida    = nrow(taxa_limpa),
  Pct_retidas   = round(100 * nrow(taxa_limpa) / nrow(taxa), 1),
  Com_addSpecies  = FALSE,
  Tempo_assign_s  = round(tempo_assign, 0),
  Nomenclatura    = "GTDB",
  Nota_uso        = "classificador_independente_harmonizado_no_script5; Species_nao_principal",
  stringsAsFactors = FALSE)
write.csv(run_metadata,
          file.path(output_path, paste0("metadata_", BANCO_NOME, ".csv")),
          row.names = FALSE, quote = TRUE)
log_msg(paste("OK: metadata_", BANCO_NOME, ".csv"), "SAVE")

cat("\n=============================================================\n")
cat("GREENGENES2 2024.09 — CONCLUIDO\n")
cat(sprintf("ASVs classificadas: %d  |  Tempo: %.1f s\n",
            nrow(taxa_limpa), tempo_assign))
cat(sprintf("Boot Genus: %.1f%%  |  Cobertura Genus: %.1f%%\n",
            boot_medio["Genus"], cob_depois["Genus"]))
cat("NOTA: Phylum sera harmonizado para NCBI/ICNP no Script 5; demais ranks GTDB sao preservados.\n")
cat("PROXIMO PASSO: 2d_BEExact.R; 02e_GSRDB_V3V4_QIIME2.R; depois 04_rblast.R\n")
cat("=============================================================\n\n")

rm(taxa, taxa_limpa, taxa_boot_obj, boots, boots_limpa, seqtab)
limpeza_mem()
log_msg("Finalizado", "FINAL")
})
