#!/usr/bin/env Rscript
# Migra arquivos soltos do output_V1 legado para as pastas por etapa.
# Padrao seguro: MIGRATION_DRY_RUN=TRUE. Nenhum arquivo e sobrescrito.

.file_args <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[
  grepl("^--file=", commandArgs(trailingOnly = FALSE))
])
.script_dir <- if (length(.file_args)) {
  dirname(normalizePath(.file_args[[1L]], mustWork = FALSE))
} else {
  getwd()
}
.lib_candidates <- unique(Filter(nzchar, c(
  Sys.getenv("PIPELINE_LIB_DIR", unset = ""),
  file.path(.script_dir, "..", "R"),
  file.path(.script_dir, "..", "..", "R"),
  file.path(getwd(), "R"),
  file.path(getwd(), "..", "R")
)))
.boot <- file.path(.lib_candidates, "pipeline_bootstrap.R")
.boot <- .boot[file.exists(.boot)]
if (!length(.boot)) {
  stop("pipeline_bootstrap.R nao localizado; defina PIPELINE_LIB_DIR.", call. = FALSE)
}
source(.boot[[1L]], local = .GlobalEnv)
rm(.file_args, .script_dir, .lib_candidates, .boot)

ler_bool <- function(nome, padrao = TRUE) {
  x <- tolower(trimws(Sys.getenv(nome, unset = if (padrao) "true" else "false")))
  if (x %in% c("true", "1", "yes", "sim", "s")) return(TRUE)
  if (x %in% c("false", "0", "no", "nao", "não", "n")) return(FALSE)
  stop(nome, " deve ser booleano.", call. = FALSE)
}

run_pipeline_script("00_migrar_saida_legada.R", "preflight", function(ctx) {
  dry_run <- ler_bool("MIGRATION_DRY_RUN", TRUE)
  remover_origem_identica <- ler_bool("MIGRATION_REMOVE_IDENTICAL", FALSE)
  raiz <- ctx$output_root
  nomes_etapas <- unname(pipeline_stage_names())

  stage_for_name <- function(nome) {
    # Diretorios historicos completos.
    mapa_dirs <- c(
      BLAST_NCBI16S = "blast", GSR_DB_V3V4 = "gsr",
      comparacao_bancos = "taxonomia", analises = "analises",
      deseq2 = "deseq2", diagnostico_inext = "inext",
      filogenia = "filogenia", faprotax = "faprotax",
      assoc_taxon_fator = "associacoes", auditoria_ASVs_nao_quimericas = "auditoria_asvs"
    )
    if (nome %in% names(mapa_dirs)) return(unname(mapa_dirs[[nome]]))

    # Arquivos canonicos e prefixos.
    if (grepl("^(seqtab|ASV_|ASVs|metadata_final|pipeline_parametros|filter_stats|merge_stats|track_|length_|chimera_|primer_|cutadapt_|sample_depth|sensitivity_|retencao_|riqueza_|REF_)", nome)) return("dada2")
    if (grepl("silva138|species_fonte_silva", nome, ignore.case = TRUE)) return("silva")
    if (grepl("rdp19|species_fonte_rdp", nome, ignore.case = TRUE)) return("rdp")
    if (grepl("gg2|greengenes", nome, ignore.case = TRUE)) return("greengenes2")
    if (grepl("taxa_bee|beexact|species_fonte_bee", nome, ignore.case = TRUE)) return("beexact")
    if (grepl("^(taxa_consenso|asvs_contaminantes|metadata_consenso|tabela_nao_identificadas)", nome)) return("taxonomia")
    if (grepl("^(phyloseq_|core9_|plus10_|all10_|comparacao_core9|sensibilidade_plus10)", nome)) return("phyloseq")
    if (grepl("preflight", nome, ignore.case = TRUE)) return("preflight")
    NA_character_
  }

  entradas <- list.files(raiz, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  entradas <- entradas[!basename(entradas) %in% nomes_etapas]
  manifest <- list(); k <- 0L

  registrar <- function(origem, destino, acao, status, detalhe = "") {
    k <<- k + 1L
    manifest[[k]] <<- data.frame(
      Origem = origem, Destino = destino, Acao = acao, Status = status,
      Detalhe = detalhe, stringsAsFactors = FALSE
    )
  }

  mover_um <- function(origem, destino) {
    if (file.exists(destino) || dir.exists(destino)) {
      if (file.exists(origem) && file.exists(destino) && !dir.exists(origem) && !dir.exists(destino)) {
        iguais <- identical(unname(tools::md5sum(origem)), unname(tools::md5sum(destino)))
        if (iguais) {
          if (!dry_run && remover_origem_identica) unlink(origem, force = TRUE)
          registrar(origem, destino, "duplicado_identico", if (dry_run) "SIMULADO" else "OK",
                    if (remover_origem_identica) "origem removivel" else "origem preservada")
          return(invisible(TRUE))
        }
      }
      registrar(origem, destino, "nao_mover", "CONFLITO", "destino ja existe e difere")
      return(invisible(FALSE))
    }
    if (dry_run) {
      registrar(origem, destino, "mover", "SIMULADO")
      return(invisible(TRUE))
    }
    dir.create(dirname(destino), recursive = TRUE, showWarnings = FALSE)
    ok <- file.rename(origem, destino)
    if (!ok && file.exists(origem) && !dir.exists(origem)) {
      ok <- file.copy(origem, destino, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
      if (ok) unlink(origem, force = TRUE)
    }
    registrar(origem, destino, "mover", if (ok) "OK" else "ERRO")
    invisible(ok)
  }

  diretorios_conteudo_direto <- c(
    "BLAST_NCBI16S", "GSR_DB_V3V4", "analises", "deseq2",
    "diagnostico_inext", "filogenia", "faprotax",
    "assoc_taxon_fator", "auditoria_ASVs_nao_quimericas"
  )

  mover_diretorio_conteudo <- function(origem, destino_root) {
    itens <- list.files(origem, full.names = TRUE, recursive = TRUE,
                        all.files = TRUE, no.. = TRUE, include.dirs = FALSE)
    if (!length(itens)) {
      registrar(origem, destino_root, "diretorio_vazio", "AVISO")
      return(invisible(TRUE))
    }
    prefixo <- paste0(normalizePath(origem, mustWork = TRUE), .Platform$file.sep)
    for (item in itens) {
      rel <- sub(prefixo, "", normalizePath(item, mustWork = TRUE), fixed = TRUE)
      mover_um(item, file.path(destino_root, rel))
    }
    if (!dry_run) {
      vazios <- list.dirs(origem, recursive = TRUE, full.names = TRUE)
      for (d in rev(vazios)) if (!length(list.files(d, all.files = TRUE, no.. = TRUE))) unlink(d, recursive = TRUE)
    }
    invisible(TRUE)
  }

  for (origem in entradas) {
    nome <- basename(origem)
    etapa <- stage_for_name(nome)
    if (is.na(etapa)) {
      registrar(origem, NA_character_, "nao_classificado", "AVISO",
                "revisao manual necessaria; arquivo preservado")
      next
    }
    stage_root <- ctx$layout$stages[[etapa]]$root
    if (dir.exists(origem)) {
      destino_root <- if (nome %in% diretorios_conteudo_direto) {
        stage_root
      } else {
        file.path(stage_root, nome)
      }
      mover_diretorio_conteudo(origem, destino_root)
    } else {
      mover_um(origem, file.path(stage_root, nome))
    }
  }

  tab <- if (length(manifest)) do.call(rbind, manifest) else data.frame(
    Origem = character(), Destino = character(), Acao = character(),
    Status = character(), Detalhe = character(), stringsAsFactors = FALSE
  )
  arq <- file.path(ctx$stage$tables, paste0(
    "manifesto_migracao_", if (dry_run) "simulacao" else "executada", ".csv"
  ))
  utils::write.csv(tab, arq, row.names = FALSE, quote = TRUE, na = "")
  cat("Modo:", if (dry_run) "SIMULACAO" else "EXECUCAO", "\n")
  cat("Itens avaliados:", nrow(tab), "\nManifesto:", arq, "\n")
  if (any(tab$Status == "CONFLITO")) {
    stop("Migracao encontrou conflitos; consulte o manifesto. Nenhum conflito foi sobrescrito.", call. = FALSE)
  }
})
