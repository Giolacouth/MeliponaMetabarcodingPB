options(stringsAsFactors = FALSE)

.pipeline_lib_dir <- local({
  source_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (is.null(source_file) || !nzchar(source_file)) {
    normalizePath(file.path(getwd(), "R"), mustWork = FALSE)
  } else {
    dirname(normalizePath(source_file, mustWork = FALSE))
  }
})

pipeline_stage_names <- function() {
  c(
    preflight = "00_preflight",
    primers = "00_primers",
    dada2 = "01_dada2",
    auditoria_asvs = "01b_auditoria_asvs",
    silva = "02a_silva",
    rdp = "02b_rdp",
    greengenes2 = "02c_greengenes2",
    beexact = "02d_beexact",
    gsr = "02e_gsr",
    blast = "04_blast",
    taxonomia = "05_taxonomia",
    phyloseq = "06_phyloseq",
    inext = "06b_inext",
    analises = "08_analises",
    deseq2 = "09_deseq2_sensibilidade",
    graficos = "10_graficos",
    faprotax = "11_faprotax"
  )
}

pipeline_layout <- function(base_path, version, output_root = NULL) {
  if (is.null(output_root) || !nzchar(output_root)) {
    output_root <- file.path(base_path, "results", paste0("output_", version))
  }

  stage_names <- pipeline_stage_names()
  stages <- lapply(stage_names, function(stage_name) {
    stage_root <- file.path(output_root, stage_name)
    list(
      root = stage_root,
      tables = file.path(stage_root, "tabelas"),
      objects = file.path(stage_root, "objetos"),
      figures = file.path(stage_root, "figuras"),
      logs = file.path(stage_root, "logs")
    )
  })

  list(output_root = output_root, stages = stages)
}

pipeline_contracts <- function(layout) {
  stage <- layout$stages
  dada2 <- stage$dada2$root
  gsr_rds <- file.path(stage$gsr$root, "rds")
  taxonomia <- stage$taxonomia$root
  phyloseq <- stage$phyloseq$root

  list(
    seqtab_global_nochim = file.path(dada2, "seqtab_global_nochim.rds"),
    seqtab_nochim = file.path(dada2, "seqtab_nochim.rds"),
    seqtab_auxiliar = file.path(dada2, "seqtab_auxiliar.rds"),
    seqtab_nochim_prefiltro = file.path(dada2, "seqtab_nochim_prefiltro.rds"),
    seqtab_pre_filtro_comprimento = file.path(dada2, "seqtab_pre_filtro_comprimento.rds"),
    seqtab_prechimera = file.path(dada2, "seqtab_prechimera.rds"),
    chimera_flags_baseline = file.path(dada2, "chimera_flags_baseline.rds"),
    asv_sequences = file.path(dada2, "ASV_sequences.tsv"),
    metadata_final = file.path(dada2, "metadata_final.tsv"),
    pipeline_parametros = file.path(dada2, "pipeline_parametros_execucao.csv"),
    taxa_silva = file.path(stage$silva$root, "taxa_silva138_limpa.rds"),
    species_silva = file.path(stage$silva$root, "species_fonte_silva138.rds"),
    taxa_rdp = file.path(stage$rdp$root, "taxa_rdp19_limpa.rds"),
    taxa_gg2 = file.path(stage$greengenes2$root, "taxa_gg2_limpa.rds"),
    taxa_beexact = file.path(stage$beexact$root, "taxa_beexact_limpa.rds"),
    taxa_gsr_disable = file.path(gsr_rds, "taxa_gsr_v3v4_confidence_disable_limpa.rds"),
    taxa_gsr_disable_alias = file.path(gsr_rds, "taxa_gsr_v3v4_limpa.rds"),
    taxa_gsr07 = file.path(gsr_rds, "taxa_gsr_v3v4_confidence_0.7_limpa.rds"),
    taxa_consenso = file.path(taxonomia, "taxa_consenso_final.rds"),
    contaminantes = file.path(taxonomia, "asvs_contaminantes_excluir.csv"),
    metadata_taxonomia = file.path(taxonomia, "metadata_consenso_taxonomico.csv"),
    taxa_gsr07_alinhada = file.path(
      taxonomia,
      "sensibilidade_gsr07",
      "taxa_gsr07_alinhada_analise.rds"
    ),
    phyloseq_core9 = file.path(phyloseq, "phyloseq_core9_primeira_run.rds"),
    phyloseq_plus10 = file.path(phyloseq, "phyloseq_plus10_com_auxiliar.rds"),
    core9_dist_bray = file.path(phyloseq, "core9_dist_bray_rel.rds"),
    core9_dist_jaccard = file.path(phyloseq, "core9_dist_jaccard_binary.rds"),
    plus10_dist_bray = file.path(phyloseq, "plus10_dist_bray_rel.rds"),
    plus10_dist_jaccard = file.path(phyloseq, "plus10_dist_jaccard_binary.rds"),
    inext_core9 = file.path(stage$inext$root, "core9_iNEXT_result.rds")
  )
}

pipeline_context <- function(stage_name) {
  stage_names <- pipeline_stage_names()
  if (!stage_name %in% names(stage_names)) {
    stop("Etapa desconhecida: ", stage_name, call. = FALSE)
  }

  default_root <- normalizePath(file.path(.pipeline_lib_dir, ".."), mustWork = FALSE)
  base_path <- Sys.getenv("PIPELINE_PROJECT_DIR", unset = default_root)
  base_path <- normalizePath(base_path, mustWork = FALSE)
  version <- Sys.getenv("PIPELINE_VERSION", unset = "V1")
  raw_path <- Sys.getenv(
    "PIPELINE_RAW_DIR",
    unset = file.path(base_path, "data", "raw")
  )
  output_root <- Sys.getenv("PIPELINE_OUTPUT_DIR", unset = "")
  layout <- pipeline_layout(base_path, version, output_root)

  list(
    base_path = base_path,
    raw_path = normalizePath(raw_path, mustWork = FALSE),
    output_root = layout$output_root,
    version = version,
    layout = layout,
    stage = layout$stages[[stage_name]],
    contracts = pipeline_contracts(layout)
  )
}

run_pipeline_script <- function(script_name, stage_name, code) {
  if (!is.function(code)) {
    stop("code deve ser uma funcao que recebe ctx.", call. = FALSE)
  }

  ctx <- pipeline_context(stage_name)
  for (path in unname(unlist(ctx$stage[c("root", "tables", "objects", "figures", "logs")]))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  message("Executando ", script_name, " [", stage_name, "]")
  code(ctx)
  invisible(ctx)
}
