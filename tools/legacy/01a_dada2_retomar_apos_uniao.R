###############################################################################
# RETOMADA SEGURA DO SCRIPT 01 APOS A FALHA NA UNIAO DAS TABELAS
#
# Use SOMENTE na mesma sessao R em que 01a_dada2_multirun.R parou em
# "=== UNINDO TABELAS ===". Este arquivo nao refaz Cutadapt, learnErrors,
# dada() ou mergePairs(); ele continua a partir das tabelas ja existentes
# em memoria.
###############################################################################

.resume_inicio <- Sys.time()
.resume_obrigatorios <- c(
  "seqtab_list", "runs_ordenados", "run_principal", "run_auxiliar",
  "output_path", "metadata", "track_list", "filter_list", "merge_list",
  "pool_used", "derepFs_aux", "derepRs_aux", "sn_aux_surv",
  "errF_principal", "errR_principal"
)
.resume_faltantes <- .resume_obrigatorios[
  !vapply(.resume_obrigatorios, exists, logical(1), envir = .GlobalEnv,
          inherits = FALSE)
]
if (length(.resume_faltantes) > 0L) {
  stop(
    "A sessao R nao preserva todos os objetos necessarios para a retomada: ",
    paste(.resume_faltantes, collapse = ", "),
    ". Nao feche a sessao apos a falha; se ela ja foi encerrada, execute o ",
    "script completo corrigido.",
    call. = FALSE
  )
}
if (!is.list(seqtab_list) || length(seqtab_list) != 2L ||
    !setequal(names(seqtab_list), runs_ordenados)) {
  stop(
    "seqtab_list nao contem exatamente as duas corridas esperadas: ",
    paste(runs_ordenados, collapse = ", "),
    call. = FALSE
  )
}
if (!all(vapply(seqtab_list[runs_ordenados], is.matrix, logical(1)))) {
  stop("Uma ou mais entradas de seqtab_list nao sao matrizes.", call. = FALSE)
}

.resume_log <- file.path(
  output_path,
  paste0("retomada_script01_", format(.resume_inicio, "%Y%m%dT%H%M%S"), ".log")
)
.resume_con <- file(.resume_log, open = "wt", encoding = "UTF-8")
sink(.resume_con, type = "output", split = TRUE)
sink(.resume_con, type = "message", append = TRUE)

.resume_status <- "ERRO"
.resume_erro <- NA_character_
tryCatch(
  withCallingHandlers({
    ###############################################################################
    # 12. UNIAO DAS TABELAS E COLAPSO DE ASVs (todas as corridas)
    ###############################################################################
    cat("\n=== UNINDO TABELAS ===\n")
    if (!setequal(names(seqtab_list), runs_ordenados) || length(seqtab_list) != 2L) {
      stop(
        "As duas corridas esperadas nao produziram tabelas de sequencias. Esperadas: ",
        paste(runs_ordenados, collapse = ", "), "; produzidas: ",
        paste(names(seqtab_list), collapse = ", ")
      )
    }
    seqtab_all <- dada2::mergeSequenceTables(
      tables = unname(seqtab_list[runs_ordenados])
    )
    cat("Tabela combinada:", nrow(seqtab_all), "amostras x", ncol(seqtab_all), "ASVs\n")

    saveRDS(seqtab_all, file.path(output_path, "seqtab_pre_filtro_comprimento.rds"))
    length_filter_audit <- data.frame(
      Sequence   = colnames(seqtab_all),
      Length     = nchar(colnames(seqtab_all)),
      TotalReads = as.numeric(colSums(seqtab_all)),
      Prevalence = as.integer(colSums(seqtab_all > 0)),
      stringsAsFactors = FALSE)
    write.table(length_filter_audit,
                file.path(output_path, "ASVs_pre_filtro_comprimento.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)


    ###############################################################################
    # 13. DISTRIBUICAO DE COMPRIMENTOS E RASTREAMENTO POR LOTE
    ###############################################################################
    cat("\n=== DISTRIBUICAO DE COMPRIMENTOS ===\n")
    seq_lengths <- nchar(getSequences(seqtab_all))

    # 1. Estatisticas Globais
    length_distribution <- table(seq_lengths)
    reads_por_comprimento <- tapply(colSums(seqtab_all), seq_lengths, sum)
    pico     <- as.integer(names(which.max(reads_por_comprimento)))
    pico_asv <- as.integer(names(which.max(length_distribution)))
    cat(sprintf("Pico Global -> Por reads: %d pb | Por ASVs: %d pb\n", pico, pico_asv))

    # 2. Aplicacao do Filtro de Janela Alvo
    keep_len <- seq_lengths >= expected_min_len & seq_lengths <= expected_max_len
    seqtab_all.filt <- seqtab_all[, keep_len, drop = FALSE]

    if (ncol(seqtab_all.filt) == 0 || sum(seqtab_all.filt) == 0)
      stop("Filtro de comprimento removeu todas as ASVs/reads. Revise expected_min_len/expected_max_len.")

    saveRDS(seqtab_all.filt, file.path(output_path, "seqtab_prechimera.rds"))
    length_removed <- data.frame(
      Sequence   = colnames(seqtab_all)[!keep_len],
      Length     = nchar(colnames(seqtab_all)[!keep_len]),
      TotalReads = as.numeric(colSums(seqtab_all)[!keep_len]),
      Prevalence = as.integer(colSums(seqtab_all > 0)[!keep_len]),
      Reason     = "Fora_da_janela_de_comprimento",
      stringsAsFactors = FALSE)
    write.table(length_removed,
                file.path(output_path, "ASVs_removidas_filtro_comprimento.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    if (nrow(length_removed) > 0L) {
      fasta_length <- Biostrings::DNAStringSet(length_removed$Sequence)
      names(fasta_length) <- paste0("LEN_REMOVED_",
                                    sprintf("%06d", seq_len(nrow(length_removed))),
                                    "|length=", length_removed$Length,
                                    "|reads=", length_removed$TotalReads)
      Biostrings::writeXStringSet(
        fasta_length,
        filepath = file.path(output_path, "ASVs_removidas_filtro_comprimento.fa"),
        format = "fasta")
    }

    # collapseNoMismatch nao aplicado: ambas as corridas usam truncLen identico
    # c(283, 278). Sem truncagem diferencial, nao ha sub-sequencias para colapsar.

    # 3. VERIFICACAO DE RETENCAO LOCAL
    # Avalia se a amostra isolada ou algum lote especifico sofreu descarte desproporcional
    reads_antes_amostra  <- rowSums(seqtab_all)
    reads_depois_amostra <- rowSums(seqtab_all.filt)
    retencao_amostra     <- round((reads_depois_amostra / reads_antes_amostra) * 100, 1)

    amostras_criticas <- names(retencao_amostra[retencao_amostra < 85])
    if (length(amostras_criticas) > 0) {
      warning("ALERTA LOCAL: As seguintes amostras perderam >15% das reads no filtro de comprimento:\n  ",
              paste(amostras_criticas, collapse = ", "))
    }

    write.csv(data.frame(SampleID = names(retencao_amostra), Retencao_Pct = retencao_amostra, Run = run_of_sample[names(retencao_amostra)]),
              file.path(output_path, "retencao_comprimento_por_amostra.csv"), row.names = FALSE)

    # 4. Visualizacao Aprimorada (Janela de Corte)
    df_plot_len <- data.frame(
      Comprimento = as.numeric(names(reads_por_comprimento)),
      Reads       = as.numeric(reads_por_comprimento)
    )

    p_len <- ggplot(df_plot_len, aes(x = Comprimento, y = Reads)) +
      geom_col(fill = "#2C7FB8", width = 1) +
      geom_vline(
        xintercept = c(expected_min_len, expected_max_len),
        color = "#D95F02",
        linetype = "dashed",
        linewidth = 1
      ) +
      annotate(
        "text",
        x = expected_min_len + 2,
        y = max(df_plot_len$Reads) * 0.92,
        label = paste0("Mínimo: ", expected_min_len, " pb"),
        color = "#D95F02",
        hjust = 0
      ) +
      annotate(
        "text",
        x = expected_max_len - 2,
        y = max(df_plot_len$Reads) * 0.92,
        label = paste0("Máximo: ", expected_max_len, " pb"),
        color = "#D95F02",
        hjust = 1
      ) +
      labs(
        title = "Distribuição de comprimento das ASVs",
        subtitle = sprintf(
          "Janela aceita: %d a %d pb",
          expected_min_len,
          expected_max_len
        ),
        x = "Comprimento (pb)",
        y = "Abundância total de leituras"
      ) +
      theme_classic()

    ggsave(file.path(plot_path, "length_distribution.pdf"), p_len, width = 8, height = 5)

    # 5. Resumo de Saida
    pct_len_reads   <- round(sum(seqtab_all.filt) / sum(seqtab_all) * 100, 1)
    pct_len_asv <- round(ncol(seqtab_all.filt) / ncol(seqtab_all) * 100, 1)

    cat(sprintf("Filtro aplicado: %d a %d pb\n", expected_min_len, expected_max_len))
    cat(sprintf("ASVs retidas globais   : %d de %d (%.1f%%)\n", ncol(seqtab_all.filt), ncol(seqtab_all), pct_len_asv))
    cat(sprintf("Reads retidas globais  : %.1f%%\n", pct_len_reads))
    if (pct_len_reads < 90) {
      cat("!! ALERTA GLOBAL: < 90% das reads totais foram retidas. Inspecione o PDF.\n")
    }

    ################################################################################
    # 14. REMOCAO E AUDITORIA DE QUIMERAS (tabela COMBINADA — 10 amostras)
    ###############################################################################
    cat("\n=== REMOCAO E AUDITORIA DE QUIMERAS (abordagem combinada) ===\n")

    # isBimeraDenovoTable: identificacao com vetor logico para auditoria
    # Nao aceita "method" — esse parametro e de removeBimeraDenovo
    bimeras_logico <- isBimeraDenovoTable(
      seqtab_all.filt,
      minSampleFraction          = 0.9,   # sequencia bimerica em >=90% das amostras
      ignoreNNegatives           = 1,
      minFoldParentOverAbundance = 1.5,   # padrao documentado
      minParentAbundance         = 2,
      allowOneOff                = FALSE,
      multithread                = FALSE,
      verbose                    = TRUE
    )

    names(bimeras_logico) <- colnames(seqtab_all.filt)
    saveRDS(bimeras_logico, file.path(output_path, "chimera_flags_baseline.rds"))
    saveRDS(seqtab_all.filt[, bimeras_logico, drop = FALSE],
            file.path(output_path, "seqtab_quimeras_removidas.rds"))
    chimera_audit_complete <- data.frame(
      Sequence   = colnames(seqtab_all.filt),
      Length     = nchar(colnames(seqtab_all.filt)),
      TotalReads = as.numeric(colSums(seqtab_all.filt)),
      Prevalence = as.integer(colSums(seqtab_all.filt > 0)),
      IsChimera  = as.logical(bimeras_logico),
      stringsAsFactors = FALSE)
    write.table(chimera_audit_complete,
                file.path(output_path, "auditoria_completa_quimeras.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    asvs_removed  <- sum(bimeras_logico)
    seqs_quimeras <- colnames(seqtab_all.filt)[bimeras_logico]
    seqtab.nochim <- seqtab_all.filt[, !bimeras_logico, drop = FALSE]

    if (ncol(seqtab.nochim) == 0 || sum(seqtab.nochim) == 0)
      stop("Remocao de quimeras removeu todas as ASVs/reads. Revise parametros.")

    asvs_before <- ncol(seqtab_all.filt)
    asvs_after  <- ncol(seqtab.nochim)
    pct_nochim  <- round(sum(seqtab.nochim) / sum(seqtab_all.filt) * 100, 1)
    # Auditoria das quimeras removidas
    if (asvs_removed > 0) {
      auditoria_quimeras <- data.frame(
        Sequence    = seqs_quimeras,
        Total_Reads = colSums(seqtab_all.filt)[seqs_quimeras],
        Prevalencia = colSums(seqtab_all.filt[, seqs_quimeras, drop = FALSE] > 0),
        stringsAsFactors = FALSE)
      write.csv(auditoria_quimeras,
                file.path(output_path, "auditoria_quimeras_removidas.csv"),
                row.names = FALSE)
      cat(sprintf("auditoria_quimeras_removidas.csv: %d sequencias.\n", asvs_removed))
    } else {
      cat("Nenhuma quimera detectada.\n")
    }

    write.csv(data.frame(ASVs_before = asvs_before, ASVs_after = asvs_after,
                         ASVs_removed = asvs_removed,
                         Pct_removed = round(100 * asvs_removed / asvs_before, 1)),
              file.path(output_path, "chimera_stats.csv"), row.names = FALSE)

    cat("Retencao pos-quimera:", pct_nochim, "% (esperado 75-95%)\n")

    chimera_check <- data.frame(
      SampleID     = rownames(seqtab.nochim),
      Reads_before = rowSums(seqtab_all.filt)[rownames(seqtab.nochim)],
      Reads_after  = rowSums(seqtab.nochim))
    chimera_check$Pct_retention <- round(100 * chimera_check$Reads_after / chimera_check$Reads_before, 1)
    write.csv(chimera_check, file.path(output_path, "chimera_retention_per_sample.csv"), row.names = FALSE)


    if (any(chimera_check$Pct_retention < 70))
      warning("Amostras com retencao pos-quimera <70%.")

    sem_quimera <- chimera_check$SampleID[chimera_check$Pct_retention == 100]
    if (length(sem_quimera) > 0) {
      cat(sprintf("ATENCAO: %d amostra(s) com 100%% retencao pos-quimera.\n",
                  length(sem_quimera)))
      cat("Causa provavel: baixa diversidade unica ou amostra com poucos candidatos\n")
      cat("a sequencias parentais. Declarar como limitacao nos Metodos.\n")
    }

    # Salvar tabela pre-filtro (combinada, 10 amostras)
    saveRDS(seqtab.nochim, file.path(output_path, "seqtab_nochim_prefiltro.rds"))
    cat("seqtab_nochim_prefiltro.rds salvo (combinado, pre-filtro — para utilitarios)\n")


    # --------------------------------------------------------------------------
    # Separacao principal / auxiliar (pos-quimera, pre-filtro de frequencia)
    # --------------------------------------------------------------------------
    amostras_principal <- metadata$SampleID[metadata$Run == run_principal]
    amostras_auxiliar  <- metadata$SampleID[metadata$Run == run_auxiliar]

    amostras_main_pres <- intersect(amostras_principal, rownames(seqtab.nochim))
    amostras_aux_pres  <- intersect(amostras_auxiliar,  rownames(seqtab.nochim))

    seqtab_main_nochim <- seqtab.nochim[amostras_main_pres, , drop = FALSE]
    seqtab_aux_nochim  <- seqtab.nochim[amostras_aux_pres,  , drop = FALSE]

    seqtab_main_nochim <- seqtab_main_nochim[
      , colSums(seqtab_main_nochim) > 0, drop = FALSE]
    seqtab_aux_nochim  <- seqtab_aux_nochim[
      , colSums(seqtab_aux_nochim)  > 0, drop = FALSE]

    cat(sprintf("\nSeparacao pos-quimera:\n"))
    cat(sprintf("  seqtab_main_nochim : %d amostras | %d ASVs\n",
                nrow(seqtab_main_nochim), ncol(seqtab_main_nochim)))
    cat(sprintf("  seqtab_aux_nochim  : %d amostra  | %d ASVs\n",
                nrow(seqtab_aux_nochim),  ncol(seqtab_aux_nochim)))

    if (nrow(seqtab_main_nochim) == 0 || ncol(seqtab_main_nochim) == 0)
      stop("Tabela principal vazia apos separacao pos-quimera.")
    # A corrida auxiliar integra o universo canônico. Não é permitido seguir com
    # tabela degenerada ou remover a amostra silenciosamente.
    aux_vazia <- (nrow(seqtab_aux_nochim) != 1L || ncol(seqtab_aux_nochim) == 0L)
    if (aux_vazia) {
      stop(
        "Tabela auxiliar inválida após a remoção de quimeras: amostras=",
        nrow(seqtab_aux_nochim), "; ASVs=", ncol(seqtab_aux_nochim),
        ". Esperado: 1 amostra e pelo menos 1 ASV.",
        call. = FALSE
      )
    }

    ###############################################################################
    # 14B. PROFUNDIDADE POR AMOSTRA (pos-quimera, pre-filtro)
    #
    # Usa seqtab.nochim:
    #   - depois da remocao de quimeras;
    #   - antes dos filtros de abundancia e prevalencia;
    #   - inclui as 10 amostras das duas corridas.
    ###############################################################################

    cat("\n=== PROFUNDIDADE POR AMOSTRA (pos-quimera, pre-filtro) ===\n")

    # --------------------------------------------------------------------------
    # 1. Calculo e validacao
    # --------------------------------------------------------------------------

    sample_depth <- rowSums(seqtab.nochim)

    if (length(sample_depth) == 0L) {
      stop(
        "Nao foi possivel calcular a profundidade: ",
        "seqtab.nochim nao contem amostras."
      )
    }

    if (is.null(names(sample_depth))) {
      stop(
        "seqtab.nochim nao possui nomes de amostras nas linhas."
      )
    }

    if (anyNA(sample_depth)) {
      stop(
        "Foram encontrados valores NA na profundidade pos-quimera."
      )
    }

    if (any(!is.finite(sample_depth))) {
      stop(
        "Foram encontrados valores nao finitos na profundidade pos-quimera."
      )
    }

    if (any(sample_depth < 0)) {
      stop(
        "Foram encontrados valores negativos na profundidade pos-quimera."
      )
    }

    amostras_sem_reads <- names(sample_depth)[sample_depth == 0]

    if (length(amostras_sem_reads) > 0L) {
      warning(
        "Amostras sem reads apos remocao de quimeras: ",
        paste(amostras_sem_reads, collapse = ", ")
      )
    }

    # --------------------------------------------------------------------------
    # 2. Tabela de profundidade
    # --------------------------------------------------------------------------

    sample_depth_df <- data.frame(
      SampleID = names(sample_depth),
      SampleLabel = unname(setNames(metadata$SampleLabel, metadata$SampleID)[names(sample_depth)]),
      Run = unname(run_of_sample[names(sample_depth)]),
      Reads_pos_quimera = as.numeric(sample_depth),
      stringsAsFactors = FALSE
    )

    if (anyNA(sample_depth_df$SampleLabel) || any(sample_depth_df$SampleLabel == "")) {
      stop("SampleLabel nao encontrado para uma ou mais amostras da tabela de profundidade.")
    }

    if (anyNA(sample_depth_df$Run)) {
      amostras_sem_run <- sample_depth_df$SampleID[
        is.na(sample_depth_df$Run)
      ]

      stop(
        "Corrida nao encontrada para a(s) amostra(s): ",
        paste(amostras_sem_run, collapse = ", ")
      )
    }

    sample_depth_df <- sample_depth_df[
      order(sample_depth_df$Reads_pos_quimera),
      ,
      drop = FALSE
    ]

    write.csv(
      sample_depth_df,
      file.path(output_path, "sample_depth_prefiltro.csv"),
      row.names = FALSE
    )

    cat("Profundidade pos-quimera por amostra:\n")
    print(sample_depth_df, row.names = FALSE)

    cat("\nResumo da profundidade pos-quimera:\n")
    print(summary(sample_depth_df$Reads_pos_quimera))

    cat(
      sprintf(
        paste0(
          "Minima: %d | Mediana: %.1f | Media: %.1f | Maxima: %d\n"
        ),
        min(sample_depth_df$Reads_pos_quimera),
        median(sample_depth_df$Reads_pos_quimera),
        mean(sample_depth_df$Reads_pos_quimera),
        max(sample_depth_df$Reads_pos_quimera)
      )
    )

    # --------------------------------------------------------------------------
    # 3. Histograma geral
    # --------------------------------------------------------------------------

    hist_file <- file.path(
      plot_path,
      "sample_depth_hist.pdf"
    )

    grDevices::pdf(
      hist_file,
      width = 8,
      height = 6
    )

    tryCatch(
      {
        hist(
          sample_depth_df$Reads_pos_quimera,
          main = "Profundidade por amostra apos remocao de quimeras",
          xlab = "Reads pos-quimera",
          ylab = "Numero de amostras",
          breaks = "Sturges"
        )

        abline(
          v = median(sample_depth_df$Reads_pos_quimera),
          lty = 2,
          lwd = 2
        )

        legend(
          "topright",
          legend = paste0(
            "Mediana = ",
            format(
              median(sample_depth_df$Reads_pos_quimera),
              big.mark = ".",
              scientific = FALSE,
              trim = TRUE
            )
          ),
          lty = 2,
          lwd = 2,
          bty = "n"
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )

    # --------------------------------------------------------------------------
    # 4. Grafico por amostra
    # --------------------------------------------------------------------------

    barplot_file <- file.path(
      plot_path,
      "sample_depth_por_amostra.pdf"
    )

    grDevices::pdf(
      barplot_file,
      width = 13,
      height = 8
    )

    tryCatch(
      {
        # Aumenta a margem esquerda para comportar os nomes completos
        par(
          mar = c(5, 14, 4, 2) + 0.1,
          mgp = c(2.5, 0.8, 0)
        )

        barplot(
          height = sample_depth_df$Reads_pos_quimera,
          names.arg = sample_depth_df$SampleLabel,
          horiz = TRUE,
          las = 1,
          cex.names = 0.75,
          xlab = "Reads pós-quimera",
          ylab = "",
          main = paste0(
            "Profundidade por amostra ",
            "(pós-quimera, pré-filtro)"
          ),
          xlim = c(
            0,
            max(sample_depth_df$Reads_pos_quimera) * 1.08
          )
        )

        abline(
          v = median(sample_depth_df$Reads_pos_quimera),
          lty = 2,
          lwd = 2
        )

        legend(
          "bottomright",
          legend = paste0(
            "Mediana = ",
            format(
              median(sample_depth_df$Reads_pos_quimera),
              big.mark = ".",
              decimal.mark = ",",
              scientific = FALSE,
              trim = TRUE
            ),
            " reads"
          ),
          lty = 2,
          lwd = 2,
          bty = "n"
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )

    cat("Arquivos gerados:\n")
    cat(" - ", hist_file, "\n", sep = "")
    cat(" - ", barplot_file, "\n", sep = "")
    cat(
      " - ",
      file.path(output_path, "sample_depth_prefiltro.csv"),
      "\n",
      sep = ""
    )
    ###############################################################################
    # 14C. ANALISE DE SENSIBILIDADE DOS FILTROS
    # Executada sobre seqtab_main_nochim (9 amostras) — os limiares sao
    # selecionados para a analise primaria. min_prev avaliado dentro do n=9.
    ###############################################################################
    cat("\n=== SENSIBILIDADE DOS FILTROS (corrida principal, n=9) ===\n")
    evaluate_threshold <- function(seqtab, min_reads, min_prev) {
      total_reads <- colSums(seqtab)
      prevalence  <- colSums(seqtab > 0)   # n de amostras com reads > 0
      keep    <- !(total_reads < min_reads & prevalence < min_prev)
      seqtab2 <- seqtab[, keep, drop = FALSE]
      if (sum(seqtab) == 0) stop("Tabela vazia na analise de sensibilidade.")
      data.frame(
        min_reads     = min_reads,
        min_prev      = min_prev,
        ASVs          = ncol(seqtab2),
        Reads         = sum(seqtab2),
        Pct_reads     = round(sum(seqtab2) / sum(seqtab) * 100, 2),
        Mean_Richness = mean(rowSums(seqtab2 > 0)),
        Mean_Shannon  = mean(vegan::diversity(seqtab2, index = "shannon")))
    }
    grid <- expand.grid(min_reads = c(1, 5, 10, 20), min_prev = c(1, 2, 3))
    sensitivity_results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
      evaluate_threshold(seqtab_main_nochim, grid$min_reads[i], grid$min_prev[i])))
    write.csv(sensitivity_results,
              file.path(output_path, "sensitivity_thresholds.csv"), row.names = FALSE)
    cat("Sensibilidade sobre corrida principal (", nrow(seqtab_main_nochim),
        " amostras):\n")
    print(sensitivity_results)

    ###############################################################################
    # 15. FILTRO DE PREVALENCIA + ABUNDANCIA
    #
    # Corrida principal (seqtab_main_nochim, n=9):
    #   logica E: descartar se (reads < min_total_reads) E (prevalencia < min_prevalence)
    #   prevalencia = numero de amostras com reads > 0
    #   min_prevalence = 2 aplicado sobre 9 amostras (~22%)
    #
    # Corrida auxiliar (seqtab_aux_nochim, n=1):
    #   filtro por leitura minima apenas (prevalencia maxima = 1; min_prev invalido)
    ###############################################################################
    cat("\n=== FILTRO DE PREVALENCIA + ABUNDANCIA ===\n")

    # --- Corrida principal ---
    cat(sprintf("Corrida principal (n=%d amostras):\n", nrow(seqtab_main_nochim)))
    cat(sprintf("  Criterio: descartar se (reads < %d) E (prevalencia < %d)\n",
                min_total_reads, min_prevalence))

    asv_reads_main <- colSums(seqtab_main_nochim)
    asv_prev_main  <- colSums(seqtab_main_nochim > 0)   # > 0, nao > 2
    descartar_main <- (asv_reads_main < min_total_reads) & (asv_prev_main < min_prevalence)
    manter_main    <- !descartar_main

    # Usa colnames() enquanto a tabela ainda contem todas as colunas (pre-filtro).
    asvs_desc_principal_seq <- data.frame(
      Sequence    = colnames(seqtab_main_nochim)[descartar_main],
      Length      = nchar(colnames(seqtab_main_nochim)[descartar_main]),
      Reads_total = as.numeric(asv_reads_main[descartar_main]),
      Prevalencia = as.integer(asv_prev_main[descartar_main]),
      Reason      = rep(paste0("Reads_total < ", min_total_reads,
                               " E Prevalencia < ", min_prevalence),
                        sum(descartar_main)),
      stringsAsFactors = FALSE)
    asvs_desc_principal_seq <- asvs_desc_principal_seq[
      order(asvs_desc_principal_seq$Reads_total, decreasing = TRUE), , drop = FALSE]
    write.table(asvs_desc_principal_seq,
                file.path(output_path, "ASVs_descartadas_principal_com_sequencia.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # --- CHECKPOINT DE FILTRO: interrupcao obrigatoria -------------------------
    if (is.na(CENARIO_FILTRO)) {
      cat("\n===================== CHECKPOINT DE FILTRO =====================\n")
      cat(sprintf("Limiares atuais: descartar se reads < %d E prevalencia < %d.\n",
                  min_total_reads, min_prevalence))
      cat(sprintf("Seriam descartadas %d de %d ASVs na corrida principal.\n",
                  sum(descartar_main), length(descartar_main)))
      cat("Auditoria: ASVs_descartadas_principal_com_sequencia.tsv\n")
      cat("Para comparar cenarios, rode 01b_auditoria_quimeras_recuperacao_ASVs.R\n")
      cat("(le seqtab_nochim_prefiltro.rds, ja salvo).\n")
      cat("Depois, edite CENARIO_FILTRO ('aplicar' ou 'sem_filtro') no topo\n")
      cat("deste script e execute-o novamente.\n")
      cat("================================================================\n")
      stop("Filtro nao confirmado. Defina CENARIO_FILTRO e re-execute.",
           call. = FALSE)
    } else if (identical(CENARIO_FILTRO, "sem_filtro")) {
      manter_main    <- rep(TRUE, length(manter_main))
      descartar_main <- !manter_main
      cat("CENARIO_FILTRO = 'sem_filtro': nenhuma ASV removida por frequencia (principal).\n")
    } else if (!identical(CENARIO_FILTRO, "aplicar")) {
      stop("CENARIO_FILTRO invalido: use 'aplicar' ou 'sem_filtro'.", call. = FALSE)
    }
    # --- fim CHECKPOINT ---

    asvs_main_antes  <- rowSums(seqtab_main_nochim > 0)
    reads_main_antes <- sum(seqtab_main_nochim)
    seqtab_main_nochim <- seqtab_main_nochim[, manter_main, drop = FALSE]
    pct_reads_main <- round(sum(seqtab_main_nochim) / reads_main_antes * 100, 1)

    cat(sprintf("  ASVs antes: %d | descarte: %d | retidas: %d\n",
                length(manter_main), sum(descartar_main), sum(manter_main)))
    cat(sprintf("  Reads retidas: %.1f%%\n", pct_reads_main))

    amostras_vazias_main <- rownames(seqtab_main_nochim)[
      rowSums(seqtab_main_nochim) == 0]
    if (length(amostras_vazias_main) > 0) {
      stop(
        "Amostra(s) da corrida principal ficaram com zero reads apos o filtro: ",
        paste(amostras_vazias_main, collapse = ", "),
        ". O pipeline nao removera amostras silenciosamente. Reveja os filtros."
      )
    }

    # --- Corrida auxiliar ---
    cat(sprintf("\nCorrida auxiliar (n=%d amostra):\n", nrow(seqtab_aux_nochim)))
    cat(sprintf("  Criterio: descartar se (reads < %d) [min_prevalence nao aplicavel]\n",
                min_total_reads))

    if (nrow(seqtab_aux_nochim) > 0 && ncol(seqtab_aux_nochim) > 0 && sum(seqtab_aux_nochim) > 0) {
      asv_reads_aux <- colSums(seqtab_aux_nochim)
      descartar_aux <- asv_reads_aux < min_total_reads
      manter_aux    <- !descartar_aux
      if (identical(CENARIO_FILTRO, "sem_filtro")) {
        manter_aux    <- rep(TRUE, length(manter_aux))
        descartar_aux <- !manter_aux
        cat("CENARIO_FILTRO = 'sem_filtro': nenhuma ASV removida por frequencia (auxiliar).\n")
      }

      asvs_desc_aux_seq <- data.frame(
        Sequence    = colnames(seqtab_aux_nochim)[descartar_aux],
        Length      = nchar(colnames(seqtab_aux_nochim)[descartar_aux]),
        Reads_total = as.numeric(asv_reads_aux[descartar_aux]),
        Prevalencia = as.integer(colSums(seqtab_aux_nochim > 0)[descartar_aux]),
        Reason      = rep(paste0("Reads_total < ", min_total_reads,
                                 " na corrida auxiliar"), sum(descartar_aux)),
        stringsAsFactors = FALSE)
      write.table(asvs_desc_aux_seq,
                  file.path(output_path, "ASVs_descartadas_auxiliar_com_sequencia.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      reads_aux_antes <- sum(seqtab_aux_nochim)
      seqtab_aux_nochim <- seqtab_aux_nochim[, manter_aux, drop = FALSE]
      pct_reads_aux <- round(sum(seqtab_aux_nochim) / reads_aux_antes * 100, 1)
    } else {
      stop(
        "Tabela auxiliar vazia antes do filtro de abundancia. ",
        "A amostra S10 nao pode ser removida silenciosamente.",
        call. = FALSE
      )
    }

    cat(sprintf("  ASVs antes: %d | descarte: %d | retidas: %d\n",
                length(manter_aux), sum(descartar_aux), sum(manter_aux)))
    cat(sprintf("  Reads retidas: %.1f%%\n", pct_reads_aux))

    amostras_vazias_aux <- rownames(seqtab_aux_nochim)[
      rowSums(seqtab_aux_nochim) == 0]
    if (length(amostras_vazias_aux) > 0) {
      stop(
        "A amostra da corrida auxiliar ficou com zero reads apos o filtro: ",
        paste(amostras_vazias_aux, collapse = ", "),
        ". O pipeline nao removera essa amostra silenciosamente; reveja o limiar.",
        call. = FALSE
      )
    }

    # Auditoria das ASVs descartadas da corrida principal
    asvs_desc_df <- data.frame(
      ASV_idx     = which(descartar_main),
      Reads_total = asv_reads_main[descartar_main],
      Prevalencia = asv_prev_main[descartar_main],
      stringsAsFactors = FALSE)
    asvs_desc_df <- asvs_desc_df[order(asvs_desc_df$Reads_total, decreasing = TRUE), ]
    limiar_alerta <- min_total_reads * 10
    asvs_alerta   <- asvs_desc_df[asvs_desc_df$Reads_total >= limiar_alerta, ]

    cat(sprintf("\n=== AUDITORIA DE ASVs DESCARTADAS (corrida principal) ===\n"))
    cat(sprintf("Total descartadas: %d | Com reads >= %d: %d\n",
                nrow(asvs_desc_df), limiar_alerta, nrow(asvs_alerta)))
    if (nrow(asvs_alerta) > 0) {
      cat(sprintf("ATENCAO: %d ASV(s) com >= %d reads — inspecionar manualmente:\n",
                  nrow(asvs_alerta), limiar_alerta))
      print(head(asvs_alerta, 20))
    }
    write.csv(asvs_desc_df,
              file.path(output_path, "asvs_descartadas_auditoria.csv"),
              row.names = FALSE)

    # Riqueza antes/depois por amostra (corrida principal)
    asvs_main_depois <- rowSums(seqtab_main_nochim > 0)
    amostras_comuns  <- intersect(names(asvs_main_antes), names(asvs_main_depois))
    riqueza_check <- data.frame(
      SampleID      = amostras_comuns,
      ASVs_antes    = asvs_main_antes[amostras_comuns],
      ASVs_depois   = asvs_main_depois[amostras_comuns])
    riqueza_check$ASVs_removidas <- riqueza_check$ASVs_antes - riqueza_check$ASVs_depois
    riqueza_check$Pct_retido <- round(
      riqueza_check$ASVs_depois / riqueza_check$ASVs_antes * 100, 1)
    print(riqueza_check, row.names = FALSE)
    write.csv(riqueza_check,
              file.path(output_path, "riqueza_por_amostra_filtro.csv"),
              row.names = FALSE)
    if (any(riqueza_check$Pct_retido < 50))
      warning("Amostra(s) perderam >50% das ASVs no filtro de frequencia.")

    # Profundidade final — corrida principal (base para rarefacao no Script 8)
    if (nrow(seqtab_main_nochim) == 0 || ncol(seqtab_main_nochim) == 0)
      stop("Tabela principal vazia apos filtro de prevalencia/abundancia.")

    depth_final_main <- rowSums(seqtab_main_nochim)
    cat(sprintf("\nProfundidade final — corrida principal (seqtab_nochim.rds):\n"))
    print(sort(depth_final_main))
    cat(sprintf("Minima: %d | Mediana: %d | Maxima: %d\n",
                min(depth_final_main), as.integer(median(depth_final_main)),
                max(depth_final_main)))
    min_idx_fin <- which.min(depth_final_main)
    cat(sprintf("Amostra limitante (rarefacao): %s (%d reads)\n",
                names(depth_final_main)[min_idx_fin], depth_final_main[min_idx_fin]))

    depth_final_aux <- rowSums(seqtab_aux_nochim)
    cat(sprintf("\nProfundidade final — corrida auxiliar (seqtab_auxiliar.rds):\n"))
    print(depth_final_aux)


    ###############################################################################
    # 15B. SENSIBILIDADE DO MODELO DE ERRO DA AUXILIAR
    #
    # Compara:
    #   1. modelo proprio da corrida auxiliar;
    #   2. modelo transferido da corrida principal.
    #
    # Os dois cenarios repetem:
    #   dada -> merge -> comprimento -> quimeras -> filtro de abundancia.
    #
    # Nenhuma saida primaria e sobrescrita.
    ###############################################################################

    cat(
      "\n=== SENSIBILIDADE COMPLETA: ",
      "MODELO PROPRIO vs. TRANSFERIDO ===\n",
      sep = ""
    )

    sens_ok <- (
      exists("derepFs_aux") &&
        exists("derepRs_aux") &&
        !is.null(errF_principal) &&
        !is.null(errR_principal) &&
        !isTRUE(aux_vazia)
    )

    processar_cenario_sensibilidade <- function(
        seqtab_aux_pre,
        nome_modelo) {

      seqtab_comb <- dada2::mergeSequenceTables(
        seqtab_list[[run_principal]],
        seqtab_aux_pre
      )

      comprimentos <- nchar(
        dada2::getSequences(seqtab_comb)
      )

      keep_len_sens <- (
        comprimentos >= expected_min_len &
          comprimentos <= expected_max_len
      )

      seqtab_len <- seqtab_comb[
        ,
        keep_len_sens,
        drop = FALSE
      ]

      if (ncol(seqtab_len) == 0L || sum(seqtab_len) == 0) {
        stop(
          "Cenario ",
          nome_modelo,
          ": tabela vazia apos filtro de comprimento."
        )
      }

      bimeras_sens <- dada2::isBimeraDenovoTable(
        seqtab_len,
        minSampleFraction = 0.9,
        ignoreNNegatives = 1,
        minFoldParentOverAbundance = 1.5,
        minParentAbundance = 2,
        allowOneOff = FALSE,
        multithread = FALSE,
        verbose = TRUE
      )

      seqtab_nochim_sens <- seqtab_len[
        ,
        !bimeras_sens,
        drop = FALSE
      ]

      main_sens <- seqtab_nochim_sens[
        intersect(
          amostras_principal,
          rownames(seqtab_nochim_sens)
        ),
        ,
        drop = FALSE
      ]

      aux_sens <- seqtab_nochim_sens[
        intersect(
          amostras_auxiliar,
          rownames(seqtab_nochim_sens)
        ),
        ,
        drop = FALSE
      ]

      main_sens <- main_sens[
        ,
        colSums(main_sens) > 0,
        drop = FALSE
      ]

      aux_sens <- aux_sens[
        ,
        colSums(aux_sens) > 0,
        drop = FALSE
      ]

      main_pre_filtro <- main_sens
      aux_pre_filtro  <- aux_sens

      if (ncol(main_sens) > 0L) {
        reads_main_sens <- colSums(main_sens)
        prev_main_sens  <- colSums(main_sens > 0)

        descartar_main_sens <- (
          reads_main_sens < min_total_reads &
            prev_main_sens < min_prevalence
        )

        main_sens <- main_sens[
          ,
          !descartar_main_sens,
          drop = FALSE
        ]
      }

      if (ncol(aux_sens) > 0L) {
        aux_sens <- aux_sens[
          ,
          colSums(aux_sens) >= min_total_reads,
          drop = FALSE
        ]
      }

      list(
        Modelo = nome_modelo,
        seqtab_combinada_prequimera = seqtab_comb,
        seqtab_main_final = main_sens,
        seqtab_aux_posquimera = aux_pre_filtro,
        seqtab_aux_final = aux_sens,
        N_ASVs_aux_prequimera = ncol(seqtab_aux_pre),
        Reads_aux_prequimera = sum(seqtab_aux_pre),
        N_ASVs_aux_posquimera = ncol(aux_pre_filtro),
        Reads_aux_posquimera = sum(aux_pre_filtro),
        N_ASVs_aux_final = ncol(aux_sens),
        Reads_aux_final = sum(aux_sens),
        N_ASVs_main_final = ncol(main_sens),
        Reads_main_final = sum(main_sens)
      )
    }

    if (sens_ok) {

      pool_aux <- unname(pool_by_run[run_auxiliar])

      dadaF_transf <- dada2::dada(
        derepFs_aux,
        err = errF_principal,
        pool = pool_aux,
        multithread = TRUE
      )

      dadaR_transf <- dada2::dada(
        derepRs_aux,
        err = errR_principal,
        pool = pool_aux,
        multithread = TRUE
      )

      if (inherits(dadaF_transf, "dada")) {
        dadaF_transf <- setNames(
          list(dadaF_transf),
          sn_aux_surv
        )
      }

      if (inherits(dadaR_transf, "dada")) {
        dadaR_transf <- setNames(
          list(dadaR_transf),
          sn_aux_surv
        )
      }

      mergers_transf <- dada2::mergePairs(
        dadaF_transf,
        derepFs_aux,
        dadaR_transf,
        derepRs_aux,
        minOverlap = minOverlap_par,
        maxMismatch = maxMismatch_par,
        justConcatenate = FALSE,
        verbose = TRUE
      )

      if (is.data.frame(mergers_transf)) {
        mergers_transf <- setNames(
          list(mergers_transf),
          sn_aux_surv
        )
      }

      seqtab_aux_proprio_pre <- seqtab_list[[run_auxiliar]]



      seqtab_aux_transf_pre <- dada2::makeSequenceTable(
        mergers_transf
      )

      cenario_proprio <- processar_cenario_sensibilidade(
        seqtab_aux_pre = seqtab_aux_proprio_pre,
        nome_modelo = "proprio"
      )

      cenario_transferido <- processar_cenario_sensibilidade(
        seqtab_aux_pre = seqtab_aux_transf_pre,
        nome_modelo = "transferido"
      )

      resumo_sens <- data.frame(
        Modelo = c(
          cenario_proprio$Modelo,
          cenario_transferido$Modelo
        ),
        ASVs_prequimera = c(
          cenario_proprio$N_ASVs_aux_prequimera,
          cenario_transferido$N_ASVs_aux_prequimera
        ),
        Reads_prequimera = c(
          cenario_proprio$Reads_aux_prequimera,
          cenario_transferido$Reads_aux_prequimera
        ),
        ASVs_posquimera = c(
          cenario_proprio$N_ASVs_aux_posquimera,
          cenario_transferido$N_ASVs_aux_posquimera
        ),
        Reads_posquimera = c(
          cenario_proprio$Reads_aux_posquimera,
          cenario_transferido$Reads_aux_posquimera
        ),
        ASVs_final = c(
          cenario_proprio$N_ASVs_aux_final,
          cenario_transferido$N_ASVs_aux_final
        ),
        Reads_final = c(
          cenario_proprio$Reads_aux_final,
          cenario_transferido$Reads_aux_final
        ),
        ASVs_main_final = c(
          cenario_proprio$N_ASVs_main_final,
          cenario_transferido$N_ASVs_main_final
        ),
        Reads_main_final = c(
          cenario_proprio$Reads_main_final,
          cenario_transferido$Reads_main_final
        ),
        stringsAsFactors = FALSE
      )

      seqs_proprio <- colnames(
        cenario_proprio$seqtab_aux_final
      )

      seqs_transferido <- colnames(
        cenario_transferido$seqtab_aux_final
      )

      seqs_uniao <- union(
        seqs_proprio,
        seqs_transferido
      )

      seqs_compartilhadas_sens <- intersect(
        seqs_proprio,
        seqs_transferido
      )

      jaccard_sens <- if (length(seqs_uniao) > 0L) {
        length(seqs_compartilhadas_sens) /
          length(seqs_uniao)
      } else {
        NA_real_
      }

      abund_proprio <- setNames(
        numeric(length(seqs_uniao)),
        seqs_uniao
      )

      abund_transferido <- setNames(
        numeric(length(seqs_uniao)),
        seqs_uniao
      )

      if (length(seqs_proprio) > 0L) {
        abund_proprio[seqs_proprio] <- colSums(
          cenario_proprio$seqtab_aux_final
        )
      }

      if (length(seqs_transferido) > 0L) {
        abund_transferido[seqs_transferido] <- colSums(
          cenario_transferido$seqtab_aux_final
        )
      }

      bray_curtis_sens <- if (
        sum(abund_proprio + abund_transferido) > 0
      ) {
        sum(abs(abund_proprio - abund_transferido)) /
          sum(abund_proprio + abund_transferido)
      } else {
        NA_real_
      }

      auditoria_sens <- data.frame(
        Sequence = seqs_uniao,
        Reads_modelo_proprio = as.numeric(
          abund_proprio[seqs_uniao]
        ),
        Reads_modelo_transferido = as.numeric(
          abund_transferido[seqs_uniao]
        ),
        Compartilhada = seqs_uniao %in%
          seqs_compartilhadas_sens,
        stringsAsFactors = FALSE
      )

      metricas_concordancia <- data.frame(
        ASVs_compartilhadas =
          length(seqs_compartilhadas_sens),
        ASVs_exclusivas_proprio =
          length(setdiff(
            seqs_proprio,
            seqs_transferido
          )),
        ASVs_exclusivas_transferido =
          length(setdiff(
            seqs_transferido,
            seqs_proprio
          )),
        Jaccard = jaccard_sens,
        Bray_Curtis = bray_curtis_sens,
        stringsAsFactors = FALSE
      )

      saveRDS(
        seqtab_aux_proprio_pre,
        file.path(
          output_path,
          "seqtab_aux_modelo_proprio_prequimera.rds"
        )
      )

      saveRDS(
        seqtab_aux_transf_pre,
        file.path(
          output_path,
          "seqtab_aux_modelo_transferido_prequimera.rds"
        )
      )

      saveRDS(
        cenario_proprio$seqtab_aux_final,
        file.path(
          output_path,
          "seqtab_aux_modelo_proprio_final.rds"
        )
      )

      saveRDS(
        cenario_transferido$seqtab_aux_final,
        file.path(
          output_path,
          "seqtab_aux_modelo_transferido_final.rds"
        )
      )

      write.csv(
        resumo_sens,
        file.path(
          output_path,
          "sensibilidade_modelo_auxiliar_resumo.csv"
        ),
        row.names = FALSE
      )

      write.csv(
        metricas_concordancia,
        file.path(
          output_path,
          "sensibilidade_modelo_auxiliar_concordancia.csv"
        ),
        row.names = FALSE
      )

      write.csv(
        auditoria_sens,
        file.path(
          output_path,
          "sensibilidade_modelo_auxiliar_ASVs.csv"
        ),
        row.names = FALSE
      )

      cat("\nResumo da sensibilidade:\n")
      print(resumo_sens, row.names = FALSE)

      cat("\nConcordancia final:\n")
      print(metricas_concordancia, row.names = FALSE)

      rm(
        dadaF_transf,
        dadaR_transf,
        mergers_transf,
        seqtab_aux_transf_pre
      )

    } else {
      warning(
        "Sensibilidade nao executada: ",
        "dereps, modelos ou tabela auxiliar indisponiveis."
      )
    }

    rm(
      list = intersect(
        c("derepFs_aux", "derepRs_aux"),
        ls()
      )
    )

    gc()
    ###############################################################################
    # 16. RASTREAMENTO GLOBAL DE READS
    # "final" = reads na saida filtrada de cada amostra (main ou aux)
    ###############################################################################
    cat("\n=== RASTREAMENTO DE READS ===\n")

    final_per_sample <- c(rowSums(seqtab_main_nochim),
                          rowSums(seqtab_aux_nochim))

    track_all <- do.call(rbind, track_list)
    # Incluir todas as amostras que passaram pelo pipeline; 0 se filtradas inteiras
    all_track_samples <- rownames(track_all)
    final_col <- setNames(rep(0L, length(all_track_samples)), all_track_samples)
    final_col[names(final_per_sample)] <- as.integer(final_per_sample)
    track_all <- cbind(track_all, final = final_col)

    raw_col   <- reads_raw_R1[rownames(track_all)]
    track_all <- cbind(raw = raw_col, track_all)

    track_pct <- round(sweep(track_all, 1, track_all[, "raw"], "/") * 100, 1)
    cat("\nAbsoluto:\n"); print(track_all)
    cat("\nPercentual (em relacao as reads brutas):\n"); print(track_pct)
    global_retention <- round(sum(track_all[, "final"]) /
                                sum(track_all[, "raw"], na.rm = TRUE) * 100, 1)
    cat("\nRetencao global (bruto->final, todas as amostras):", global_retention, "%\n")

    filter_stats <- do.call(rbind, filter_list)
    merge_stats  <- do.call(rbind, merge_list)

    # Retencao sequencial: cada etapa sobre a IMEDIATAMENTE anterior.
    # Usa track_all (denominadores encadeados), nao a subtracao de percentuais
    # com bases distintas.
    retencao_seq <- data.frame(
      SampleID            = rownames(track_all),
      Run                 = run_of_sample[rownames(track_all)],
      bruto_para_input    = round(100 * track_all[, "input"]     / track_all[, "raw"],       1),
      input_para_filtrado = round(100 * track_all[, "filtradas"] / track_all[, "input"],     1),
      filtrado_para_denoised = round(100 * track_all[, "denoised_R1"] / track_all[, "filtradas"], 1),
      denoised_para_merge = round(100 * track_all[, "merged"]    / track_all[, "denoised_R1"], 1),
      merge_para_final    = round(100 * track_all[, "final"]     / track_all[, "merged"],     1),
      bruto_para_final    = round(100 * track_all[, "final"]     / track_all[, "raw"],        1),
      row.names = NULL
    )
    write.csv(retencao_seq,
              file.path(output_path, "retencao_sequencial.csv"), row.names = FALSE)

    # Alerta sobre a etapa de merge com denominador denoised -> merge.
    quedas_merge <- retencao_seq$SampleID[
      !is.na(retencao_seq$denoised_para_merge) & retencao_seq$denoised_para_merge < 70]
    if (length(quedas_merge) > 0)
      warning("Amostras com retencao denoised->merge < 70%: ",
              paste(quedas_merge, collapse = ", "))

    ###############################################################################
    # 17. CONTAMINANTES CONHECIDOS
    ###############################################################################
    contaminant_genera_ref <- c(
      "Ralstonia", "Pseudomonas", "Bradyrhizobium", "Cutibacterium",
      "Propionibacterium", "Acinetobacter", "Sphingomonas", "Stenotrophomonas",
      "Delftia", "Comamonas", "Methylobacterium", "Herbaspirillum",
      "Burkholderia", "Escherichia", "Shigella", "Staphylococcus")
    writeLines(contaminant_genera_ref,
               file.path(output_path, "REF_generos_contaminacao_conhecida.txt"))
    cat("\n[Sem controle negativo] Lista de contaminantes exportada.\n")

    ###############################################################################
    # 18. DICIONARIO GLOBAL DE IDs, FASTA UNIFICADO E TABELAS
    #
    # DESIGN:
    #   O universo de ASVs e particionado em tres subconjuntos:
    #     (a) principal          : colnames(seqtab_main_nochim)  → IDs ASV_1…N
    #     (b) compartilhadas     : intersecao entre principal e auxiliar
    #     (c) exclusivas_aux     : colnames(seqtab_aux_nochim) \ principal → IDs ASV_AUX_1…M
    #
    #   O FASTA global (ASVs.fa) contem (a) + (c), garantindo que o Script 4
    #   (rBLAST) processe a totalidade das sequencias de ambas as corridas.
    #   O prefixo ASV_AUX_ distingue as exclusivas da auxiliar sem alterar os
    #   IDs ASV_N da corrida principal. Os scripts a jusante preservam a unidade ASV.
    #
    #   Scripts 2a-2e e 4: operam sobre seqtab_global_nochim.rds (main + aux),
    #   classificando TODAS as ASVs em uma unica passagem. O particionamento por
    #   origem (principal vs. auxiliar_exclusiva) e preservado em ASV_sequences.tsv
    #   e pode ser auditado por todos os scripts a jusante.
    ###############################################################################
    cat("\n=== DICIONARIO GLOBAL DE IDs, FASTA E TABELAS ===\n")

    asv_seqs_main    <- colnames(seqtab_main_nochim)
    asv_seqs_aux_all <- colnames(seqtab_aux_nochim)

    # Particionamento
    asv_seqs_exclusivas_aux  <- setdiff(asv_seqs_aux_all, asv_seqs_main)
    asv_seqs_compartilhadas  <- intersect(asv_seqs_aux_all, asv_seqs_main)

    cat(sprintf("ASVs corrida principal             : %d\n", length(asv_seqs_main)))
    cat(sprintf("ASVs corrida auxiliar (total)      : %d\n", length(asv_seqs_aux_all)))
    cat(sprintf("  compartilhadas com principal     : %d\n", length(asv_seqs_compartilhadas)))
    cat(sprintf("  exclusivas da auxiliar           : %d\n", length(asv_seqs_exclusivas_aux)))

    # Esquema de IDs:
    #   ASV_1 … ASV_N      : corrida principal (retrocompativel com Scripts 2-8)
    #   ASV_AUX_1 … ASV_AUX_M : exclusivas da auxiliar
    asv_ids_main         <- paste0("ASV_", seq_along(asv_seqs_main))
    asv_ids_excl_aux     <- if (length(asv_seqs_exclusivas_aux) > 0)
      paste0("ASV_AUX_", seq_along(asv_seqs_exclusivas_aux)) else character(0)

    # Dicionario global: sequencia -> ID
    seq2id_global <- c(
      setNames(asv_ids_main,     asv_seqs_main),
      setNames(asv_ids_excl_aux, asv_seqs_exclusivas_aux))

    # Renomear colunas dos seqtabs com IDs consistentes
    seqtab_main_nochim_ids           <- seqtab_main_nochim
    colnames(seqtab_main_nochim_ids) <- asv_ids_main

    seqtab_aux_nochim_ids            <- seqtab_aux_nochim
    colnames(seqtab_aux_nochim_ids)  <- seq2id_global[asv_seqs_aux_all]
    # Blindagem: todo colname da auxiliar deve ter ID resolvido. NA aqui faria
    # pivot_longer descartar a coluna silenciosamente, perdendo ASVs da tabela long.
    if (ncol(seqtab_aux_nochim_ids) > 0 && anyNA(colnames(seqtab_aux_nochim_ids)))
      stop("IDs nao resolvidos para ", sum(is.na(colnames(seqtab_aux_nochim_ids))),
           " ASV(s) da auxiliar. Verifique seq2id_global.")

    # --------------------------------------------------------------------------
    # FASTA GLOBAL E LOCAL via writeFasta (dada2 p.63)
    # --------------------------------------------------------------------------
    all_seqs_global <- c(asv_seqs_main, asv_seqs_exclusivas_aux)
    all_ids_global  <- c(asv_ids_main,  asv_ids_excl_aux)

    # Fasta Global (Entrada Script 4 - rBLAST)
    writeLines(c(rbind(paste0(">", all_ids_global), all_seqs_global)),
               file.path(output_path, "ASVs.fa"))

    # FASTA restrito corrida principal
    writeLines(c(rbind(paste0(">", asv_ids_main), asv_seqs_main)),
               file.path(output_path, "ASVs_principal.fa"))

    # FASTA restrito corrida auxiliar (se houver exclusivas)
    if (length(asv_seqs_exclusivas_aux) > 0)
      writeLines(c(rbind(paste0(">", asv_ids_excl_aux), asv_seqs_exclusivas_aux)),
                 file.path(output_path, "ASVs_auxiliar_exclusivas.fa"))

    # --------------------------------------------------------------------------
    # Tabela de sequencias global com anotacao de origem
    # --------------------------------------------------------------------------
    asv_headers_global <- data.frame(
      ASV_ID   = all_ids_global,
      Sequence = all_seqs_global,
      Length   = nchar(all_seqs_global),
      Origem   = c(rep("principal",        length(asv_seqs_main)),
                   rep("auxiliar_exclusiva", length(asv_seqs_exclusivas_aux))),
      stringsAsFactors = FALSE)
    write.table(asv_headers_global,
                file.path(output_path, "ASV_sequences.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # --------------------------------------------------------------------------
    # Tabela de contagens de ASVs e tabela long — corrida principal
    # --------------------------------------------------------------------------
    write.table(seqtab_main_nochim_ids,
                file.path(output_path, "ASV_table.tsv"),
                sep = "\t", quote = FALSE, col.names = NA)

    asv_long <- as.data.frame(seqtab_main_nochim_ids)
    asv_long$SampleID <- rownames(asv_long)
    asv_long <- tidyr::pivot_longer(asv_long, cols = starts_with("ASV_"),
                                    names_to = "ASV_ID", values_to = "Reads")
    asv_long <- asv_long[asv_long$Reads > 0, ]
    write.table(asv_long, file.path(output_path, "ASV_table_long.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Tabela long da corrida auxiliar com IDs consistentes
    if (nrow(seqtab_aux_nochim_ids) > 0) {
      asv_long_aux <- as.data.frame(seqtab_aux_nochim_ids)
      asv_long_aux$SampleID <- rownames(asv_long_aux)
      asv_long_aux <- tidyr::pivot_longer(
        asv_long_aux,
        cols = dplyr::starts_with("ASV_"),
        names_to = "ASV_ID", values_to = "Reads")
      asv_long_aux <- asv_long_aux[asv_long_aux$Reads > 0, ]
      write.table(asv_long_aux,
                  file.path(output_path, "ASV_table_auxiliar_long.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      cat(sprintf("ASV_table_auxiliar_long.tsv: %d linhas (ASVs presentes na auxiliar)\n",
                  nrow(asv_long_aux)))
    }

    # --------------------------------------------------------------------------
    # seqtab_global_nochim: merge pos-filtro de main + aux
    # Uso: entrada canonica dos Scripts 2a-2e e 4.
    # Classifica TODAS as ASVs (main + exclusivas_aux) em uma unica passagem.
    # O particionamento por origem e feito no Script 7 via ASV_sequences.tsv.
    # --------------------------------------------------------------------------
    # A junção global é obrigatória: 9 amostras principais + 1 auxiliar.
    if (isTRUE(aux_vazia) || nrow(seqtab_aux_nochim) != 1L ||
        ncol(seqtab_aux_nochim) == 0L) {
      stop(
        "A tabela auxiliar final não satisfaz o contrato de 1 amostra com ASVs. ",
        "O seqtab global não será criado incompleto.",
        call. = FALSE
      )
    }
    seqtab_global_nochim <- mergeSequenceTables(
      seqtab_main_nochim,
      seqtab_aux_nochim
    )
    # A tabela global com IDs e derivada da tabela global por sequencia, evitando
    # usar mergeSequenceTables() sobre rotulos que nao sao sequencias biologicas.
    ids_globais <- unname(seq2id_global[colnames(seqtab_global_nochim)])
    if (anyNA(ids_globais) || anyDuplicated(ids_globais) > 0L) {
      stop("Falha ao mapear as sequencias globais para ASV_ID unicos.", call. = FALSE)
    }
    seqtab_global_nochim_ids <- seqtab_global_nochim
    colnames(seqtab_global_nochim_ids) <- ids_globais

    amostras_esperadas <- metadata$SampleID
    if (nrow(seqtab_global_nochim) != length(amostras_esperadas) ||
        !setequal(rownames(seqtab_global_nochim), amostras_esperadas)) {
      stop(
        "seqtab_global_nochim não preservou o universo amostral: esperadas=",
        length(amostras_esperadas), "; encontradas=", nrow(seqtab_global_nochim),
        call. = FALSE
      )
    }
    seqtab_global_nochim <- seqtab_global_nochim[
      amostras_esperadas, , drop = FALSE
    ]
    seqtab_global_nochim_ids <- seqtab_global_nochim_ids[
      amostras_esperadas, , drop = FALSE
    ]
    saveRDS(seqtab_global_nochim,     file.path(output_path, "seqtab_global_nochim.rds"))
    saveRDS(seqtab_global_nochim_ids, file.path(output_path, "seqtab_global_nochim_ids.rds"))

    cat(sprintf("seqtab_global_nochim.rds: %d amostras | %d ASVs (colunas = sequencias)\n",
                nrow(seqtab_global_nochim), ncol(seqtab_global_nochim)))
    cat(sprintf("seqtab_global_nochim_ids.rds: %d amostras | %d ASVs (colunas = IDs)\n",
                nrow(seqtab_global_nochim_ids), ncol(seqtab_global_nochim_ids)))
    # --------------------------------------------------------------------------
    # Saidas primarias
    # --------------------------------------------------------------------------
    saveRDS(seqtab_main_nochim, file.path(output_path, "seqtab_nochim.rds"))
    saveRDS(seqtab_aux_nochim,  file.path(output_path, "seqtab_auxiliar.rds"))

    cat(sprintf("\nseqtab_nochim.rds   : %d amostras | %d ASVs (principal — auditoria e analise core9)\n",
                nrow(seqtab_main_nochim), ncol(seqtab_main_nochim)))
    cat(sprintf("seqtab_auxiliar.rds : %d amostra  | %d ASVs (auxiliar  — analise suplementar)\n",
                nrow(seqtab_aux_nochim),  ncol(seqtab_aux_nochim)))
    cat("AVISO: seqtab_auxiliar.rds usa modelo de erros PROPRIO da corrida auxiliar\n")
    cat("       (amostra unica; ajuste inspecionado em modelo_erros_*_proprio.pdf).\n")
    cat("       Para classificacao taxonomica completa da auxiliar (incluindo\n")
    cat("       ASVs exclusivas), usar seqtab_global_nochim.rds nos Scripts 2a-2e e 4.\n")
    cat("       Declarar ambas as limitacoes na secao de Material e Metodos.\n")

    ###############################################################################
    # 19. ESTATISTICAS E RELATORIO FINAL
    ###############################################################################
    write.csv(filter_stats, file.path(output_path, "filter_stats.csv"), row.names = FALSE)
    write.csv(merge_stats,  file.path(output_path, "merge_stats.csv"),  row.names = FALSE)
    write.csv(track_all,    file.path(output_path, "track_reads.csv"))
    write.csv(track_pct,    file.path(output_path, "track_reads_percentual.csv"))
    write.csv(as.data.frame(length_distribution),
              file.path(output_path, "length_distribution.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(log_path, "sessionInfo.txt"))

    report_file <- file.path(output_path,
                             paste0("pipeline_report_", pipeline_version, ".txt"))
    sink(report_file)
    tryCatch({
      cat(
        "RELATORIO FINAL — PIPELINE DADA2 (multi-run) ",pipeline_version, "\n\n", sep = ""
      )
      cat("dada2:", as.character(packageVersion("dada2")), "\n\n")
      cat("Corridas processadas:\n"); print(table(metadata$Run))
      cat(sprintf("\nCorrida principal : %s (%d amostras) — learnErrors proprio\n",
                  run_principal, run_n_amostras[run_principal]))
      cat(sprintf(  "Corrida auxiliar  : %s (%d amostra) — learnErrors proprio\n",
                  run_auxiliar, run_n_amostras[run_auxiliar]))
      cat("\nPolitica de pooling por corrida:\n")
      print(pool_by_run)
      cat("Fallback de pooling:", pool_fallback, "\n")
      cat("Pool por corrida:\n")
      for (r in names(pool_used)) cat(" ", r, ":", pool_used[[r]], "\n")
      cat(sprintf("\nAmostras finais (principal): %d\n", nrow(seqtab_main_nochim)))
      cat(sprintf("ASVs finais     (principal): %d  [IDs ASV_1…ASV_%d]\n",
                  ncol(seqtab_main_nochim), ncol(seqtab_main_nochim)))
      cat(sprintf("Amostras finais (auxiliar) : %d\n", nrow(seqtab_aux_nochim)))
      cat(sprintf("ASVs finais     (auxiliar) : %d  [compartilhadas: %d | exclusivas: %d]\n",
                  ncol(seqtab_aux_nochim),
                  length(asv_seqs_compartilhadas),
                  length(asv_seqs_exclusivas_aux)))
      cat(sprintf("ASVs.fa (FASTA global)     : %d sequencias\n", length(all_ids_global)))
      cat("Reads brutas (pre-cutadapt):",
          format(sum(reads_raw_R1, na.rm = TRUE), big.mark = ","), "\n")
      cat("Retencao global:", global_retention, "%\n")
      cat("Retencao pos-quimera (combinado):", pct_nochim, "%\n")
      cat("Filtro comprimento:", expected_min_len, "-", expected_max_len, "pb\n")
      cat("Filtro principal: reads <", min_total_reads,
          "E prevalencia <", min_prevalence, "/ 9 amostras (logica E)\n")
      cat("Filtro auxiliar: reads <", min_total_reads, "apenas\n")
      cat(sprintf("Profundidade minima principal (rarefacao): %d reads — %s\n",
                  min(depth_final_main), names(which.min(depth_final_main))))
      cat("Controle negativo: AUSENTE.\n")
      cat("\nLIMITACAO METODOLOGICA:\n")
      cat("  O modelo de erros da corrida auxiliar (", run_auxiliar, ") foi estimado\n")
      cat("  a partir de uma unica amostra. Embora a retencao de reads tenha sido\n")
      cat("  elevada, a menor diversidade de amostras pode reduzir a estabilidade\n")
      cat("  do ajuste em faixas de qualidade com poucas observacoes. A amostra foi\n")
      cat("  adicionalmente processada com o modelo da corrida principal como analise\n")
      cat("  de sensibilidade (ver sensibilidade_modelo_auxiliar_resumo.csv). A amostra\n")
      cat("  auxiliar nao foi incluida nos testes estatisticos principais entre grupos.\n")
      print(sessionInfo())
    }, finally = {
      sink()
    })

    ###############################################################################
    # 20. FINALIZACAO
    ###############################################################################
    cat("\n=== PIPELINE DADA2 (multi-run) FINALIZADO ===\n")
    cat("Saidas em:", output_path, "\n")
    cat("-- OBJETOS CANONICOS E ANALISE PRINCIPAL --\n")
    cat(" - seqtab_nochim.rds            (principal, pos-filtro — colunas = sequencias)\n")
    cat(" - ASVs.fa                      (principal + exclusivas_aux — FASTA para BLAST)\n")
    cat(" - ASV_table.tsv / ASV_table_long.tsv (contagens por ASV; apenas corrida principal)\n")
    cat(" - ASV_sequences.tsv            (todos os IDs com coluna Origem)\n")
    cat("-- ANALISE SUPLEMENTAR (auxiliar) --\n")
    cat(" - seqtab_auxiliar.rds          (auxiliar, pos-filtro; modelo de erros proprio; modelo principal apenas em sensibilidade)\n")
    cat(" - seqtab_global_nochim.rds     (main + aux — colunas = sequencias; entrada Scripts 2-3)\n")
    cat(" - seqtab_global_nochim_ids.rds (main + aux — colunas = IDs ASV/ASV_AUX)\n")
    cat(" - ASVs_principal.fa            (FASTA so da corrida principal)\n")
    cat(" - ASVs_auxiliar_exclusivas.fa  (FASTA das ASVs exclusivas da auxiliar)\n")
    cat(" - ASV_table_auxiliar_long.tsv\n")
    cat("-- UTILITARIOS --\n")
    cat(" - seqtab_nochim_prefiltro.rds  (combinado, pre-filtro — auditoria/sensibilidade; nao e entrada do Script 9)\n")
    cat(" - track_reads.csv              (raw->input->filtradas->denoised->merged->final)\n")
    cat(" - sensitivity_thresholds.csv   (avaliado sobre corrida principal, n=9)\n")
    cat(" - cutadapt_retention.csv / primer_residual_check.csv\n")
    cat("Proximo passo: 01b_auditoria_quimeras_recuperacao_ASVs.R (validar limiares de filtro)\n")

    .resume_status <- "SUCESSO"
  }, warning = function(w) {
    cat("[WARN] ", conditionMessage(w), "\n", sep = "", file = stderr())
  }),
  error = function(e) {
    .resume_erro <<- conditionMessage(e)
    stop(e)
  },
  finally = {
    .resume_fim <- Sys.time()
    cat("\nSTATUS_RETOMADA=", .resume_status, "\n", sep = "")
    cat("INICIO=", format(.resume_inicio, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("FIM=", format(.resume_fim, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("DURACAO_SEGUNDOS=",
        round(as.numeric(difftime(.resume_fim, .resume_inicio, units = "secs")), 3),
        "\n", sep = "")
    if (!is.na(.resume_erro)) cat("ERRO=", .resume_erro, "\n", sep = "")
    while (sink.number(type = "message") > 0L) sink(type = "message")
    while (sink.number(type = "output") > 0L) sink(type = "output")
    close(.resume_con)
  }
)

cat("Retomada concluida. Log: ", .resume_log, "\n", sep = "")
