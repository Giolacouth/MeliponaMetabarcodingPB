###############################################################################
# FUNCOES ESTATISTICAS EXATAS COMPARTILHADAS
#
#
###############################################################################

FUNCOES_EXATAS_VERSAO <- "1.3.0"

.numero_alocacoes_fixadas <- function(tamanhos) {
  tamanhos <- as.integer(tamanhos)
  if (!length(tamanhos) || anyNA(tamanhos) || any(tamanhos < 1L)) return(NA_real_)
  valor <- exp(lfactorial(sum(tamanhos)) - sum(lfactorial(tamanhos)))
  if (!is.finite(valor)) return(valor)
  round(valor)
}

matriz_permutacoes_rotulos_unicos <- function(
    grupo,
    incluir_observada = FALSE,
    max_alocacoes = 100000L
) {
  g <- droplevels(as.factor(grupo))
  if (anyNA(g)) stop("grupo contem NA; remova-os antes da enumeracao.", call. = FALSE)
  if (nlevels(g) < 2L) stop("Sao necessarios ao menos dois grupos.", call. = FALSE)

  niveis <- levels(g)
  tamanhos <- as.integer(table(g))
  n <- length(g)
  n_aloc <- .numero_alocacoes_fixadas(tamanhos)
  if (!is.finite(n_aloc) || n_aloc > as.numeric(max_alocacoes)) {
    stop(
      "Numero de alocacoes rotuladas excede o limite: ", n_aloc,
      " > ", max_alocacoes, ".", call. = FALSE
    )
  }

  n_linhas <- as.integer(n_aloc - if (isTRUE(incluir_observada)) 0L else 1L)
  if (n_linhas < 1L) stop("Nenhuma permutacao alternativa disponivel.", call. = FALSE)

  out <- matrix(NA_integer_, nrow = n_linhas, ncol = n)
  colnames(out) <- if (!is.null(names(grupo))) names(grupo) else as.character(seq_len(n))

  indices_fonte <- lapply(niveis, function(nv) which(g == nv))
  grupo_obs_idx <- match(g, niveis)
  atribuicao <- integer(n)
  linha <- 0L
  k <- length(tamanhos)

  emitir <- function() {
    if (!isTRUE(incluir_observada) && all(atribuicao == grupo_obs_idx)) {
      return(invisible(NULL))
    }
    perm <- integer(n)
    for (j in seq_len(k)) {
      pos_destino <- which(atribuicao == j)
      perm[pos_destino] <- indices_fonte[[j]]
    }
    linha <<- linha + 1L
    out[linha, ] <<- perm
    invisible(NULL)
  }

  enumerar <- function(indices_livres, pos_grupo) {
    if (pos_grupo == k) {
      atribuicao[indices_livres] <<- k
      emitir()
      return(invisible(NULL))
    }
    escolhas <- combn(indices_livres, tamanhos[pos_grupo], simplify = FALSE)
    for (sel in escolhas) {
      atribuicao[sel] <<- pos_grupo
      enumerar(setdiff(indices_livres, sel), pos_grupo + 1L)
    }
    invisible(NULL)
  }

  enumerar(seq_len(n), 1L)
  if (linha != n_linhas || anyNA(out)) {
    stop(
      "Falha interna ao gerar permutacoes unicas: esperado ", n_linhas,
      ", obtido ", linha, ".", call. = FALSE
    )
  }

  # Validacao defensiva da matriz entregue a vegan/indicspecies.
  indice_esperado <- seq_len(n)
  linhas_validas <- apply(
    out, 1L,
    function(idx) identical(sort(as.integer(idx)), indice_esperado)
  )
  if (!all(linhas_validas)) {
    stop("Matriz customizada contem linha que nao e permutacao completa.", call. = FALSE)
  }
  chaves_rotulos <- apply(
    out, 1L,
    function(idx) paste(as.integer(g[as.integer(idx)]), collapse = ",")
  )
  if (anyDuplicated(chaves_rotulos) > 0L) {
    stop("Matriz customizada contem alocacoes de rotulos duplicadas.", call. = FALSE)
  }
  chave_observada <- paste(as.integer(g), collapse = ",")
  if (!isTRUE(incluir_observada) && chave_observada %in% chaves_rotulos) {
    stop("Ordenacao observada foi incluida indevidamente na matriz customizada.", call. = FALSE)
  }

  attr(out, "n_alocacoes_totais") <- as.numeric(n_aloc)
  attr(out, "n_permutacoes_alternativas") <- as.numeric(n_linhas)
  attr(out, "p_min_teorico") <- 1 / as.numeric(n_aloc)
  attr(out, "desenho") <- paste(niveis, tamanhos, sep = "=", collapse = "; ")
  attr(out, "inclui_observada") <- isTRUE(incluir_observada)
  attr(out, "metodo") <- "Enumeracao completa das alocacoes unicas dos rotulos com tamanhos fixos"
  out
}

resumo_permutacoes_rotulos <- function(matriz, conjunto, fator) {
  if (!is.matrix(matriz)) stop("matriz deve ser uma matriz de permutacoes.", call. = FALSE)
  data.frame(
    Conjunto = conjunto,
    Fator = fator,
    Desenho = attr(matriz, "desenho"),
    N_alocacoes_rotulos_totais = attr(matriz, "n_alocacoes_totais"),
    N_permutacoes_alternativas = nrow(matriz),
    Inclui_ordenacao_observada = isTRUE(attr(matriz, "inclui_observada")),
    p_min_teorico = attr(matriz, "p_min_teorico"),
    Metodo = attr(matriz, "metodo"),
    stringsAsFactors = FALSE
  )
}

.preparar_teste_postos <- function(y, grupo) {
  ok <- is.finite(y) & !is.na(grupo)
  y <- as.numeric(y[ok])
  g <- droplevels(as.factor(grupo[ok]))

  if (length(y) < 2L) stop("Menos de duas observacoes validas.", call. = FALSE)
  if (nlevels(g) < 2L) stop("Menos de dois grupos validos.", call. = FALSE)
  if (any(table(g) < 1L)) stop("Grupo vazio apos remocao de NA.", call. = FALSE)

  ranks <- rank(y, ties.method = "average")
  tab_empates <- table(y)
  n <- length(y)
  correcao_empates <- if (n > 1L) {
    1 - sum(tab_empates^3 - tab_empates) / (n^3 - n)
  } else {
    NA_real_
  }

  list(
    y = y,
    g = g,
    ranks = ranks,
    tamanhos = as.integer(table(g)),
    niveis = levels(g),
    correcao_empates = as.numeric(correcao_empates),
    todos_iguais = length(unique(y)) == 1L
  )
}

.calcular_h_kw_por_indices <- function(ranks, grupo_idx, tamanhos, correcao_empates) {
  if (!is.finite(correcao_empates) || correcao_empates <= 0) return(0)
  n <- length(ranks)
  k <- length(tamanhos)
  somas <- numeric(k)
  for (j in seq_len(k)) somas[j] <- sum(ranks[grupo_idx == j])
  h_sem_correcao <- 12 / (n * (n + 1)) * sum((somas^2) / tamanhos) - 3 * (n + 1)
  max(0, as.numeric(h_sem_correcao / correcao_empates))
}

# Teste exato por enumeracao completa para dois ou mais grupos. A estatistica e
# o H de Kruskal-Wallis com a mesma correcao de empates do teste padrao. Para
# dois grupos, equivale a um teste exato de permutacao baseado em postos.
teste_postos_exato_grupos_fixos <- function(
    y,
    grupo,
    tamanhos_esperados = NULL,
    max_alocacoes = 1000000L
) {
  erro <- NA_character_
  prep <- tryCatch(
    .preparar_teste_postos(y, grupo),
    error = function(e) {
      erro <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(prep)) {
    return(data.frame(
      H = NA_real_, df = NA_integer_, p_exato = NA_real_,
      p_assintotico_diagnostico = NA_real_, N_permutacoes_exatas = NA_real_,
      p_min_teorico = NA_real_, Epsilon2 = NA_real_,
      Correcao_empates = NA_real_, Todos_valores_iguais = NA,
      Desenho = NA_character_, Metodo = "Teste exato por enumeracao",
      Nota = "Teste nao executado.", Erro = erro,
      Versao_funcao = FUNCOES_EXATAS_VERSAO, stringsAsFactors = FALSE
    ))
  }

  tamanhos <- prep$tamanhos
  desenho <- paste(prep$niveis, tamanhos, sep = "=", collapse = "; ")

  if (!is.null(tamanhos_esperados)) {
    esperados <- sort(as.integer(tamanhos_esperados))
    if (!identical(sort(tamanhos), esperados)) {
      return(data.frame(
        H = NA_real_, df = length(tamanhos) - 1L, p_exato = NA_real_,
        p_assintotico_diagnostico = NA_real_, N_permutacoes_exatas = NA_real_,
        p_min_teorico = NA_real_, Epsilon2 = NA_real_,
        Correcao_empates = prep$correcao_empates,
        Todos_valores_iguais = prep$todos_iguais,
        Desenho = desenho, Metodo = "Teste exato por enumeracao",
        Nota = paste0(
          "Desenho observado diferente do esperado: ", desenho,
          "; esperado=", paste(esperados, collapse = "/"), "."
        ),
        Erro = NA_character_, Versao_funcao = FUNCOES_EXATAS_VERSAO,
        stringsAsFactors = FALSE
      ))
    }
  }

  n_aloc <- .numero_alocacoes_fixadas(tamanhos)
  if (!is.finite(n_aloc) || n_aloc > as.numeric(max_alocacoes)) {
    return(data.frame(
      H = NA_real_, df = length(tamanhos) - 1L, p_exato = NA_real_,
      p_assintotico_diagnostico = NA_real_, N_permutacoes_exatas = n_aloc,
      p_min_teorico = if (is.finite(n_aloc) && n_aloc > 0) 1 / n_aloc else NA_real_,
      Epsilon2 = NA_real_, Correcao_empates = prep$correcao_empates,
      Todos_valores_iguais = prep$todos_iguais, Desenho = desenho,
      Metodo = "Teste exato por enumeracao",
      Nota = paste0("Numero de alocacoes excede o limite: ", n_aloc, "."),
      Erro = NA_character_, Versao_funcao = FUNCOES_EXATAS_VERSAO,
      stringsAsFactors = FALSE
    ))
  }

  n <- length(prep$y)
  k <- length(tamanhos)
  grupo_obs_idx <- match(prep$g, prep$niveis)
  h_obs <- .calcular_h_kw_por_indices(
    prep$ranks, grupo_obs_idx, tamanhos, prep$correcao_empates
  )

  if (isTRUE(prep$todos_iguais)) {
    return(data.frame(
      H = 0, df = k - 1L, p_exato = 1,
      p_assintotico_diagnostico = 1,
      N_permutacoes_exatas = as.numeric(n_aloc),
      p_min_teorico = 1 / as.numeric(n_aloc), Epsilon2 = 0,
      Correcao_empates = prep$correcao_empates,
      Todos_valores_iguais = TRUE, Desenho = desenho,
      Metodo = "Teste exato por enumeracao completa dos rotulos",
      Nota = "Todos os valores sao iguais; H=0 e p=1 por definicao operacional.",
      Erro = NA_character_, Versao_funcao = FUNCOES_EXATAS_VERSAO,
      stringsAsFactors = FALSE
    ))
  }

  total <- 0L
  maiores_iguais <- 0L
  atribuicao <- integer(n)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(h_obs))

  enumerar <- function(indices_livres, pos_grupo) {
    if (pos_grupo == k) {
      atribuicao[indices_livres] <<- k
      h_perm <- .calcular_h_kw_por_indices(
        prep$ranks, atribuicao, tamanhos, prep$correcao_empates
      )
      total <<- total + 1L
      if (h_perm >= h_obs - tol) maiores_iguais <<- maiores_iguais + 1L
      return(invisible(NULL))
    }

    escolhas <- combn(indices_livres, tamanhos[pos_grupo], simplify = FALSE)
    for (sel in escolhas) {
      atribuicao[sel] <<- pos_grupo
      enumerar(setdiff(indices_livres, sel), pos_grupo + 1L)
    }
    invisible(NULL)
  }

  enumerar(seq_len(n), 1L)

  if (!isTRUE(all.equal(as.numeric(total), as.numeric(n_aloc), tolerance = 0))) {
    stop(
      "Falha interna na enumeracao: esperado ", n_aloc,
      " e obtido ", total, ".", call. = FALSE
    )
  }

  p_exato <- maiores_iguais / total
  df <- k - 1L
  p_assint <- stats::pchisq(h_obs, df = df, lower.tail = FALSE)
  denom <- n - k
  eps2 <- if (denom > 0L) max(0, min(1, (h_obs - k + 1) / denom)) else NA_real_

  data.frame(
    H = h_obs,
    df = df,
    p_exato = p_exato,
    p_assintotico_diagnostico = p_assint,
    N_permutacoes_exatas = as.numeric(total),
    p_min_teorico = 1 / as.numeric(total),
    Epsilon2 = eps2,
    Correcao_empates = prep$correcao_empates,
    Todos_valores_iguais = FALSE,
    Desenho = desenho,
    Metodo = "Teste de postos exato por enumeracao completa dos rotulos",
    Nota = paste0(
      "Todos os ", total,
      " arranjos rotulados com tamanhos fixos foram enumerados; ",
      "o p assintotico e apenas diagnostico."
    ),
    Erro = NA_character_,
    Versao_funcao = FUNCOES_EXATAS_VERSAO,
    stringsAsFactors = FALSE
  )
}

# Interface estrita usada no core9 por BeeSpecies.
teste_kw_exato_234 <- function(y, grupo) {
  teste_postos_exato_grupos_fixos(
    y = y,
    grupo = grupo,
    tamanhos_esperados = c(2L, 3L, 4L),
    max_alocacoes = 1260L
  )
}

# Comparacoes pareadas exatas apos um teste global. Cada par enumera todas as
# realocacoes mantendo seus dois tamanhos. Ajuste multiplo deve ser feito fora
# ou por esta funcao com p.adjust.method.
posthoc_postos_exato_pares <- function(
    y,
    grupo,
    p.adjust.method = "holm",
    max_alocacoes = 1000000L
) {
  ok <- is.finite(y) & !is.na(grupo)
  y <- as.numeric(y[ok])
  g <- droplevels(as.factor(grupo[ok]))
  pares <- combn(levels(g), 2L, simplify = FALSE)

  out <- lapply(pares, function(par) {
    sel <- g %in% par
    z <- teste_postos_exato_grupos_fixos(
      y[sel], droplevels(g[sel]), max_alocacoes = max_alocacoes
    )
    data.frame(
      Grupo1 = par[1L],
      Grupo2 = par[2L],
      Comparison = paste(par, collapse = " versus "),
      N1 = sum(g[sel] == par[1L]),
      N2 = sum(g[sel] == par[2L]),
      z,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$p_ajustado <- stats::p.adjust(out$p_exato, method = p.adjust.method)
  out$Metodo_ajuste <- p.adjust.method
  rownames(out) <- NULL
  out
}

###############################################################################
# PERMANOVA EXATA DE UM FATOR POR ENUMERACAO DAS ALOCACOES
###############################################################################

.preparar_permanova_exata <- function(dist_obj, grupo) {
  if (!inherits(dist_obj, "dist")) stop("dist_obj deve ser da classe dist.", call. = FALSE)
  labs <- attr(dist_obj, "Labels")
  if (is.null(labs) || anyNA(labs) || any(!nzchar(labs)) || anyDuplicated(labs)) {
    stop("A distancia deve possuir Labels validos e unicos.", call. = FALSE)
  }
  if (is.null(names(grupo))) {
    if (length(grupo) != length(labs)) stop("grupo sem nomes e comprimento divergente.", call. = FALSE)
    names(grupo) <- labs
  }
  if (anyDuplicated(names(grupo))) stop("grupo possui nomes duplicados.", call. = FALSE)
  g <- grupo[labs]
  if (anyNA(g)) stop("Existem amostras da distancia sem grupo.", call. = FALSE)
  g <- droplevels(as.factor(g))
  if (nlevels(g) < 2L) stop("Menos de dois grupos.", call. = FALSE)
  tamanhos <- as.integer(table(g))
  if (any(tamanhos < 2L)) stop("PERMANOVA requer ao menos duas unidades por grupo.", call. = FALSE)

  d <- as.matrix(dist_obj)
  if (anyNA(d) || any(!is.finite(d)) || any(d < 0) ||
      max(abs(d - t(d))) > 1e-10 || max(abs(diag(d))) > 1e-10) {
    stop("Matriz de distancia invalida.", call. = FALSE)
  }
  n <- nrow(d)
  a <- -0.5 * (d^2)
  j <- diag(n) - matrix(1 / n, n, n)
  gower <- j %*% a %*% j
  list(
    labels = labs, grupo = g, niveis = levels(g), tamanhos = tamanhos,
    n = n, k = nlevels(g), gower = gower,
    total_ss = sum(diag(gower))
  )
}

.pseudo_f_permanova <- function(gower, grupo_idx, k, total_ss) {
  n <- nrow(gower)
  x <- stats::model.matrix(~ 0 + factor(grupo_idx, levels = seq_len(k)))
  h_full <- x %*% solve(crossprod(x), t(x))
  h0 <- matrix(1 / n, n, n)
  ss_model <- sum(diag((h_full - h0) %*% gower))
  ss_res <- total_ss - ss_model
  df_model <- k - 1L
  df_res <- n - k
  f <- if (df_res > 0L && ss_res > 0) {
    (ss_model / df_model) / (ss_res / df_res)
  } else if (ss_model > 0 && ss_res <= sqrt(.Machine$double.eps)) {
    Inf
  } else {
    0
  }
  c(F = as.numeric(f), R2 = if (total_ss > 0) ss_model / total_ss else NA_real_,
    SS_model = ss_model, SS_residual = ss_res)
}

permanova_exata_grupos_fixos <- function(
    dist_obj,
    grupo,
    tamanhos_esperados = NULL,
    max_alocacoes = 100000L
) {
  erro <- NA_character_
  prep <- tryCatch(
    .preparar_permanova_exata(dist_obj, grupo),
    error = function(e) { erro <<- conditionMessage(e); NULL }
  )
  if (is.null(prep)) {
    return(data.frame(
      F = NA_real_, R2 = NA_real_, df_model = NA_integer_, df_residual = NA_integer_,
      p_exato = NA_real_, N_alocacoes_exatas = NA_real_, p_min_teorico = NA_real_,
      Desenho = NA_character_, Metodo = "PERMANOVA exata por enumeracao",
      Nota = "Teste nao executado.", Erro = erro,
      Versao_funcao = FUNCOES_EXATAS_VERSAO, stringsAsFactors = FALSE
    ))
  }
  desenho <- paste(prep$niveis, prep$tamanhos, sep = "=", collapse = "; ")
  if (!is.null(tamanhos_esperados) &&
      !identical(sort(prep$tamanhos), sort(as.integer(tamanhos_esperados)))) {
    return(data.frame(
      F = NA_real_, R2 = NA_real_, df_model = prep$k - 1L,
      df_residual = prep$n - prep$k, p_exato = NA_real_,
      N_alocacoes_exatas = .numero_alocacoes_fixadas(prep$tamanhos),
      p_min_teorico = NA_real_, Desenho = desenho,
      Metodo = "PERMANOVA exata por enumeracao",
      Nota = "Desenho observado difere do esperado.", Erro = NA_character_,
      Versao_funcao = FUNCOES_EXATAS_VERSAO, stringsAsFactors = FALSE
    ))
  }
  n_aloc <- .numero_alocacoes_fixadas(prep$tamanhos)
  if (!is.finite(n_aloc) || n_aloc > max_alocacoes) {
    return(data.frame(
      F = NA_real_, R2 = NA_real_, df_model = prep$k - 1L,
      df_residual = prep$n - prep$k, p_exato = NA_real_,
      N_alocacoes_exatas = n_aloc,
      p_min_teorico = if (is.finite(n_aloc)) 1 / n_aloc else NA_real_,
      Desenho = desenho, Metodo = "PERMANOVA exata por enumeracao",
      Nota = paste0("Numero de alocacoes excede o limite: ", n_aloc, "."),
      Erro = NA_character_, Versao_funcao = FUNCOES_EXATAS_VERSAO,
      stringsAsFactors = FALSE
    ))
  }

  obs_idx <- match(prep$grupo, prep$niveis)
  obs <- .pseudo_f_permanova(prep$gower, obs_idx, prep$k, prep$total_ss)
  tol <- sqrt(.Machine$double.eps) * max(1, abs(obs[["F"]]))
  atribuicao <- integer(prep$n)
  total <- 0L
  maiores_iguais <- 0L

  enumerar <- function(indices_livres, pos_grupo) {
    if (pos_grupo == prep$k) {
      atribuicao[indices_livres] <<- prep$k
      z <- .pseudo_f_permanova(prep$gower, atribuicao, prep$k, prep$total_ss)
      total <<- total + 1L
      if (z[["F"]] >= obs[["F"]] - tol) maiores_iguais <<- maiores_iguais + 1L
      return(invisible(NULL))
    }
    escolhas <- combn(indices_livres, prep$tamanhos[pos_grupo], simplify = FALSE)
    for (sel in escolhas) {
      atribuicao[sel] <<- pos_grupo
      enumerar(setdiff(indices_livres, sel), pos_grupo + 1L)
    }
    invisible(NULL)
  }
  enumerar(seq_len(prep$n), 1L)
  if (total != n_aloc) stop("Falha interna na enumeracao PERMANOVA.", call. = FALSE)

  data.frame(
    F = unname(obs[["F"]]), R2 = unname(obs[["R2"]]),
    df_model = prep$k - 1L, df_residual = prep$n - prep$k,
    p_exato = maiores_iguais / total,
    N_alocacoes_exatas = as.numeric(total), p_min_teorico = 1 / total,
    Desenho = desenho,
    Metodo = "PERMANOVA exata de um fator por enumeracao completa das alocacoes",
    Nota = paste0(
      "Foram enumeradas ", total,
      " alocacoes unicas com tamanhos fixos. O teste nao corrige dependencia, ",
      "confundimento ou heterogeneidade de dispersao."
    ),
    Erro = NA_character_, Versao_funcao = FUNCOES_EXATAS_VERSAO,
    stringsAsFactors = FALSE
  )
}

permanova_exata_pares <- function(dist_obj, grupo, p.adjust.method = "holm",
                                  max_alocacoes = 100000L) {
  labs <- attr(dist_obj, "Labels")
  if (is.null(names(grupo))) names(grupo) <- labs
  g <- droplevels(as.factor(grupo[labs]))
  pares <- combn(levels(g), 2L, simplify = FALSE)
  dmat <- as.matrix(dist_obj)
  out <- lapply(pares, function(par) {
    ids <- labs[g %in% par]
    dsub <- stats::as.dist(dmat[ids, ids, drop = FALSE])
    gsub <- droplevels(g[g %in% par])
    names(gsub) <- ids
    z <- permanova_exata_grupos_fixos(dsub, gsub, max_alocacoes = max_alocacoes)
    data.frame(
      Grupo1 = par[[1L]], Grupo2 = par[[2L]],
      Comparison = paste(par, collapse = " versus "),
      N1 = sum(gsub == par[[1L]]), N2 = sum(gsub == par[[2L]]),
      z, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$p_ajustado <- stats::p.adjust(out$p_exato, method = p.adjust.method)
  out$Metodo_ajuste <- p.adjust.method
  rownames(out) <- NULL
  out
}

FUNCOES_EXATAS_VERSAO <- "1.3.0"
