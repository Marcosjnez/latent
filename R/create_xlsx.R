create_xlsx <- function(..., file, sheet_names = NULL, digits = 3L,
                        breaks = c(-0.50, -0.20, 0.20, 0.50),
                        colors = c("#FA8072", "#EA9999", NA,
                                   "#90EE90", "#6AA84F")) {

  #### Check arguments ####

  matrices <- list(...)
  expressions <- as.list(substitute(list(...)))[-1L]

  if(length(matrices) == 0L) {
    stop("At least one matrix or data.frame must be provided")
  }

  if(!all(vapply(matrices, is.matrix, logical(1L)) |
          vapply(matrices, is.data.frame, logical(1L)))) {
    stop("All objects supplied through ... must be matrices or data.frames")
  }

  if(length(digits) != 1L || !is.numeric(digits) || is.na(digits) ||
     digits < 0L || digits != as.integer(digits)) {
    stop("digits must be a non-negative integer")
  }

  digits <- as.integer(digits)

  if(length(breaks) < 1L || !is.numeric(breaks) ||
     any(!is.finite(breaks)) ||
     is.unsorted(breaks, strictly = TRUE)) {
    stop("breaks must contain increasing numeric values")
  }

  if(length(colors) != length(breaks)+1L) {
    stop("colors must contain one more value than breaks")
  }

  #### Normalize colors ####

  normalize_color <- function(x) {

    if(length(x) != 1L) {
      stop("Each color must contain a single value")
    }

    if(is.na(x)) {

      #### Result ####

      return(NA_character_)

    }

    x <- as.character(x)

    rgb <- tryCatch(
      grDevices::col2rgb(x),
      error = function(e) {
        stop("Invalid color: ", x)
      }
    )

    result <- grDevices::rgb(
      red = rgb[1L, 1L],
      green = rgb[2L, 1L],
      blue = rgb[3L, 1L],
      maxColorValue = 255
    )

    #### Result ####

    return(result)

  }

  colors <- vapply(colors, normalize_color, character(1L))

  #### Sheet names ####

  if(is.null(sheet_names)) {

    sheet_names <- names(matrices)

    if(is.null(sheet_names)) {
      sheet_names <- rep("", length(matrices))
    }

    unnamed <- which(is.na(sheet_names) | sheet_names == "")

    for(i in unnamed) {
      sheet_names[i] <- deparse(expressions[[i]], nlines = 1L)
    }

  } else {

    if(!is.character(sheet_names) ||
       length(sheet_names) != length(matrices)) {
      stop("sheet_names must contain one name for each matrix or data.frame")
    }

  }

  if(any(is.na(sheet_names)) || any(sheet_names == "")) {
    stop("All worksheets must have a valid name")
  }

  sheet_names <- gsub("[\\[\\]:*?/\\\\]", "_", sheet_names)
  sheet_names <- substr(sheet_names, 1L, 31L)
  sheet_names <- make.unique(sheet_names)

  #### Create workbook ####

  wb <- openxlsx::createWorkbook()

  styles <- lapply(colors, FUN = \(x) {

    if(is.na(x)) {

      #### Result ####

      return(NULL)

    }

    result <- openxlsx::createStyle(bgFill = x)

    #### Result ####

    return(result)

  })

  number_format <- if(digits == 0L) {
    "0"
  } else {
    paste0("0.", paste(rep("0", digits), collapse = ""))
  }

  number_style <- openxlsx::createStyle(numFmt = number_format)

  breaks_text <- format(
    breaks,
    digits = 17L,
    trim = TRUE,
    scientific = FALSE,
    decimal.mark = "."
  )

  #### Add matrices and data.frames ####

  for(i in seq_along(matrices)) {

    x <- matrices[[i]]
    sheet <- sheet_names[i]
    matrix_input <- is.matrix(x)

    openxlsx::addWorksheet(wb, sheet)

    openxlsx::writeData(
      wb, sheet, x,
      rowNames = TRUE
    )

    if(nrow(x) == 0L || ncol(x) == 0L) {
      next
    }

    for(j in seq_len(ncol(x))) {

      #### Identify numeric cells ####

      column <- if(matrix_input) {
        x[, j]
      } else {
        x[[j]]
      }

      convert_numeric <- matrix_input && is.character(column)

      if(convert_numeric) {

        values <- suppressWarnings(as.numeric(column))

      } else if(is.numeric(column)) {

        values <- column

      } else {

        next

      }

      numeric_rows <- which(is.finite(values))

      if(length(numeric_rows) == 0L) {
        next
      }

      rows <- numeric_rows+1L
      cols <- j+1L

      #### Write numeric matrix entries as Excel numbers ####

      if(convert_numeric) {

        runs <- split(
          numeric_rows,
          cumsum(c(TRUE, diff(numeric_rows) != 1L))
        )

        for(run in runs) {

          openxlsx::writeData(
            wb, sheet,
            unname(values[run]),
            startRow = run[1L]+1L,
            startCol = cols,
            colNames = FALSE,
            rowNames = FALSE
          )

        }

      }

      #### Number format ####

      openxlsx::addStyle(
        wb, sheet,
        style = number_style,
        rows = rows,
        cols = cols,
        gridExpand = TRUE,
        stack = TRUE
      )

      #### Conditional formatting ####

      cell <- paste0(
        openxlsx::int2col(cols),
        min(rows)
      )

      rules <- character(length(breaks)+1L)

      # First interval: x <= first break
      rules[1L] <- paste0(
        cell, "<=", breaks_text[1L]
      )

      # Intermediate intervals:
      # previous break < x <= current break
      if(length(breaks) > 1L) {

        for(k in 2:length(breaks)) {

          rules[k] <- paste0(
            cell, ">", breaks_text[k-1L],
            ",", cell, "<=", breaks_text[k]
          )

        }

      }

      # Last interval: x > last break
      rules[length(rules)] <- paste0(
        cell, ">",
        breaks_text[length(breaks)]
      )

      rules <- paste0(
        "IFERROR(AND(ISNUMBER(",
        cell, "),",
        rules,
        "),FALSE)"
      )

      for(k in seq_along(styles)) {

        if(!is.null(styles[[k]])) {

          openxlsx::conditionalFormatting(
            wb, sheet,
            cols = cols,
            rows = rows,
            type = "expression",
            rule = rules[k],
            style = styles[[k]]
          )

        }

      }

    }

  }

  #### Save workbook ####

  openxlsx::saveWorkbook(
    wb,
    file = file,
    overwrite = TRUE
  )

  #### Result ####

  return(invisible(file))

}
