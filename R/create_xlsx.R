create_xlsx <- function(..., file, sheet_names = NULL, digits = 3L,
                        breaks = c(-0.50, -0.20, 0.20, 0.50),
                        colors = c("#FA8072", "#EA9999", NA,
                                   "#90EE90", "#6AA84F")) {
  
  #### Check arguments ####
  
  matrices <- list(...)
  expressions <- as.list(substitute(list(...)))[-1L]
  
  if(length(matrices) == 0L) {
    stop("At least one matrix must be provided")
  }
  
  if(!all(vapply(matrices, is.matrix, logical(1L)))) {
    stop("All objects supplied through ... must be matrices")
  }
  
  if(!all(vapply(matrices, is.numeric, logical(1L)))) {
    stop("All matrices supplied through ... must be numeric")
  }
  
  if(length(digits) != 1L || !is.numeric(digits) || is.na(digits) ||
     digits < 0L || digits != as.integer(digits)) {
    stop("digits must be a non-negative integer")
  }
  
  digits <- as.integer(digits)
  
  if(length(breaks) != 4L || !is.numeric(breaks) ||
     any(!is.finite(breaks)) ||
     is.unsorted(breaks, strictly = TRUE)) {
    stop("breaks must contain four increasing numeric values")
  }
  
  if(length(colors) != 5L) {
    stop("colors must contain five color values")
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
      stop("sheet_names must contain one name for each matrix")
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
  
  #### Add matrices ####
  
  for(i in seq_along(matrices)) {
    
    x <- matrices[[i]]
    sheet <- sheet_names[i]
    
    openxlsx::addWorksheet(wb, sheet)
    
    openxlsx::writeData(
      wb, sheet, x,
      rowNames = TRUE
    )
    
    if(nrow(x) == 0L || ncol(x) == 0L) {
      next
    }
    
    rows <- 2:(nrow(x)+1L)
    cols <- 2:(ncol(x)+1L)
    
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
    
    if(!is.null(styles[[1]])) {
      openxlsx::conditionalFormatting(
        wb, sheet,
        cols = cols,
        rows = rows,
        rule = paste0("<=", breaks[1]),
        style = styles[[1]]
      )
    }
    
    if(!is.null(styles[[2]])) {
      openxlsx::conditionalFormatting(
        wb, sheet,
        cols = cols,
        rows = rows,
        rule = paste0(
          "AND(B2>", breaks[1],
          ",B2<", breaks[2], ")"
        ),
        style = styles[[2]]
      )
    }
    
    if(!is.null(styles[[3]])) {
      openxlsx::conditionalFormatting(
        wb, sheet,
        cols = cols,
        rows = rows,
        rule = paste0(
          "AND(B2>=", breaks[2],
          ",B2<=", breaks[3], ")"
        ),
        style = styles[[3]]
      )
    }
    
    if(!is.null(styles[[4]])) {
      openxlsx::conditionalFormatting(
        wb, sheet,
        cols = cols,
        rows = rows,
        rule = paste0(
          "AND(B2>", breaks[3],
          ",B2<", breaks[4], ")"
        ),
        style = styles[[4]]
      )
    }
    
    if(!is.null(styles[[5]])) {
      openxlsx::conditionalFormatting(
        wb, sheet,
        cols = cols,
        rows = rows,
        rule = paste0(">=", breaks[4]),
        style = styles[[5]]
      )
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