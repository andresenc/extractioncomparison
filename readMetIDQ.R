library(xlsx)
library(dplyr)

## read metabolite names and classes, raw data and LOD/LOQ/Valid flag ##
################################################################################
readMetIDQ <- function(file, sampleName="Sample Name", colSampleName=FALSE,
                       filterIndicators=TRUE, groupAnnotation=FALSE,
                       verbose=TRUE) {
  if (verbose) {
    print("------------------------------")
    print("Loading full sheet...")
  }
  ## read the full Excel sheet ##
  fullSheet <- read.xlsx(file, sheetIndex = 1)
  colnames(fullSheet) <- NULL
  fullSheet <- as.matrix(fullSheet)
  fullSheet[fullSheet == "NA"] <- NA
  
  ## remove only NA columns and rows at the end of the sheet
  fullNArows <- rowSums(is.na(fullSheet)) == ncol(fullSheet)
  for (i in length(fullNArows):1) {
    if (isFALSE(fullNArows[i])) {
      keepRows <- 1:i
      break
    }
  }
  fullNAcols <- colSums(is.na(fullSheet)) == nrow(fullSheet)
  for (i in length(fullNAcols):1) {
    if (isFALSE(fullNAcols[i])) {
      keepCols <- 1:i
      break
    }
  }
  fullSheet <- fullSheet[keepRows, keepCols]
  
  ## get range of data
  colClass <- ceiling(which(fullSheet == "Class") / nrow(fullSheet))
  colsData <- (colClass+1):ncol(fullSheet)
  colSampleType <- ceiling(which(fullSheet == "Sample Type") / nrow(fullSheet))
  rowsData <- which(fullSheet[, colSampleType] == "Sample")
  
  ## metabolite names and classes
  rowClass <- which(fullSheet == "Class") %% nrow(fullSheet)
  metabolites <- fullSheet[rowClass-1, colsData]
  classes <- fullSheet[rowClass, colsData]
  names(metabolites) <- classes
  if (verbose) {
    print(paste("Number of metabolites:", length(metabolites)))
    print(paste("Number of classes:", length(unique(names(metabolites)))))
  }
  
  ## Sample name
  if (isFALSE(colSampleName)) {
    colSampleName <- ceiling(which(fullSheet == sampleName) / nrow(fullSheet))
  }
  sampleNames <- fullSheet[rowsData, colSampleName]
  
  ## raw data
  dataRaw <- as.data.frame(fullSheet[rowsData, colsData])
  rownames(dataRaw) <- sampleNames
  colnames(dataRaw) <- metabolites
  dataRaw[] <- sapply(dataRaw[], as.numeric)
  
  ## get the background color
  if (verbose) {
    print("------------------------------")
    print("Reading background color...")
  }
  wb <- loadWorkbook(file)
  sheet1 <- getSheets(wb)[[1]]
  rows <- getRows(sheet1, rowIndex = rowsData+1) # first row was the header at the data read in
  cells <- getCells(rows, colIndex = colsData)
  styles <- sapply(cells, getCellStyle)
  bg <- sapply(styles, function(style) {
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) "")
    rgb <- paste(rgb, collapse = "")
    return(rgb)
  })
  rowIndex <- as.numeric(unlist(lapply(strsplit(names(bg),split = "\\."),
                                       "[", 1)))
  rowIndex <- rowIndex - min(rowIndex) + 1 # start at 1
  colIndex <- as.numeric(unlist(lapply(strsplit(names(bg), split = "\\."),
                                       "[", 2)))
  colIndex <- colIndex - min(colIndex) + 1 # start at 1
  matBG <- matrix("", ncol = max(colIndex), nrow = max(rowIndex))
  for (i in 1:length(bg)) {
    matBG[rowIndex[i], colIndex[i]] <- bg[i]
  }
  dataFlag <- as.data.frame(matBG)
  dataFlag <- dataFlag[rowSums(dataFlag == "") != ncol(dataFlag), ]
  rownames(dataFlag) <- sampleNames
  colnames(dataFlag) <- metabolites
  
  # dataFlag <- as.matrix(dataFlag)
  dataFlag[dataFlag == "6a5acd"] <- "LOD"
  dataFlag[dataFlag == "87ceeb"] <- "LOQ"
  dataFlag[dataFlag == "00cd66"] <- "Valid"
  dataFlag[is.na(dataRaw)] <- NA
  if (verbose) {
    sumLOD <- sum(dataFlag == "LOD", na.rm = TRUE)
    sumLOQ <- sum(dataFlag == "LOQ", na.rm = TRUE)
    sumValid <- sum(dataFlag == "Valid", na.rm = TRUE)
    total <- nrow(dataFlag)*ncol(dataFlag)
    print(paste0("Number of LODs: ", sumLOD,
                 " (", round(sumLOD/total*100, 2), "%)"))
    print(paste0("Number of LOQs: ", sumLOQ,
                 " (", round(sumLOQ/total*100, 2),"%)"))
    print(paste0("Number of Valids: ", sumValid,
                 " (", round(sumValid/total*100, 2),"%)"))
  }
  
  ## filter matabolism indicators
  if (isTRUE(filterIndicators) & "Metabolism Indicators" %in% names(metabolites)) {
    notIndicators <- which(names(metabolites) != "Metabolism Indicators")
    metabolites <- metabolites[notIndicators]
    if (verbose) {
      print("------------------------------")
      print("Filtering metabolism indicators...")
      print(paste("Number of remaining metabolites:", length(metabolites)))
      print(paste("Number of remaining classes:",
                  length(unique(names(metabolites)))))
    }
    
    dataRaw <- dataRaw[, metabolites]
    dataFlag <- dataFlag[, metabolites]
  }

  ## annotate with groups
  if (!isFALSE(groupAnnotation)) {
    if (verbose) {
      print("------------------------------")
      print("Annotating data...")
    }
    colName <- groupAnnotation[["Group Column"]]
    colGroup <- ceiling(which(fullSheet == colName) / nrow(fullSheet))
    groupEntries <- fullSheet[rowsData, colGroup]
    groupNameOpt <- names(groupAnnotation)[2:length(groupAnnotation)]
    for (groupName in rev(groupNameOpt)) {
      options <- groupAnnotation[[groupName]]
      if (grepl(options[1], paste(groupEntries, collapse = ""))) {
        pattern <- options[1]
        other <- options[2]
      } else {
        pattern <- options[2]
        other <- options[1]
      }
      groups <- sapply(groupEntries, function(GroupEntry) {
        if_else(grepl(pattern, GroupEntry, ignore.case = TRUE),
                pattern,
                other)
      })
      dataRaw <- mutate(dataRaw,
                        tmpName = groups,
                        .before = 1)
      colnames(dataRaw)[1] <- groupName
      dataFlag <- mutate(dataFlag,
                         tmpName = groups,
                         .before = 1)
      colnames(dataFlag)[1] <- groupName
    }
    if (verbose) {
      print(paste("Annotated colums:",
                  paste(names(groupAnnotation)[2:length(groupAnnotation)],
                        collapse = ", ")))
    }
    dataRaw <- arrange(dataRaw, !!sym(groupNameOpt[1]), !!sym(groupNameOpt[2]))
    dataFlag <- arrange(dataFlag, !!sym(groupNameOpt[1]), !!sym(groupNameOpt[2]))
  }
  
  outList <- list(dataRaw, dataFlag, metabolites)
  return(outList)
}


## rename rownames to tissue name
################################################################################
renameRownames <- function(inputList, patterns, newNames, startwith = TRUE) {
  dataRaw <- inputList[[1]]
  dataFlag <- inputList[[2]]
  metabolites <- inputList[[3]]
  
  for (i in seq_along(patterns)) {
    if (startwith) {
      pat <- paste0("^", patterns[i])
    } else {
      pat <- patterns[i]
    }
    newRownames <- str_replace_all(rownames(dataRaw), pat, newNames[i])
    rownames(dataRaw) <- newRownames
    rownames(dataFlag) <- newRownames
  }
  outList <- list(dataRaw, dataFlag, metabolites)  
}


## Shiny input ##
################################################################################
toInputFormat <- function(dataList) {
  dataRaw <- dataList[[1]]
  dataFlag <- dataList[[2]]
  metabolites <- dataList[[3]]
  tissue <- str_split(rownames(dataRaw)[1], "_")[[1]][2]
  methods <- str_remove_all(rownames(dataRaw), "\\.[1-3].+")
  shinyInput <- lapply(unique(methods), function(m) {
    tmp_df <- dataRaw[methods == m, ]
    means <- apply(tmp_df, 2, function(x) mean(x, na.rm=TRUE))
    sds <- apply(tmp_df, 2, function(x) sd(x, na.rm=TRUE))
    cvs <- sds / means
    tmp_flag <- dataFlag[methods == m, ]
    lods <- apply(tmp_flag, 2, function(x) if_else(sum(x != "LOD", na.rm = TRUE) >= 2, 1, 0))
    # lods <- apply(tmp_flag, 2, function(x) mean(x != "LOD", na.rm = TRUE))
    dfOut <- data.frame(Mean = means,
                        CV = cvs,
                        Method = factor(names(method_vector[method_vector == m]),
                                        levels = names(method_vector)), 
                        Metabolite = factor(metabolites, levels = metabolites),
                        Class = factor(names(metabolites), levels = unique(names(metabolites))),
                        SizeA = sum(tmp_df[1, ], na.rm = TRUE),
                        SizeB = sum(tmp_df[2, ], na.rm = TRUE),
                        SizeC = ifelse(nrow(tmp_df) == 3, sum(tmp_df[3, ], na.rm = TRUE), NA),
                        LOD = lods,
                        Tissue = factor(tissue, levels = tissues),
                        Methods = factor(m, levels = method_vector))
  }) %>%
    bind_rows()
  
  ## Imputation and transformation
  validMetabos <- shinyInput %>%
    group_by(Metabolite) %>%
    summarise(Valid = sum(LOD) > 0) %>%
    data.frame() %>%
    select(Valid) %>%
    unlist()
  dataImputed <- apply(dataRaw[, validMetabos], 2, function(vec) {
    vec[vec == 0] <- NA
    vec[is.na(vec)] <- min(vec, na.rm = TRUE) / 5
    return(vec)
  })
  dataTransformed <- log2(dataImputed)
  
  ## Anova; Groups
  dfGroups <- lapply(metabolites, function(metabolite) {
    if (metabolite %in% colnames(dataTransformed)) {
      anova <- aov(Conc ~ Method, data = data.frame(Method = methods,
                                                    Conc = dataTransformed[, metabolite]))
      postHoc <- HSD.test(anova, "Method", group = TRUE)
      groups <- toupper(postHoc$groups$groups)
      dfOut <- data.frame(Metabolite = factor(metabolite, levels = metabolites),
                          Class = names(metabolites[which(metabolites == metabolite)]),
                          Method = factor(str_replace(rownames(postHoc$groups), "\\.", " "),
                                          levels = method_vector),
                          Group = groups)
    } else {
      dfOut <- data.frame(Metabolite = factor(metabolite, levels = metabolites),
                          Class = names(metabolites[which(metabolites == metabolite)]),
                          Method = factor(method_vector,
                                          levels = method_vector),
                          Group = NA)
    }
    return(dfOut)
  }) %>%
    bind_rows() %>%
    arrange(Method, Class, Metabolite)
  
  shinyInput$Group <- dfGroups$Group
  
  ## add filtered metabolites as a column with NA
  dataTransformedOut <- sapply(metabolites, function(metabolite) {
    if (metabolite %in% colnames(dataTransformed)) {
      dataTransformed[, metabolite]
    } else {
      rep(NA, nrow(dataTransformed))
    }
  })
  rownames(dataTransformedOut) <- rownames(dataRaw)
  colnames(dataTransformedOut) <- colnames(dataRaw)
  dataTransformedOut <- as.data.frame(dataTransformedOut)
  
  return(list(shinyInput, dataTransformedOut))
}