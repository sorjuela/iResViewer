#' Interactive visualization of differential gene expression results
#'
#' @param geneInfo Data frame with gene annotation information. Should have at
#'   least two columns: gene (the gene ID) and symbol (the gene symbol).
#' @param abundances Named list with data frames in "long" format (one row per
#'   gene/sample combination), containing abundance estimates. These will be
#'   used to illustrate the abundance pattern for selected genes. Each data
#'   frame must have at least four columns: sample (the sample ID), gene (the
#'   gene ID), value (the abundance) and group (a sample group label, used to
#'   order and color the points in the plot).
#' @param geneInfo2 Data frame with gene annotation information. Should have at
#'   least two columns: gene (the gene ID) and symbol (the gene symbol).
#' @param abundances2 Named list with data frames in "long" format (one row per
#'   gene/sample combination), containing abundance estimates. These will be
#'   used to illustrate the abundance pattern for selected genes. Each data
#'   frame must have at least four columns: sample (the sample ID), gene (the
#'   gene ID), value (the abundance) and group (a sample group label, used to
#'   order and color the points in the plot).
#'
#' @param appTitle App title
#' @param ... Additional arguments (currently not used)
#'
#' @author Charlotte Soneson
#'
#' @import shiny ggplot2
#' @importFrom utils write.csv
#' @importFrom dplyr %>% filter arrange
#' @importFrom ggrepel geom_label_repel
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'   dashboardBody tabBox
#' @export
#'
iResViewer <- function(geneInfo = NULL,
                       abundances = list(),
                       geneInfo2 = NULL,
                       abundances2 = list(),
                       appTitle = "iResViewer", ...) {
  options(ucscChromosomeNames = FALSE)
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  pLayout <- function() {
    shinydashboard::dashboardPage(
      skin = "green",

      shinydashboard::dashboardHeader(title = appTitle,
                                      titleWidth = nchar(appTitle) * 20),

      shinydashboard::dashboardSidebar(
        shiny::textInput(inputId = "sel.gene", label = "Selected gene")
      ),

      shinydashboard::dashboardBody(shiny::fluidRow(
        do.call(shinydashboard::tabBox,
                c(
                  width = 12,

            
                  ## ======================================================== ##
                  ## Tabs with abundances
                  ## ======================================================== ##
                  lapply(names(abundances), function(w)
                    shiny::tabPanel(
                      w,
                      shiny::div(style = "position:relative",
                                 shiny::uiOutput(paste0(w, "_abundance_ui")),
                                 shiny::uiOutput(paste0(w, "_abundance_hover_info"))))
                  ),
                  
                  ## ======================================================== ##
                  ## Tabs with abundances (Exon)
                  ## ======================================================== ##
                  lapply(names(abundances2), function(w)
                    shiny::tabPanel(
                      w,
                      shiny::div(style = "position:relative",
                                 shiny::uiOutput(paste0(w, "_abundance_ui")),
                                 shiny::uiOutput(paste0(w, "_abundance_hover_info"))))
                )
        )
      ))
    ))
  }

  server_function <- function(input, output, session) {
    options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

    ## ====================================================================== ##
    ## Abundance plots
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_abundance_ui")]] <- shiny::renderUI({
        shiny::plotOutput(paste0(nm, "_abundance"), hover = paste0(nm, "_abundance_hover"),
                          width = "90%", height = "800px")
      })
    }, names(abundances))

    Map(function(nm) {
      output[[paste0(nm, "_abundance")]] <- shiny::renderPlot({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo$gene)),
                  tolower(geneInfo$symbol)))) return(NULL)
          id <- (geneInfo %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            df <- abundances[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            ggplot(df, aes(x = sample, y = value, group = gene, col = group)) +
              geom_line(col = "black") + geom_point(size = 3) +
              theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                 legend.position = "bottom",
                                 axis.text.y = element_text(size = 14),
                                 axis.title.y = element_text(size = 16)) +
              coord_cartesian(ylim = c(2,15), expand = FALSE) + 
              guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + xlab("") +
              ylab("Intensity")
          }
        }
      })
    }, names(abundances))

    Map(function(nm) {
      output[[paste0(nm, "_abundance_hover_info")]] <- shiny::renderUI({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo$gene)),
                  tolower(geneInfo$symbol)))) return(NULL)
          id <- (geneInfo %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            hover <- input[[paste0(nm, "_abundance_hover")]]
            df <- abundances[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            res2 <- shiny::nearPoints(df, hover, threshold = 5, maxpoints = 1)
            res2$genename <- geneInfo$symbol[match(res2$gene, geneInfo$gene)]
            left_pct <- (hover$x - hover$domain$left)/(hover$domain$right - hover$domain$left)
            top_pct <- (hover$domain$top - hover$y)/(hover$domain$top - hover$domain$bottom)
            left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
            top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
            style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                            "left:", left_px + 2, "px; top:", top_px + 2, "px;")
            shiny::wellPanel(
              style = style,
              shiny::p(shiny::HTML(paste0("<b> Sample: </b>", res2$sample, "<br/>",
                                          "<b> ID: </b>", res2$gene, "<br/>",
                                          "<b> Gene: </b>", res2$genename, "<br/>",
                                          "<b> abundance: </b>", round(res2$value, 3), "<br/>")))
            )
          }
        }
      })
    }, names(abundances))
    
    ## ====================================================================== ##
    ## Abundance plots (exons)
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_abundance_ui")]] <- shiny::renderUI({
        shiny::plotOutput(paste0(nm, "_abundance"), hover = paste0(nm, "_abundance_hover"),
                          width = "90%", height = "800px")
      })
    }, names(abundances2))
    
    Map(function(nm) {
      output[[paste0(nm, "_abundance")]] <- shiny::renderPlot({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo2$gene)),
                  tolower(geneInfo2$symbol)))) return(NULL)
          id <- (geneInfo2 %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            df <- abundances2[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            ggplot(df, aes(x = sample, y = value, group = gene, col = group)) +
              geom_line(col = "black") + geom_point(size = 3) +
              theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
                                 legend.position = "bottom",
                                 axis.text.y = element_text(size = 14),
                                 axis.title.y = element_text(size = 16)) +
              coord_cartesian(ylim = c(2,15), expand = FALSE) + 
              guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + xlab("") +
              ylab("Intensity")
          }
        }
      })
    }, names(abundances2))
    
    Map(function(nm) {
      output[[paste0(nm, "_abundance_hover_info")]] <- shiny::renderUI({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo2$gene)),
                  tolower(geneInfo2$symbol)))) return(NULL)
          id <- (geneInfo2 %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            hover <- input[[paste0(nm, "_abundance_hover")]]
            df <- abundances2[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            res2 <- shiny::nearPoints(df, hover, threshold = 5, maxpoints = 1)
            res2$genename <- geneInfo2$symbol[match(res2$gene, geneInfo2$gene)]
            left_pct <- (hover$x - hover$domain$left)/(hover$domain$right - hover$domain$left)
            top_pct <- (hover$domain$top - hover$y)/(hover$domain$top - hover$domain$bottom)
            left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
            top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
            style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                            "left:", left_px + 2, "px; top:", top_px + 2, "px;")
            shiny::wellPanel(
              style = style,
              shiny::p(shiny::HTML(paste0("<b> Sample: </b>", res2$sample, "<br/>",
                                          "<b> ID: </b>", res2$gene, "<br/>",
                                          "<b> Gene: </b>", res2$genename, "<br/>",
                                          "<b> abundance: </b>", round(res2$value, 3), "<br/>")))
            )
          }
        }
      })
    }, names(abundances2))

  }

  shinyApp(ui = pLayout, server = server_function)
}
