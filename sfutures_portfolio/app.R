# Packages-----
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
library(shiny)
library(shinydashboard)
library(ggplot2)
library(sensiPhy)
library(dplyr)
library(magrittr)
library(stringr)
library(glue)
library(fundiversity)
library(picante)
library(ggalt)     # for geom_encircle
library(ape)       # read.tree / drop.tip, base phylogeny plotting

# Functions needed --------------------------------------------------------
remove_species_by_traits <- function(tree, traits, trait1, weight1, trait2 = NULL, weight2 = NULL,
                                     percentage_remove, uncertainty = 0) {
  if (!trait1 %in% names(traits)) stop("The specified trait1 is not in the traits table.")
  if (!is.null(trait2) && !trait2 %in% names(traits)) stop("The specified trait2 is not in the traits table.")
  
  tree_species <- tree$tip.label
  matched_traits <- data.frame(species = rownames(traits), traits)
  matched_traits <- matched_traits[match(tree_species, matched_traits$species), ]
  
  if (any(is.na(matched_traits))) stop("Not all species in the tree have corresponding traits.")
  
  normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  matched_traits$trait1_normalized <- normalize(matched_traits[[trait1]]) * weight1
  if (!is.null(trait2) && !is.null(weight2)) {
    matched_traits$trait2_normalized <- normalize(matched_traits[[trait2]]) * weight2
    matched_traits$score <- matched_traits$trait1_normalized + matched_traits$trait2_normalized
  } else {
    matched_traits$score <- matched_traits$trait1_normalized
  }
  
  matched_traits$score <- matched_traits$score + rnorm(nrow(matched_traits), mean = 0, sd = uncertainty)
  
  num_remove <- ceiling(nrow(matched_traits) * percentage_remove / 100)
  species_to_remove <- head(matched_traits[order(-matched_traits$score), "species"], num_remove)
  
  new_tree <- drop.tip(tree, species_to_remove)
  return(list("new_tree" = new_tree, "removed_species" = species_to_remove))
}

plot_phylogeny_grouped <- function(tree,
                                   degraded_species,
                                   added_species,
                                   col_degraded = "red",
                                   col_added = "green",
                                   col_background = gray(0.85),
                                   layout = "circular",
                                   size_tip = 1) {
  edge_colors <- rep("grey50", nrow(tree$edge))
  
  get_path_edges <- function(tree, tip) {
    tip_idx <- which(tree$tip.label == tip)
    node <- tip_idx
    edges <- integer()
    while (node != (Ntip(tree) + 1)) {
      parent_edge <- which(tree$edge[, 2] == node)
      if (length(parent_edge) == 0) break
      edges <- c(edges, parent_edge)
      node <- tree$edge[parent_edge, 1]
    }
    edges
  }
  
  for (sp in degraded_species) {
    edge_colors[get_path_edges(tree, sp)] <- col_degraded
  }
  for (sp in added_species) {
    edge_colors[get_path_edges(tree, sp)] <- col_added
  }
  
  type <- if (layout == "circular") "fan" else "phylogram"
  
  plot.phylo(tree, type = type, show.tip.label = FALSE, edge.color = edge_colors)
  
  tip_colors <- rep(col_background, length(tree$tip.label))
  tip_colors[tree$tip.label %in% degraded_species] <- col_degraded
  tip_colors[tree$tip.label %in% added_species]    <- col_added
  
  tiplabels(pch = 19, col = tip_colors, cex = size_tip)
}

# Load data ---------------------------------------------------------------
trees   <- read.tree("phylogenetic_tree_allsp2020_S3.tre")
simucom <- read.csv("simulated_communities_N100000.csv")
traits  <- read.csv("functional_traits_species_imputation.csv") %>%
  filter(alien == 0) %>%
  mutate(
    wood_density.y = wd_g_cm3.y,
    max_height.y = max_height_m.y,
    sla.y = sla_mm_mg.y
  )
rownames(traits) <- traits$tip_name

# Process data ------------------------------------------------------------
colmaxFD   = "#a6d854"
colmaxPD   = "#8da0cb"
colmaxDIV  = '#e78ac3'

Nsimu <- nrow(simucom)
dd <- simucom[sample(seq_len(nrow(simucom)), Nsimu), ]
topfd <- quantile(simucom$zfd, probs = 0.95)
toppd <- quantile(simucom$zpd, probs = 0.95)
max_pd <- max(simucom$pd)
max_fd <- max(simucom$fd)

comp <- match_dataphy(formula = sla.y ~ wood_density.y + max_height.y, data = traits, phy = trees)

totsp <- 51
deg_prob = 0.15

removed_species <- remove_species_by_traits(tree = comp$phy,
                                            traits = comp$data %>% dplyr::select(wood_density.y, sla.y, max_height.y),
                                            trait1 = "max_height.y",
                                            trait2 = "wood_density.y",
                                            weight1 = 0.5,
                                            weight2 = 0.5,
                                            percentage_remove = 100 * (1 - deg_prob))
lost_sp <- removed_species$removed_species
all_sp <- comp$phy$tip.label
deg_sp <- base::setdiff(all_sp, lost_sp)

# Plot general pattern to avoid compouting again
# p1------
# pp1 <- ggplot() +
#   geom_point(data = dd, aes(x = zpd, y = zfd), alpha = .5, color = "gray") 

# Process about page from README.
about_content <-  includeMarkdown("README.md")

  # UI ----------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Community Restoration Dashboard"),
  dashboardSidebar(
    width = 320,
    div(class = "sidebar-info",
        p(class = "small",
          "This interactive app simulates community restoration scenarios. ",
          "Use the selector to explore different strategies and click ",
          em("Select another community"),
          " to sample a new simulated community.")
    ),
    
    sidebarMenu(
      menuItem("Simulator", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # controls block below About
    div(class = "sidebar-controls",
        tags$h4(style = "margin-top:0;", "Controls"),
        selectInput("type", "Select restoration approach:",
                    choices = c("Random" = "rand",
                                "Maximize PD" = "pd",
                                "Maximize FD" = "fd",
                                "Portfolio (high PD & FD)" = "both",
                                "Trade-off [high FD | low PD]" = "tpd",
                                "Trade-off [high PD | low FD]" = "tfd"),
                    selected = "rand"),
        actionButton("generate_btn", "Select another community", icon = icon("refresh"))
        # note: no width="100%" here
    )
  ) ,
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .sidebar-info { padding: 10px 15px 0 15px; }
        .sidebar-info p { white-space: normal; word-wrap: break-word; line-height: 1.3; }
        .content-wrapper, .right-side { overflow-x: hidden; }
        .box .caption { display:block; font-size: 12px; color:#555; margin-top: 6px; }
        @media (max-width: 1400px) { .bottom-charts .col-sm-6 { width: 100%; } }
        .about-box .box-body { font-size: 15px; line-height: 1.5; }
      "))
    ),
    tabItems(
      # --- Dashboard tab ---
      tabItem(
        tabName = "dashboard",
        fluidRow(
          valueBoxOutput("vb_fd_max", width = 3),
          valueBoxOutput("vb_fd_ref", width = 3),
          valueBoxOutput("vb_pd_max", width = 3),
          valueBoxOutput("vb_pd_ref", width = 3)
        ),
        fluidRow(
          column(
            width = 6,
            shinydashboard::box(
              title = "Portfolio of communities",
              width = 12,
              plotOutput("plot1", height = "520px"),
              tags$small(class = "caption",
                         "Green point = simulated restored community. Heatmap shows density of simulations in standardized PD–FD space. ",
                         "Dashed lines mark the 95th percentile thresholds; arrow indicates the distance to the selected community from 95th percentile of FD and PD."
              )
            )
          ),
          column(
            width = 6,
            shinydashboard::box(
              title = "Phylogenetic Tree & Trait Space",
              width = 12,
              splitLayout(
                cellWidths = c("50%", "50%"),
                plotOutput("plot2", height = "520px"),
                plotOutput("plot3", height = "520px")
              ),
              tags$small(class = "caption",
                         strong("Left:"),
                         " Circular phylogeny with ",
                         span(style = "color:tomato;font-weight:600;", "red"),
                         " = species present in the degraded ecosystem and ",
                         span(style = "color:green;font-weight:600;", "green"),
                         " = species added to the restored community.",
                         strong("Right:"),
                         " Trait space (wood density × max height); green points/hull = restored community, red = degraded ecosystem. "
              )
            )
          )
        ),
        fluidRow(
          class = "bottom-charts",
          column(
            width = 6,
            fluidRow(
              column(
                6,
                shinydashboard::box(
                  title = "Phylogenetic Diversity",
                  width = 12,
                  plotOutput("plot4", height = "340px"),
                  tags$small(class = "caption",
                             "Bars show PD as proportion of reference (100%). Dashed lines at 25%, 50%, 75%."
                  )
                )
              ),
              column(
                6,
                shinydashboard::box(
                  title = "Functional Diversity",
                  width = 12,
                  plotOutput("plot5", height = "340px"),
                  tags$small(class = "caption",
                             "Bars show FD as proportion of reference (100%). Dashed lines at 25%, 50%, 75%."
                  )
                )
              )
            )
          # ),
          # column(width = 6,
          #        box(
          #          title = "Controls",
          #          width = 12, status = "primary", solidHeader = TRUE,
          #          selectInput("type", "Select restoration approach:",
          #                      choices = c("Random" = "rand",
          #                                  "Maximize PD" = "pd", 
          #                                  "Maximize FD" = "fd", 
          #                                  "Portfolio (high PD & FD)" = "both",
          #                                  "Trade-off [high FD | low PD]" = "tpd",
          #                                  "Trade-off [high PD | low FD]" = "tfd"),
          #                      selected = "rand"),
          #          actionButton("generate_btn", "Select another community", icon = icon("refresh"))
          #        )
          )
        )
      ),
      # --- About tab ---
      tabItem(
        tabName = "about",
        fluidPage(
          fluidRow(
            column(
              width = 10,
              shinydashboard::box(
                title = "About",
                width = 12, status = "primary", solidHeader = TRUE, class = "about-box",
                about_content
              )
            )
          )
        )
      )
    )
  )
)

# Server ------------------------------------------------------------------
server <- function(input, output, session) {
  plots <- reactiveValues()
  metrics <- reactiveValues(fd_max = NA, fd_ref = NA, pd_max = NA, pd_ref = NA)
  
  generatePlots <- function() {
    crop_pd   <- simucom %>% dplyr::filter(pd == max_pd)
    crop_fd   <- simucom %>% dplyr::filter(fd == max_fd)
    crop_both <- simucom %>% dplyr::filter(zfd > topfd & zpd > toppd)
    crop_tfd  <- simucom %>% dplyr::filter(zfd < -1 & zpd > toppd)
    crop_tpd  <- simucom %>% dplyr::filter(zfd > topfd & zpd < -1)
    
    df <- switch(input$type,
                 rand = simucom,
                 pd   = crop_pd,
                 fd   = crop_fd,
                 tpd  = crop_tpd,
                 tfd  = crop_tfd,
                 both = crop_both)
    
    restored_simu <- sample.int(nrow(df), 1)
    restored_data <- df[restored_simu, , drop = FALSE]
    restored_sp <- str_split(df[restored_simu, ]$composition, pattern = "-")[[1]]
    restored_list <- c(deg_sp, restored_sp)
    
    trait_community <- traits %>% dplyr::filter(tip_name %in% restored_list)
    
    deg_data <- comp$data %>% dplyr::filter(tip_name %in% deg_sp) %>% dplyr::select(max_height.y, wood_density.y)
    res_data <- comp$data %>% dplyr::filter(tip_name %in% restored_list) %>% dplyr::select(max_height.y, wood_density.y)
    
    fd_ref <- as.numeric(fd_fric(comp$data[, c("max_height.y", "wood_density.y")])[2])
    fd_deg <- as.numeric(fd_fric(deg_data)[2]) / fd_ref
    fd_res <- as.numeric(fd_fric(res_data)[2]) / fd_ref
    fd_ref <- fd_ref / fd_ref
    
    all_presence = (all_sp %in% all_sp) * 1
    res_presence = (all_sp %in% restored_list) * 1
    deg_presence = (all_sp %in% deg_sp) * 1
    
    species_presence <- t(data.frame(all_presence, res_presence, deg_presence))
    colnames(species_presence) <- all_sp
    
    pd_comm <- pd(samp = species_presence, tree = comp$phy)
    pd_ref <- pd_comm$PD[1]
    pd_deg <- pd_comm$PD[3] / pd_ref
    pd_res <- pd_comm$PD[2] / pd_ref
    pd_ref <- pd_ref / pd_ref
    
    datlevel <- data.frame(
      pd = c(pd_ref, pd_deg, pd_res),
      fd = c(fd_ref, fd_deg, fd_res),
      sr = c(51, 7, 15),
      type = c("reference", "degraded", "restored")
    )
    
    pd_max <- round((restored_data$fd / max_fd) * 100) # (kept as in your workflow if intended)
    fd_max <- round((restored_data$fd / max_fd) * 100)
    
    metrics$fd_max <- fd_max
    metrics$fd_ref <- round((datlevel$fd[3]) * 100)
    metrics$pd_max <- round((restored_data$pd / max_pd) * 100)
    metrics$pd_ref <- round((datlevel$pd[3]) * 100)
    
    # p1------
    p1 <- ggplot() +
      geom_bin_2d(data = dd, aes(x = zpd, y = zfd), alpha = .5) +
    #p1 <- pp1 +
      scale_fill_gradient(low = "white", high = "black") +
      geom_point(data = dd %>% dplyr::filter(zpd >= toppd & zfd >= topfd), aes(x = zpd, y = zfd),
                 alpha = .4, size = 1, color = colmaxDIV) +
      geom_point(data = restored_data, aes(x = zpd, y = zfd), color = "green3", size = 7) +
      geom_hline(yintercept = topfd, linetype = 2) +
      geom_vline(xintercept = toppd, linetype = 2) +
      annotate("segment",
               x = toppd, y = topfd,
               xend = restored_data$zpd, yend = restored_data$zfd,
               arrow = arrow(length = grid::unit(0.3, "cm")),
               colour = "darkgreen", size = 1.5) +
      scale_y_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4)) +
      scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4, 4)) +
      theme_classic(base_size = 14) +
      labs(y = "Standardized Functional diversity", 
           x = "Standardized Phylogenetic diversity",
           fill = "") +
      annotate("text", x = 3.5, y = 3.5, label = "Portfolio\nHigh PD&FD", hjust = 1, vjust = 1,
               size = 6, color = colmaxDIV, fontface = "bold") +
      annotate("text", x = 3, y = -3, label = "Top 5%\nPD", hjust = 1, vjust = 0,
               size = 6, color = colmaxPD, fontface = "bold") +
      annotate("text", x = -3, y = 3, label = "Top 5%\nFD", hjust = 0.5, vjust = 1,
               size = 6, color = colmaxFD, fontface = "bold") +
      annotate("text", x = -3.5, y = -3.5, label = "Low\nPD & FD", hjust = 0.5, vjust = 0,
               size = 6, color = "black", fontface = "bold") +
      ggtitle(glue("FD = {metrics$fd_ref}% vs. reference | {metrics$fd_max}% vs. max\n",
                   "PD = {metrics$pd_ref}% vs. reference | {metrics$pd_max}% vs. max")) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        plot.margin = margin(15, 20, 20, 20)
      )
    
    p2 <- function() {
      plot_phylogeny_grouped(tree = trees,
                             degraded_species = deg_sp,
                             added_species = restored_sp,
                             col_degraded = "tomato",
                             col_added = "green3",
                             col_background = gray(0.85),
                             layout = "circular",
                             size_tip = 1)
    }
    
    p3 <- ggplot() +
      geom_point(data = comp$data, aes(y = max_height.y, x = wood_density.y)) +
      geom_encircle(data = comp$data, aes(y = max_height.y, x = wood_density.y),
                    size = 2, expand = 0, s_shape = 1) +
      geom_encircle(data = trait_community, aes(y = max_height.y, x = wood_density.y),
                    alpha = 0.2, fill = "green1", expand = 0, s_shape = 1) +
      geom_encircle(data = trait_community %>% dplyr::filter(tip_name %in% deg_sp),
                    aes(y = max_height.y, x = wood_density.y),
                    alpha = 0.3, fill = "tomato", expand = 0, s_shape = 1) +
      geom_point(data = trait_community, aes(y = max_height.y, x = wood_density.y),
                 color = "green3", size = 5) +
      geom_point(data = trait_community %>% dplyr::filter(tip_name %in% deg_sp),
                 aes(y = max_height.y, x = wood_density.y),
                 color = "tomato", size = 5) +
      theme_classic() +
      labs(y = "Maximum height (m)", x = "Wood density (g/cm³)") +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        plot.margin = margin(15, 20, 20, 20)
      )
    
    p4 <- ggplot(datlevel, aes(y = pd, x = type, fill = type)) +
      geom_col(alpha = 0.8, show.legend = FALSE) +
      geom_text(aes(label = paste0(round(pd * 100), "%")), 
                vjust = -0.6, size = 5) +
      scale_fill_manual(values = c("reference" = "grey70", "degraded" = "tomato", "restored" = "green3")) +
      scale_x_discrete(limits = c("reference", "degraded", "restored")) +
      geom_hline(yintercept = c(.25, .5, .75), linetype = 2) +
      coord_cartesian(ylim = c(0, 1.05), clip = "off") +
      theme_classic() +
      labs(y = "PD (vs. reference)", x = "") +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.margin = margin(20, 20, 30, 20)
      )
    
    p5 <- ggplot(datlevel, aes(y = fd, x = type, fill = type)) +
      geom_col(alpha = 0.8, show.legend = FALSE) +
      geom_text(aes(label = paste0(round(fd * 100), "%")), 
                vjust = -0.6, size = 5) +
      scale_fill_manual(values = c("reference" = "grey70", "degraded" = "tomato", "restored" = "green3")) +
      scale_x_discrete(limits = c("reference", "degraded", "restored")) +
      geom_hline(yintercept = c(.25, .5, .75), linetype = 2) +
      coord_cartesian(ylim = c(0, 1.05), clip = "off") +
      theme_classic() +
      labs(y = "FD (vs. reference)", x = "") +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.margin = margin(20, 20, 30, 20)
      )
    

# Plot general relationship between PD, TD, FD ----------------------------

    
    list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
  }
  
  observeEvent(list(input$generate_btn, input$type), {
    newPlots <- generatePlots()
    plots$p1 <- newPlots$p1
    plots$p2 <- newPlots$p2
    plots$p3 <- newPlots$p3
    plots$p4 <- newPlots$p4
    plots$p5 <- newPlots$p5
  })
  
  observe({
    req(is.null(plots$p1))
    newPlots <- generatePlots()
    plots$p1 <- newPlots$p1
    plots$p2 <- newPlots$p2
    plots$p3 <- newPlots$p3
    plots$p4 <- newPlots$p4
    plots$p5 <- newPlots$p5
  })
  
  output$plot1 <- renderPlot({ plots$p1 })
  output$plot2 <- renderPlot({ plots$p2() })
  output$plot3 <- renderPlot({ plots$p3 })
  output$plot4 <- renderPlot({ plots$p4 })
  output$plot5 <- renderPlot({ plots$p5 })
  
  output$vb_fd_max <- renderValueBox({
    valueBox(paste0(metrics$fd_max, "%"), "FD vs. max potential", icon = icon("chart-line"), color = "light-blue")
  })
  output$vb_fd_ref <- renderValueBox({
    valueBox(paste0(metrics$fd_ref, "%"), "FD vs. reference", icon = icon("chart-bar"), color = "light-blue")
  })
  output$vb_pd_max <- renderValueBox({
    valueBox(paste0(metrics$pd_max, "%"), "PD vs. max potential", icon = icon("chart-line"), color = "purple")
  })
  output$vb_pd_ref <- renderValueBox({
    valueBox(paste0(metrics$pd_ref, "%"), "PD vs. reference", icon = icon("chart-bar"), color = "purple")
  })
}

shinyApp(ui = ui, server = server)
