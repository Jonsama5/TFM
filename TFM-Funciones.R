# =========================================================
# FUNCIONES PARA EL ANÁLISIS DE LA MICROBIOTA CON PHYLOSEQ
# =========================================================

#############                               #############
###        FUNCIONES DE RESUMEN DE PARÁMETROS         ### 
#############                               #############
phy_statistics <- function(x) {
  # Calcular las estadísticas globales
  
  # Total de lecturas para todas las muestras
  total_reads <- sum(sample_sums(x))
  
  # Número total de OTUs en el dataset
  n_otu <- nrow(otu_table(x))
  
  # Media y mediana de lecturas (en este caso por todas las muestras)
  mean_reads <- mean(sample_sums(x))
  median_reads <- median(sample_sums(x))
  
  # Desviación estándar de lecturas (para todas las muestras)
  devs <- sd(sample_sums(x))
  
  # Mínimo y máximo de lecturas entre todas las muestras
  min_reads <- min(sample_sums(x))
  max_reads <- max(sample_sums(x))
  
  # Calcular los índices de diversidad globales (Shannon y Simpson)
  diversity <- estimate_richness(x, measures = c("Shannon", "Simpson"))
  
  # Calcular la distancia beta
  beta_b <- distance(x, method = c("bray"), weighted = TRUE)
  beta_j <- distance(x, method = c("jaccard"), weighted = TRUE)
  
  # Crear un tibble con los resultados globales
  stats_tibble <- tibble(
    TotalReads = total_reads,
    NumOTUs = n_otu,
    MeanReads = mean_reads,
    MedianReads = median_reads,
    SDReads = devs,
    MinReads = min_reads,
    MaxReads = max_reads,
    ShannonDiversity = mean(diversity$Shannon),  # Promedio global de Shannon
    SimpsonDiversity = mean(diversity$Simpson),   # Promedio global de Simpson
    BrayDiversity = mean(beta_b), # Promedio global de Bray-Curtis
    JaccardDiversity = mean(beta_j) # Promedio global de Jaccard
  )
  
  # Devolver el tibble con las estadísticas globales
  return(stats_tibble)
}

#############                               #############
### FUNCIONES DE CREACIÓN Y VISUALIZACIÓN DE GRÁFICOS ### 
#############                               #############

#Creador de gráficos para cada rango taxonómico
plot_bar_by_ranks <- function(physeq, paleta = NULL) {
  library(phyloseq)
  library(ggplot2)
  library(scales)
  
  #Lista donde se guardaron los gráficos generados
  plot_bar_list <- list()
  
  #Iterar sobre cada rango taxonómico del objeto physeq
  for (rank in rank_names(physeq)) {
    taxa <- get_taxa_unique(physeq, rank)
    #Si se define una paleta de colores y es válida (puede colorear cada taxón),
    #entonces se utiliza. Si no se define (NULL) o tiene insuficiente cantidad
    #de colores, se genera automáticamente (much to my dismay) una válida.
    if (is.null(paleta) || length(paleta) < length(taxa)) { #Si "paleta" es NuLL o menor que el número de taxones
      autopaleta <- hue_pal()(length(taxa)) #Hue_pal genera colores (por hue, chroma, luminance) con la longitud del número de taxones
      names(autopaleta) <- taxa #El nombre de los colores (código) se asocia a cada taxón
      pal <- autopaleta
    } else {
      pal <- paleta
    }
    
    #Creación del gráfico
    plot_bar_list[[rank]] <- plot_bar(physeq, fill = rank) +
      labs(title = paste("Composición por", rank)) +
      scale_fill_manual(values = pal)
  }
  
  return(plot_bar_list)
}
paleta <-  c("#FFFFFF", "#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000",
             "#000080", "#800080", "#008080", "#C0C0C0", "#808080", "#FFA500", "#FFFF99", "#FFD700", "#FF69B4",
             "#87CEEB", "#FF6347", "#FF7F50", "#20B2AA", "#ADFF2F", "#7FFFD4", "#40E0D0", "#FF4500", "#DA70D6", "#FFDAB9",
             "#DB7093", "#FFDEAD", "#00FA9A", "#B0E0E6", "#FFC0CB", "#E6E6FA", "#00FF7F", "#4682B4", "#F08080", "#DAA520",
             "#FFA07A", "#228B22", "#FF8C00", "#FF69B4", "#BA55D3", "#7FFF00", "#6A5ACD", "#48D1CC", "#CD5C5C", "#F0E68C",
             "#00CED1", "#20B2AA", "#9932CC", "#8B0000", "#4B0082", "#7CFC00", "#D2691E", "#BC8F8F", "#8A2BE2", "#FA8072",
             "#2E8B57", "#FF1493", "#1E90FF", "#FFD700", "#B22222", "#FF6347", "#40E0D0", "#808000", "#FF8C00", "#FF69B4",
             "#DB7093", "#20B2AA", "#00FF7F", "#7FFF00", "#4682B4", "#F08080", "#DAA520", "#FFA07A", "#9932CC", "#8B0000",
             "#7CFC00", "#D2691E", "#BC8F8F", "#8A2BE2", "#FA8072", "#2E8B57", "#FF1493", "#1E90FF", "#B22222", "#FF6347",
             "#000000", "#FFFFFF", "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400", "#8B008B", "#556B2F", "#FF8C00",
             "#9932CC", "#8B0000", "#E9967A", "#8FBC8F", "#483D8B", "#2F4F4F", "#00CED1", "#9400D3", "#FF1493", "#00BFFF",
             "#696969", "#1E90FF", "#B22222", "#FFFAF0", "#778795", "#FF00FF", "#DCDCDC", "#F8F8FF", "#FFD700", "#DAA520",
             "#ADFF2F", "#FF69B4", "#CD5C5C", "#8B4513", "#F0E68C", "#20B2AA", "#6A5ACD","#9370DB", "#00FF7F", "#9ACD32", 
             "#BDB76B",     "#778899", "#FF8C00", "#BA55D3", "#4169E1", "#F08080", "#20B2AA", "#FF4500", "#ADFF2F", "#8A2BE2", 
             "#87CEEB", "#800000")

paleta_genus <- c("#000000", "#A52A2A", "#FF69B4", "#6A5ACD", "#1E90FF", "#7B68EE", "#808080", "#FFD700", "#556B2F", 
                  "#00FA9A", "#FF8C00", "#DA70D6", "#F8F8FF", "#000080", "#DC143C", "#9932CC", "#2F4F4F", "#FF1493", 
                  "#B0C4DE", "#FFDAB9", "#ADFF2F")


#Visualizador general de gráficos
plotwatcher <- function(plotlist) {
  cat("Welcome to plotwatcher, write a number in the terminal for the following:") 
  option <- readline(prompt = "Select an option: \n(1) Visualize plot by plot\n(2) Visualize plots together\n(3) Exit\n")
  i <- 1 
  
  while (option != "3") {
    if (option == "1") {
      # Bucle para visualizar los gráficos uno por uno
      repeat {
        cat("Visualizing plot by plot\n")
        print(plotlist[[i]])  # Mostrar el gráfico actual
        plotindex <- readline(prompt = "Choose an option: \n(1) Next plot\n(2) Previous plot\n(3) Exit\n")
        
        if (plotindex == "1" && i < length(plotlist)) {
          i <- i + 1  # Avanzar al siguiente gráfico
        } else if (plotindex == "2" && i > 1) {
          i <- i - 1  # Retroceder al gráfico anterior
        } else if (plotindex == "3") {
          cat("Exiting plot by plot visualizer\n") 
          break  # Salir si elige "3"
        } else {
          cat("Invalid option or end of plot list\n")
          break #Vuelve al menú principal
        }
      }
    } else if (option == "2") {
      cat("Visualizing all plots together\n")
      library(patchwork)
      combined <- wrap_plots(lapply(plotlist, function(p) p + theme(legend.position = "none")))
      print(combined)
      
      break  # Salir después de mostrar todos los gráficos, no vuelvo al menú porque va a tardar y no quiero sobrecargar
    } else {
      cat("Invalid option. Please try again\n")
    }
    # Pedir de nuevo una opción si no es la opción 3
    option <- readline(prompt = "Select an option: \n(1) Visualize plot by plot\n(2) Visualize plots together\n(3) Exit\n")
  }
  cat("Exiting the function.\n")
}

# Función para seleccionar los N taxones más abundantes
get_top_taxa <- function(physeq_object, top_n) {
  taxa_sums_all <- taxa_sums(physeq_object)
  top_taxa <- names(sort(taxa_sums_all, decreasing = TRUE)[1:top_n])
  pruned_object <- prune_taxa(top_taxa, physeq_object)
  return(pruned_object)
}

#############                               #############
##############   FUNCIONES DE DIVERSIDAD   ############## 
#############                               #############

#Creador de gráficos de riqueza alfa para cada rango taxonómico
plot_richness_by_variables <- function(physeq) {
  library(phyloseq)
  library(ggplot2)
  
  #Lista donde se guardaron los gráficos generados
  plot_bar_list <- list()
  
  #Iterar sobre cada rango taxonómico del objeto physeq
  for (variable in colnames(sample_data(physeq))) {
    plot_bar_list[[variable]] <- plot_richness(physeq, x = variable, measures = c("Shannon", "Simpson"), title = paste("Diversidad por", variable)) + 
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      geom_boxplot(alpha = 0.6) +
      theme_bw() +
      theme(panel.grid = element_line(color = "grey", size = 0.5), 
            plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
            axis.text.x = element_text(angle = 90))
    ylab(label = "Índice alfa diversidad")
  }
  
  return(plot_bar_list)
}

#Test estadístico Kruskal Wallis
kw_by_variables <- function(physeq) {
  df_alpha <- estimate_richness(physeq, measures = c("Shannon", "Simpson")) %>%
    cbind(sample_data(physeq))
  
  results <- list()
  
  for (variable in sample_variables(physeq)) {
    # Saltar si la variable tiene solo un valor único
    if (length(unique(df_alpha[[variable]])) < 2) next
    
    formula_shannon <- as.formula(paste("Shannon ~", variable))
    formula_simpson <- as.formula(paste("Simpson ~", variable))
    
    test_shannon <- kruskal.test(formula_shannon, data = df_alpha)
    test_simpson <- kruskal.test(formula_simpson, data = df_alpha)
    
    results[[variable]] <- tibble(
      Variable = variable,
      Index = c("Shannon", "Simpson"),
      Statistic = c(test_shannon$statistic, test_simpson$statistic),
      Df = c(test_shannon$parameter, test_simpson$parameter),
      P_value = c(test_shannon$p.value, test_simpson$p.value)
    )
  }
  
  bind_rows(results)
} 

#Test estadístico dunn posterior a Kruskall-Wallis
dunn_by_variables <- function(physeq) {
  df_alpha <- estimate_richness(physeq, measures = c("Shannon", "Simpson")) %>%
    cbind(sample_data(physeq))
  
  results <- list()
  
  for (variable in sample_variables(physeq)) {
    # Saltar variables con < 3 niveles únicos
    if (length(unique(df_alpha[[variable]])) < 3) next
    
    for (index in c("Shannon", "Simpson")) {
      formula <- as.formula(paste(index, "~", variable))
      
      dunn_res <- dunnTest(formula, data = df_alpha, method = "bh")
      
      dunn_df <- dunn_res$res %>%
        mutate(
          Variable = variable,
          Index = index
        ) %>%
        rename(
          Comparison = Comparison,
          Z = Z,
          P_uncorrected = P.unadj,
          P_adjusted = P.adj
        ) %>%
        select(Variable, Index, Comparison, Z, P_uncorrected, P_adjusted)
      
      results[[paste(variable, index, sep = "_")]] <- dunn_df
    }
  }
  
  bind_rows(results)
} 

#Creador de gráficos de riqueza beta (Bray-Curtis y Jaccard) para cada rango taxonómico
plot_bray_richness_by_variables <- function(physeq){
  plot_bar_list <- list()
  dist_bray <- distance(physeq, method = "bray")
  ordination <- ordinate(physeq, method = "PCoA", distance = dist_bray)
  for (variable in colnames(sample_data(physeq))) {
    if (is.numeric(sample_data(physeq)[[variable]])) {
      scale <- scale_color_gradient(low = "blue", high = "red")
    } else {
      scale <- scale_color_manual(values = paleta)
    }
    plot_bar_list[[variable]] <- plot_ordination(phy_b, ordination, color = variable) +
      scale +
      geom_point(size = 3) +
      theme_minimal() +
      ggtitle("Diversidad beta- Bray Curtis")
  }
  return(plot_bar_list)
} 
plot_jaccard_richness_by_variables <- function(physeq){
  plot_bar_list <- list()
  dist_jaccard <- distance(physeq, method = "jaccard")
  ordination <- ordinate(physeq, method = "PCoA", distance = dist_jaccard)
  for (variable in colnames(sample_data(physeq))) {
    if (is.numeric(sample_data(physeq)[[variable]])) {
      scale <- scale_color_gradient(low = "blue", high = "red")
    } else {
      scale <- scale_color_manual(values = paleta)
    }
    plot_bar_list[[variable]] <- plot_ordination(phy_b, ordination, color = variable) +
      scale +
      geom_point(size = 3) +
      theme_minimal() +
      ggtitle("Diversidad beta- Jaccard")
  }
  return(plot_bar_list)
} 
              
        