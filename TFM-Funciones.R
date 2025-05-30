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
  #USO: Obtener datos de distribución de los datos y métricas de diversidad de
      #un objeto phyloseq
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
  #USO: Obtener una lista de gráficos composicionales de un objeto phyloseq en
      #función del rango taxonómico para su visualización con plotwatcher()
plot_bar_by_vars <- function(physeq, paleta = NULL) {
  library(phyloseq)
  library(ggplot2)
  library(scales)
  
  #Lista donde se guardaron los gráficos generados
  plot_bar_list <- list()
  
  #Iterar sobre cada rango taxonómico del objeto physeq
  for (var in sample_variables(physeq)) {
    #Creación del gráfico
    plot_bar_list[[var]] <- plot_bar(physeq, fill = "Genus") +
      facet_wrap(as.formula(paste("~", var)), scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste0("Composición por Género agrupada por ", var)) +
      scale_fill_manual(values = paleta)
  }
  
  return(plot_bar_list)
}
  #USO: Obtener una lista de gráficos composicionales de un objeto phyloseq en
      #función de las variables para su visualización con plotwatcher()
paleta <-  c("#FFFFFF", "#000000", "#FF0000", "#00FF00", "#0000FF", "#", "#FF00FF", "#00FFFF", "#800000", "#008000",
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
  #USO: Introduce el resultado de los creadores de gráficos para ver cada gráfico
      #de la lista (opción 1), ver en conjunto (opción 2) y salir (opción 3)
# Función para seleccionar los N taxones más abundantes
get_top_taxa <- function(physeq_object, top_n) {
  taxa_sums_all <- taxa_sums(physeq_object)
  top_taxa <- names(sort(taxa_sums_all, decreasing = TRUE)[1:top_n])
  pruned_object <- prune_taxa(top_taxa, physeq_object)
  return(pruned_object)
}
  #USO: Proporciona un objeto phyloseq y un número n de taxones para filtrar
      #n taxones más abundantes
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
  #USO: Obtener una lista de gráficos de diversidad alfa de un objeto phyloseq en
      #función de las variables para su visualización con plotwatcher()
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
  #USO: Realizar un test Kruskal-Wallis para la diversidad alfa de un objeto phyloseq
      #en función de las variables
#Test estadístico dunn posterior a Kruskall-Wallis
dunn_by_variables <- function(physeq) {
  df_alpha <- estimate_richness(physeq, measures = c("Shannon", "Simpson")) %>%
    cbind(sample_data(physeq))
  
  results <- list()
  
  for (variable in sample_variables(physeq)) {
    # Saltar variables con < 3 niveles únicos
    if (length(unique(df_alpha[[variable]])) < 3) next
    
    # Convertir la variable a factor para evitar errores
    df_alpha[[variable]] <- as.factor(df_alpha[[variable]])
    
    for (index in c("Shannon", "Simpson")) {
      formula <- as.formula(paste(index, "~", variable))
      
      # Manejar errores con tryCatch
      dunn_res <- tryCatch({
        dunnTest(formula, data = df_alpha, method = "bh")
      }, error = function(e) {
        message(paste("Error con variable:", variable, "y índice:", index))
        return(NULL)
      })
      
      if (is.null(dunn_res)) next
      
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
  #USO: Realizar un test post hoc Dunn para la diversidad alfa de un objeto phyloseq
      #en función de las variables tras aplicar "kw_by_variables"
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
  #USO: Obtener una lista de gráficos de diversidad beta de un objeto phyloseq en
      #función de las variables para su visualización con plotwatcher()             
#Creador de gráfico Deseq2 para el análisis diferencial de abundancias en función de variables
create_MAplot_from_phy <- function(physeq, condition_var, title = NULL) {
  library(DESeq2)
  library(phyloseq)
  library(ggplot2)
  library(ggrepel)
  # Asegurarse que la variable existe en el metadata
  if (!(condition_var %in% colnames(sample_data(physeq)))) {
    stop(paste("La variable", condition_var, "no está en los metadatos del phyloseq"))
  }
  # Guardar la variable título si necesario
  if (is.null(title)) {
    title = condition_var
  }
  # Crear objeto DESeq2 desde physeq usando la variable de diseño
  dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", condition_var)))
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Preparar dataframe para MA plot
  ma_data <- data.frame(
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    padj = res$padj,
    OTU = rownames(res)
  )
  
  # Extraer taxonomía para poder etiquetarla
  tax <- as.data.frame(tax_table(physeq))
  ma_data$Genus <- tax[match(ma_data$OTU, rownames(tax)), "Genus"]
  
  # Filtrar para etiquetas (padj < 0.05) y valores inválidos
  label_data <- subset(ma_data, padj < 0.05)
  label_data$Label <- ifelse(is.na(label_data$Genus), label_data$OTU, label_data$Genus)
  ma_data <- subset(ma_data, is.finite(baseMean) & is.finite(log2FoldChange) & !is.na(padj))
  
  # Crear plot
  plot <- ggplot(ma_data, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), size = 3, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
    geom_text_repel(data = label_data, aes(label = Label), size = 4, max.overlaps = 20) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 17, face = "bold"), 
      axis.text.x = element_text(angle = 70, vjust = 1, hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 8), 
      axis.title.y.left = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 14), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
    ) +
    labs(
      title = title,
      x = "Log10(Abundancia promedio + 1)",
      y = paste("Cambio en Log2 (por", condition_var, ")"),
      color = "Significativo (padj < 0.05)"
    )
  return(plot)
}
  #USO: Obtener un gráfico de abundancia diferencial de un objeto phyloseq proporcionado
      #en función de la variable dependiente (categórica) proporcionada
