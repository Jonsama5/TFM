setwd("C:/Users/jonco/Desktop/Bioinformática/UNIR - 2º Cuatrimestre/TFM/Datos/")

library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(microbiome) #Microbiome package
library(phyloseq) #Main microbiome package
library(patchwork) #wrap plots
library(Rtsne) #t-sne package
library(vegan) #Microbiome package; used for PERMANOVA
library(FSA) #dunn test
library(gt) #Plot tibbles

# ============================
# PREPROCESAMIENTO
# ============================
#Cargo los datos
bacteria_csv <- read.csv("BACTERIA_SPP_TABLE.csv")
virus_csv <- read.csv("VIRUSES_SPP_TABLE.csv")
metadata <- read.csv("Muestras mosquito Cameroon - R metadata.csv")
metadata[,1] <- sub("-.*", "", metadata[,1])

#Hay que hacer arreglos porque hay un par de duplicados. Además, separo otu y taxa.

#Arreglo y separo en bacterias
colnames(bacteria_csv) <- sub("\\..*", "", colnames(bacteria_csv))
bacteria_csv <- subset(bacteria_csv, select = -AN0007)
bacteria_csv <- subset(bacteria_csv, select = -AN0425.1)

bacteria_otu <- bacteria_csv[,1:221]
bacteria_tax <- bacteria_csv[,222:length(bacteria_csv)]

#Arreglo y separo en virus
colnames(virus_csv) <- sub("\\..*", "", colnames(virus_csv))
virus_csv <- subset(virus_csv, select = -AN0007)
virus_csv <- subset(virus_csv, select = -AN0425.1)

virus_otu <- virus_csv[,1:221]
virus_tax <- virus_csv[,222:length(virus_csv)]
#Voy a crear un dataset global
global_csv <- rbind(bacteria_csv, virus_csv)
global_otu <- global_csv[,1:221]
global_tax <- global_csv[,222:length(global_csv)]
# ============================
#     PREPARACIÓN PHYLOSEQ
# ============================
#Pongo rownames en común entre la tabla OTU y TAX
new_row_names_b <- sprintf("otu%03d", seq_len(nrow(bacteria_otu)))
new_row_names_v <- sprintf("otu%03d", seq_len(nrow(virus_otu)))
new_row_names_global <- sprintf("otu%03d", seq_len(nrow(global_otu)))
rownames(bacteria_otu) <- new_row_names_b
rownames(bacteria_tax) <- new_row_names_b
rownames(virus_otu) <- new_row_names_v
rownames(virus_tax) <- new_row_names_v
rownames(global_otu) <- new_row_names_global
rownames(global_tax) <- new_row_names_global
metadata <- metadata %>%
  mutate(Zone = case_when(
    Location %in% c("Doulougou", "Massila", "Daiguene", "Moussourtouk", "Makabay (Djarengo)", 
                    "Moulva", "Laf", "Badjawa", "Lougol", "Mayo Lebride", "Lamoudan", "Lainde Mbana", "Mayo Boki") ~ "Norte",
    Location %in% c("Tibati", "Palama", "Carrefour Poli", "Mgbandji", "Nkolondom", "Djaba", 
                    "Balda Bouri", "Banda", "Teckel", "Mabarangal'L", "Beka Goto", "Tekel") ~ "Central",
    Location %in% c("Bamendi", "Manda", "Mfelap", "Manchoutvi") ~ "Oeste",
    Location %in% c("Oitibili", "Ahala", "Obala", "Essos") ~ "Suroeste",
    Location %in% c("Gado Badzere", "Zembe Borongo", "Mayos", "Mbalmayo", "Lougol", "Mayo Dafan", 
                    "Avebe", "Nlozok", "DombÃ©") ~ "Sureste",
    Location %in% c("Afan-Essokye", "Foulassi I") ~ "Sur",
    TRUE ~ NA_character_  # En caso de que no coincida con ninguna localidad, asigna NA
  ))
samples <- metadata
samples <- samples %>% tibble::column_to_rownames("Sample.ID")

# ============================
#          PHYLOSEQ
# ============================  
#Bacterias
bacteria_OTU <- otu_table(bacteria_otu, taxa_are_rows = TRUE)
bacteria_taxa <- as.matrix(bacteria_tax)
bacteria_TAX <- tax_table(bacteria_taxa)

SAMPLES <- sample_data(samples)

phy_b <- phyloseq(bacteria_OTU, bacteria_TAX, SAMPLES)

#Virus
virus_OTU <- otu_table(virus_otu, taxa_are_rows = TRUE)
virus_taxa <- as.matrix(virus_tax)
virus_TAX <- tax_table(virus_taxa)

phy_v <- phyloseq(virus_OTU, virus_TAX, SAMPLES)
#Global
global_OTU <- otu_table(global_otu, taxa_are_rows = TRUE)
global_taxa <- as.matrix(global_tax)
global_TAX <- tax_table(global_taxa)

phy_g <- phyloseq(global_OTU, global_TAX, SAMPLES)

# ============================
#    ANÁLISIS PRELIMINAR
# ============================
phy_g
phy_b
phy_v
#Lecturas
sum(readcount(phy_g))
sum(readcount(phy_b))
sum(readcount(phy_v))
#Características
summarize_phyloseq(phy_g)
summarize_phyloseq(phy_b)
summarize_phyloseq(phy_v)
#Rarefacción
total = median(sample_sums(phy_g))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_g_a <- phy_g
phy_g_a = transform_sample_counts(phy_g, standf)

total = median(sample_sums(phy_b))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_b_a <- phy_b
phy_b_a = transform_sample_counts(phy_b, standf)

total = median(sample_sums(phy_v))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_v_a <- phy_v
phy_v_a = transform_sample_counts(phy_v, standf)
#Número de muestras por punto de muestreo (para hacer el mapa, no voy a sobrecargarlo
#con sitios que solo tienen 1, de cara a agruparlas por zonas)
ggplot(SAMPLES, aes(x = Location)) +
  geom_bar() +
  labs(title = "Número de muestras por punto de muestreo",
       x = "Punto de muestreo",
       y = "Frecuencia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

# ============================
# ANÁLISIS EXPLORATORIO (ML)
# ============================

#Creo un DF para usarlo en estos gráficos
t_df <- data.frame(t(global_otu))
t_df$ID <- rownames(t(global_otu))
colnames(t_df)[colnames(t_df) == "ID"] <- "SampleID"
colnames(metadata)[colnames(metadata) == "Sample.ID"] <- "SampleID"
rownames(t_df) <- NULL
df <- merge.data.frame(metadata, t_df, by = "SampleID")

  #Creo una paleta de colores distinguible
# Crear una función de gradiente
gradiente_func <- colorRampPalette(c("red", "blue", "green", "purple", "orange", "black", "cyan", "pink"))
# Generar 256 colores
gradiente_colores <- gradiente_func(length(unique(df$Location)) * 20)  # 20 veces más colores
# Escoger aleatoriamente
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
             "#ADFF2F", "#FF69B4", "#CD5C5C", "#8B4513", "#F0E68C", "#20B2AA", "#6A5ACD","#9370DB", "#00FF7F", "#9ACD32", "#BDB76B",     "#778899", "#FF8C00", "#BA55D3", "#4169E1", "#F08080", "#20B2AA", "#FF4500", "#ADFF2F", "#8A2BE2", "#87CEEB", "#800000")

################ PCA ################ 
#PCA --> No utilizable, los CP explican muy poca variabilidad
phy_pca <- transformSampleCounts(phy_g, function(x) x/sum(x)) #Proporciones relativas
otu_mat_g <- as(otu_table(phy_pca), "matrix") #Pasamos a una matriz de abundancias (lo extraemos del objeto clase phyloseq)
otu_mat_g <- t(otu_mat_g) #Las muestras deben ser filas
pca_res_g <- prcomp(otu_mat_g, center = TRUE, scale. = TRUE)
#Para plotearlo:
pca_scores <- as.data.frame(pca_res_g$x) #Coordenadas
sample_data_df <- as.data.frame(sample_data(phy_g)) #Extraigo metadatos y los guardo en un DF
sample_data_df$SampleID <- rownames(sample_data_df) #Creo columna de muestra en metadatos para hacer un bind
pca_plot_df <- cbind(pca_scores, sample_data_df) #Hago cbind pero tengo que eliminar la colum

ggplot(pca_plot_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", round(summary(pca_res_g)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res_g)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()
#Voy a mirar los autovalores
library(stat)
library(factoextra)
head(get_eigenvalue(pca_res_g))

################ t-SNE ################ 
library(Rtsne)
#Primero hay que fijar una semilla
set.seed(1995)
str(DF)
#Preparamos set de datos
data.tsne <-sapply(df[,9:length(df)], as.numeric)
#Aplicamos el método al DF y extraemos distancias
data.tsne.unique <- unique(data.tsne)
tsne <- Rtsne(X = data.tsne.unique, dims = 2)
tsne_result <- data.frame(tsne$Y) #Y[,1:2] Para primeras 2 dimensiones
# Graficamos
#Colores aleatorios
library(randomcoloR)
# Crear tantos colores como ubicaciones únicas tengas
palette_location <- distinctColorPalette(k = length(unique(df$Location)))
palette_location <- as.vector(palette_location)
class(palette_location)

t_sne_location <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Location)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Location", x = "X1", y = "X2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) + 
  scale_color_manual(values = paleta)

t_sne_year <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Year)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Year", x = "X1", y = "X2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) + 
  scale_color_gradient(low = "blue", high = "red")

t_sne_month <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Month)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Month", x = "X1", y = "X2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) + 
  scale_color_gradient(low = "blue", high = "red")

#He reformateado a numérico porque si no no podría representar una escala continua
df$Latitude <- gsub(",", ".", df$Latitude)
df$Latitude <- as.numeric(df$Latitude)
df$Longitude <- gsub(",", ".", df$Longitude)
df$Longitude <- as.numeric(df$Longitude)
df$Month <- as.numeric(df$Month)
t_sne_latitude <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Latitude)) + 
  geom_point(size = 3) + 
  labs(title = "t-SNE Latitude", x = "X1", y = "X2") + 
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradient(low = "blue", high = "red")

t_sne_longitude <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Longitude)) + 
  geom_point(size = 3) + 
  labs(title = "t-SNE Longitude", x = "X1", y = "X2") + 
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_gradient(low = "blue", high = "red")

t_sne_sex <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Sex)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Sex", x = "X1", y = "X2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) + 
  scale_color_viridis_d()

list_tsne <- mget(ls(pattern = "^t_sne_"))
wrap_plots(list_tsne, ncol = 2, nrow = 3)

################ MDS ################ 

data.mds <-sapply(df[,9:length(df)], as.numeric)
# Utilizamos la funcion dist para calcular la matriz de distancias euclideas
# matriz NxN de distancias entre todos los puntos
distances <- dist(data.mds, method = 'euclidean') #Con esto creamos la matriz de distancias necesaria para el cmdscale()
# Utilizamos la función cmdscale para realizar el MSD
mds.results <- cmdscale(distances, eig=TRUE, k=2, x.ret=TRUE)
# Calculamos la varianza explicada
varianza.explicada <- mds.results$eig/sum(mds.results$eig) * 100
# Sacamos en un dataframe los puntos del mds
mds.df <- data.frame(mds.results$points)
# Grafico
mds_location <- ggplot(mds.df, aes(x=X1, y=X2, color=df$Location)) +
  geom_point(size=3) + 
  scale_color_manual(values= paleta) +
  labs(title="MDS Location", x="Dimension 1 (X1)", y="Dimension 2 (X2)") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

mds_year <- ggplot(tsne_result, aes(x = X1, y = X2, color = df$Year)) +
  geom_point(size = 3) +
  labs(title = "MDS Year", x = "X1", y = "X2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) + 
  scale_color_gradient(low = "blue", high = "red")

mds_month <- ggplot(mds.df, aes(x=X1, y=X2, color=df$Month)) +
  geom_point(size=3) + 
  scale_color_gradient(low = "blue", high = "red") +
  labs(title="MDS Month", x="Dimension 1 (X1)", y="Dimension 2 (X2)") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))

mds_latitude <- ggplot(mds.df, aes(x=X1, y=X2, color=df$Latitude)) +
  geom_point(size=3) + 
  scale_color_gradient(low = "blue", high = "red") +
  labs(title="MDS Latitude", x="Dimension 1 (X1)", y="Dimension 2 (X2)") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) 

mds_longitude <- ggplot(mds.df, aes(x=X1, y=X2, color=df$Longitude)) +
  geom_point(size=3) + 
  scale_color_gradient(low = "blue", high = "red") +
  labs(title="MDS Longitude", x="Dimension 1 (X1)", y="Dimension 2 (X2)") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) 

mds_sex <- ggplot(mds.df, aes(x=X1, y=X2, color=df$Sex)) +
  geom_point(size=3) + 
  scale_color_viridis_d() +
  labs(title="MDS Sex", x="Dimension 1 (X1)", y="Dimension 2 (X2)") +
  theme_classic() + 
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5)) 
list_mds <- mget(ls(pattern = "^mds_"))
wrap_plots(list_mds, ncol = 2, nrow = 3)

# ============================
# ANÁLISIS COMPOSICIONAL
# ============================

#### Gráfico de sectores de la composición por reino ####
# Obtener la tabla de taxonomía
tax_table_data <- tax_table(phy_g)
# Asegúrate de que la columna correspondiente al reino esté presente (normalmente es la primera columna)
kingdom_data <- tax_table_data[, "Superkingdom"]
# Contar la cantidad de muestras por cada reino
kingdom_counts <- table(kingdom_data)
# Convertir a un data.frame para usarlo en ggplot
kingdom_df <- as.data.frame(kingdom_counts)
colnames(kingdom_df) <- c("Superkingdom", "Count")
# Crear el gráfico de "quesito"
Plot_sector <- ggplot(kingdom_df, aes(x = "", y = Count, fill = Superkingdom)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Composición de Reinos en el Microbioma") +
  theme_void() +  # Elimina el fondo para que sea un gráfico de "quesito"
  theme(legend.title = element_blank())   # Opcional: quitar el título de la leyenda

########### Gráficos composicionales ################

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

#Visualizador general de gráficos
plotwatcher <- function(plotlist) {
  cat("Welcome to plotwatcher, write a number in the terminal for the following:") 
  option <- readline(prompt = "Select an option: \n(1) Visualize plot by plot\n(2) Visualize plots together\n(3) Exit\n")
  plotlist <- plotlist
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


#Objetos
  #Composición
global_por_rango <- plot_bar_by_ranks(phy_g_a)
bacteria_composition <- plot_bar_by_ranks(phy_b_a)
virus_composition <- plot_bar_by_ranks(phy_v_a)

  #Visualizador
Plot_sector
plotwatcher(global_por_rango)
plotwatcher(bacteria_composition)
plotwatcher(virus_composition)

################ RESULTADOS DEL APARTADO ################ 
plotwatcher(global_por_rango)
plotwatcher(bacteria_composition)
plotwatcher(virus_composition)
#Quedaría elegir gráfico definitivo, ver si lo hago con pruning (porque como está
#Se visualiza bien), y solucionar los colores para fijarlos

# ============================
# ANÁLISIS DE DIVERSIDAD
# ============================

################ Riqueza general ################ 

Riqueza_bacteriana <- plot_richness(phy_b, x = "Country", measures = c("Shannon", "Simpson"), title = "Diversidad bacteriana") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral <- plot_richness(phy_v, x = "Country", measures = c("Shannon", "Simpson"), title = "Diversidad viral") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

################ Riqueza por variable ################ 

#ALPHA
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
  }##Función de riqueza alfa por variables
  
    #Objetos de riqueza
    Riqueza_g <- plot_richness_by_variables(phy_g)
    Riqueza_b <- plot_richness_by_variables(phy_b)
    Riqueza_v <- plot_richness_by_variables(phy_v)
    #Visualizador
    plotwatcher(Riqueza_g)
    plotwatcher(Riqueza_b)
    plotwatcher(Riqueza_v)

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
  } ##Función de Kruskal Wallis por variables
  
    #Resultados (Visualización)
    bacteria_alpha_kw <- kw_by_variables(phy_b_a)
    bacteria_alpha_kw <- bacteria_alpha_kw %>%
      gt() %>%
      fmt_number(columns = where(is.numeric), decimals = 3) %>%
      tab_header(title = "Resultados Kruskal-Wallis para diversidad alfa por variable")
    
    virus_alpha_kw <- kw_by_variables(phy_v_a)
    virus_alpha_kw <- virus_alpha_kw %>%
      gt() %>%
      fmt_number(columns = where(is.numeric), decimals = 3) %>%
      tab_header(title = "Resultados Kruskal-Wallis para diversidad alfa por variable")
    
  #Test estadístico dunn para analizar qué grupos son diferentes significativamente
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
    } ##Función de dunn por variables
    bacteria_dunn_res <- dunn_by_variables(phy_b_a)
    virus_dunn_res <- dunn_by_variables(phy_v_a)
    
    #Resultados (visualización)
     bacteria_dunn_res <- bacteria_dunn_res %>%
      filter(P_adjusted < 0.05) %>%
      gt() %>%
      fmt_number(columns = where(is.numeric), decimals = 3) %>%
      tab_header(title = "Resultados de pruebas de Dunn post-hoc (Bacteria)")
     virus_dunn_res <- virus_dunn_res %>%
       filter(P_adjusted < 0.05) %>%
       gt() %>%
       fmt_number(columns = where(is.numeric), decimals = 3) %>%
       tab_header(title = "Resultados de pruebas de Dunn post-hoc (Virus)")

#BETA
     #Creador de gráficos de riqueza beta para cada rango taxonómico
     plot_bray_richness_by_variables <- function(physeq){
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
     } ##Función de diversidad beta Bray-Curtis por variables
     plot_jaccard_richness_by_variables <- function(physeq){
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
     } ##Función de diversidad beta Jaccard por variables
     
       #Objetos de riqueza
       bacteria_bray <- plot_bray_richness_by_variables(phy_b)
       bacteria_jaccard <- plot_jaccard_richness_by_variables(phy_b)
       virus_bray <- plot_bray_richness_by_variables(phy_v)
       virus_jaccard <- plot_jaccard_richness_by_variables(phy_v)
       
       #Visualizador
       plotwatcher(bacteria_bray)
       plotwatcher(bacteria_jaccard)
       plotwatcher(virus_bray)
       plotwatcher(virus_jaccard)
     
       #Test estadístico PERMANOVA para analizar qué variables producen efectos significativos sobre la diversidad
        #Preparación (matriz de distancias y metadatos)
       dist_bray_b <- distance(phy_b_a, method = "bray")
       dist_bray_v <- distance(phy_v_a, method = "bray")
       dist_jaccard_b <- distance(phy_b_a, method = "jaccard")
       dist_jaccard_v <- distance(phy_v_a, method = "jaccard")
       sample_df_b <- as(sample_data(phy_b_a), "data.frame") #Extraer metadatos
       sample_df_v <- as(sample_data(phy_v_a), "data.frame") #Extraer metadatos
        
        #PERMANOVA
       permanova_b_res_location <- adonis2(dist_bray_b ~ Location, data = sample_df_b) #Permanova por localidad
       permanova_b_res_year <- adonis2(dist_bray_b ~ Year, data = sample_df_b) #Permanova por año
       permanova_b_res_latitude <- adonis2(dist_bray_b ~ Latitude, data = sample_df_b) #Permanova por latitud
       permanova_b_res_longitude <- adonis2(dist_bray_b ~ Longitude, data = sample_df_b) #Permanova por longitud
       permanova_b_res_sex <- adonis2(dist_bray_b ~ Sex, data = sample_df_b) #Permanova por sexo
       
       permanova_v_res_location <- adonis2(dist_bray_v ~ Location, data = sample_df) #Permanova por localidad
       permanova_v_res_year <- adonis2(dist_bray_v ~ Year, data = sample_df) #Permanova por año
       permanova_v_res_latitude <- adonis2(dist_bray_v ~ Latitude, data = sample_df) #Permanova por latitud
       permanova_v_res_longitude <- adonis2(dist_bray_v ~ Longitude, data = sample_df) #Permanova por longitud
       permanova_v_res_sex <- adonis2(dist_bray_v ~ Sex, data = sample_df) #Permanova por sexo

       permanova_b_res_location2 <- adonis2(dist_jaccard_b ~ Location, data = sample_df_b) #Permanova por localidad
       permanova_b_res_year2 <- adonis2(dist_jaccard_b ~ Year, data = sample_df_b) #Permanova por año
       permanova_b_res_latitude2 <- adonis2(dist_jaccard_b ~ Latitude, data = sample_df_b) #Permanova por latitud
       permanova_b_res_longitude2 <- adonis2(dist_jaccard_b ~ Longitude, data = sample_df_b) #Permanova por longitud
       permanova_b_res_sex2 <- adonis2(dist_jaccard_b ~ Sex, data = sample_df_b) #Permanova por sexo
       
       permanova_v_res_location2 <- adonis2(dist_jaccard_v ~ Location, data = sample_df) #Permanova por localidad
       permanova_v_res_year2 <- adonis2(dist_jaccard_v ~ Year, data = sample_df) #Permanova por año
       permanova_v_res_latitude2 <- adonis2(dist_jaccard_v ~ Latitude, data = sample_df) #Permanova por latitud
       permanova_v_res_longitude2 <- adonis2(dist_jaccard_v ~ Longitude, data = sample_df) #Permanova por longitud
       permanova_v_res_sex2 <- adonis2(dist_jaccard_v ~ Sex, data = sample_df) #Permanova por sexo
       
        #visualización PERMANOVA (Bray)
         # Convertir resultados bacterias
         permanova_tbl_b <- bind_rows(
           Location  = as_tibble(permanova_b_res_location, rownames = "Term"),
           Year      = as_tibble(permanova_b_res_year, rownames = "Term"),
           Latitude  = as_tibble(permanova_b_res_latitude, rownames = "Term"),
           Longitude = as_tibble(permanova_b_res_longitude, rownames = "Term"),
           Sex       = as_tibble(permanova_b_res_sex, rownames = "Term"),
           .id = "Variable"
         ) %>% mutate(Grupo = "Bacterias")
         
         # Convertir resultados virus
         permanova_tbl_v <- bind_rows(
           Location  = as_tibble(permanova_v_res_location, rownames = "Term"),
           Year      = as_tibble(permanova_v_res_year, rownames = "Term"),
           Latitude  = as_tibble(permanova_v_res_latitude, rownames = "Term"),
           Longitude = as_tibble(permanova_v_res_longitude, rownames = "Term"),
           Sex       = as_tibble(permanova_v_res_sex, rownames = "Term"),
           .id = "Variable"
         ) %>% mutate(Grupo = "Virus")
         
         permanova_tbl_bray <- bind_rows(permanova_tbl_b, permanova_tbl_v)
         # Usamos dplyr para dejar en blanco los valores repetidos de la columna Variable
         permanova_tbl_bray <- permanova_bray %>%
           arrange(Grupo, Variable, Term) %>%
           mutate(Variable = if_else(duplicated(Variable) & duplicated(paste(Grupo, Variable)), "", Variable))
         Permanova_res_bray <-permanova_tbl_bray %>%
           select(Grupo, Variable, Term, everything()) %>%
           gt() %>%
           fmt_number(columns = where(is.numeric), decimals = 3) %>%
           tab_header(title = "Resultados de PERMANOVA (Bray-Curtis)")

         #visualización PERMANOVA (Jaccard)
           # Convertir resultados bacterias
           permanova_tbl_b2 <- bind_rows(
             Location  = as_tibble(permanova_b_res_location2, rownames = "Term"),
             Year      = as_tibble(permanova_b_res_year2, rownames = "Term"),
             Latitude  = as_tibble(permanova_b_res_latitude2, rownames = "Term"),
             Longitude = as_tibble(permanova_b_res_longitude2, rownames = "Term"),
             Sex       = as_tibble(permanova_b_res_sex2, rownames = "Term"),
             .id = "Variable"
           ) %>% mutate(Grupo = "Bacterias")
           
           # Convertir resultados virus
           permanova_tbl_v2 <- bind_rows(
             Location  = as_tibble(permanova_v_res_location2, rownames = "Term"),
             Year      = as_tibble(permanova_v_res_year2, rownames = "Term"),
             Latitude  = as_tibble(permanova_v_res_latitude2, rownames = "Term"),
             Longitude = as_tibble(permanova_v_res_longitude2, rownames = "Term"),
             Sex       = as_tibble(permanova_v_res_sex2, rownames = "Term"),
             .id = "Variable"
           ) %>% mutate(Grupo = "Virus")
           
           permanova_tbl2 <- bind_rows(permanova_tbl_b2, permanova_tbl_v2)
           
           # Usamos dplyr para dejar en blanco los valores repetidos de la columna Variable
           permanova_tbl_clean2 <- permanova_tbl2 %>%
             arrange(Grupo, Variable, Term) %>%
             mutate(Variable = if_else(duplicated(Variable) & duplicated(paste(Grupo, Variable)), "", Variable))
           
           Permanova_res_jaccard <- permanova_tbl_clean2 %>%
             select(Grupo, Variable, Term, everything()) %>%
             gt() %>%
             fmt_number(columns = where(is.numeric), decimals = 3) %>%
             tab_header(title = "Resultados de PERMANOVA (Jaccard)")

################ RESULTADOS DEL APARTADO ################ 
  #Riqueza alfa general
  Riqueza_bacteriana
  Riqueza_viral
  #Riqueza alfa por variable
  plotwatcher(Riqueza_g)
  plotwatcher(Riqueza_b)
  plotwatcher(Riqueza_v)
  #Kruskal wallis para riqueza alfa
  bacteria_alpha_kw
  virus_alpha_kw
  #Test dunn para kruskal wallis pairwise (sorteado por significativo)
  bacteria_dunn_res
  virus_dunn_res
  #Riqueza beta por variable
  plotwatcher(bacteria_bray)
  plotwatcher(bacteria_jaccard)
  plotwatcher(virus_bray)
  plotwatcher(virus_jaccard)
  #Test estadístico PERMANOVA por variable de diversidad beta
  Permanova_res_bray
  Permanova_res_jaccard
