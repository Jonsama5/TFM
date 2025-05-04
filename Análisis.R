setwd("C:/Users/jonco/Desktop/Bioinformática/UNIR - 2º Cuatrimestre/TFM/Datos/")

library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(microbiome)
library(phyloseq)
library(patchwork)
library(Rtsne)
library(vegan)
library(FSA)

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
SAMPLES$SampleID <- rownames(SAMPLES)
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
#Tengo dos métodos:
#Método A
total = median(sample_sums(phy_g))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_g_a <- phy_g
phy_g_a = transform_sample_counts(phy_g, standf)
#Método B
total = median(sample_sums(phy_g))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_g_b <- phy_g
phy_g_b = transform_sample_counts(phy_g, standf)
#Estos dos hacen un escalado proporcional y no eliminan OTUs, el siguiente
#si que los elimina
phy_g_rarefied <- rarefy_even_depth(phy_g, sample.size = min(sample_sums(phy_g)), rngseed = 42)
phy_g_rarefied
summarize_phyloseq(phy_g_rarefied)

#Probando graficos con rarefacciones A, B, C + quitando singletons
  #Eliminar singletos
phy_g_filtered <- prune_taxa(taxa_sums(phy_g) > 1, phy_g)
phy_g_c_filtered <- prune_taxa(taxa_sums(phy_g_rarefied) > 1, phy_g_rarefied)
phy_g_filtered <- prune_taxa(taxa_sums(phy_g) > 1, phy_g)

A <- plot_bar(phy_g, fill = "Phylum", title = "Sin rarefaccion") +
  theme(legend.position = "none") #Sin rarefaccion
B <- plot_bar(phy_g_filtered, fill = "Phylum", title = "Sin rarefaccion + singletons") +
  theme(legend.position = "none") #Sin rarefaccion + singletons
C <- plot_bar(phy_g_a, fill = "Phylum", title = "Rarefaccion A") +
  theme(legend.position = "none") #Rarefaccion A
D <- plot_bar(phy_g_b, fill = "Phylum", title = "Rarefaccion B") +
  theme(legend.position = "none") #Rarefaccion B
E <- plot_bar(phy_g_rarefied, fill = "Phylum", title = "Rarefaccion C") +
  theme(legend.position = "none") #Rarefaccion C
G <- plot_bar(phy_g_c_filtered, fill = "Phylum", title = "Rarefaccion C + singletons") +
  theme(legend.position = "none") #Rarefaccion C + singletons

plot_list1 <- list(A, B, C)
plot_list2 <- list(D, E, G)
wrap_plots(plot_list1, ncol = 3)
wrap_plots(plot_list2, ncol = 3)
wrap_plots(list(D, E))

  #Miro quitando el que tiene más reads (20k)
readcount(phy_g)
phy_g <- subset_samples(physeq = phy_g, sample_names(phy_g) != "AN0582")
#Manualmente miro los singletons
taxa_sums <- taxa_sums(phy_g)
singletons <- taxa_sums[taxa_sums == 1]



plot_composition(phy_g)
plot_heatmap(phy_g)
plot_landscape(phy_g, transformation = "compositional")
plot_core(phy_g)

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


get_pca_var(pca_res_g)
fviz_pca_var(pca_res_g, col.var = "cos2", alpha.var = 0.5, gradient.cols = c("blue", "yellow", "red"),
              labelsize = 1)


### t-SNE ###
library(Rtsne)
#PrimeroRtsne#Primero hay que fijar una semilla
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

#Más colores
# Crear una función de gradiente
gradiente_func <- colorRampPalette(c("red", "blue", "green", "purple", "orange", "black", "cyan", "pink"))
# Generar 256 colores
gradiente_colores <- gradiente_func(length(unique(df$Location)) * 20)  # 20 veces más colores
# Escoger aleatoriamente
paleta <-  c("#FFFFBB", "#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000",
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
                      "#BDB76B", "#778899", "#FF8C00", "#BA55D3", "#4169E1", "#F08080", "#20B2AA", "#FF4500", "#ADFF2F", "#8A2BE2", 
                      "#87CEEB", "#800000", "#bb9e5a", "#a552d8", "#1c92d7", "#84d72e", "#c81d90", "#d75325", "#5e2eb7", "#3f6a43", 
                      "#b8c124", "#c62b38", "#d8a7df", "#51a1b1","#e5b43c", "#4cd3d4", "#cc4fe2", "#2a312b", "#675c58", "#57da5b", 
                      "#d1e232", "#c2307f", "#b624b0", "#9fdbc4", "#6c4b17", "#8a8cda", "#e480ad", "#3c367a", "#dc68ba", "#4c3541", 
                      "#29a948", "#6783f1", "#D4AC0D", "#1F618D", "#45B39D", "#B03A2E", "#6C3483", "#117864", "#CA6F1E", "#2E86C1", 
                      "#A93226", "#229954", "#884EA0", "#F5B041", "#2C3E50", "#AF601A", "#641E16", "#2471A3", "#D68910", "#0E6251", 
                      "#B2BABB", "#76448A", "#1D8348", "#D35400", "#566573", "#CB4335", "#7FB3D5", "#2980B9", "#A3E4D7", "#7D6608", 
                      "#196F3D", "#78281F", "#2F4F4F", "#FF1493", "#5F9EA0", "#FFD700", "#8B0000", "#7CFC00", "#D2691E", "#00CED1",
                        "#DC143C", "#556B2F", "#8B008B", "#8FBC8F", "#483D8B", "#00BFFF", "#A52A2A", "#DAA520",
                        "#800080", "#FF4500", "#9ACD32", "#20B2AA", "#FF00FF", "#00008B", "#CD5C5C", "#F4A460",
                        "#191970", "#FFA07A", "#808000", "#F5DEB3", "#C71585", "#FF6347", "#66CDAA", "#DB7093",
                        "#9932CC", "#B8860B", "#708090", "#ADFF2F", "#E9967A", "#BA55D3", "#B0C4DE", "#DDA0DD",
                        "#B22222", "#3CB371", "#F0E68C", "#BC8F8F", "#4682B4", "#FFB6C1", "#D8BFD8", "#008B8B",
                        "#9400D3", "#CD853F", "#5F5F5F", "#A9A9A9", "#D3D3D3", "#9B30FF", "#800000", "#FF69B4", "#8B7D7B", "#F4A300",
             "#7D9EC0", "#D8A7C5", "#B03060", "#A52A2A", "#708090", "#4B0082", "#FF4500", "#00FA9A",
             "#B22222", "#FFD700", "#3CB371", "#9932CC", "#8B4513", "#32CD32", "#FF8C00", "#FF7F50",
             "#F5F5F5", "#C0C0C0", "#4682B4", "#D2B48C", "#000080", "#B8860B", "#A52A2A", "#3E8E41",
             "#D2691E", "#F08080", "#FF6347", "#DDA0DD", "#00008B", "#D3D3D3", "#2E8B57", "#FF1493",
             "#ADFF2F", "#8B008B", "#8A2BE2", "#E6E6FA", "#FF00FF", "#7FFFD4", "#006400", "#E9967A",
             "#F4A460", "#FFD700", "#C71585", "#DC143C", "#C0FF3E", "#FFD700", "#008B8B", "#6A5ACD", "#A9A9A9", "#B0E0E6", "#8A2BE2", "#4B0082", "#7CFC00", "#FF1493", "#E6E6FA", "#D3D3D3",
             "#90EE90", "#D8BFD8", "#C71585", "#FFE4B5", "#ADFF2F", "#32CD32", "#00FA9A", "#FF6347",
             "#FF8C00", "#FF4500", "#7B68EE", "#D2B48C", "#8B4513", "#BC8F8F", "#D2691E", "#C71585",
             "#F4A460", "#8B008B", "#BDB76B", "#F0E68C", "#6A5ACD", "#708090", "#DDA0DD", "#FF7F50",
             "#20B2AA", "#FF00FF", "#4682B4", "#B0C4DE", "#DC143C", "#FF6347", "#4E9F3D", "#F0E68C",
             "#C0C0C0", "#B8860B", "#A52A2A", "#1E90FF", "#7FFFD4", "#8A2BE2", "#6A5ACD", "#D2691E",
             "#00CED1", "#F5F5F5", "#2E8B57", "#40E0D0", "#D3D3D3", "#7CFC00", "#FFD700", "#BC8F8F")
length(paleta)

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
### MDS ###

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






rank_names(phy_g)
rank_names(phy_g)
summarize_phyloseq(phy_b)
meta(phy_b)
tax_table(phy_b)
abundances(phy_b)
abundances(phy_b, "compositional")
readcount(phy_b)[1:5]
reads <- data.frame(reads = readcount(phy_b))
DF <- psmelt(phy_b)
summary(DF)
DF <- rownames_to_column(DF)
table(DF)
sample_sums(phy_b)
sample_names(phy_b)
sample_variables(phy_b)
ntaxa(phy_b)
topx <- top_taxa(phy_b, n = 15)
taxa(phy_b)
#Rarefacción
total = median(sample_sums(phy_g))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_g = transform_sample_counts(phy_g, standf)
phy_g


#Taxones más abundantes
phy15 <- prune_taxa(taxa = topx, phy_b)
phy15



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
ggplot(kingdom_df, aes(x = "", y = Count, fill = Superkingdom)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Composición de Reinos en el Microbioma") +
  theme_void() +  # Elimina el fondo para que sea un gráfico de "quesito"
  theme(legend.title = element_blank())   # Opcional: quitar el título de la leyenda

########### Gráficos composicionales ################
  #Comenzaré visualizando el global, más que nada porque espero muchos
  #microorganismos no clasificados
plot_bar(phy_g_a, fill = "Phylum", title = "Rarefaccion A") +
    scale_fill_manual(values = paleta) #Usando rarefacción A
  #Viendo que no se observa bien del todo, voy a hacer un pruning, primero veo
  #las abundancias y luego las voy a cortar
abundance_g <- taxa_sums(phy_g)
tax_table_data <- tax_table(phy_g)
abundance_sorted <- sort(abundance_data, decreasing = TRUE)
total_abundance <- sum(abundance_data)
half_abundance <- sum(total_abundance*0.5)
cumulative_abundance <- cumsum(abundance_sorted)
taxa_to_keep <- names(cumulative_abundance[cumulative_abundance <= half_abundance | cumulative_abundance == half_abundance])

# Realizar el prune_taxa para conservar solo estos taxones
top_taxa(phy_g, n = 20)
phy_top_rarefied <- prune_taxa(top_taxa(phy_g_rarefied, n=20), phy_g_rarefied)
composition_rarefied <- plot_bar(phy_top_rarefied, fill="Genus") + scale_fill_manual(values=paleta)
phy_top_a <- prune_taxa(top_taxa(phy_g_a, n=20), phy_g_a)
composition_a <- plot_bar(phy_top_a, fill="Genus") + scale_fill_manual(values=paleta)
    #He graficado la composicion de los mas abundante de rareficiado (unica opcion, si no queda fatal)
    #Comparo los metodos proporcionales con absoluto
comparacion_composicion_rareficada <- wrap_plots(list(composition_a, composition_rarefied),nrow = 2)
print(comparacion_composicion_rareficada)
#Los unk son porque no se ha elegido taxonomia en la BD, pero si que se tiene guardada
#como especie potencial. Eso habria que comentarlo. Tengo que mirar con solo bacterias 
#y solo virus porque a pesar de que los virus son muchos menos, al graficarlos globalmente
#se representan bastante.
top_taxa(phy_b, n = 20)
phy_top_b <- prune_taxa(top_taxa(phy_b, n=20), phy_b)
composition_b <- plot_bar(phy_top_b, fill="Phylum", title = "Composition b") + scale_fill_manual(values=paleta)
print(composition_b)


#Ordenar los gráficos por cada conjunto y rarefacción:
  #Rarefacciones a bacterias y virus solo + top taxones
  
#Método A
total = median(sample_sums(phy_b))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_b_a <- phy_b
phy_b_a = transform_sample_counts(phy_b, standf)
phy_b_a <- prune_taxa(top_taxa(phy_b_a, n = 20), phy_b_a)

total = median(sample_sums(phy_v))
standf = function(x, t=total) round(t*(x/sum(x)))  
phy_v_a <- phy_v
phy_v_a = transform_sample_counts(phy_v, standf)
phy_v_a <- prune_taxa(top_taxa(phy_v_a, n=20), phy_v_a)
#Metodo C
phy_b_rarefied <- rarefy_even_depth(phy_b, sample.size = min(sample_sums(phy_b)), rngseed = 42)
phy_b_rarefied <- prune_taxa(top_taxa(phy_b_rarefied, n=20), phy_b_rarefied)
phy_v_rarefied <- rarefy_even_depth(phy_v, sample.size = min(sample_sums(phy_v)), rngseed = 42)
phy_v_rarefied <- prune_taxa(top_taxa(phy_v_rarefied, n=20), phy_v_rarefied)
  #Nuestros objetos son: phy_b/v_a y phy_b/v:rarefied.
  #Ahora grafico la composicion con ambas rarefacciones y sin ella
composition_b_a <- plot_bar(phy_b_a, fill="Genus", title = "Composition b-a") + scale_fill_manual(values=paleta)
composition_b_rarefied <- plot_bar(phy_b_rarefied, fill="Genus", title = "Composition b-rarefied") + scale_fill_manual(values=paleta)
composition_b <- plot_bar(phy_b, fill="Family", title = "Composition b") + scale_fill_manual(values=paleta)

composition_v_a <- plot_bar(phy_v_a, fill="Genus", title = "Composition v-a") + scale_fill_manual(values=paleta)
composition_v_rarefied <- plot_bar(phy_v_rarefied, fill="Genus", title = "Composition v-rarefied") + scale_fill_manual(values=paleta)
composition_v <- plot_bar(phy_v, fill="Genus", title = "Composition v") + scale_fill_manual(values=paleta)

comparacion_composicion_b <- wrap_plots(list(composition_b, composition_b_a, composition_b_rarefied), nrow = 3)
comparacion_composicion_v <- wrap_plots(list(composition_v, composition_v_a, composition_v_rarefied), nrow = 3)
    #Todo esto es a nivel de filo, pero habría que ir cambiandolo a diferentes


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
global_por_rango <- plot_bar_by_ranks(phy_g_a)

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
plotwatcher(global_por_rango)


# ============================
# ANÁLISIS DE DIVERSIDAD
# ============================

#En este apartado encontrarás:
  #1 Diversidad alfa y beta general (todas las muestras) para dos índices en cada tipo
  #2 Diversidad alfa y beta para los diferentes tipos de variables con los mismos dos índices
  #3 Test esadísticos para evauar diferencias entre variables

#1#
#!!! MUY IMPORTANTE, EN ESTE CASO EL SIMPSON ES SIMPSON INVERSO, OSEA 0 BAJA DIVERSIDAD 1 ALTA !!!
  #Diversidad global en cada conjunto (ALFA)
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

#2#
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
Riqueza_g <- plot_richness_by_variables(phy_g)
Riqueza_b <- plot_richness_by_variables(phy_b)
Riqueza_v <- plot_richness_by_variables(phy_v)
plotwatcher(Riqueza_b)
plotwatcher(Riqueza_v)
plotwatcher(Riqueza_g)

#3#
  #Kruskal wallis
df_alpha <- estimate_richness(phy_b, measures = c("Shannon", "Simpson")) %>%
  cbind(sample_data(phy_b))
kw_shannon <- kruskal.test(Shannon ~ Location, data = df_alpha)
kw_simpson <- kruskal.test(Simpson ~ Location, data = df_alpha)
library(FSA) #El dunn test es para ver entre que grupos existen estas difs
dunnTest(Shannon ~ Location, data = df_alpha, method = "bh") #El metodo es por false discovery rate(investigar)
# Convertir los resultados a tibble
alpha_stats <- tibble(
  Index = c("Shannon", "Simpson"),
  Statistic = c(kw_shannon$statistic, kw_simpson$statistic),
  Df = c(kw_shannon$parameter, kw_simpson$parameter),
  P_value = c(kw_shannon$p.value, kw_simpson$p.value)
)

# Ver resultado
print(alpha_stats)

alpha_stats %>%
  gt() %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_header(title = "Resultados Kruskal-Wallis para diversidad alfa")


#1#
#Diversidad global en cada conjunto (BETA), #Aqui solo se aprecia la variable sex
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
}
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
}

bacteria_bray <- plot_bray_richness_by_variables(phy_b)
plotwatcher(bacteria_bray)
bacteria_jaccard <- plot_jaccard_richness_by_variables(phy_b)
plotwatcher(bacteria_jaccard)
virus_bray <- plot_bray_richness_by_variables(phy_v)
plotwatcher(virus_bray)
virus_jaccard <- plot_jaccard_richness_by_variables(phy_v)
plotwatcher(virus_jaccard)
#3#
  #PERMANOVA
library(vegan)
sample_df <- as(sample_data(phy_b), "data.frame")
permanova_res <- adonis2(dist_bray ~ Location, data = sample_df)
permanova_tbl <- as.data.frame(permanova_res) %>%
  rownames_to_column(var = "Term") %>%
  as_tibble()
print(permanova_tbl)
library(gt)

permanova_tbl %>%
  gt() %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_header(title = "Resultados de PERMANOVA (Bray-Curtis)")


###### Gráficos ###########

plot_bar(phy15, fill = "Family", title = "Título") 

plot_heatmap(phy15, method = "NMDS", distance = "bray")

plot_richness(phy15, measures = c("Shannon", "Simpson"))

phy.ord <- ordinate(phy, "NMDS", "bray")
plot_ordination(phy, phy.ord, type = "biplot") + geom_point(size=3) 
plot_ordination(phy, phy.ord, type="split", color="Class", 
                title="biplot") +  
  geom_point(size=3)

plot_net(phy15, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="Class", point_label="Genus")
net <- make_network(phy15, max.dist = 2)
plot_network(net, phy15, type = "taxa")
plot_tree(phy)






