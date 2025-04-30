setwd("C:/Users/jonco/Desktop/Bioinformática/UNIR - 2º Cuatrimestre/TFM/Datos/")

library(ggplot2)
library(readxl)
library(dplyr)
library(tibble)
library(microbiome)
library(phyloseq)
library(patchwork)
library(Rtsne)

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

### t-SNE ###
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
#Viendo que no se observa bien del todo, voy a hacer un pruning
#Realizar el prune_taxa para conservar solo estos taxones
top_taxa(phy_g, n = 20)
phy_top_rarefied <- prune_taxa(top_taxa(phy_g_rarefied, n=20), phy_g_rarefied)
composition_rarefied <- plot_bar(phy_top_rarefied, fill="Genus") + scale_fill_manual(values=paleta)
phy_top_a <- prune_taxa(top_taxa(phy_g_a, n=20), phy_g_a)
composition_a <- plot_bar(phy_top_a, fill="Genus") + scale_fill_manual(values=paleta)

#Los unk son porque no se ha elegido taxonomia en la BD, pero si que se tiene guardada
#como especie potencial. Eso habria que comentarlo. Tengo que mirar con solo bacterias 
#y solo virus porque a pesar de que los virus son muchos menos, al graficarlos globalmente
#se representan bastante.
top_taxa(phy_b, n = 20)
phy_top_b <- prune_taxa(top_taxa(phy_b, n=20), phy_b)
composition_b <- plot_bar(phy_top_b, fill="Phylum", title = "Composition b") + scale_fill_manual(values=paleta)
print(composition_b)