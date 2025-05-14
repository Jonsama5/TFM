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
library(tidyr)
source("TFM-Funciones.R")

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
  # Reemplaza los NA o strings vacíos por "Unidentified genus"
bacteria_tax$Genus[is.na(bacteria_tax$Genus) | bacteria_tax$Genus == ""] <- "Unidentified genus"
bacteria_tax$Genus[grepl("^unk_", bacteria_tax$Genus)] <- "Unclassified"

#Arreglo y separo en virus
colnames(virus_csv) <- sub("\\..*", "", colnames(virus_csv))
virus_csv <- subset(virus_csv, select = -AN0007)
virus_csv <- subset(virus_csv, select = -AN0425.1)

virus_otu <- virus_csv[,1:221]
virus_tax <- virus_csv[,222:length(virus_csv)]
virus_tax$Genus[is.na(virus_tax$Genus) | virus_tax$Genus == ""] <- "Unidentified genus"
virus_tax$Genus[grepl("^unk_", virus_tax$Genus)] <- "Unclassified"

#Voy a crear un dataset global
global_csv <- rbind(bacteria_csv, virus_csv)
global_otu <- global_csv[,1:221]
global_tax <- global_csv[,222:length(global_csv)]
global_tax$Genus[is.na(global_tax$Genus) | global_tax$Genus == ""] <- "Unidentified genus"
virus_tax$Genus[grepl("^unk_", virus_tax$Genus)] <- "Unclassified"

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
phy_g_a = transform_sample_counts(phy_g,function(x) 100* x/sum(x))
phy_b_a = transform_sample_counts(phy_b,function(x) 100* x/sum(x))
phy_v_a = transform_sample_counts(phy_v,function(x) 100* x/sum(x))
#Número de muestras por punto de muestreo (para hacer el mapa, no voy a sobrecargarlo
#con sitios que solo tienen 1, de cara a agruparlas por zonas)
sample_per_location <- ggplot(SAMPLES, aes(x = Location)) +
  geom_bar() +
  labs(title = "Número de muestras por punto de muestreo",
       x = "Punto de muestreo",
       y = "Frecuencia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))
#Resumen de parámetros globales, bacterianos y virales
phy_b_statistics <- phy_statistics(phy_b)
phy_v_statistics <- phy_statistics(phy_v)
phy_g_statistics <- phy_statistics(phy_g)
# Crear el tibble final con los estadísticos
final_statistics <- tibble(
  Type = c("Bacteria", "Virus", "Global"),
  TotalReads = c(phy_b_statistics$TotalReads, phy_v_statistics$TotalReads, phy_g_statistics$TotalReads),
  NumOTUs = c(phy_b_statistics$NumOTUs, phy_v_statistics$NumOTUs, phy_g_statistics$NumOTUs),
  MeanReads = c(phy_b_statistics$MeanReads, phy_v_statistics$MeanReads, phy_g_statistics$MeanReads),
  MedianReads = c(phy_b_statistics$MedianReads, phy_v_statistics$MedianReads, phy_g_statistics$MedianReads),
  SDReads = c(phy_b_statistics$SDReads, phy_v_statistics$SDReads, phy_g_statistics$SDReads),
  MinReads = c(phy_b_statistics$MinReads, phy_v_statistics$MinReads, phy_g_statistics$MinReads),
  MaxReads = c(phy_b_statistics$MaxReads, phy_v_statistics$MaxReads, phy_g_statistics$MaxReads),
  ShannonDiversity = c(phy_b_statistics$ShannonDiversity, phy_v_statistics$ShannonDiversity, phy_g_statistics$ShannonDiversity),
  SimpsonDiversity = c(phy_b_statistics$SimpsonDiversity, phy_v_statistics$SimpsonDiversity, phy_g_statistics$SimpsonDiversity),
  BrayDiversity = c(phy_b_statistics$BrayDiversity, phy_v_statistics$BrayDiversity, phy_g_statistics$BrayDiversity),
  JaccardDiversity = c(phy_b_statistics$JaccardDiversity, phy_v_statistics$JaccardDiversity, phy_g_statistics$JaccardDiversity)
)
colnames(final_statistics) <- c(
  "Type",
  "Total Reads",
  "Number of OTUs",
  "Mean Reads",
  "Median Reads",
  "Standard Deviation of Reads",
  "Minimum Reads",
  "Maximum Reads",
  "Shannon Diversity",
  "Simpson Diversity",
  "Bray-Curtis Diversity",
  "Jaccard Diversity"
)
# Usamos pivot_wider para convertir las métricas en filas y los tipos en columnas
final_statistics_wide <- final_statistics %>%
  pivot_longer(cols = -Type, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Type, values_from = Value)

# Crear la tabla con gt sin anotación científica
gt_table <- final_statistics_wide %>%
  gt() %>%
  tab_header(
    title = "Resumen de Estadísticos"
  ) %>%
  cols_label(
    Metric = "Parámetro",
    `Bacteria` = "Bacteria",
    `Virus` = "Virus",
    `Global` = "Global"
  ) %>%
  tab_options(
    table.width = pct(80)
  ) %>%
  fmt_number(decimals = 2)

# Mostrar la tabla
gt_table

# ============================
# ANÁLISIS COMPOSICIONAL
# ============================

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

#Proporción de lecturas por reino
plot_reads <- plot_bar(phy_g_a, fill = "Superkingdom") +
  labs(
    title = "Abundancia relativa total de Bacteria y Virus",
    x = "Muestras",
    y = "Abundancia relativa"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Dominio") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
#Filtrar por los 50 taxones más abundantes
phy_b_top <- get_top_taxa(phy_b_a, top_n = 50)
phy_v_top <- get_top_taxa(phy_v_a, top_n = 50)

#Composición
global_por_rango <- plot_bar_by_ranks(phy_g_a)
bacteria_composition <- plot_bar_by_ranks(phy_b_a)
virus_composition <- plot_bar_by_ranks(phy_v_a)

bacteria_composition_top <- plot_bar_by_ranks(phy_b_top)
virus_composition_top <- plot_bar_by_ranks(phy_v_top)

  #Composición a nivel de género
 bacteria_genera <- plot_bar(phy_b_top, fill = "Genus") +
  labs(title = paste("Composición de la microbiota bacteriana por género")) +
  theme_minimal(base_size = 12) +
  labs(x = "Muestras", fill = "Género") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
  ) +
  scale_y_continuous(limits = c(0, 100)) +  # Limita el eje Y hasta 100
  scale_fill_manual(values = paleta_genus)

 #Virus
virus_genera <- plot_bar(phy_v_top, fill = "Genus") +
  labs(title = paste("Composición de la microbiota vírica por género")) +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = paleta_genus) +
  labs(x = "Muestras", fill = "Género") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
  ) +
  scale_y_continuous(limits = c(0, 100))  # Limita el eje Y hasta 100

tax_df_v_top <- as.data.frame(tax_table(phy_v_top))
otu_sums_v_top <- taxa_sums(phy_v_top)
tax_df_v_top$Abundances <- otu_sums_v_top

write.csv(tax_df_b_top, file = "Abundancias de los 50 taxones bacterianos más abundantes.csv")
write.csv(tax_df_v_top, file = "Abundancias de los 50 taxones víricos más abundantes.csv")

#Visualizador
plotwatcher(global_por_rango)
plotwatcher(bacteria_composition)
plotwatcher(virus_composition)
plotwatcher(bacteria_composition_top)
plotwatcher(virus_composition_top)

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

################ Riqueza alfa global, bacteriana y viral ################ 
Riqueza_g <- plot_richness_by_variables(phy_g)
Riqueza_b <- plot_richness_by_variables(phy_b)
Riqueza_v <- plot_richness_by_variables(phy_v)
#Visualizador
plotwatcher(Riqueza_g)
plotwatcher(Riqueza_b)
plotwatcher(Riqueza_v)
#Resultados con test estadístico Kruskal-Wallis (Visualización)
bacteria_alpha_kw <- kw_by_variables(phy_b)
bacteria_alpha_kw <- bacteria_alpha_kw %>%
  gt() %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_header(title = "Resultados Kruskal-Wallis para diversidad alfa por variable")

virus_alpha_kw <- kw_by_variables(phy_v)
virus_alpha_kw <- virus_alpha_kw %>%
  gt() %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_header(title = "Resultados Kruskal-Wallis para diversidad alfa por variable")

#Prueba posterior de dunn
bacteria_dunn_res <- dunn_by_variables(phy_b)
virus_dunn_res <- dunn_by_variables(phy_v)

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

################ Riqueza beta bacteriana y viral ################ 
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

#Gráficos de diversidad beta por variable
riqueza_b_location <- plot_richness(phy_b, x = "Location", measures = c("Shannon", "Simpson"), title = "Diversidad bacteriana por localización") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_b_year <- plot_richness(phy_b, x = "Year", measures = c("Shannon", "Simpson"), title = "Diversidad bacteriana por año") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_b_zona <- plot_richness(phy_b, x = "Zone", measures = c("Shannon", "Simpson"), title = "Diversidad bacteriana por zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_b_sex <- plot_richness(phy_b, x = "Sex", measures = c("Shannon", "Simpson"), title = "Diversidad bacteriana por sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_v_location <- plot_richness(phy_v, x = "Location", measures = c("Shannon", "Simpson"), title = "Diversidad vírica por localización") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_v_year <- plot_richness(phy_v, x = "Year", measures = c("Shannon", "Simpson"), title = "Diversidad vírica por año") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_v_zona <- plot_richness(phy_v, x = "Zone", measures = c("Shannon", "Simpson"), title = "Diversidad vírica por zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

riqueza_v_sex <- plot_richness(phy_v, x = "Sex", measures = c("Shannon", "Simpson"), title = "Diversidad vírica por sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_boxplot(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_text(angle = 90))
ylab(label = "Índice alfa diversidad")

dist_bj <- distance(phy_b, method = "jaccard")
dist_vj <- distance(phy_v, method = "jaccard")
ordination_bj <- ordinate(phy_b, method = "PCoA", distance = dist_bj)
ordination_vj <- ordinate(phy_b, method = "PCoA", distance = dist_vj)

dist_bb <- distance(phy_b, method = "bray")
dist_vb <- distance(phy_v, method = "bray")
ordination_bb <- ordinate(phy_b, method = "PCoA", distance = dist_bb)
ordination_vb <- ordinate(phy_b, method = "PCoA", distance = dist_vb)


riqueza_b_jaccard_location <- plot_ordination(phy_b, ordination_bj, color = "Location") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_b_jaccard_year <- plot_ordination(phy_b, ordination_bj, color = "Year") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_gradient(low = "blue", high = "red")

riqueza_b_jaccard_zona <- plot_ordination(phy_b, ordination_bj, color = "Zona") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_b_jaccard_sex <- plot_ordination(phy_b, ordination_bj, color = "Sex") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_b_bray_location <- plot_ordination(phy_b, ordination_bb, color = "Location") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

riqueza_b_bray_year <- plot_ordination(phy_b, ordination_bb, color = "Year") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Bray-Curtis") + 
  scale_color_gradient(low = "blue", high = "red")

riqueza_b_bray_zona <- plot_ordination(phy_b, ordination_bb, color = "Zona") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

riqueza_b_bray_sex <- plot_ordination(phy_b, ordination_bb, color = "Sex") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

###

riqueza_v_jaccard_location <- plot_ordination(phy_v, ordination_vj, color = "Location") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_v_jaccard_year <- plot_ordination(phy_v, ordination_vj, color = "Year") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Jaccard") + 
  scale_color_gradient(low = "blue", high = "red")

riqueza_v_jaccard_zona <- plot_ordination(phy_v, ordination_vj, color = "Zona") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_v_jaccard_sex <- plot_ordination(phy_v, ordination_vj, color = "Sex") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Jaccard") + 
  scale_color_manual(values =  paleta)

riqueza_v_bray_location <- plot_ordination(phy_v, ordination_vb, color = "Location") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

riqueza_v_bray_year <- plot_ordination(phy_v, ordination_vb, color = "Year") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Bray-Curtis") + 
  scale_color_gradient(low = "blue", high = "red")

riqueza_v_bray_zona <- plot_ordination(phy_v, ordination_vb, color = "Zona") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

riqueza_v_bray_sex <- plot_ordination(phy_v, ordination_vb, color = "Sex") +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Bray-Curtis") + 
  scale_color_manual(values =  paleta)

#Test estadístico PERMANOVA para analizar qué variables producen efectos significativos sobre la diversidad
#Preparación (matriz de distancias y metadatos)
dist_bray_b <- distance(phy_b, method = "bray")
dist_bray_v <- distance(phy_v, method = "bray")
dist_jaccard_b <- distance(phy_b, method = "jaccard")
dist_jaccard_v <- distance(phy_v, method = "jaccard")
sample_df_b <- as(sample_data(phy_b), "data.frame") #Extraer metadatos
sample_df_v <- as(sample_data(phy_v), "data.frame") #Extraer metadatos

#PERMANOVA
permanova_b_res_location <- adonis2(dist_bray_b ~ Location, data = sample_df_b) #Permanova por localidad
permanova_b_res_year <- adonis2(dist_bray_b ~ Year, data = sample_df_b) #Permanova por año
permanova_b_res_latitude <- adonis2(dist_bray_b ~ Latitude, data = sample_df_b) #Permanova por latitud
permanova_b_res_longitude <- adonis2(dist_bray_b ~ Longitude, data = sample_df_b) #Permanova por longitud
permanova_b_res_sex <- adonis2(dist_bray_b ~ Sex, data = sample_df_b) #Permanova por sexo

permanova_v_res_location <- adonis2(dist_bray_v ~ Location, data = sample_df_v) #Permanova por localidad
permanova_v_res_year <- adonis2(dist_bray_v ~ Year, data = sample_df_v) #Permanova por año
permanova_v_res_latitude <- adonis2(dist_bray_v ~ Latitude, data = sample_df_v) #Permanova por latitud
permanova_v_res_longitude <- adonis2(dist_bray_v ~ Longitude, data = sample_df_v) #Permanova por longitud
permanova_v_res_sex <- adonis2(dist_bray_v ~ Sex, data = sample_df_v) #Permanova por sexo

permanova_b_res_location2 <- adonis2(dist_jaccard_b ~ Location, data = sample_df_b) #Permanova por localidad
permanova_b_res_year2 <- adonis2(dist_jaccard_b ~ Year, data = sample_df_b) #Permanova por año
permanova_b_res_latitude2 <- adonis2(dist_jaccard_b ~ Latitude, data = sample_df_b) #Permanova por latitud
permanova_b_res_longitude2 <- adonis2(dist_jaccard_b ~ Longitude, data = sample_df_b) #Permanova por longitud
permanova_b_res_sex2 <- adonis2(dist_jaccard_b ~ Sex, data = sample_df_b) #Permanova por sexo

permanova_v_res_location2 <- adonis2(dist_jaccard_v ~ Location, data = sample_df_v) #Permanova por localidad
permanova_v_res_year2 <- adonis2(dist_jaccard_v ~ Year, data = sample_df_v) #Permanova por año
permanova_v_res_latitude2 <- adonis2(dist_jaccard_v ~ Latitude, data = sample_df_v) #Permanova por latitud
permanova_v_res_longitude2 <- adonis2(dist_jaccard_v ~ Longitude, data = sample_df_v) #Permanova por longitud
permanova_v_res_sex2 <- adonis2(dist_jaccard_v ~ Sex, data = sample_df_v) #Permanova por sexo

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
  permanova_tbl_bray <- permanova_tbl_bray %>%
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
  
# ============================
# RESULTADOS DEL SCRIPT
# ============================

### Análisis preliminar ###
sample_per_location #Histograma para conocer la distribución de muestras
gt_table            #Tabla con parámetros de los datos

### Composición del microbioma ###
plot_reads  #Gráfico composicional del número de lecturas bacterianas y vírica

bacteria_genera #Gráfico de los géneros
virus_genera    #Gráfico de los géneros

### Análisis de diversidad ###
Riqueza_bacteriana #Gráfico de diversidad alfa bacteriana general
Riqueza_viral      #Gráfico de diversidad alfa viral general

#Riqueza por variables (relevantes)
riqueza_b_location
riqueza_b_zona
riqueza_b_year
riqueza_b_sex
riqueza_v_location
riqueza_v_zona
riqueza_v_year
riqueza_v_sex

riqueza_b_bray_location
riqueza_b_bray_zona
riqueza_b_bray_year

riqueza_b_bray_zona
riqueza_b_jaccard_location
riqueza_b_jaccard_zona
riqueza_b_jaccard_year
riqueza_b_jaccard_sex

riqueza_v_bray_location
riqueza_v_bray_zona
riqueza_v_bray_year
riqueza_v_bray_sex
riqueza_v_bray_zona
riqueza_v_jaccard_location
riqueza_v_jaccard_zona
riqueza_v_jaccard_year
riqueza_v_jaccard_sex
#Test estadísticos
#Kruskal wallis para alfa diversidad con dunn test para ver los grupos significativos
bacteria_alpha_kw
virus_alpha_kw

bacteria_dunn_res
virus_dunn_res

#PERMANOVA para beta diversidad
Permanova_res_bray
Permanova_res_jaccard
