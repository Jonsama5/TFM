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
library(ggpubr) #Arrange plots
source("TFM-Funciones.R")

# ============================
# PREPROCESAMIENTO
# ============================
#Cargo los datos
bacteria_csv <- read.csv("BACTERIA_SPP_TABLE.csv")
virus_csv <- read.csv("VIRUSES_SPP_TABLE.csv")
plas_csv <- read.csv("plasmodio_SPP_table.csv")
metadata <- read.csv("Muestras mosquito Cameroon - R metadata.csv")
metadata[,1] <- sub("-.*", "", metadata[,1])

#Hay que hacer arreglos porque hay un par de duplicados. Además, separo otu y taxa.

#Arreglo y separo en bacterias
colnames(bacteria_csv) <- sub("\\..*", "", colnames(bacteria_csv))
bacteria_csv <- subset(bacteria_csv, select = -AN0007)
bacteria_csv <- subset(bacteria_csv, select = -AN0425.1)

bacteria_otu <- bacteria_csv[,1:221]
bacteria_tax <- bacteria_csv[,222:length(bacteria_csv)]
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
#Arreglo y separo en plasmodio
colnames(plas_csv)[1:223] <- sub("\\..*", "", colnames(plas_csv)[1:223])
colnames(plas_csv)[colnames(plas_csv) == "AN007"] <- "AN0070"
plas_csv <- subset(plas_csv, select = -AN0070)
plas_csv <- subset(plas_csv, select = -AN0455)

plas_otu <- plas_csv[,1:221]
plas_tax <- plas_csv[,222:length(plas_csv)]
#Voy a crear un dataset global
global_csv <- rbind(bacteria_csv, virus_csv, plas_csv)
global_otu <- global_csv[,1:221]
global_tax <- global_csv[,222:length(global_csv)]
global_tax$Genus[is.na(global_tax$Genus) | global_tax$Genus == ""] <- "Unidentified genus"
global_tax$Genus[grepl("^unk_", global_tax$Genus)] <- "Unclassified"

# ============================
#     PREPARACIÓN PHYLOSEQ
# ============================
#Pongo rownames en común entre la tabla OTU y TAX
new_row_names_b <- sprintf("otu%03d", seq_len(nrow(bacteria_otu)))
new_row_names_v <- sprintf("otu%03d", seq_len(nrow(virus_otu)))
new_row_names_p <- sprintf("otu%03d", seq_len(nrow(plas_otu)))
new_row_names_global <- sprintf("otu%03d", seq_len(nrow(global_otu)))
rownames(bacteria_otu) <- new_row_names_b
rownames(bacteria_tax) <- new_row_names_b
rownames(virus_otu) <- new_row_names_v
rownames(virus_tax) <- new_row_names_v
rownames(plas_otu) <- new_row_names_p
rownames(plas_tax) <- new_row_names_p
rownames(global_otu) <- new_row_names_global
rownames(global_tax) <- new_row_names_global
#Introduzco metadatos
  #Localidades de los puntos de muestreo
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
  #Abundancias en los metadatos (para futuras regresiones o contrastes)
plas_otu #Abundancias
plas_tax #Clasificación
plas_by_sample <- colSums(plas_otu)
metadata$Plasmodium_abundance <- plas_by_sample
  #Mosquitos infectados
metadata <- metadata %>%
  mutate(Infected = if_else(Plasmodium_abundance > 0, "Yes", "No"))
metadata <- metadata %>%
  mutate(across(c(Country, Location, Sex, Zone, Infected), as.factor))

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

#Plasmodio
plas_OTU <- otu_table(plas_otu, taxa_are_rows = TRUE)
plas_taxa <- as.matrix(plas_tax)
plas_TAX <- tax_table(plas_taxa)

phy_p <- phyloseq(plas_OTU, plas_TAX, SAMPLES)

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
phy_p
#Lecturas
sum(readcount(phy_g))
sum(readcount(phy_b))
sum(readcount(phy_v))
sum(readcount(phy_p))
#Características
summarize_phyloseq(phy_g)
summarize_phyloseq(phy_b)
summarize_phyloseq(phy_v)
summarize_phyloseq(phy_p)
#Rarefacción
phy_g_a = transform_sample_counts(phy_g,function(x) 100* x/sum(x))
phy_b_a = transform_sample_counts(phy_b,function(x) 100* x/sum(x))
phy_v_a = transform_sample_counts(phy_v,function(x) 100* x/sum(x))
phy_p_a = transform_sample_counts(phy_p,function(x) 100* x/sum(x))
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
phy_p_statistics <- phy_statistics(phy_p)
phy_g_statistics <- phy_statistics(phy_g)
# Crear el tibble final con los estadísticos
final_statistics <- tibble(
  Type = c("Bacteria", "Virus", "Plasmodium", "Global"),
  TotalReads = c(phy_b_statistics$TotalReads, phy_v_statistics$TotalReads, phy_p_statistics$TotalReads, phy_g_statistics$TotalReads),
  NumOTUs = c(phy_b_statistics$NumOTUs, phy_v_statistics$NumOTUs, phy_p_statistics$NumOTUs, phy_g_statistics$NumOTUs),
  MeanReads = c(phy_b_statistics$MeanReads, phy_v_statistics$MeanReads, phy_p_statistics$MeanReads, phy_g_statistics$MeanReads),
  MedianReads = c(phy_b_statistics$MedianReads, phy_v_statistics$MedianReads, phy_p_statistics$MedianReads, phy_g_statistics$MedianReads),
  SDReads = c(phy_b_statistics$SDReads, phy_v_statistics$SDReads, phy_p_statistics$SDReads, phy_g_statistics$SDReads),
  MinReads = c(phy_b_statistics$MinReads, phy_v_statistics$MinReads, phy_p_statistics$MinReads, phy_g_statistics$MinReads),
  MaxReads = c(phy_b_statistics$MaxReads, phy_v_statistics$MaxReads, phy_p_statistics$MaxReads,phy_g_statistics$MaxReads),
  ShannonDiversity = c(phy_b_statistics$ShannonDiversity, phy_v_statistics$ShannonDiversity, phy_p_statistics$ShannonDiversity, phy_g_statistics$ShannonDiversity),
  SimpsonDiversity = c(phy_b_statistics$SimpsonDiversity, phy_v_statistics$SimpsonDiversity, phy_p_statistics$SimpsonDiversity, phy_g_statistics$SimpsonDiversity),
  BrayDiversity = c(phy_b_statistics$BrayDiversity, phy_v_statistics$BrayDiversity, phy_p_statistics$BrayDiversity, phy_g_statistics$BrayDiversity),
  JaccardDiversity = c(phy_b_statistics$JaccardDiversity, phy_v_statistics$JaccardDiversity, phy_p_statistics$JaccardDiversity, phy_g_statistics$JaccardDiversity)
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
Microbiome_data_summary <- final_statistics_wide %>%
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
Microbiome_data_summary

# ============================
# ANÁLISIS COMPOSICIONAL
# ============================

#Filtrar por los 50 taxones más abundantes
phy_b_top <- get_top_taxa(phy_b_a, top_n = 50)
phy_v_top <- get_top_taxa(phy_v_a, top_n = 50)

#Composición
global_por_rango <- plot_bar_by_ranks(phy_g_a)
bacteria_composition <- plot_bar_by_ranks(phy_b_a)
virus_composition <- plot_bar_by_ranks(phy_v_a)

bacteria_composition_top <- plot_bar_by_ranks(phy_b_top)
virus_composition_top <- plot_bar_by_ranks(phy_v_top)

#Proporción de lecturas por reino
plot_reads <- plot_bar(phy_g_a, fill = "Superkingdom") +
  labs(title = paste("A")) +
  theme_minimal(base_size = 12) +
  labs(x = element_blank(), fill = "Género", y = "Porcentaje de lecturas total") +
  theme(
    plot.title = element_text(hjust = 0, size = 17, face = "bold"), 
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    axis.title.y.left = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
  ) +
  scale_y_continuous(limits = c(0, 100)) +  # Limita el eje Y hasta 100
  scale_fill_manual(values = c("#8A2BE2","#FFFF00", "#FF0000"))
  #Composición a nivel de género
 bacteria_genera <- plot_bar(phy_b_top, fill = "Genus") +
  labs(title = paste("B")) +
  theme_minimal(base_size = 12) +
  labs(x = element_blank(), fill = "Género", y = "Abundancia relativa bacteriana") +
  theme(
    plot.title = element_text(hjust = 0, size = 17, face = "bold"), 
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 0.5, size = 6),
    axis.text.y = element_text(size = 8),
    axis.title.y.left = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
  ) +
  scale_y_continuous(limits = c(0, 100)) +  # Limita el eje Y hasta 100
  scale_fill_manual(values = paleta_genus)

  #Virus
  virus_genera <- plot_bar(phy_v_top, fill = "Genus")  +
    labs(title = paste("C")) +
    theme_minimal(base_size = 12) +
    labs(x = element_blank(), fill = "Género", y = "Abundancia relativa vírica") +
    theme(
      plot.title = element_text(hjust = 0, size = 17, face = "bold"), 
      axis.text.x = element_text(angle = 70, vjust = 1, hjust = 0.5, size = 6),
      axis.text.y = element_text(size = 8), 
      axis.title.y.left = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 14),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.grid.major.y = element_line(color = "grey", size = 0.5)  # Añade la cuadrícula para el eje Y
    ) +
    scale_y_continuous(limits = c(0, 100)) +  # Limita el eje Y hasta 100
    scale_fill_manual(values = paleta_genus)

Composition_plots <-  ggarrange(plot_reads, bacteria_genera, virus_genera, 
            ncol = 1, 
            nrow = 3, 
            align = "v")
  
  
tax_df_b_top <- as.data.frame(tax_table(phy_b_top))
otu_sums_b_top <- taxa_sums(phy_b_top)
tax_df_b_top$Abundances <- otu_sums_b_top

tax_df_v_top <- as.data.frame(tax_table(phy_v_top))
otu_sums_v_top <- taxa_sums(phy_v_top)
tax_df_v_top$Abundances <- otu_sums_v_top

tax_df_p <- as.data.frame(tax_table(phy_p))
otu_sums_p <- taxa_sums(phy_p)
tax_df_p$Abundances <- otu_sums_p

write.csv(tax_df_b_top, file = "Abundancias de los 50 taxones bacterianos más abundantes.csv")
write.csv(tax_df_v_top, file = "Abundancias de los 50 taxones víricos más abundantes.csv")

  #Composición por variables
  bacteria_by_vars <- plot_bar_by_vars(phy_b_top, paleta)
  virus_by_vars <- plot_bar_by_vars(phy_v_top, paleta)

  #Análisis de abundancia diferencial en función de variables
  AD_infected <- create_MAplot_from_phy(phy_g, "Infected")
  phy_g <- subset_samples(phy_g, !is.na(Zone))
  AD_zone <- create_MAplot_from_phy(phy_g, "Zone")
  AD_sex <- create_MAplot_from_phy(phy_g, "Sex")

MA_plots <-  ggarrange(AD_infected, AD_zone, AD_sex, legend = "right",common.legend = TRUE, align = "v", ncol = 1, nrow = 3)
  
#Visualizador
plotwatcher(global_por_rango)
plotwatcher(bacteria_composition)
plotwatcher(virus_composition)
plotwatcher(bacteria_composition_top)
plotwatcher(virus_composition_top)
plotwatcher(bacteria_by_vars)
plotwatcher(virus_by_vars)

# ============================
# ANÁLISIS DE DIVERSIDAD
# ============================

################ Riqueza general ################ 

Riqueza_general <- plot_richness(phy_g, x = "Country", measures = c("Shannon", "Simpson"), title = "Diversidad alfa general") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana <- plot_richness(phy_b, x = "Country", measures = c("Shannon", "Simpson"), title = "Diversidad alfa bacteriana") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral <- plot_richness(phy_v, x = "Country", measures = c("Shannon", "Simpson"), title = "Diversidad alfa viral") + 
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

### Riqueza alfa por variables en conjunto global, bacteriana y viral ###
phy_g_filtrado <- subset_samples(phy_g, Zone != "NA")
phy_g_filtrado <- subset_samples(phy_g_filtrado, Sex != "Undefined")

Riqueza_general_sex <- plot_richness(phy_g_filtrado, x = "Country", measures = c("Shannon"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_zone <- plot_richness(phy_g_filtrado, x = "Country", measures = c("Shannon"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_infected <- plot_richness(phy_g, x = "Country", measures = c("Shannon"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_variables_SH <- wrap_plots(Riqueza_general_sex, Riqueza_general_zone, Riqueza_general_infected)

Riqueza_general_sex <- plot_richness(phy_g_filtrado, x = "Country", measures = c("Simpson"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_zone <- plot_richness(phy_g_filtrado, x = "Country", measures = c("Simpson"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin() +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_infected <- plot_richness(phy_g, x = "Country", measures = c("Simpson"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_general_variables_SS <- wrap_plots(Riqueza_general_sex, Riqueza_general_zone, Riqueza_general_infected)

phy_b_filtrado <- subset_samples(phy_b, Zone != "NA")
phy_b_filtrado <- subset_samples(phy_b_filtrado, Sex != "Undefined")

Riqueza_bacteriana_sex <- plot_richness(phy_b_filtrado, x = "Country", measures = c("Shannon"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_zone <- plot_richness(phy_b_filtrado, x = "Country", measures = c("Shannon"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_infected <- plot_richness(phy_b, x = "Country", measures = c("Shannon"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_variables_SH <- wrap_plots(Riqueza_bacteriana_sex, Riqueza_bacteriana_zone, Riqueza_bacteriana_infected)

Riqueza_bacteriana_sex <- plot_richness(phy_b_filtrado, x = "Country", measures = c("Simpson"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_zone <- plot_richness(phy_b_filtrado, x = "Country", measures = c("Simpson"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_infected <- plot_richness(phy_b, x = "Country", measures = c("Simpson"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_bacteriana_variables_SS <- wrap_plots(Riqueza_bacteriana_sex, Riqueza_bacteriana_zone, Riqueza_bacteriana_infected)

phy_v_filtrado <- subset_samples(phy_v, Zone != "NA")
phy_v_filtrado <- subset_samples(phy_v_filtrado, Sex != "Undefined")

Riqueza_viral_sex <- plot_richness(phy_v_filtrado, x = "Country", measures = c("Shannon"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_zone <- plot_richness(phy_v_filtrado, x = "Country", measures = c("Shannon"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_infected <- plot_richness(phy_v, x = "Country", measures = c("Shannon"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_variables_SH <- wrap_plots(Riqueza_viral_sex, Riqueza_viral_zone, Riqueza_viral_infected)

Riqueza_viral_sex <- plot_richness(phy_v_filtrado, x = "Country", measures = c("Simpson"), title = "Sexo") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Sex, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_zone <- plot_richness(phy_v_filtrado, x = "Country", measures = c("Simpson"), title = "Zona") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Zone, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_infected <- plot_richness(phy_v, x = "Country", measures = c("Simpson"), title = "Estado de infección") + 
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  geom_violin(alpha = 0.6) +
  facet_wrap(~ Infected, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black", face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x.bottom =  element_blank()) +
  ylab(label = "Índice alfa diversidad")

Riqueza_viral_variables_SS <- wrap_plots(Riqueza_viral_sex, Riqueza_viral_zone, Riqueza_viral_infected)

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

dist_bj <- distance(phy_b, method = "jaccard")
dist_vj <- distance(phy_v, method = "jaccard")
ordination_bj <- ordinate(phy_b, method = "PCoA", distance = dist_bj)
ordination_vj <- ordinate(phy_b, method = "PCoA", distance = dist_vj)

dist_bb <- distance(phy_b, method = "bray")
dist_vb <- distance(phy_v, method = "bray")
ordination_bb <- ordinate(phy_b, method = "PCoA", distance = dist_bb)
ordination_vb <- ordinate(phy_b, method = "PCoA", distance = dist_vb)

dist_g <- distance(phy_g, method = "jaccard")
ordination_g <- ordinate(phy_g, method = "PCoA", distance = dist_g)
plot_ordination(phy_g, ordination_g) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_manual(values =  paleta)

Riqueza_bacteriana_jaccard <- plot_ordination(phy_b, ordination_bj) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Jaccard") + 
  scale_color_manual(values =  paleta)

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

Riqueza_bacteriana_bray <- plot_ordination(phy_b, ordination_bb) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta bacteriana - Bray-Curtis") + 
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

Riqueza_viral_jaccard <- plot_ordination(phy_v, ordination_vj) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Jaccard") + 
  scale_color_manual(values =  paleta)

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

Riqueza_viral_bray <- plot_ordination(phy_v, ordination_vb) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("Diversidad beta vírica - Bray-Curtis") + 
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

Lista_beta <- list(Riqueza_bacteriana_bray, Riqueza_bacteriana_jaccard, Riqueza_viral_bray, Riqueza_viral_jaccard)
Riqueza_beta <- wrap_plots(Lista_beta)

#Test estadístico PERMANOVA para analizar qué variables producen efectos significativos sobre la diversidad
#Preparación (matriz de distancias y metadatos)
dist_bray_b <- distance(phy_b, method = "bray")
dist_bray_v <- distance(phy_v, method = "bray")
dist_jaccard_b <- distance(phy_b, method = "jaccard")
dist_jaccard_v <- distance(phy_v, method = "jaccard")
sample_df_b <- as(sample_data(phy_b), "data.frame") #Extraer metadatos
sample_df_v <- as(sample_data(phy_v), "data.frame") #Extraer metadatos

#Comprobación de homocedasticidad
dist_homo <- phyloseq::distance(phy_g, method = "bray")

sexo <- sample_data(phy_g)$Sex
zona <- sample_data(phy_g)$Zone
infected <- sample_data(phy_g)$Infected
year <- sample_data(phy_g)$Year
latitude <- sample_data(phy_g)$Latitude
longitude <- sample_data(phy_g)$Longitude

disp_sexo <- betadisper(dist_homo, sexo)
disp_zona <- betadisper(dist_homo, zona)
disp_infected <- betadisper(dist_homo, infected)
disp_year <- betadisper(dist_homo, year)
disp_latitude <- betadisper(dist_homo, latitude)
disp_longitude <- betadisper(dist_homo, longitude)

anova(disp_sexo) #Pvalor significativo --> comprobar con Tukey que grupos difieren
anova(disp_zona) #Pvalor significativo --> comprobar con Tukey que grupos difieren
anova(disp_infected) #Pvalor no significativo --> existe homocedasticidad
anova(disp_year) #Pvalor significativo --> comprobar con Tukey que grupos difieren
anova(disp_latitude) #Pvalor significativo --> comprobar con Tukey que grupos difieren
anova(disp_longitude) #Pvalor significativo --> comprobar con Tukey que grupos difieren

TukeyHSD(disp_sexo) #En efecto, no hay homocedasticidad entre machos y hembras, sin embargo esto 
#quiere decir que tienen varianzas diferentes pero ya me esperaba esto porque las hembras tienen microbiota unica.
TukeyHSD(disp_zona) #Solo ocurre entre Suroeste-Central, no influye en el análisis (no hubo resultados de ellos)
TukeyHSD(disp_infected) #Existe homocedasticidad
TukeyHSD(disp_year) #Ocurre con todas las comparaciones concernientes a 2005
TukeyHSD(disp_latitude)$group %>% #Solo ocurre en una categoría --> aceptamos homocedasticidad
  as.data.frame() %>%
  tibble::rownames_to_column("comparison") %>%
  filter(`p adj` < 0.05)
TukeyHSD(disp_longitude)$longitude %>% #Ocurre poquísimo --> aceptamos homocedasticidad
  as.data.frame() %>%
  tibble::rownames_to_column("comparison") %>%
  filter(`p adj` < 0.05)

#PERMANOVA
permanova_b_res_location <- adonis2(dist_bray_b ~ Location, data = sample_df_b) #Permanova por localidad
permanova_b_res_year <- adonis2(dist_bray_b ~ Year, data = sample_df_b) #Permanova por año
permanova_b_res_latitude <- adonis2(dist_bray_b ~ Latitude, data = sample_df_b) #Permanova por latitud
permanova_b_res_longitude <- adonis2(dist_bray_b ~ Longitude, data = sample_df_b) #Permanova por longitud
permanova_b_res_sex <- adonis2(dist_bray_b ~ Sex, data = sample_df_b) #Permanova por sexo
permanova_b_res_infected <- adonis2(dist_bray_b ~ Infected, data = sample_df_b) #Permanova por infeccion

permanova_v_res_location <- adonis2(dist_bray_v ~ Location, data = sample_df_v) #Permanova por localidad
permanova_v_res_year <- adonis2(dist_bray_v ~ Year, data = sample_df_v) #Permanova por año
permanova_v_res_latitude <- adonis2(dist_bray_v ~ Latitude, data = sample_df_v) #Permanova por latitud
permanova_v_res_longitude <- adonis2(dist_bray_v ~ Longitude, data = sample_df_v) #Permanova por longitud
permanova_v_res_sex <- adonis2(dist_bray_v ~ Sex, data = sample_df_v) #Permanova por sexo
permanova_v_res_infected <- adonis2(dist_bray_v ~ Infected, data = sample_df_v) #Permanova por infeccion

permanova_b_res_location2 <- adonis2(dist_jaccard_b ~ Location, data = sample_df_b) #Permanova por localidad
permanova_b_res_year2 <- adonis2(dist_jaccard_b ~ Year, data = sample_df_b) #Permanova por año
permanova_b_res_latitude2 <- adonis2(dist_jaccard_b ~ Latitude, data = sample_df_b) #Permanova por latitud
permanova_b_res_longitude2 <- adonis2(dist_jaccard_b ~ Longitude, data = sample_df_b) #Permanova por longitud
permanova_b_res_sex2 <- adonis2(dist_jaccard_b ~ Sex, data = sample_df_b) #Permanova por sexo
permanova_b_res_infected2 <- adonis2(dist_jaccard_b ~ Infected, data = sample_df_b) #Permanova por infeccion

permanova_v_res_location2 <- adonis2(dist_jaccard_v ~ Location, data = sample_df_v) #Permanova por localidad
permanova_v_res_year2 <- adonis2(dist_jaccard_v ~ Year, data = sample_df_v) #Permanova por año
permanova_v_res_latitude2 <- adonis2(dist_jaccard_v ~ Latitude, data = sample_df_v) #Permanova por latitud
permanova_v_res_longitude2 <- adonis2(dist_jaccard_v ~ Longitude, data = sample_df_v) #Permanova por longitud
permanova_v_res_infected2 <- adonis2(dist_jaccard_v ~ Infected, data = sample_df_v) #Permanova por infeccion

#visualización PERMANOVA (Bray)
  # Convertir resultados bacterias
  permanova_tbl_b <- bind_rows(
    Location  = as_tibble(permanova_b_res_location, rownames = "Term"),
    Year      = as_tibble(permanova_b_res_year, rownames = "Term"),
    Latitude  = as_tibble(permanova_b_res_latitude, rownames = "Term"),
    Longitude = as_tibble(permanova_b_res_longitude, rownames = "Term"),
    Sex       = as_tibble(permanova_b_res_sex, rownames = "Term"),
    Infected  = as_tibble(permanova_b_res_infected, rownames = "Term"),
    .id = "Variable"
  ) %>% mutate(Grupo = "Bacterias")
  
  # Convertir resultados virus
  permanova_tbl_v <- bind_rows(
    Location  = as_tibble(permanova_v_res_location, rownames = "Term"),
    Year      = as_tibble(permanova_v_res_year, rownames = "Term"),
    Latitude  = as_tibble(permanova_v_res_latitude, rownames = "Term"),
    Longitude = as_tibble(permanova_v_res_longitude, rownames = "Term"),
    Sex       = as_tibble(permanova_v_res_sex, rownames = "Term"),
    Infected       = as_tibble(permanova_v_res_infected, rownames = "Term"),
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
    Infected       = as_tibble(permanova_b_res_infected2, rownames = "Term"),
    .id = "Variable"
  ) %>% mutate(Grupo = "Bacterias")
  
  # Convertir resultados virus
  permanova_tbl_v2 <- bind_rows(
    Location  = as_tibble(permanova_v_res_location2, rownames = "Term"),
    Year      = as_tibble(permanova_v_res_year2, rownames = "Term"),
    Latitude  = as_tibble(permanova_v_res_latitude2, rownames = "Term"),
    Longitude = as_tibble(permanova_v_res_longitude2, rownames = "Term"),
    Sex       = as_tibble(permanova_v_res_sex2, rownames = "Term"),
    Infected       = as_tibble(permanova_v_res_infected2, rownames = "Term"),
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
# ANÁLISIS DEL PARÁSITO
# ============================  
#Metadata contiene las abundancias de parásito
total_plas_reads <- sum(metadata$Plasmodium_abundance)  
print(plas_csv) #Especies detectadas

#Correlación entre taxones y abundancia del parásito

#Correlación entre taxones y presencia/ausencia del parásito (hecho con DEseq2
#en el apartado de análisis de abundancia diferencial por estado de infección)
# -- CONFIGURACIÓN -- #
phy <- phy_g  #objeto phyloseq
variable_numerica <- "Infected"#variable dependiente

#Preprocesado
  #Extraer abundancias y metadatos
otu <- as.data.frame(otu_table(phy))
if (taxa_are_rows(phy)) otu <- t(otu)  #muestras como filas

meta <- as.data.frame(sample_data(phy))

#Verifica que la variable exista para evitar errores
if (!(variable_numerica %in% colnames(meta))) {
  stop(paste("La variable", variable_numerica, "no está en sample_data"))
}

#Extraer variable numérica
num_var <- meta[[variable_numerica]]

#Verifica que sea numérica para evitar errores
if (!is.numeric(num_var)) {
  stop(paste("La variable", variable_numerica, "no es numérica"))
}

#Obtener abundancias por fila normalizadas
otu_rel <- otu / rowSums(otu)

#Corelaciones
cor_results <- apply(otu_rel, 2, function(x) {
  suppressWarnings(cor.test(x, num_var, method = "spearman"))
})

#Ordenar resultados
res_df <- data.frame(
  OTU = names(cor_results),
  Spearman_rho = sapply(cor_results, function(x) x$estimate),
  p_value = sapply(cor_results, function(x) x$p.value),
  stringsAsFactors = FALSE
)

res_df$padj <- p.adjust(res_df$p_value, method = "BH")

#Añadir etiquetas taxonómicas
if (!is.null(tax_table(phy))) {
  tax <- as.data.frame(tax_table(phy))
  res_df <- cbind(res_df, tax[match(res_df$OTU, rownames(tax)), ])
}
#Exportar y guardar como objeto
write.csv(res_df, file = "correlacion_OTUs_sexo.csv", row.names = FALSE)
correlacion_OTUs_vs_parasito <- read.csv("correlacion_OTUs_vs_parasito.csv")

#Correlación entre variables y abundancia del parásito
  #Variables cuantitativas
correlacion_plas_vs_año <- cor.test(metadata$Plasmodium_abundance, metadata$Year, method = "spearman")
  #Variables categóricas (K < 2)
kruskal.test(Plasmodium_abundance ~ Zone, data = metadata)
correlacion_plas_vs_zona <- pairwise.wilcox.test(metadata$Plasmodium_abundance, metadata$Zone, p.adjust.method = "fdr")

kruskal.test(Plasmodium_abundance ~ Sex, data = metadata)
correlacion_plas_vs_sexo <- pairwise.wilcox.test(metadata$Plasmodium_abundance, metadata$Sex, p.adjust.method = "fdr")

# ============================
# RESULTADOS DEL SCRIPT
# ============================

### Análisis preliminar ###
sample_per_location #Histograma para conocer la distribución de muestras
Microbiome_data_summary #Tabla con parámetros de los datos

### Composición del microbioma ###
Composition_figure #Lecturas por reino + 50 géneros más abundantes de bacteria y virus

MA_plots #Taxones difernecialmente abundantes entre variables categóricas

### Análisis de diversidad ###
Riqueza_general
Riqueza_bacteriana #Gráfico de diversidad alfa bacteriana general
Riqueza_viral      #Gráfico de diversidad alfa viral general
wrap_plots(Riqueza_general, Riqueza_bacteriana, Riqueza_viral) #Juntos
Riqueza_beta       #Gráfico de diversidad beta general

#Riqueza por variables (relevantes)
  #Alfa
Riqueza_general_variables_SH
Riqueza_bacteriana_variables_SH
Riqueza_viral_variables_SH

Riqueza_general_variables_SS
Riqueza_bacteriana_variables_SS
Riqueza_viral_variables_SS

riqueza_b_location #Cada variable sola
riqueza_b_zona
riqueza_b_year
riqueza_b_sex
riqueza_v_location
riqueza_v_zona
riqueza_v_year
riqueza_v_sex

  #Beta
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

#Análisis del parásito
plas_csv #Datos para el nº de infectados y abundancia
correlacion_OTUs_vs_parasito 
correlacion_plas_vs_año
correlacion_plas_vs_zona
correlacion_plas_vs_sexo
#Otros resultados 
  #Datos numéricos de las abundancias de los gráficos
  df_b_genus <- tax_df_b_top %>%
    group_by(Genus) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
  df_b_family <- tax_df_b_top %>%
    group_by(Family) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
  df_b_species <- tax_df_b_top %>%
    group_by(Species) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
  
  df_v_genus <- tax_df_v_top %>%
    group_by(Genus) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
  df_v_family <- tax_df_v_top %>%
    group_by(Family) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
  df_v_species <- tax_df_v_top %>%
    group_by(Species) %>%
    summarise(Abundances = sum(Abundances, na.rm = TRUE)) %>%
    arrange(desc(Abundances))
