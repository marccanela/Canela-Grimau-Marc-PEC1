# Data exploration script

# Packages
library(dplyr)
library(ggplot2)

# Data import
setwd(dir = "~/Desktop/Canela-Grimau-Marc-PEC1")
load(file = "data_and_metadata.rda")

# Filter numeric variables
numeric_vars <- DataInfo_S013 %>%
  filter(varTpe %in% c("numeric", "integer")) %>%
  pull(VarName)
DataValues_S013_numeric <- DataValues_S013 %>%
  select(all_of(numeric_vars))

# =============== PCA =====================

time_vars <- numeric_vars[grepl("_T[0-9]+$", numeric_vars)]
data_pca <- DataValues_S013_numeric %>%
  select(all_of(c(time_vars, "Group")))
data_pca_clean <- data_pca %>%
  select(where(~ !any(is.na(.)) && sd(.) != 0))
pca_result <- prcomp(data_pca_clean[, -ncol(data_pca_clean)], scale=TRUE)

summary(pca_result)
plot(pca_result)

loading_scores <- pca_result$rotation

top_pc1 <- loading_scores %>%
  as.data.frame() %>%
  mutate(variable = rownames(loading_scores)) %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:5)
top_pc1['PC1']

top_pc2 <- loading_scores %>%
  as.data.frame() %>%
  mutate(variable = rownames(loading_scores)) %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:5)
top_pc2['PC2']

pc1 <- sub("_T0$", "", row.names(top_pc1['PC1']))
pc1 <- gsub("\\.", " ", pc1)
pc1 <- sub(" (?!.* )", ":", pc1, perl = TRUE)
pc1_info <- AAInformation_S006 %>%
  filter(Metabolite.abbreviation %in% pc1)
pc1_info

pc2 <- sub("_T0$", "", row.names(top_pc2['PC2']))
pc2 <- gsub("\\.", " ", pc2)
pc2 <- sub(" (?!.* )", ":", pc2, perl = TRUE)
pc2_info <- AAInformation_S006 %>%
  filter(Metabolite.abbreviation %in% pc2)
pc2_info

plot_pca_data <- as.data.frame(pca_result$x)
plot_pca_data$group <- as.character(data_pca_clean$Group)

pca_plot <- ggplot(plot_pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = "PCA of Metabolomic Data",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  theme(legend.position = "top")
pca_plot

var.test(plot_pca_data$PC1[plot_pca_data$group == '1'],
         plot_pca_data$PC1[plot_pca_data$group == '2'])
t.test(PC1 ~ group, data = plot_pca_data, var.equal = TRUE)

plot_pca_data$group <- as.character(DataValues_S013$GENDER)

pca_plot <- ggplot(plot_pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = "PCA of Metabolomic Data",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  theme(legend.position = "top")
pca_plot

var.test(plot_pca_data$PC1[plot_pca_data$group == 'F'],
         plot_pca_data$PC1[plot_pca_data$group == 'M'])
t.test(PC1 ~ group, data = plot_pca_data, var.equal = TRUE)

plot_pca_data$group <- as.character(DataValues_S013$SURGERY)

pca_plot <- ggplot(plot_pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(title = "PCA of Metabolomic Data",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  theme(legend.position = "top")
pca_plot

var.test(plot_pca_data$PC1[plot_pca_data$group == 'by pass'],
         plot_pca_data$PC1[plot_pca_data$group == 'tubular'])
t.test(PC1 ~ group, data = plot_pca_data, var.equal = TRUE)

# =============== Dendrogram =====================

data_pca_clean$Group <- NULL
hierarchical_data <- scale(data_pca_clean)
matdis <- dist(hierarchical_data)
hierarchical_data.hc <- hclust(matdis, method="average")
plot(hierarchical_data.hc)
cor(matdis, cophenetic(hierarchical_data.hc))

