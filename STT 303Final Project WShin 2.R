library(tidyverse)
library(dplyr)
#loading Metadata
metadata <- read_tsv("https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
head(metadata)

#loading Gene Expression
getwd()
setwd("~/Desktop/School/2025 STT303")
getwd()
list.files()
gct_file <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"

gene_expression <- read_tsv(
  gct_file,
  skip = 2,
  col_types = cols(
    Name        = col_character(),
    Description = col_character(),
    .default    = col_double()
  )
)
head(gene_expression)

#clean metadata

selected_tissues <- c("Lung", "Liver", "Heart - Left Ventricle")
clean_metadata <- metadata %>%
  select(SAMPID, SMTSD) %>%    
  filter(!is.na(SMTSD)) %>%    
  filter(SMTSD %in% selected_tissues) %>%  
  distinct(SAMPID, .keep_all = TRUE) 

clean_metadata %>% 
  count(SMTSD)

#clean gene expression
sample_ids <- clean_metadata$SAMPID

gene_expression_cols <- colnames(gene_expression)
common_ids <- intersect(sample_ids, gene_expression_cols)
missing_ids <- setdiff(sample_ids, gene_expression_cols)

common_ids
missing_ids
length(missing_ids) 
head(missing_ids, 10) 

gene_expression_sub <- gene_expression %>%
  select(Name, all_of(common_ids))

#rows and columns
gene_expression_mat <- gene_expression_sub %>% 
  column_to_rownames(var = "Name") %>% 
  as.matrix()

gene_expression_t   <- t(gene_expression_mat)

gene_expression_df  <- gene_expression_t %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SAMPID")

#merge
merged_df <- clean_metadata %>% 
  select(SAMPID, SMTSD) %>%  
  left_join(gene_expression_df, by = "SAMPID")

#clean merge data
merged_df_clean <- merged_df[complete.cases(merged_df), ]
anyNA(merged_df_clean) 
table(merged_df_clean$SMTSD)

gene_cols   <- setdiff(names(merged_df_clean), c("SAMPID","SMTSD"))
expr_mat    <- as.matrix( merged_df_clean[gene_cols] )
avg_expr    <- rowMeans(expr_mat, na.rm = TRUE)
total_expr  <- rowSums(expr_mat, na.rm = TRUE)

#calculate total and average expression
sample_stats <- merged_df_clean %>%
  select(SAMPID, SMTSD) %>%
  mutate(
    AverageExpr = avg_expr,
    TotalExpr   = total_expr
  )

#visualize total and average expression by tissue
ggplot(sample_stats, aes(x = SMTSD, y = AverageExpr, fill = SMTSD)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Average Gene Expression by Tissue",
    x     = "Tissue",
    y     = "Mean TPM"
  )

ggplot(sample_stats, aes(x = SMTSD, y = TotalExpr, fill = SMTSD)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Total Gene Expression by Tissue",
    x     = "Tissue",
    y     = "Mean TPM"
  )

tissue_summary <- sample_stats %>%
  group_by(SMTSD) %>%
  summarise(
    N             = n(),
    Mean_AvgExpr  = mean(AverageExpr),
    Median_AvgExpr= median(AverageExpr),
    Mean_TotalExpr= mean(TotalExpr),
    Median_TotalExpr= median(TotalExpr)
  )

tissue_summary

#average each gene
gene_means <- merged_df_clean %>%
  pivot_longer(
    cols      = all_of(gene_cols),
    names_to  = "Gene",
    values_to = "Expression"
  ) %>%
  group_by(SMTSD, Gene) %>%
  summarise(
    MeanExpr = mean(Expression, na.rm = TRUE),
    .groups  = "drop"
  )

top10 <- gene_means %>%
  group_by(Gene) %>%
  summarise(overall_var = var(MeanExpr)) %>%
  arrange(desc(overall_var)) %>%
  slice_head(n = 10) %>%
  pull(Gene)

top50 <- gene_means %>%
  group_by(Gene) %>%
  summarise(overall_var = var(MeanExpr)) %>%
  arrange(desc(overall_var)) %>%
  slice_head(n = 50) %>%
  pull(Gene)

#top 50 graph
gene_means %>%
  filter(Gene %in% top50) %>%
  ggplot(aes(x = Gene, y = MeanExpr, fill = SMTSD)) +
  geom_col(position = "dodge") +
  labs(
    title = "Average Expression of Top 50 Variable Genes by Tissue",
    x     = "Gene",
    y     = "Mean TPM"
  )

#top 10 graph
gene_means %>%
  filter(Gene %in% top10) %>%
  ggplot(aes(x = Gene, y = MeanExpr, fill = SMTSD)) +
  geom_col(position = "dodge") +
  labs(
    title = "Average Expression of Top 10 Variable Genes by Tissue",
    x     = "Gene",
    y     = "Mean TPM"
  )

#matrix of mean 
heat_mat <- gene_means %>%
  filter(Gene %in% top50) %>%
  pivot_wider(
    names_from  = SMTSD,
    values_from = MeanExpr
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

heat_df <- as.data.frame(heat_mat) %>%
  rownames_to_column("Gene") %>%       # gene names as a column
  pivot_longer(
    cols      = -Gene,                # all columns except "Gene"
    names_to  = "Tissue",             # column names become values in "Tissue"
    values_to = "MeanExpr"            # cell values go into "MeanExpr"
  )

ggplot(heat_df, aes(x = Tissue, y = Gene, fill = MeanExpr)) +
  geom_tile() +
  scale_fill_gradient(
    low  = "white",      # color for low values
    high = "firebrick"   # color for high values
  )

#PCA
pca_res <- prcomp(
  merged_df_clean[gene_cols])

pca_res

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
pc1_pct <- round(var_explained[1], 1)
pc2_pct <- round(var_explained[2], 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  setNames(c("PC1", "PC2")) %>%
  mutate(Tissue = merged_df_clean$SMTSD)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "PCA of GTEx Tissue Samples",
    x     = paste0("PC1 (", pc1_pct, "%)"),
    y     = paste0("PC2 (", pc2_pct, "%)")
  )

loadings <- as.data.frame(pca_res$rotation)

top_pc1 <- loadings %>%
  select(PC1) %>%
  mutate(Gene = rownames(.)) %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n = 10)
top_pc1

view(top10)

#building a model
#splitting data to test and sample data
set.seed(2025)

model_df <- merged_df_clean %>%
  select(SAMPID, SMTSD, all_of(top50)) %>%
  mutate(SMTSD = factor(SMTSD))

train_df <- model_df %>%
  group_by(SMTSD) %>%
  slice_sample(prop = 0.8)

test_df <- anti_join(model_df, train_df, by = "SAMPID")

train_df %>% count(SMTSD) %>% rename(Train = n)
test_df  %>% count(SMTSD) %>% rename(Test  = n)

#training
train_mat    <- as.matrix(train_df[ , top50])
train_mean   <- colMeans(train_mat)
train_sd     <- apply(train_mat, 2, sd)

#scaling
train_scaled_mat <- scale(train_mat, 
                          center = train_mean, 
                          scale  = train_sd)
train_scaled_df  <- as.data.frame(train_scaled_mat) %>%
  mutate(SMTSD = train_df$SMTSD)

test_mat        <- as.matrix(test_df[ , top50])
test_scaled_mat <- scale(test_mat, 
                         center = train_mean, 
                         scale  = train_sd)
test_scaled_df  <- as.data.frame(test_scaled_mat) %>%
  mutate(SMTSD = test_df$SMTSD)

#logistic model
install.packages("nnet")
library(nnet)
logistic_model <- multinom(
  SMTSD ~ .,
  data  = train_scaled_df,
  trace = FALSE    # suppress iteration log
)
summary(logistic_model)

#predicting
pred   <- predict(logistic_model, newdata = test_scaled_df)
pred_result <- table(
  Actual    = test_scaled_df$SMTSD,
  Predicted = pred
)
print(pred_result)

#accuracy
accuracy <- sum(diag(pred_result)) / sum(pred_result)
cat("Overall accuracy:", round(accuracy * 100, 1), "%\n")

#regression
reg_df <- merged_df_clean %>%
  mutate(AverageExpr = rowMeans(select(., all_of(gene_cols))))
lm_avg <- lm(AverageExpr ~ SMTSD, data = reg_df)
summary(lm_avg)

#clustering
set.seed(2025)
km <- kmeans(train_scaled_mat, centers = 3, nstart = 25)
table(Cluster = km$cluster, Tissue = train_scaled_df$SMTSD)


#precision f1
precision <- diag(pred_result) / colSums(pred_result)
recall    <- diag(pred_result) / rowSums(pred_result)
f1        <- 2 * precision * recall / (precision + recall)

class_metrics <- data.frame(
  Tissue    = names(precision),
  Precision = round(precision, 3),
  Recall    = round(recall,    3),
  F1        = round(f1,        3)
)

print(class_metrics)
