if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")

library(ggplot2)
library(dplyr)
library(reshape2)

set.seed(1)  # reproducible random numbers
genes <- c("CDH1","EPCAM","VIM","ZEB1","SNAI1","MMP9")   # mix epithelial & mesenchymal genes
samples <- paste0("S", 1:6)                             # S1..S6 samples
condition <- c("Normal","Normal","Normal","Cancer","Cancer","Cancer")  # sample labels
genes
samples
condition
expr_mat <- matrix(nrow = length(genes), ncol = length(samples))
expr_mat
rownames(expr_mat) <- genes
colnames(expr_mat) <- samples

for (g in 1:length(genes)) {
  gene_name <- genes[g]
  
  if (gene_name %in% c("CDH1","EPCAM")) {
    expr_mat[g, 1:3] <- round(rnorm(3, mean = 8, sd = 0.4), 1)   # Normal: higher
    expr_mat[g, 4:6] <- round(rnorm(3, mean = 5, sd = 0.4), 1)   # Cancer: lower
  } else {
    expr_mat[g, 1:3] <- round(rnorm(3, mean = 5, sd = 0.4), 1)   # Normal: lower
    expr_mat[g, 4:6] <- round(rnorm(3, mean = 8, sd = 0.4), 1)   # Cancer: higher
  }
}

# View the matrix
expr_mat



# Convert to data frame for easy handling
expr_df <- as.data.frame(expr_mat)
expr_df$Gene <- rownames(expr_df)   # add gene column  
expr_df
expr_long <- melt(expr_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
expr_long$Condition <- rep(condition, each = length(genes))
expr_long
# 4. Compute summary statistics per Gene per Condition (mean ± sd)
summary_table <- expr_long %>%
  group_by(Gene, Condition) %>%
  summarise(mean_expr = mean(Expression), sd_expr = sd(Expression), n = n(), .groups = "drop")

summary_table

# 6. Plot: mean expression per gene, side-by-side (Normal vs Cancer) with error bars (sd)
p <- ggplot(summary_table, aes(x = Gene, y = mean_expr, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
                position = position_dodge(width = 0.7), width = 0.2) +
  labs(title = "Mean gene expression (Normal vs Cancer)",
       x = "Gene", y = "Mean expression (a.u.)") +
  theme_minimal() +
  scale_fill_manual(values = c("Normal" = "#1f77b4", "Cancer" = "#d62728"))

print(p)
ggsave("gene_expression_summary_plot.png", plot = p, width = 8, height = 4)
write.csv(summary_table, "gene_expression_summary_table.csv", row.names = FALSE)
