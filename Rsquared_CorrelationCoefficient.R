

# Load data
data("mtcars")
df <- mtcars
df$cyl <- as.factor(df$cyl)


sp <- ggscatter(cor, x = "Prot.logFC", y = "RNA.logFC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0, label.y = 5)


ggscatter(cor, x = "Prot.logFC", y = "RNA.logFC", add = "reg.line") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0, label.y = 5
  )


# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::
sp <- ggscatter(df, x = "wt", y = "mpg",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 30)

# Use R2 instead of R
ggscatter(df, x = "wt", y = "mpg", add = "reg.line") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3
  )


summary(lm(cor$Prot.logFC~cor$RNA.logFC))$adj.r.squared

cor.test(cor$Prot.logFC,cor$RNA.logFC, method = "pearson")







