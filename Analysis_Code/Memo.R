# Setting
install.packages("readxl")
install.packages("tidyverse")
install.packages("dplyr")


# qq-Plots
library(ggplot2)
set.seed(123) # 再現性を確保するためにシードを設定

# 正規分布に従う理想的なデータを生成
normal_data <- data.frame(value = rnorm(100, mean = 0, sd = 1))

# Q-Qプロットを作成
ggplot(normal_data, aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Ideal Q-Q Plot (Normal Distribution)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()

# 偏った（歪んだ）データを生成（例：指数分布）
skewed_data <- data.frame(value = rexp(100))

# Q-Qプロットを作成
ggplot(skewed_data, aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Non-Normal Q-Q Plot (Skewed Distribution)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()
