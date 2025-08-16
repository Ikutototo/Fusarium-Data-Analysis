# Setting ---------------------------------
setwd("/Users/ikutosaito/Documents/Positron/Fusarium-R-Analysis/Figure")

rm(list = ls()) # Clear workspace

# library
library(dplyr)
library(purrr)
library(readxl)
library(ggplot2)


# Antimicrobial Susceptibility Testing (AST)_Data
# Data imports Antimicrobial Susceptibility Testing (AST)
Data <- readxl::read_excel(
    # path = "~/Documents/Positron/Fusarium-R-Analysis/DiseaseIndex_v2.xlsx",
    path = "~/Documents/Positron/Fusarium-R-Analysis/Input_Data/250816_Data.xlsx",
    col_names = TRUE
)

colnames(Data)


# Disease-Index ----------------------------------------------------------
Data |>
    filter(date == "250724") |>
    group_by(Treatments, TreatsConc) |>
    summarise(
        Sum_DiseaseIndex = sum(DiseaseIndex, na.rm = TRUE),
        Parameter = n() * 4,
        n = n(),
        .groups = "drop"
    ) |>
    mutate(
        Incidence = (Sum_DiseaseIndex / Parameter) * 100
    ) |>
    ggplot(aes(x = Treatments, y = Incidence, fill = factor(TreatsConc))) + # TreatsConcで色分け
    geom_bar(width = 0.5, stat = "identity", position = "dodge", alpha = 0.7) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(
            size = 20,
            angle = 20,
            hjust = 1,
            colour = "black"
        ),
        axis.text.y = element_text(size = 20, hjust = 1, colour = "black")
    )


# FreshWeight
FrechWeight_Plot <-
    Data |>
    filter(date == "250724") |>
    group_by(Treatments, TreatsConc, Spores_Order) |>
    summarise(
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
    ) |>
    ggplot(aes(
        x = Treatments,
        y = FreshWeight_mean,
        fill = factor(TreatsConc)
    )) +
    geom_bar(width = 0.5, stat = "identity", position = "dodge", alpha = 0.7) +
    # 有意性マークを追加
    geom_text(
        data = ttest_results,
        aes(
            x = Treatments,
            y = y_position,
            label = significance,
            group = TreatsConc
        ),
        position = position_dodge(width = 0.5),
        vjust = 0,
        size = 6,
        color = "red",
        fontface = "bold"
    ) +
    geom_errorbar(
        aes(
            ymin = FreshWeight_mean - FreshWeight_sd,
            ymax = FreshWeight_mean + FreshWeight_sd
        ),
        position = position_dodge(width = 0.5),
        width = 0.3,
        linewidth = 0.4,
        linetype = 1,
        alpha = 0.8,
        color = "black"
    ) +
    xlab("") +
    ylab("") +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "plain", color = "black"),
        strip.text = element_text(size = 20, face = "bold", color = "black"),
        strip.background = element_rect(fill = "lightgray", color = "gray50"),
        axis.title.x = element_text(
            size = 23,
            colour = "black",
            face = "bold",
            vjust = -0.5
        ),
        axis.title.y = element_text(size = 23, colour = "black", face = "bold"),
        axis.text.x = element_text(
            size = 20,
            color = "black",
            angle = 90,
            face = "bold"
        ),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90")
    ) +
    geom_jitter(
        data = Data |> filter(date == "250724"),
        aes(
            x = Treatments,
            y = FreshWeight,
            group = TreatsConc
        ),
        size = 2.3,
        alpha = 0.65,
        show.legend = FALSE,
        color = "#666666",
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5) # geom_bar position = dodgeと合わせること
    ) +
    scale_fill_manual(
        values = c(
            "0.2" = "#466ee5e2",
            "0.5" = "#4bd22aff",
            "1" = "#f681ddff",
            "0" = "#2acad2ce"
        ),
        labels = c(
            "No Treatment",
            "ImE 0.2(%)",
            "ImE 0.5(%)",
            "ImE 1.0(%)"
        )
    )


Control_Data <- Data |>
    filter(date == "250724", Treatments == "NC", TreatsConc == 0)

summary_stats <- Data |>
    filter(date == "250724") |>
    group_by(Treatments, TreatsConc, Spores_Order) |>
    summarise(
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        Parameter = n() * 4,
        .groups = "drop"
    )

ttest_results <- Data |>
    filter(date == 250724 & !(Treatments == "NC" & TreatsConc == 0)) |>
    group_by(Treatments, TreatsConc) |>
    summarise(
        ttest_result = list(t.test(
            FreshWeight,
            Control_Data$FreshWeight
        )),
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        Parameter = n() * 4,
        .groups = "drop"
    ) |>
    mutate(
        p_value = map_dbl(ttest_result, ~ .x$p.value),
        t_statistic = map_dbl(ttest_result, ~ .x$statistic),
        significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            p_value < 0.1 ~ ".",
            TRUE ~ "ns"
        ),
        y_position = FreshWeight_mean + FreshWeight_sd + 1
    )


# 250816 -----------------------------------------------------------------
Control_Data <- Data |>
    filter(Treatments == "NC")

summary_stats <- Data |>
    group_by(Treatments, TreatsConc, Spores_Order) |>
    summarise(
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        Parameter = n() * 4,
        .groups = "drop"
    )

# 正規性検定
NormalDistribution <- Data |>
    group_by(Treatments, Spores_Order) |>
    summarise(
        shapiro_p = shapiro.test(FreshWeight)$p.value,
        n = n(),
        .groups = "drop"
    )


ttest_results <- Data |>
    filter(Treatments != "NC") |>
    group_by(Treatments, TreatsConc, Spores_Order) |>
    summarise(
        ttest_result = list(t.test(
            FreshWeight,
            Control_Data$FreshWeight
        )),
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        Parameter = n() * 4,
        .groups = "drop"
    ) |>
    mutate(
        p_value = map_dbl(ttest_result, ~ .x$p.value),
        t_statistic = map_dbl(ttest_result, ~ .x$statistic),
        significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            p_value < 0.1 ~ ".",
            TRUE ~ "ns"
        ),
        y_position = FreshWeight_mean + FreshWeight_sd + 1 # 1の値は調整用
    )

# FreshWeight
FrechWeight_Plot <-
    Data |>
    group_by(Treatments, TreatsConc, Spores_Order) |>
    summarise(
        FreshWeight_mean = mean(FreshWeight, na.rm = TRUE),
        FreshWeight_sd = sd(FreshWeight, na.rm = TRUE),
        FreshWeight_se = sd(FreshWeight, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
    ) |>
    ggplot(aes(
        x = factor(Treatments, levels = c("NC", "PC", "ImE")),
        y = FreshWeight_mean,
        fill = factor(Spores_Order)
    )) +
    geom_bar(width = 0.5, stat = "identity", position = "dodge", alpha = 0.7) +
    # 有意性マークを追加
    geom_text(
        data = ttest_results,
        aes(
            x = Treatments,
            y = y_position,
            label = significance,
            group = Spores_Order
        ),
        position = position_dodge(width = 0.5),
        vjust = 0,
        size = 6,
        color = "red",
        fontface = "bold",
        show.legend = FALSE,
    ) +
    geom_errorbar(
        aes(
            ymin = pmax(0, FreshWeight_mean - FreshWeight_sd),
            ymax = FreshWeight_mean + FreshWeight_sd
        ),
        position = position_dodge(width = 0.5), # geom_barと合わせること
        width = 0.13,
        linewidth = 0.5,
        linetype = 1,
        alpha = 0.7,
        show.legend = FALSE,
        color = "black"
    ) +
    xlab("") +
    ylab("") +
    scale_y_continuous(
        breaks = c(0, 3, 6, 9, 12, 15, 18),
        labels = c(
            "0g",
            "3g",
            "6g",
            "9g",
            "12g",
            "15g",
            "18g"
        )
    ) +
    theme(
        legend.position = c(0.5, 0.80), # 位置を調整
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "plain", color = "black"),
        strip.text = element_text(size = 20, face = "bold", color = "black"),
        strip.background = element_rect(fill = "lightgray", color = "gray50"),
        axis.title.x = element_text(
            size = 23,
            colour = "black",
            face = "bold",
            vjust = -0.5
        ),
        axis.title.y = element_text(size = 23, colour = "black", face = "bold"),
        # axis.text.x = element_text(
        #     size = 20,
        #     color = "black",
        #     angle = 90,
        #     face = "bold",
        #     vjust = 0.5
        # ),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30, color = "black", face = "bold"),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90")
    ) +
    geom_jitter(
        data = Data,
        aes(x = Treatments, y = FreshWeight),
        size = 2.3,
        alpha = 0.65,
        show.legend = FALSE,
        color = "#666666",
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5) # geom_bar position = dodgeと合わせること
    ) +
    guides(
        fill = guide_legend(
            override.aes = list(
                shape = 22,
                size = 6
            )
        )
    ) +
    scale_fill_manual(
        values = c(
            "6" = "#325A9B",
            "0" = "#1C8356"
        ),
        labels = c(
            "No Spores",
            expression(1 %*% 10^6 ~ "Spores/mL")
        )
    )
