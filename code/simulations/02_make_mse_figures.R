# ============================================================
# Figures: Monte Carlo out-of-sample MSE
# Paper: Forecasting Dynamic Panel Models with Shrinkage of Fixed Effects
# Authors: Magdalena Cornejo and Walter Sosa-Escudero
#
# Input:  output/tables/table3_all_scenarios.csv
# Output: output/figures/*.pdf and *.eps
# ============================================================

library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)

# -----------------------------
# 1) Paths and data
# -----------------------------
input_file <- file.path("data","processed","table3_all_scenarios.csv")
output_dir <- file.path("output", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

full_data <- read.csv(input_file)

# -----------------------------
# 2) Plot settings
# -----------------------------
colours_models <- c(
  "OLS"         = "#e31a1c",  # Rojo fuerte
  "LSDV"        = "#fb9a99",  # Rojo claro (LSDV)
  "AH"          = "#00441b",  # Verde oscuro
  "GMM1"        = "#238b45",  # Verde medio
  "GMM2"        = "#74c476",  # Verde claro
  "EBMLE"       = "#08306b",  # Azul Negro (Muy oscuro)
  "URE"         = "#08519c",  # Azul Real
  "LASSO"       = "#2171b5",  # Azul Denim
  "RIDGE"       = "#6baed6",  # Azul Cielo Fuerte
  "ELASTIC_NET" = "#9ecae1"   # Celeste 
)

linetypes_grupos <- c(
  "Shrinkage-type" = "solid", 
  "OLS-type"  = "dashed", 
  "IV-type"   = "dotted"
)

niveles_modelos <- names(colours_models)

# -----------------------------
# 3) Main MSE figure function
# -----------------------------
crear_plot_mse <- function(df_input, n_val, t_val, g_val, titulo_expr) {
  
  data_plot <- df_input %>%
    filter(N == n_val, T0 == t_val, gamma == g_val, Measure == "MSE") %>%
    pivot_longer(cols = all_of(c("OLS", "LSDV", "AH", "GMM1", "GMM2", "EBMLE", "URE", "LASSO", "RIDGE", "ELASTIC_NET")),
                 names_to = "Model", values_to = "MSE") %>%
    mutate(
      Group = case_when(
        Model %in% c("EBMLE", "URE", "LASSO", "RIDGE", "ELASTIC_NET") ~ "Shrinkage-type",
        Model %in% c("OLS", "LSDV") ~ "OLS-type",
        Model %in% c("AH", "GMM1", "GMM2") ~ "IV-type"
      ),
      Model = factor(Model, levels = niveles_modelos),
      Group = factor(Group, levels = names(linetypes_grupos))
    )
  
  graf <- ggplot(data_plot, aes(x = sparse, y = MSE, 
                                color = Model, linetype = Group, group = Model)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = colours_models) +
    scale_linetype_manual(values = linetypes_grupos) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")) +
    labs(title = titulo_expr, y = "MSE", x = expression(paste("Sparsity fraction of ", gamma))) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  legend_text <- ggplot() +
    theme_void() +
    # OLS-type
    annotate("text", x = 0.05, y = 10, label = "OLS-type (dashed)", fontface = "bold", hjust = 0, size = 3.8) +
    annotate("text", x = 0.05, y = 9.4, label = "  ├─ OLS", color = colours_models["OLS"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 8.8, label = "  └─ LSDV",  color = colours_models["LSDV"], hjust = 0, size = 3.4) +
    
    # IV-type
    annotate("text", x = 0.05, y = 7.8, label = "IV-type (dotted)", fontface = "bold", hjust = 0, size = 3.8) +
    annotate("text", x = 0.05, y = 7.2, label = "  ├─ AH",   color = colours_models["AH"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 6.6, label = "  ├─ GMM1", color = colours_models["GMM1"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 6.0, label = "  └─ GMM2", color = colours_models["GMM2"], hjust = 0, size = 3.4) +
    
    # Shrinkage-type
    annotate("text", x = 0.05, y = 5.0, label = "Shrinkage-type (solid)", fontface = "bold", hjust = 0, size = 3.8) +
    annotate("text", x = 0.05, y = 4.4, label = "  ├─ EBMLE", color = colours_models["EBMLE"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 3.8, label = "  ├─ URE",   color = colours_models["URE"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 3.2, label = "  ├─ LASSO", color = colours_models["LASSO"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 2.6, label = "  ├─ RIDGE", color = colours_models["RIDGE"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 2.0, label = "  └─ ENET",  color = colours_models["ELASTIC_NET"], hjust = 0, size = 3.4) +
    
    coord_cartesian(xlim = c(0, 1), ylim = c(1.5, 10.5), clip = "off")
  plot_grid(graf, ggplotGrob(legend_text), ncol = 2, rel_widths = c(3.5, 1.5), align = "h")
  }

# -----------------------------
# 4) Main MSE panels
# Change T0/T*/H here when producing alternative figures.
# ----------------------------- 
p1 <- crear_plot_mse(full_data, 100, 10, 0.2, expression(gamma == 0.2 ~ ", N=100, T*" == 8 ~ ","~ H == 2))
p2 <- crear_plot_mse(full_data, 100, 10, 0.8, expression(gamma == 0.8 ~ ", N=100, T*" == 8 ~ ","~ H == 2))
p3 <- crear_plot_mse(full_data, 20, 10, 0.2,  expression(gamma == 0.2 ~ ", N=20, T*" == 8 ~ ","~ H == 2))
p4 <- crear_plot_mse(full_data, 20, 10, 0.8,  expression(gamma == 0.8 ~ ", N=20, T*" == 8 ~ ","~ H == 2))

final_grid <- plot_grid(p1, p2, p3, p4, ncol = 2)
final_grid

# -----------------------------
# 5) Save main MSE figure
# -----------------------------
ggsave(file.path(output_dir, "mse_T10.pdf"), final_grid, width = 14, height = 8)
ggsave(file.path(output_dir, "mse_T10.eps"), final_grid, device = "eps", width = 14, height = 8)

# -----------------------------
# 6) MSE figure: shrinkage estimators only
# -----------------------------
crear_plot_shrinkage <- function(df_input, n_val, t_val, g_val, titulo_expr) {
  
  modelos_shrinkage <- c("EBMLE", "URE", "LASSO", "RIDGE", "ELASTIC_NET")
  
  data_plot <- df_input %>%
    filter(N == n_val, T0 == t_val, gamma == g_val, Measure == "MSE") %>%
    pivot_longer(cols = all_of(modelos_shrinkage),
                 names_to = "Model", values_to = "MSE") %>%
    mutate(Model = factor(Model, levels = niveles_modelos))
  
  graf <- ggplot(data_plot, aes(x = sparse, y = MSE, color = Model, group = Model)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.8) + 
    scale_color_manual(values = colours_models) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
    labs(title = titulo_expr, y = "MSE", x = expression(paste("Sparsity fraction of ", gamma))) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9)
    )
  legend_text <- ggplot() +
    theme_void() +
    annotate("text", x = 0.05, y = 6, label = "Shrinkage-type", fontface = "bold", hjust = 0, size = 3.8) +
    annotate("text", x = 0.05, y = 5, label = "  ├─ EBMLE", color = colours_models["EBMLE"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 4, label = "  ├─ URE",   color = colours_models["URE"],   hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 3, label = "  ├─ LASSO", color = colours_models["LASSO"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 2, label = "  ├─ RIDGE", color = colours_models["RIDGE"], hjust = 0, size = 3.4) +
    annotate("text", x = 0.05, y = 1, label = "  └─ ENET",  color = colours_models["ELASTIC_NET"], hjust = 0, size = 3.4) +
    coord_cartesian(xlim = c(0, 1.2), ylim = c(0.5, 6.5), clip = "off")
  
  plot_grid(graf, ggplotGrob(legend_text), ncol = 2, rel_widths = c(3.5, 1.5), align = "h")
}

ps1 <- crear_plot_shrinkage(full_data, 100, 5, 0.2, expression(gamma == 0.2 ~ "," ~ N == 100 ~ "," ~ T^"*" == 4 ~ "," ~ H == 1))
ps2 <- crear_plot_shrinkage(full_data, 100, 5, 0.8, expression(gamma == 0.8 ~ "," ~ N == 100 ~ "," ~ T^"*" == 4 ~ "," ~ H == 1))
ps3 <- crear_plot_shrinkage(full_data, 20, 5, 0.2, expression(gamma == 0.2 ~ "," ~ N == 20 ~ "," ~ T^"*" == 4 ~ "," ~ H == 1))
ps4 <- crear_plot_shrinkage(full_data, 20, 5, 0.8, expression(gamma == 0.8 ~ "," ~ N == 20 ~ "," ~ T^"*" == 4 ~ "," ~ H == 1))

final_grid_shrinkage <- plot_grid(ps1, ps2, ps3, ps4, ncol = 2)

final_grid_shrinkage

ggsave(file.path(output_dir, "mse_shrinkage_T5.pdf"), final_grid_shrinkage, width = 14, height = 9)
ggsave(file.path(output_dir, "mse_shrinkage_T5.eps"), final_grid_shrinkage, device = "eps", width = 14, height = 9)
