rm(list = ls())

# Load libraries
library(sf)
library(ggplot2)
library(readxl)

repos_path <- "/Users/js6249/repos/postdoc/metapop_covid19_kor"

# Date to plot
start_date_str <- c("February 18 2020", "August 3 2020", "October 12 2020", "February 1 2020")
end_date_str <- c("May 5 2020", "October 11 2020", "February 25 2021", "February 25 2021")
start_date_plot <- as.Date(start_date_str, format = "%B %d %Y")
end_date_plot <- as.Date(end_date_str, format = "%B %d %Y")

# Read shape file
korea_adm <- st_read(sprintf("%s/data/geographic/gadm41_KOR_shp/gadm41_KOR_1.shp", repos_path))

# Read case data
colname <- as.character(read_excel(sprintf("%s/data/incidence/COVID-19_confirmed_230303.xlsx", repos_path), 
                                   sheet = "시도별(17개시도+검역)", n_max = 1, col_names = FALSE))
case_df <- read_excel(sprintf("%s/data/incidence/COVID-19_confirmed_230303.xlsx", repos_path),
                      sheet = "시도별(17개시도+검역)", skip = 6, col_names = colname)

# Process case data
case_df <- case_df[c(1, 3:19)]# Remove useless column
for (i in 1:17)
{
  case_df[colname[i + 2]] <- lapply(case_df[colname[i + 2]], function(x) as.numeric(sub("-", "0", x)))
}

# Get date index
date_start <- as.Date("2020/01/20")
start_date_idx <- as.numeric(start_date_plot - date_start) + 1
end_date_idx <- as.numeric(end_date_plot - date_start) + 1

# Extract case data to plot
case_plot = data.frame(row.names = colname[3:19])
for (i in 1:4)
{
  date_idx_span = start_date_idx[i]:end_date_idx[i]
  case_plot[sprintf("wave%d", i)] = as.data.frame(colSums(case_df[date_idx_span, 2:18]))
}

colname_eng <- c("Seoul", "Busan", "Daegu", "Incheon", "Gwangju", "Daejeon",
                 "Ulsan", "Sejong", "Gyeonggi-do", "Gangwon-do", "Chungcheongbuk-do",
                 "Chungcheongnam-do", "Jeollabuk-do", "Jeollanam-do", "Gyeongsangbuk-do",
                 "Gyeongsangnam-do", "Jeju")# Replace column names to match with adm data
case_plot$NAME_1 <- colname_eng

# Merge with adm data
df_merge <- merge(korea_adm, case_plot, by = "NAME_1") 

# Plot heatmap
blue <- rgb(0, 114, 189, maxColorValue = 255)
orange <- rgb(217, 83, 25, maxColorValue = 255)
yellow <- rgb(237, 177, 32, maxColorValue = 255)
violet <- rgb(126, 47, 142, maxColorValue = 255)
green <- rgb(119, 172, 48, maxColorValue = 255)
skyblue <- rgb(77, 190, 238, maxColorValue = 255)
red <- rgb(162, 20, 47, maxColorValue = 255)

blue0 <- rgb(255 * 0.85, 147 + (255 - 147) * 0.85, 189 + (255 - 189) * 0.85, maxColorValue = 255)
orange0 <- rgb(217 + (255 - 217) * 0.85, 83 + (255 - 83) * 0.85, 25 + (255 - 25) * 0.85, maxColorValue = 255)
yellow0 <- rgb(237 + (255 - 237) * 0.85, 177 + (255 - 177) * 0.85, 32 + (255 - 32) * 0.85, maxColorValue = 255)
violet0 <- rgb(126 + (255 - 126) * 0.85, 47 + (255 - 47) * 0.85, 142 + (255 - 142) * 0.85, maxColorValue = 255)
green0 <- rgb(119 + (255 - 119) * 0.85, 172 + (255 - 172) * 0.85, 48 + (255 - 48) * 0.85, maxColorValue = 255)

h1 = ggplot(data = df_merge) +
  geom_sf(aes(fill = wave1), color = "grey") +
  scale_fill_gradient(low = blue0, high = blue, name = "Cumulative cases") +
  labs(title = "1st wave", subtitle = sprintf("%s ~ %s", start_date_str[1], end_date_str[1])) +
  theme_void() +
  theme(legend.text = element_text(size = 6, family = "Helvetica"), 
        legend.title = element_text(size = 6, family = "Helvetica"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "Helvetica-Bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6, family = "Helvetica"))
print(h1)
ggsave(sprintf("%s/result/plot_heatmap_case/wave1.png", repos_path), h1, width = 8.5, height = 6, units = 'cm')

h2 = ggplot(data = df_merge) +
  geom_sf(aes(fill = wave2), color = "grey") +
  scale_fill_gradient(low = blue0, high = blue, name = "Cumulative cases") +
  labs(title = "2nd wave", subtitle = sprintf("%s ~ %s", start_date_str[2], end_date_str[2])) +
  theme_void() +
  theme(legend.text = element_text(size = 6, family = "Helvetica"),
        legend.title = element_text(size = 6, family = "Helvetica"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "Helvetica-Bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6, family = "Helvetica"))
print(h2)
ggsave(sprintf("%s/result/plot_heatmap_case/wave2.png", repos_path), h2, width = 8.5, height = 6, units = 'cm')

h3 = ggplot(data = df_merge) +
  geom_sf(aes(fill = wave3), color = "grey") +
  scale_fill_gradient(low = blue0, high = blue, name = "Cumulative cases") +
  labs(title = "3rd wave", subtitle = sprintf("%s ~ %s", start_date_str[3], end_date_str[3])) +
  theme_void() +
  theme(legend.text = element_text(size = 6, family = "Helvetica"), 
        legend.title = element_text(size = 6, family = "Helvetica"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "Helvetica-Bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6, family = "Helvetica"))
print(h3)
ggsave(sprintf("%s/result/plot_heatmap_case/wave3.png", repos_path), h3, width = 8.5, height = 6, units = 'cm')

h4 = ggplot(data = df_merge) +
  geom_sf(aes(fill = wave4), color = "grey") +
  scale_fill_gradient(low = green0, high = green, name = "Cumulative cases") +
  labs(title = "Total", subtitle = sprintf("%s ~ %s", start_date_str[4], end_date_str[4])) +
  theme_void() +
  theme(legend.text = element_text(size = 9, family = "Helvetica"), 
        legend.title = element_text(size = 9, family = "Helvetica"),
        plot.title = element_text(hjust = 0.5, size = 10, family = "Helvetica-Bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, family = "Helvetica"))
print(h4)
ggsave(sprintf("%s/result/plot_heatmap_case/total.png", repos_path), h4, width = 8.5, height = 6, units = 'cm')
