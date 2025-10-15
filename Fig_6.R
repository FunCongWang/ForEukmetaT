
# Extract latent variables (first two MCOA axes)
load("data/MCOA_3.Rdata")
colnames(df5)[1:2]<-c("MCOA1", "MCOA2")
#df5<-df5[rownames(df5) != "LS4",]

env.sub<-metadata[,c("MAT","MAP","Temp.season","Perc.season",
                     "pH","TP","C_P","N_P",
                     "AP","AS","AK","Ca","Mg","Al","Mn")]

response <- df5[, c("MCOA1", "MCOA2")]

# Ensure environmental variables are clean
predictors <- env.sub

library(corrplot)
library(caret)

# Compute Spearman correlation matrix
cor_matrix <- cor(predictors, method = "spearman", use = "complete.obs")

# Plot to visualize
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# Identify and remove highly correlated variables
high_cor_pairs <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)
predictors_filtered <- predictors[, !(names(predictors) %in% high_cor_pairs)]

library(VSURF)

#SynVar1
set.seed(123)
vsurf_result <- VSURF(x = predictors_filtered, y = response$MCOA1)

# Selected variables
selected_vars <- names(predictors_filtered)[vsurf_result$varselect.interp]
selected_data <- predictors_filtered[, selected_vars]

# Load packages
library(terra)           # raster handling
library(randomForest)    # random forest model

# -------------------------
# 1. Load field dataset
# -------------------------
# Columns: MCOA1, MAT, TempSeason, PrecSeason
df <- cbind(df5[1],env.sub[c(1,3,4,9)])

# -------------------------
# 2. Fit Random Forest model
# -------------------------
set.seed(123) # reproducibility
rf_mod <- randomForest(
  MCOA1 ~ MAT + Temp.season + Perc.season + AP,
  data = df,
  ntree = 500,      # number of trees
  importance = TRUE # to check variable importance
)

print(rf_mod)
importance(rf_mod)

# -------------------------
# 3. Load climate rasters
# -------------------------
mat <- rast("data/soil_climate_tif/BIO1.tif")
tseason <- rast("data/soil_climate_tif/BIO4.tif")
pseason <- rast("data/soil_climate_tif/BIO15.tif")
availp<-rast("data/soil_climate_tif/Soil_P.tif")

# Align all rasters to the template (mat)
tseason_r <- resample(tseason, mat, method = "bilinear")
pseason_r <- resample(pseason, mat, method = "bilinear")
availp_r  <- resample(availp, mat, method = "bilinear")

# Now stack
clim_stack <- c(mat, tseason_r, pseason_r, availp_r)
names(clim_stack) <- c("MAT", "Temp.season", "Perc.season","AP")

# -------------------------
# 4. Predict fungal growth
# -------------------------
fungal_growth_map <- predict(clim_stack, rf_mod, na.rm = TRUE)

# -------------------------
# 5. Mask to forest areas
# -------------------------
forest_mask <- rast("data/predict/fc_CN.tif") # original mask
forest_mask_aligned <- resample(forest_mask, fungal_growth_map, method = "near")

fungal_growth_forest <- terra::mask(x=fungal_growth_map, mask=forest_mask_aligned)

# -------------------------
# 6. Save output
# -------------------------
#writeRaster(fungal_growth_forest, "fungal_growth_rate_RF.tif", overwrite = TRUE)

library(terra)

# Load China boundary (downloaded shapefile or from rnaturalearth)
china <- vect("data/site_map/china/china.shp")  # replace with your file

# Crop prediction to China extent
fungal_growth_china <- crop(fungal_growth_forest, china)

# Mask to China boundary
fungal_growth_china <- terra::mask(fungal_growth_china, china)

# Plot
plot(fungal_growth_china, main="Predicted Fungal Growth Rate in Chinese Forests")
lines(china, lwd=1.5, col="lightgrey")


#plot with ggplot2
library(terra)
library(ggplot2)
library(colorRamps)
library(sf)

# Convert raster to data frame
r_df <- as.data.frame(fungal_growth_china, xy = TRUE, na.rm = TRUE)
names(r_df)[3] <- "MCOA1"

# Convert China boundary to data frame
china_df <- as.data.frame(china, geom="WKT") # or use fortify(sf_object) if sf

# Plot with ggplot
mapCN<-ggplot() +
  geom_tile(data = r_df, aes(x = x, y = y, fill = MCOA1)) +
  geom_sf(data = st_as_sf(china), fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option="plasma") +
  #scale_fill_gradient(low = "blue", high = "red")+
  coord_sf(xlim = c(73, 138), ylim = c(15, 55), expand = FALSE) +
  scale_y_continuous(breaks = seq(20,50,by=10))+
  theme_bw() +
  labs(title="Prediction of MCOA1 in Chinese Forests",
       fill="MCOA1",
       x= "Longitude",
       y= "Latitude") +
  #guides(fill=guide_legend(position = "inside"))+
  theme(plot.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        #legend.justification.inside = c(0.5,0.5),
        #legend.box = "horizontal",
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
mapCN

mapCN_nineline<-ggplot() +
  geom_tile(data = r_df, aes(x = x, y = y, fill = MCOA1)) +
  geom_sf(data = st_as_sf(china), fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option="plasma") +
  coord_sf(ylim = c(3,22),xlim = c(107,120))+
  scale_x_continuous(breaks = seq(108,120,4))+
  theme_bw()+theme(aspect.ratio = 1.5,plot.margin=unit(c(0,0,0,0),"mm"),
                   legend.position = "none",axis.title = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1),
                   panel.grid = element_blank())
mapCN_nineline

library(cowplot)
ggdraw()+
  draw_plot(mapCN)+
  draw_plot(mapCN_nineline,x = 0.66, y = 0.15, width = 0.25, height = 0.25)

#ggsave("plot2/predict_mcoa1_2025.09.11.pdf",height = 6,width = 8)
