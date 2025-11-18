#install.packages("openxlsx")
#install.packages("rgplates")
#install.packages("geojsonsf")
library(openxlsx)
library(rgplates)
library(geojsonsf)
library(tidyverse)

#set working directory
wd = "~/ts/labs/lab03/_Gplates/"
setwd(wd)

fn = "12_regions_points.xlsx"
xls = openxlsx::read.xlsx(xlsxFile=fn, sheet="present")
xls

# Plot your polygons on a map:
world_extent = NULL
world_extent$long = c(-180,180,180,-180)
world_extent$lat = c(90, 90,-90,-90)
df <- data.frame(lon=world_extent$long, lat=world_extent$lat)
world_polygon <- df %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
summarise(geometry = st_combine(geometry)) %>%
st_cast("POLYGON")
#plot(world_polygon, col="lightblue")
plot(world_polygon, col="lightblue", reset=FALSE)

list_of_polygons = list()
area_names = unique(xls$area)
cols = rainbow(n=length(area_names))
i = 1
for (i in 1:length(area_names))
	{
	rowTF = xls$area == area_names[i]
	
	tmpxls = xls[rowTF,]
	df <- data.frame(lon=tmpxls$long, lat=tmpxls$lat)
	
	tmp_polygon <- df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
	tmp_polygon
	plot(tmp_polygon, col=cols[i], add=TRUE)
	}

# doesn't quite look right!

# Try again:
world_extent = NULL
world_extent$long = c(-180,180,180,-180)
world_extent$lat = c(90,90,-90,-90)
df <- data.frame(lon=world_extent$long, lat=world_extent$lat)
world_polygon <- df %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
summarise(geometry = st_combine(geometry)) %>%
st_cast("POLYGON")
#plot(world_polygon, col="lightblue")
plot(world_polygon, col="lightblue", reset=FALSE)

list_of_polygons = list()
area_names = unique(xls$area)
cols = rainbow(n=length(area_names))
i = 8
for (i in 1:length(area_names))
	{
	rowTF = xls$area == area_names[i]
	tmpxls = xls[rowTF,]
	# edit points just east of -180, for asia
	if (area_names[i] == "P_Palearctic")
		{
		tmpxls$long[tmpxls$long < -169] = 180
		}
	
	df <- data.frame(lon=tmpxls$long, lat=tmpxls$lat)
	tmp_polygon <- df %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
	tmp_polygon
	plot(tmp_polygon, col=cols[i], add=TRUE)
	}




#######################################################
# Reconstruct those points back in time
#######################################################

ages = seq(from=1, to=100, by=1)
reconstructed_points_table = xls

runslow = FALSE
if (runslow)
	{
	cat("\nReconstructing paleolongitudes (x) & paleolatitudes (y) for ", length(ages), " ages, ending at ", max(ages), "...", sep="")
	i=1
	for (i in 1:length(ages))
		{
		# Nice to print each step so you can see timing
		cat(ages[i], ",", sep="")
		
		tmpresult = rgplates::reconstruct(x=xls[,c("long","lat")], age=ages[i], model="MERDITH2021")
		
		# Change the names of the columns to specify age
		colnames(tmpresult) = paste0(c("long_","lat_"), ages[i])
		tmpresult
		
		# Add tmpresult to table
		reconstructed_points_table = cbind(reconstructed_points_table, tmpresult)
		}
	# Save reconstructed_points_table to Rdata
	save(reconstructed_points_table, file="reconstructed_points_table.Rdata")
	} else {
	# loads to 'reconstructed_points_table'
	load(file="reconstructed_points_table.Rdata")
	}

reconstructed_points_table

head(reconstructed_points_table[,1:7])
tail(reconstructed_points_table[,1:7])

ages = c(0,ages)


pdffn = "polygons_through_time.pdf"
pdf(file=pdffn, width=6, height=6)

# Try again:
world_extent = NULL
world_extent$long = c(-180,180,180,-180)
world_extent$lat = c(90,90,-90,-90)
df <- data.frame(lon=world_extent$long, lat=world_extent$lat)
world_polygon <- df %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
summarise(geometry = st_combine(geometry)) %>%
st_cast("POLYGON")
#plot(world_polygon, col="lightblue")
#plot(world_polygon, col="lightblue", reset=FALSE)

timeval = 1+2*(40)
list_of_polygons = list()
area_names = unique(xls$area)
cols = rainbow(n=length(area_names))
i = 8
cat("Plotting polygons through time...a=")
reconstruction_col = 0
for (ai in 1:length((ages)))
	{
	a = ages[ai]
	cat(a,",", sep="")
	timeval = ages[ai]
	reconstruction_col = 2*ai
	plot(world_polygon, col="lightblue", reset=FALSE, xlim=c(-180, 180), ylim=c(-90,90))
	title(paste0(a, " million years ago"))

	for (i in 1:length(area_names))
		{
		rowTF = xls$area == area_names[i]
		tmpxls = reconstructed_points_table[rowTF,c(reconstruction_col:(reconstruction_col+1))]
		naTF1 = is.na(tmpxls[,1])
		naTF2 = is.na(tmpxls[,2])
		naTF = (naTF1 + naTF2) == 0
		tmpxls = tmpxls[naTF,]
		
		# edit points just east of -180, for asia
		if (area_names[i] != "N_Nearctic")
			{
			tmpxls$long[tmpxls$long < -169] = 179
			}
		
		df <- data.frame(lon=tmpxls$long, lat=tmpxls$lat)
		tmp_polygon <- df %>%
		st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
		summarise(geometry = st_combine(geometry)) %>%
		st_cast("POLYGON")
		tmp_polygon
		plot(tmp_polygon, col=cols[i], add=TRUE)
		}
	} # END for (a in ages)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)








#######################################################
# Distances between polygons
#######################################################

#######################################################
# Distances between polygons
#######################################################
xcol = 2
min_dist = array(data=NA, dim=c(length(area_names), length(area_names), ((ncol(reconstructed_points_table)-1)/2)))
centroid_dist = array(data=NA, dim=c(length(area_names), length(area_names), ((ncol(reconstructed_points_table)-1)/2)))
i=1
j=1
k=2
cat("\nCalculating distances for time slice i=\n")
for (i in 1:((ncol(reconstructed_points_table)-1)/2))
	{
	cat(i, ",", sep="")
	for (j in 1:length(area_names))
		{
		for (k in 1:length(area_names))
			{
			rowTF = xls$area == area_names[j]
			tmpxls = xls[rowTF,]
			
			rowTF2 = xls$area == area_names[k]
			tmpxls2 = xls[rowTF2,]
			
			# Edit to remove NA points - 
			# Gplates can LOSE points back in time
			# e.g. if the tectonic plate "disappears"
			long1 = reconstructed_points_table[rowTF,xcol]
			lat1 = reconstructed_points_table[rowTF,xcol+1]
			naTF1 = is.na(long1)
			naTF2 = is.na(lat1)
			naTF = (naTF1 + naTF2) > 0
			long1 = long1[naTF == FALSE]
			lat1 = lat1[naTF == FALSE]

			long2 = reconstructed_points_table[rowTF2,xcol]
			lat2 = reconstructed_points_table[rowTF2,xcol+1]
			naTF1 = is.na(long2)
			naTF2 = is.na(lat2)
			naTF = (naTF1 + naTF2) > 0
			long2 = long2[naTF == FALSE]
			lat2 = lat2[naTF == FALSE]
			
			area1 = cbind(long1,lat1)
			colnames(area1) = c("x","y")
			area2 = cbind(long2,lat2)
			colnames(area2) = c("x","y")

			ensure_closure = FALSE
			if (ensure_closure == TRUE)
				{
				# Make sure the last coordinate matches the first coordinate
				# (this is usually the case, but coordinates can get dropped)
				if ( (area1[,"x"][1] != area1[,"x"][length(area1[,"x"])]) || (area1[,"y"][1] != area1[,"y"][length(area1[,"y"])]) )
					{
					area1 = rbind(area1, area1[1,])
					}
				if ( (area2[,"x"][1] != area2[,"x"][length(area2[,"x"])]) || (area2[,"y"][1] != area2[,"y"][length(area2[,"y"])]) )
					{
					area2 = rbind(area2, area2[1,])
					}
				} # END if (ensure_closure == TRUE)
			
			make_points_unique = FALSE
			if (make_points_unique == TRUE)
				{
				pts_txt = apply(X=area1, MARGIN=1, FUN=paste0, collapse="_")
				unique(pts_txt)
				nums_to_keep = match(x=unique(pts_txt), table=pts_txt)
				area1 = area1[nums_to_keep,]

				pts_txt = apply(X=area2, MARGIN=1, FUN=paste0, collapse="_")
				unique(pts_txt)
				nums_to_keep = match(x=unique(pts_txt), table=pts_txt)
				area2 = area2[nums_to_keep,]
				}
			
			#area1_alphahull = alphahull::ahull(x=area1, alpha=1)
			#area2_alphahull = alphahull::ahull(x=area2, alpha=1)
			
			# doesn't like non-convex points
			#area1_convexhull = spatstat.geom::convexhull.xy(area1)
			#area2_convexhull = spatstat.geom::convexhull.xy(area2)
			
			#area1_ahull_mat = area1_alphahull$xahull
			#area2_ahull_mat = area2_alphahull$xahull
			#colnames(area1_ahull_mat) = c("x","y")
			#colnames(area2_ahull_mat) = c("x","y")
			#area1_ahull_mat = rbind(area1_ahull_mat, area1_ahull_mat[1,])
			#area2_ahull_mat = rbind(area2_ahull_mat, area2_ahull_mat[1,])
			
			# CRS: Coordinate Reference System
			# WGS84 (EPSG:4326)
			#poly1 <- st_polygon(list(area1_ahull_mat)) %>% st_sfc(crs = 4326)
			#poly2 <- st_polygon(list(area2_ahull_mat)) %>% st_sfc(crs = 4326)
			
			area1_mp = st_multipoint(x=area1, dim="XY") %>% st_sfc(crs = 4326)
			area2_mp = st_multipoint(x=area2, dim="XY") %>% st_sfc(crs = 4326)
			st_centroid(area1_mp)
			#poly1 <- st_concave_hull(x=st_multipoint(x=area1, dim="XY"), ratio=0.5, allow_holes=FALSE) %>% st_sfc(crs = 4326)
			#poly2 <- st_concave_hull(st_multipoint(x=area2, dim="XY"), ratio=0.5, allow_holes=FALSE) %>% st_sfc(crs = 4326)
			

			#z = suppressWarnings(try(min(spDists(x=area1, y=area2, longlat=TRUE, diagonal=FALSE), na.rm=TRUE)))
			#z = suppressWarnings(try(st_distance(x=poly1, y=poly2)))
			z = try(st_distance(x=area1_mp, y=area2_mp))
			z2 = try(st_distance(x=st_centroid(area1_mp), y=st_centroid(area2_mp)))
			
			
			
			if (class(z) == "try-error")
				{
				min_dist[j,k,i] = NA
				} else {
				min_dist[j,k,i] = z
				}
			if (class(z2) == "try-error")
				{
				centroid_dist[j,k,i] = NA
				} else {
				centroid_dist[j,k,i] = z2
				}
			}
		}
	xcol = xcol+2
	}


set_diag_to_NA <- function(ma)
	{
	diag(ma) = NA
	}

set_diag_to_NA_array <- function(ma)
	{
	TFmat = matrix(data=FALSE, nrow=dim(ma)[1], ncol=dim(ma)[2])
	diag(TFmat) = TRUE
	num_slices = dim(ma)[3]
	for (n in 1:num_slices)
		{
		diag(ma[,,n]) = NA
		}
	return(ma)
	}


#######################################################
# Distances tend to be correlated
#######################################################
xvals = c(set_diag_to_NA_array(centroid_dist))
yvals = c(set_diag_to_NA_array(min_dist))
naTF1 = !is.na(xvals)
naTF2 = !is.na(yvals)
naTF = (naTF1 + naTF2) == 2
xvals = xvals[naTF] / 1000
yvals = yvals[naTF] / 1000

library(BioGeoBEARS)
linear_regression_plot(x=xvals, y=yvals, xlabel="centroid distance (km)", ylabel="minimum distance (km)", tmppch=".", pointscol="blue")
title("Centroid distance vs. minimum distance, 0-100 Ma")





