###############
## reconstruction of a tree as polygons from a point cloud using a voxel approach
## for each voxel a best fitting plane is used to project the 3d points to the plane 
## and use the convex hull to form a polygon. If a voxel contains less than 3 points
## a random rotated circle is added
## this process is separately done for leafs and wood
## by Julian Frey and Zoe Schindler
## 2023-09-19
###############

# Load required libraries
library(rgl)
library(lidR)
library(CspStandSegmentation)
library(doParallel)

###############
## Settings
###############

filename <- "V:/RayLeaf/tree_summer_1cm.las" # las filename
vox_size <- 0.05 # voxel size in m for polygon fitting
buffer <- sqrt(2*(0.01^2))/2 # buffer polygons by this size to avoid artificial gaps
min_pts_vox <- 3 # minimum number of points per plain fitting voxel 
outdir <- "V:/RayLeaf/out/" # output directory (with tailing /)
num_cores <- 20 # number of cores for parallel processing
singletons_factor <- 1/3 # multiplication factor for voxel size to create circles from single points in voxels (if <3pts in voxel)

###############
## Functions
###############

# function to fit a planar convex hull polygon to a set of rows of3d points "mat"
# it extends the convex hull by the buffer value
# if plot_vox is set to TRUE it will plot the 3d points and the convex hull
plane_3d <- function(mat, buffer = 0.005, plot_vox = F){
  if(nrow(mat) < 3) {
    warning("Less than 3 points in voxel, return NA")
    return(NA)
    } # if less than 3 points return NA (no plane possible)
  # calculate center and center mat 
  mean_X <- mean(mat[,1])
  mean_Y <- mean(mat[,2])
  mean_Z <- mean(mat[,3])
  mat[,1] <- mat[,1] - mean_X
  mat[,2] <- mat[,2] - mean_Y
  mat[,3] <- mat[,3] - mean_Z
  
  # add plot
  if(plot_vox){
    plot3d(mat[,1],mat[,2],mat[,3], aspect = F, type = "p", col = grey(0.5, 0.5), decorate = T, axes = F, xlab = "X", ylab = "Y", zlab = "Z"); box3d()
  }
  # # retreive the orientation of the plot3d
  # if(!exists(view)){
  #   view <- par3d("userMatrix")
  # } else {
  #   rgl.viewpoint(userMatrix = view)
  # }
  # 
  # save the current view as svg
  #rgl::rgl.postscript("Z:/3_projects/A1/administration/Publications/blender/process/points.svg", fmt = "svg")
  if(nrow(mat) == 3) {
    convex_hull <- mat
    } else {
    # Perform Singular Value Decomposition (SVD)
    svd_result <- svd(mat)
    
    # Extract the normal vector of the best-fitting plane (from the last column of V)
    normal_vector <- svd_result$v[, 3]
    
    # add the plane normal to the plot
    if(plot_vox){
      #normal_vector_sc <- normal_vector * 0.01
      #arrow3d(c(0,0,0),c(normal_vector_sc[1],normal_vector_sc[2],normal_vector_sc[3]), col = "red")
      planes3d(normal_vector, alpha = 0.5, color = rgb(0,170,212,maxColorValue = 255))
    }
    #rgl::rgl.postscript("Z:/3_projects/A1/administration/Publications/blender/process/plane.svg", fmt = "svg")
    
    # Project 3D points onto the plane
    distances <- as.matrix(mat) %*% normal_vector
    points_2d <- mat - distances %*% t(normal_vector)
    
    # add the projected points to the plot with lines to the original points
    if(plot_vox){
      for(i in 1:nrow(mat)){
        segments3d(c(mat[i,1],points_2d[i,1]),c(mat[i,2],points_2d[i,2]),c(mat[i,3],points_2d[i,3]), col = "darkblue")
      }
      points3d(points_2d[,1],points_2d[,2],points_2d[,3], col = rgb(0,170,212,maxColorValue = 255))
    }
    #rgl::rgl.postscript("Z:/3_projects/A1/administration/Publications/blender/process/proj_pts_no_plane.svg", fmt = "svg")
    
    # Compute the Convex Hull from the two bigger components of the normal vector
    abs_normal_order <- order(abs(normal_vector))
    chull_idx <- chull(points_2d[,c(abs_normal_order[3], abs_normal_order[2])])
    convex_hull <- points_2d[c(chull_idx, chull_idx[1]),]
  }
  # add the convex hull to the plot
  if(plot_vox){
    polygon3d(convex_hull[,1],convex_hull[,2],convex_hull[,3], col = rgb(0,170,212,maxColorValue = 255))
  }
  #rgl::rgl.postscript("Z:/3_projects/A1/administration/Publications/blender/process/conv_hull.svg", fmt = "svg")
  
  # buffer the convex hull
  # Generate centroid for buffering
  centroid <- matrix(rep(0,3), nrow = 1)
  
  # Calculate the direction vectors from the centroid to each edge
  directions <- convex_hull - as.numeric(centroid)
  
  # Normalize the direction vectors to have a length of 'distance'
  norms <- sqrt(rowSums(directions^2))
  normalized_directions <- directions / matrix(rep(norms, each = 3), ncol = 3, byrow = TRUE) * buffer
  
  # Expand the polygon by adding the normalized directions to the original points
  convex_hull_ext <- convex_hull + normalized_directions
  
  # add the convex hull to the plot
  if(plot_vox){
    lines3d(convex_hull_ext[,1],convex_hull_ext[,2],convex_hull_ext[,3], col = "darkblue")
  }
  #rgl::rgl.postscript("Z:/3_projects/A1/administration/Publications/blender/process/conv_hull_ext.svg", fmt = "svg")
  return(data.frame(X = convex_hull_ext[,1] + mean_X,Y = convex_hull_ext[,2] + mean_Y,Z = convex_hull_ext[,3] + mean_Z))
}

# creates a planar circle in 3d space with center 'point', radius 'buffer'. 
# random_rot rotates the circle randomly in 3d and segments is the number of
# points to describe the circle
create_circle <- function(point, buffer = 0.005, random_rot = T, segments = 6){
  circle <- cbind(conicfit::calculateCircle(0,0,buffer, steps = segments+1),0) # create circle and add third dimension
  if(random_rot) {circle <- circle %*% random_rotation_matrix()}
  poly <- data.frame(X= circle[,1] + point[1],Y = circle[,2] + point[2], Z = circle[,3] +  point[3])
  return(poly)
}

# triangulation of polygons based on rgl
# if it fails on the first try it drops single points and 
# retries to avoid failures if three points for a straight line
create_shape <- function(mat, col = "limegreen") {
  poly_i <- NULL
  if(any(mat[1,] != mat[nrow(mat),])) {mat <- mat[c(1:nrow(mat),1),]}
  t <- 0
  while(is.null(poly_i) & t <= nrow(mat)){
    if(t == 0) { # try to fit triangles to all points
      i <- mat
    } else { # point removal if it does not work
      i <- mat[-t,]
    }
    if (is.null(poly_i)) {
      try({poly_i <- rgl::polygon3d(
        x = i[,1],
        y = i[,2],
        z = i[,3],
        plot = FALSE, coords = c(1,2), random = F)}, silent = T)
      if (is.null(poly_i)) {
        try({poly_i <- rgl::polygon3d(
          x = i[,1],
          y = i[,2],
          z = i[,3],
          plot = FALSE, coords = c(2,3), random = F)}, silent = T)
        if (is.null(poly_i)) {
          try({poly_i <- rgl::polygon3d(
            x = i[,1],
            y = i[,2],
            z = i[,3],
            plot = FALSE, coords = c(1,3), random = F)}, silent = T)
        }
      }
    }
    
    t <- t+1
  }
  if (!is.null(poly_i)) {
    poly_i$material$color <- col
    poly_i
  }
}

# function to get a class index from the most abundant class in a set of points
# used to choose if a voxel is treated as leafs or wood
most_abundant <- function(x) as.integer(names(sort(table(x), decreasing = T))[1])

# random rotation matrix generation stolen from soobench package
random_rotation_matrix <- function(d=3) {
  simple_rotation_matrix <- function(d, i, j, alpha) {
    R <- diag(d)
    R[i, i] <- cos(alpha)
    R[i, j] <- sin(alpha)
    R[j, i] <- -sin(alpha)
    R[j, j] <- cos(alpha)
    R
  }
  
  R <- diag(d)
  for (i in 2:d)
    R <- R %*% simple_rotation_matrix(d, 1, i, runif(1, -pi/4, pi/4))
  
  if (d > 2) {
    for (i in 2:(d-1))
      R <- R %*% simple_rotation_matrix(d, i, d, runif(1, -pi/4, pi/4))
  }
  R
}

# pol: a single mesh3d object from rgl or a list of them
# file: the filepath
# mtl_name: character containing the name for the Material specified in a separate mtl file 
write_obj <- function(pol, file, mtl_name = "Material1"){
  # generate header information
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()), paste("usemtl", mtl_name), 'o object1')
  #check if its a list or just one polygon
  if(is.list(pol)){
    l <- 0 # index for vertex numbers already used
    lines <- character() # obj text lines 
    p_i <- 1 # polygon index
    for(p in pol){
      #g <- paste0("o polygon", p_i) # create a group for evry polygon
      v = apply(round(p$vb,5), 2, function(x) paste("v", paste(x[1:3], collapse = " "))) # add the coordinates round(p$vb,5)
      f = apply(p$it+l, 2, function(x) paste("f", paste(x[1:3], collapse = " "))) # add the triangle coordinate order
      f_rev = apply(p$it+l, 2, function(x) paste("f", paste(x[3:1], collapse = " "))) # add the reverse order for backsides
      l <- l+ncol(p$vb) # increase the vertex index 
      lines <- c(lines, c(v,f,f_rev)) # add the lines
      p_i <- p_i + 1 # increase the polygon index
    }
    writeLines(c(obj_head,lines), file) # write to file
  } else {
    v = apply(pol$vb, 2, function(x) paste("v", paste(x[1:3], collapse = " ")))
    f = apply(pol$it, 2, function(x) paste("f", paste(x[1:3], collapse = " ")))
    writeLines(c(obj_head,v,f), file)
  }
}

# parallel version of the function above
# pol: a single mesh3d object from rgl or a list of them
# file: the filepath
# mtl_name: character containing the name for the Material specified in a separate mtl file 
write_obj_par <- function(pol, file, mtl_name = "Material1", cores = 5){
  # generate header information
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()), paste("usemtl", mtl_name), 'o object1')
  #check if its a list or just one polygon
  if(is.list(pol)){
    l <- 0 # index for vertex numbers already used
    cl <- makeCluster(cores)
    lines <- character()
    for(p in pol){
      v = parApply(cl,p$vb, 2, function(x) paste("v", paste(round(x[1:3],5), collapse = " "))) # add the coordinates round(p$vb,5)
      f = parApply(cl,p$it+l, 2, function(x) paste("f", paste(x[1:3], collapse = " "))) # add the triangle coordinate order
      f_rev = parApply(cl,p$it+l, 2, function(x) paste("f", paste(x[3:1], collapse = " "))) # add the reverse order for backsides
      l <- l+ncol(p$vb) # increase the vertex index 
      lines <- c(lines, c(v,f,f_rev))# add the lines
    }
    stopCluster(cl)
    writeLines(c(obj_head,lines), file) # write to file
  } else {
    v = apply(pol$vb, 2, function(x) paste("v", paste(x[1:3], collapse = " ")))
    f = apply(pol$it, 2, function(x) paste("f", paste(x[1:3], collapse = " ")))
    writeLines(c(obj_head,v,f), file)
  }
}

# function to write polygons without triangulation by rgl
# pol: a list of data.frames with 3d points for each polygon
# file: the filepath
# mtl_name: character containing the name for the Material specified in a separate mtl file 
write_obj_df <- function(pol, file, mtl_name = "Material1"){
  # generate header information
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()), paste("usemtl", mtl_name), 'o object1')
  write(obj_head, file)
  #check if its a list or just one polygon
  if(is.list(pol)){
    l <- 0 # index for vertex numbers already used
    lines <- character()
    for(p in pol){
      v = apply(p, 1, function(x) paste("v", paste(round(x,5), collapse = " "))) # add the coordinates round(p$vb,5)
      f = paste("f", paste(1:nrow(p)+l, collapse = " ")) # add the triangle coordinate order
      f_rev = paste("f", paste(rev(1:nrow(p))+l, collapse = " ")) # add the reverse order for backsides
      l <- l+nrow(p) # increase the vertex index 
      #lines <- c(lines, c(v,f,f_rev))# add the lines
      write(c("\n",v,f,f_rev), file, append = T)
    }
    #writeLines(c(obj_head,lines), file) # write to file
  } else {
    warning("pol needs to be a list with data.frames")
    stop()
  }
}


###############
## read and preprocess data 
###############

las <- readLAS(filename)

# apply offset for to avoid numerical problems
means <- apply(las@data,2,mean)
las@data$X <- las@data$X - means[1]
las@data$Y <- las@data$Y - means[2]
las@data$Z <- las@data$Z - means[3]


# classify leaf/wood (Reflectance < -5.4)
hist(las$Reflectance, breaks = 100)
abline(v = -5.4)
las$Classification[las$Reflectance < -5.4] <- 4L

# use 5cm voxel to assign the most abundant class to all points in voxel
class_vox <- las |> voxel_metrics(most_abundant(Classification), res = vox_size)
las <- las |> add_voxel_coordinates(res = vox_size)
las@data <- merge(las@data, class_vox, by.x = c("x_vox","y_vox","z_vox"), by.y = c("X","Y","Z"))

###############
## Fit a polygon to every voxel
###############

# filter leaves and calculate voxel coordinates for loop and
# create a fitting plane convex hull for every 5cm leave voxel
vegetation <- las |>  filter_poi(V1 == 4)
vegetation_vox <- vegetation |> voxel_metrics(length(Z), vox_size)
vegetation_vox_singles <- vegetation_vox[vegetation_vox$V1 < min_pts_vox,]
vegetation_vox <- vegetation_vox[vegetation_vox$V1 >= min_pts_vox,]

# Register a parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# parallel process for voxel extraction and convex hull fitting
t1 <- system.time({ # measure execution time
  polys <- foreach(v = 1:nrow(vegetation_vox), .noexport = c("las")) %dopar% {
    tryCatch({
      vox_bb <- vegetation_vox[v,]
      vox <- vegetation |> lidR::filter_poi(x_vox == vox_bb$X[1], y_vox == vox_bb$Y[1], z_vox == vox_bb$Z[1])
      plane_3d(as.matrix(vox@data[, 1:3]), buffer)
    }, error = function(e) {return(v) })
  }
  fails <- sapply(polys, is.integer)
  if(any(fails)){ 
    polys[fails] <- foreach(v = simplify2array(polys[fails]), .noexport = c("las")) %dopar% { 
      tryCatch({
        vox_bb <- vegetation_vox[v,]
        vox <- vegetation |> lidR::filter_poi(x_vox == vox_bb$X[1], y_vox == vox_bb$Y[1], z_vox == vox_bb$Z[1])
        plane_3d(as.matrix(vox@data[, 1:3]), buffer)
      }, error = function(e) {return(v) })
    }
  }
})

# Stop the parallel back end
stopCluster(cl)
print(paste("Processing time for leafs:", t1[3]/60, "min"))

# add circles to single points in voxels
vegetation_singles <- merge(las@data, vegetation_vox_singles, by.x = c("x_vox", "y_vox", "z_vox"), by.y = c("X","Y","Z"))[,c("X","Y","Z")]
vegetation_singles_polys <- list()
for(i in 1:nrow(vegetation_singles)){
  vegetation_singles_polys[[i]] <- create_circle(as.numeric(vegetation_singles[i,]), buffer = vox_size * singletons_factor)
}

# combine convex hulls and circles
veg_polys <- c(polys, vegetation_singles_polys)

# create a fitting plane convex hull for every 5cm wood voxel
bark_polys <- list()
bark <- las |>  filter_poi(V1 == 0)
bark_vox <- bark |> voxel_metrics(length(Intensity), vox_size)
bark_vox_singles <- bark_vox[bark_vox$V1 < min_pts_vox,]
bark_vox <- bark_vox[bark_vox$V1 >= min_pts_vox,]

# Register a parallel back end
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# parallel process for voxel extraction and convex hull fitting
t2 <- system.time({ # measure execution time
  bark_polys <- foreach(v = 1:nrow(bark_vox), .noexport = c("las")) %dopar% { # for some reason this fails randomly
    tryCatch({
      vox_bb <- bark_vox[v,]
      vox <- bark |> lidR::filter_poi(x_vox == vox_bb$X[1], y_vox == vox_bb$Y[1], z_vox == vox_bb$Z[1])
      plane_3d(as.matrix(vox@data[, 1:3]), buffer)
    }, error = function(e) {return(v) })
  }
  fails <- sapply(bark_polys, is.integer)
  if(any(fails)){ # therefore we try again
    bark_polys[fails] <- foreach(v = simplify2array(bark_polys[fails]), .noexport = c("las")) %dopar% { 
      tryCatch({
        vox_bb <- bark_vox[v,]
        vox <- bark |> lidR::filter_poi(x_vox == vox_bb$X[1], y_vox == vox_bb$Y[1], z_vox == vox_bb$Z[1])
        plane_3d(as.matrix(vox@data[, 1:3]), buffer)
      }, error = function(e) {return(v) })
    }
  }
})

stopCluster(cl)

# create polygons for single points
bark_singles <- merge(las@data, bark_vox_singles, by.x = c("x_vox", "y_vox", "z_vox"), by.y = c("X","Y","Z"))[,c("X","Y","Z")]
bark_singles_polys <- list()
for(i in 1:nrow(bark_singles)){
  bark_singles_polys[[i]] <- create_circle(as.numeric(bark_singles[i,]), buffer = vox_size * 0.25)
}

# combine convex hulls and circles
bark_polys <- c(bark_polys, bark_singles_polys)

# # write dtm as stl to disk
# dtm_tin <- rasterize_terrain(filter_poi(las, Classification == 2), res = 0.1)
# plot_dtm3d(dtm_tin, clear_artifacts = F)
# writeSTL(paste0(outdir, "dtm2.stl"))
# close3d()

# # triangulate polygons using rgl
# cl <- makeCluster(num_cores)
# item_rgl <- parLapply(cl = cl, X = veg_polys, fun = create_shape)
# stopCluster(cl)
# 
# # this failes for single triangle-polygons with duplicated points
# # which get excluded
# failed_rgl <- which(sapply(item_rgl,is.null))
# length(failed_rgl)
# item_rgl <-  item_rgl[!sapply(item_rgl,is.null)]
# 
# bark_item_rgl  <- lapply(bark_polys, create_shape, col = "chocolate4")
# bark_item_rgl <-  bark_item_rgl[!sapply(bark_item_rgl,is.null) & (sapply(bark_item_rgl,typeof) == "list")]

# # export bark and leaf polygons separately
# close3d()
# with(las@data, plot3d(X,Y,Z, aspect = F, type = "n", col = grey(0.5, 0.5), decorate = F))
# rgl::shade3d(shapelist3d(item_rgl, plot = FALSE), lit = F)
# writeSTL(paste0(outdir, "leaves_new.stl"))
# close3d()
# with(las@data, plot3d(X,Y,Z, aspect = F, type = "n", col = grey(0.5, 0.5), decorate = F), lwd = NULL)
# rgl::shade3d(shapelist3d(bark_item_rgl, plot = FALSE), lit = T)
# writeSTL(paste0(outdir, "bark_new.stl"))
# close3d()
write_obj_df(veg_polys, paste0(outdir, "leaves.obj"))
write_obj_df(bark_polys, paste0(outdir, "bark.obj"))

# ###################
# ## calculate bbox for blender images
# ###################
# 
# img_dims <- c(2048,1600)
# ortho_scale <- 29.2
# px_size <- ortho_scale/img_dims[1] 
# 
# ul <- means[1:2] + (img_dims/2) * px_size * c(-1,1)# upper left 
# lr  <- means[1:2] + (img_dims/2) * px_size * c(1,-1) # lower right
# 
# ###################
# ## write csv with metadata
# ###################
# 
# data.frame(
#   data = c(
#     "filename",
#     "vox_size",
#     "X_offset",
#     "Y_offset",
#     "Z_offset",
#     "px_size",
#     "ul_X",
#     "ul_Y",
#     "lr_X",
#     "lr_Y"
#   ),
#   value = c(
#     filename,
#     vox_size,
#     means[1],
#     means[2],
#     means[3],
#     px_size,
#     ul[1],
#     ul[2],
#     lr[1],
#     lr[2]
#   )
# ) |> write.csv2(paste0(outdir,"metadata.csv"))
# 
# 
# # code snipped to add worldfiles to png's created by blender
# 
# worldfile <- paste(px_size, 0,0,-px_size, ul[1], ul[2], sep = "\n")
# 
# png_dir <- "V:/RayLeaf/out/shade_sec/"
# files <- list.files(png_dir, pattern = "*.png")
# for(f in files) writeLines(worldfile, paste0( png_dir, strsplit(f, "[.]")[[1]][1], ".pgw"))
