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
library(ggplot2)
library(terra)

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

# less settings
n_bands <- 60 # number of bands spectral bands

###############
## Functions
###############

# function to fit a planar convex hull polygon to a set of rows of3d points "mat"
# it extends the convex hull by the buffer value
# @param mat a matrix with 3d points
# @param buffer the buffer value to extend the convex hull
# @param plot_vox if TRUE it will plot the 3d points and the convex hull
# @param vox_size the voxel size to fit polygons
# @param singletons_factor multiplication factor for voxel size to create circles from single points in voxels (if <3pts in voxel)
# @param segments number of segments for the single point circle

plane_3d <- function(mat, buffer = 0.005, plot_vox = F, vox_size = 0.05, singletons_factor = 1/3, segments = 6){
  if(nrow(mat) < 3) {
    return(create_circle(colMeans(mat), buffer = vox_size * singletons_factor, segments = segments))
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
# @param pol a list of data.frames with 3d points for each polygon
# @param file the filepath
# @param mtl_name character containing the name for the Material specified in a separate mtl file 
# @param fast if set to TRUE it will use less accurate rounding to speed up the process and create smaller files
# @param backfaces if set to TRUE it will add the reverse order of the triangles to create backfaces
# @param mode if set to "normal" it will write the polygons as they are, if set to "dart" or "less" it will change the order of the dimensions to be used in rtm
# @param dart_x_off if mode id "dart" subtract each x-value from this value to invert the x-axis
write_obj_df <- function(pol, file, mtl_name = "Material1",mtl_lib = NA, fast = T, backfaces = T, mode = "normal", dart_x_off = 0){
  # generate header information
  if(is.na(mtl_lib)){
    mtl_lib <- mtl_name
  }
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()),paste0("mtllib ", mtl_lib, ".mtl") ,paste("usemtl", mtl_name), 'o object1')
  write(obj_head, file)
  #check if its a list or just one polygon
  if(is.list(pol)){
    l <- 0 # index for vertex numbers already used
    lines <- character()
    if(fast){
      options(scipen = 10)
      idx <- 0
      out <- list()
      less_mat <- matrix(c(0,-1,0,0,
                           1,0,0,0,
                           0,0,1,0,
                           0,0,0,1), nrow = 4 , byrow = T)
      #out <- character()
          for(p in pol){
            p <- round(p,5)
            if(any(is.na(p))) next
            if(mode == "dart"){
              p <- p[,c(2,3,1)]
              p[,3] <- - dart_x_off - p[,3]
              p[,1] <- - p[,1]
            }
            if(mode == "less"){
              #p <- p[,c(2,3,1)]
              p <- ((as.matrix(cbind(p,1)) %*% less_mat) |> as.data.frame())[,c(2,3,1)]
              #p[,3] <- - p[,3]
              #p[,1] <- - p[,1]
            }
            
            idx <- idx + 1 # increment polygon index
            v = apply(p, 1, function(x) paste("v", paste(x, collapse = " "))) # add the coordinates round(p$vb,5)
            f = paste("f", paste(round(1:nrow(p)+l,5), collapse = " ")) # add the triangle coordinate order
            if(backfaces){
              f_rev = paste("f", paste(rev(round(1:nrow(p)+l, )), collapse = " ")) # add the reverse order for backsides
              #write(c(v,f,f_rev), file, append = T)
              out <- c(out, list(paste(v, collapse = "\n"),f,f_rev))

            } else {
              out <- c(out, list(paste(v, collapse = "\n"),f))
              #out <- c(out, v,f)
              #write(c(v,f), file, append = T)
            }
            if(idx > 500){
              #write(out, file, append = T)
              data.table::fwrite(out, file, append = T, quote = F, sep = "\n",sep2 = c(""," ",""), row.names = F, col.names = F, nThread = 20)
              out <- list()
              #out <- character()
              idx <- 0
            }
            l <- l+nrow(p) # increment the vertex index 
          }
          if(length(out) > 0){
            #write(out, file, append = T)
            data.table::fwrite(out, file, append = T, quote = F, sep = "\n",sep2 = c(""," ",""), row.names = F, col.names = F, nThread = 20)
          }
          #writeLines(c(obj_head,lines), file) # write to file
      } else {
        for(p in pol){
          v = apply(p, 1, function(x) paste("v", paste(format(x,scientific = F), collapse = " "))) # add the coordinates round(p$vb,5)
          f = paste("f", paste(format(1:nrow(p)+l, scientific = F), collapse = " ")) # add the triangle coordinate order
          if(backfaces){
            f_rev = paste("f", paste(rev(format(1:nrow(p)+l, scientific = F)), collapse = " ")) # add the reverse order for backsides
            #lines <- c(lines, c(v,f,f_rev))# add the lines
            write(c(v,f,f_rev), file, append = T)
          } else {
            l <- l+nrow(p) # increase the vertex index 
            #lines <- c(lines, c(v,f))# add the lines
            write(c(v,f), file, append = T)
          }
          l <- l+nrow(p) # increase the vertex index 
        }
        #writeLines(c(obj_head,lines), file) # write to file
    }
    } else {
      warning("pol needs to be a list with data.frames")
      stop()
    }
}

# function to write polygons without triangulation by rgl
# @param vox a data.frame with center coordinates for each voxel
# @param vox_res the voxel edge length
# @param file the filepath
# @param mtl_name character containing the name for the Material specified in a separate mtl file 
# @param fast if set to TRUE it will use less accurate rounding to speed up the process and create smaller files
# @param backfaces if set to TRUE it will add the reverse order of the triangles to create backfaces
# @param mode if set to "normal" it will write the polygons as they are, if set to "dart" or "less" it will change the order of the dimensions to be used in rtm
# @param dart_x_off if mode id "dart" subtract each x-value from this value to invert the x-axis
write_obj_vox <- function(vox, vox_res = 0.25, file, mtl_name = "Material1",mtl_lib = NA, fast = T, backfaces = T, mode = "normal", dart_x_off = 0){
  # generate header information
  if(is.na(mtl_lib)){
    mtl_lib <- mtl_name
  }
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()),paste0("mtllib ", mtl_lib, ".mtl") ,paste("usemtl", mtl_name), 'o object1')
  write(obj_head, file)
  #check if its a list or just one polygon
  if(is.data.frame(vox)){
    #vox <- vox * vox_res # convert to meters
    # generate the corners of a null cube for the voxels
    vox_corners <- cbind(expand.grid(X = vox_res/2 * c(-1,1), Y = vox_res/2 * c(-1,1), Z = vox_res/2 * c(-1,1)), n = 1:8)
    vox_sides <- rbind(which(vox_corners$X == -vox_res/2), which(vox_corners$X == vox_res/2), which(vox_corners$Y == -vox_res/2), which(vox_corners$Y == vox_res/2), which(vox_corners$Z == -vox_res/2), which(vox_corners$Z == vox_res/2))
    vox_planes <- rep(1:3, each = 2)
    
    l <- 0 # index for vertex numbers already used
      options(scipen = 10)
      idx <- 0
      out <- list()
      less_mat <- matrix(c(0,-1,0,0,
                           1,0,0,0,
                           0,0,1,0,
                           0,0,0,1), nrow = 4 , byrow = T)
      #out <- character()
      for(v in 1:nrow(vox)){
        pol <- list()
        #generate planes for the cube
        for(s in 1:nrow(vox_sides)){
          chull <- chull(vox_corners[vox_sides[s,],-vox_planes[s]])
          pol[[s]] <- data.frame(X = vox_corners[vox_sides[s,chull],1] + as.numeric(vox[v,1]), Y = vox_corners[vox_sides[s,chull],2] + as.numeric(vox[v,2]), Z = vox_corners[vox_sides[s,chull],3] + as.numeric(vox[v,3]))
        }
        out <- list(paste("g voxel",v, sep = "_"))
        idx <- 0
        for(p in pol){
          p <- round(p,5)
          if(any(is.na(p))) next
          if(mode == "dart"){
            p <- p[,c(2,3,1)]
            p[,3] <- - dart_x_off - p[,3]
            p[,1] <- - p[,1]
          }
          if(mode == "less"){
            #p <- p[,c(2,3,1)]
            p <- ((as.matrix(cbind(p,1)) %*% less_mat) |> as.data.frame())[,c(2,3,1)]
            #p[,3] <- - p[,3]
            #p[,1] <- - p[,1]
          }
          
          idx <- idx + 1 # increment polygon index
          v = apply(p, 1, function(x) paste("v", paste(x, collapse = " "))) # add the coordinates round(p$vb,5)
          f = paste("f", paste(round(1:nrow(p)+l,5), collapse = " ")) # add the triangle coordinate order
          if(backfaces){
            f_rev = paste("f", paste(rev(round(1:nrow(p)+l, )), collapse = " ")) # add the reverse order for backsides
            #write(c(v,f,f_rev), file, append = T)
            out <- c(out, list(paste(v, collapse = "\n"),f,f_rev))
            
          } else {
            out <- c(out, list(paste(v, collapse = "\n"),f))
            #out <- c(out, v,f)
            #write(c(v,f), file, append = T)
          }
          
          
          # if(idx > 500){
          #   #
          #   data.table::fwrite(out, file, append = T, quote = F, sep = "\n",sep2 = c(""," ",""), row.names = F, col.names = F, nThread = 20)
          #   out <- list()
          #   #out <- character()
          #   idx <- 0
          # }
          l <- l+nrow(p) # increment the vertex index 
        }
        data.table::fwrite(out, file, append = T, quote = F, sep = "\n",sep2 = c(""," ",""), row.names = F, col.names = F, nThread = 20)
      }
      # if(length(out) > 0){
      #   #write(out, file, append = T)
      #   data.table::fwrite(out, file, append = T, quote = F, sep = "\n",sep2 = c(""," ",""), row.names = F, col.names = F, nThread = 20)
      # }
      #writeLines(c(obj_head,lines), file) # write to file
  } else {
    warning("vox needs to be a data.frames with three collumns")
    stop()
  }
}


#' Generate a .mtl file from a color palette
#'
#' This function generates a .mtl file from a color palette. The .mtl file is a file format that describes the material properties of objects in a 3D scene. The .mtl file is used in conjunction with .obj files to describe the appearance of 3D objects. The colors from the palette will be numbered as c1 to cn_cols in the .mtl file.
#'
#' @param palette A color palette to use for the materials
#' @param filename The name of the .mtl file to generate
#' @param n_cols The number of colors to use from the palette
#'
#' @return A .mtl file with the material properties of the colors in the palette
#' @export generate_mtl_file
#'
#' @examples
#' generate_mtl_file(filename = "clipboard", n_cols = 1)
#' print(readClipboard())
generate_mtl_file <- function(palette = sample(grey.colors(n_cols)), filename = "colors.mtl", n_cols = 10) {
  col2mtl <- function(col) {
    paste(paste(c("Ka","Kd","Ks"),paste(col2rgb(col)/255, collapse = " ")) , collapse = "\n")
  }
  
  cols <- palette[1:n_cols]
  col_df <- cbind(paste0("newmtls c",1:n_cols,"\n"), sapply(cols, col2mtl))
  col_lines <- apply(col_df, 1, function(x) paste(x, collapse = " "))
  
  writeLines(col_lines, filename)
}

#Function to calculate obj mesh objects out of a tree point cloud
#
# @param f filepath to the las file 
# @param class_thresh threshold for the classification of leaves (default = -5.4)
# @param num_cores number of cores for parallel processing (default = detectCores()/2-1)
# @param vox_size voxel size in m for polygon fitting (default = 0.05)
# @param buffer buffer polygons by this size to avoid artificial gaps (default = 0.005)
# @param min_pts_vox minimum number of points per plain fitting voxel (default = 3)
# @param singletons_factor multiplication factor for voxel size to create circles from single points in voxels (if <3pts in voxel) (default = 1/3)
# @param target_dir directory to save the obj files (default = basename(f)), if NA returns a list with the polygons
# @param backfaces if TRUE it will add the reverse order of the triangles to create backfaces (default = T)
# @param mtl_prefix prefix for the material name in the mtl file usually the tree species (default = NA)
# @param mtl_lib name of the mtl library file (default = NA) becomes the same as the mtl_name (prefix + _ + leaf or bark)
# @param write_mode if set to "normal" it will write the polygons as they are, if set to "dart" or "less" it will change the order of the dimensions to be used in rtm
# @param segments number of segments for the single point circle (default = 3)
# @param dart_x_off if write_mode "dart" subtract each x-value from this value to invert the x-axis (default = 0)
tree_mesh <- function(f, class_thresh = -5.4, num_cores = parallel::detectCores()/2-1, vox_size = 0.05, buffer = 0.005, min_pts_vox = 3, singletons_factor = 1/3, segments = 3, target_dir = basename(f), backfaces = T, mtl_prefix = NA, mtl_lib = NA, write_mode = "normal", dart_x_off = 0){
  
  if(lidR::is(f, "LAS")){
    las <- f
    f <- "tree.las"
    } else {
    if(!file.exists(f)){
      stop("file does not exist")
    }
    head <- rlas::read.lasheader(f)
    # check for offset
    if(any(head[c(grep("Max",names(head)), grep("Min",names(head)))] > 9000)) {
      error("pointcloud has a large offset which might lead to numerical problems.")
    }
    las <- lidR::readTLSLAS(f, select = "0")
  }
  # classify leaf/wood (Reflectance < -5.4)
  if(is.numeric(class_thresh)){
    las@data$Classification[las$Reflectance >= class_thresh] <- 0L # wood
    las@data$Classification[las$Reflectance < class_thresh] <- 4L # leaf
  } else if(is.character(class_thresh) && class_thresh %in% names(las@data)) {
    # use FSCT classes
    # 1: Terrain, 2: Vegetation, 3: Coarse woody debris, 4: Stems/branches.
    las@data$Classification[las@data[[class_thresh]] %in% c(3,4)] <- 0L # wood
    las@data$Classification[las@data[[class_thresh]] %in% c(1,2)] <- 4L # leaf
  } else {
    stop("class_thresh needs to be a numeric threshold for Reflectance or a valid collumn in las@data with classes")
  }
  
  
  # use 5cm voxel to assign the most abundant class to all points in voxel
  class_vox <- las |> voxel_metrics(most_abundant(Classification), res = vox_size)
  las <- las |> add_voxel_coordinates(res = vox_size)
  las@data <- merge(las@data, class_vox, by.x = c("x_vox","y_vox","z_vox"), by.y = c("X","Y","Z"))
  
  # filter leaves and calculate voxel coordinates for loop and
  # create a fitting plane convex hull for every 5cm leave voxel
  vegetation <- las |>  filter_poi(V1 == 4)
  bark <- las |>  filter_poi(V1 == 0)
  if(!is.empty(vegetation)){
    vegetation_vox <- vegetation |> voxel_metrics(length(Z), vox_size)
  }
  if(!is.empty(bark)){
    bark_vox <- bark |> voxel_metrics(length(Z), vox_size)
  }
  # Register a parallel backend
  cluster_available <- tryCatch({
    unlist(foreach(i = 1) %dopar% {return(TRUE)})
  }, error = function(e) {
    F
  })
  
  if(!cluster_available){
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
  }
  # parallel process for voxel extraction and convex hull fitting
  if(!is.empty(vegetation)){
    t1 <- system.time({ # measure execution time
      veg_polys <- foreach(z = unique(vegetation_vox$Z), .noexport = c("las"), .export = c('plane_3d','create_circle','random_rotation_matrix'), .combine = "c") %dopar% {
          vox_slice <- vegetation_vox[vegetation_vox$Z == z,]
          veg_slice <- vegetation@data[vegetation@data$z_vox == z,c("X","Y","Z","x_vox","y_vox","z_vox")]
          polys <- list()
          for(v in 1:nrow(vox_slice)){
            vox_bb <- vox_slice[v,]
            vox <- veg_slice[veg_slice$z_vox == vox_bb$Z[1],]
            vox <- vox[vox$x_vox == vox_bb$X[1],]
            vox <- vox[vox$y_vox == vox_bb$Y[1],]
            polys <- c(polys,list(plane_3d(as.matrix(vox[, 1:3]), buffer, F, vox_size, singletons_factor, segments = segments)))
          }
          return(polys)
      }
    })
  }
  #veg_polys <- unlist(veg_polys, recursive = F)
  
  # create a fitting plane convex hull for every 5cm wood voxel
  bark_polys <- list()
  if(!is.empty(bark)){
    
    t2 <- system.time({ # measure execution time
      bark_polys <- foreach(z = unique(bark_vox$Z), .noexport = c("las"), .export = c('plane_3d','create_circle','random_rotation_matrix'), .combine = "c") %dopar% {
        vox_slice <- bark_vox[bark_vox$Z == z,]
        bark_slice <- bark@data[bark@data$z_vox == z,c("X","Y","Z","x_vox","y_vox","z_vox")]
        polys <- list()
        for(v in 1:nrow(vox_slice)){
          vox_bb <- vox_slice[v,]
          vox <- bark_slice[bark_slice$z_vox == vox_bb$Z[1],]
          vox <- vox[vox$x_vox == vox_bb$X[1],]
          vox <- vox[vox$y_vox == vox_bb$Y[1],]
          polys <- c(polys,list(plane_3d(as.matrix(vox[, 1:3]), buffer, F, vox_size, singletons_factor, segments = segments)))
        }
        return(polys)
      }
    })
  }
  # parallel process for voxel extraction and convex hull fitting
  # t2 <- system.time({ # measure execution time
  #   bark_polys <- foreach(v = 1:nrow(bark_vox), .noexport = c("las"), .export = c('plane_3d','create_circle', 'random_rotation_matrix')) %dopar% { # for some reason this fails randomly
  #       vox_bb <- bark_vox[v,]
  #       vox <- bark@data[bark@data$z_vox == vox_bb$Z[1],c("X","Y","Z","x_vox","y_vox","z_vox")]
  #       vox <- vox[vox$x_vox == vox_bb$X[1],]
  #       vox <- vox[vox$y_vox == vox_bb$Y[1],]
  #       plane_3d(as.matrix(vox[, 1:3]), buffer, F, vox_size, singletons_factor, segments = segments)
  #   }
  # })
  # stop the parallel backend
  if(!cluster_available){
    stopCluster(cl)
  }
  
  # generate material string
  if(is.na(mtl_prefix)){
    leaf_mtl <- "leaf"
    bark_mtl <- "bark"
  } else {
    leaf_mtl <- paste0(mtl_prefix, "_leaf")
    bark_mtl <- paste0(mtl_prefix, "_bark")
  }
  
  # check if target dir exists and write obj files
  if(is.na(target_dir)) {
    return(list(leaves = veg_polys, bark = bark_polys))
  } else {
    if(!dir.exists(target_dir)){
      dir.create(target_dir)
    }
    if(!is.empty(vegetation)) write_obj_df(veg_polys, paste0(target_dir,tools::file_path_sans_ext(basename(f)), "_leaves.obj"), mtl_name = leaf_mtl,mtl_lib = mtl_lib, backfaces = backfaces, mode = write_mode, dart_x_off = dart_x_off)
    if(!is.empty(bark)) write_obj_df(bark_polys, paste0(target_dir,tools::file_path_sans_ext(basename(f)), "_bark.obj"), mtl_name = bark_mtl,mtl_lib = mtl_lib, backfaces = backfaces, mode = write_mode, dart_x_off = dart_x_off)
  }
}

#Function to calculate a time series of sun vectors in less standard
#
# @param start/end start and end time of the time series
# @param increment time increment in seconds
# @param lat/lon latitude and longitude of the location
sun_vectors <- function(start = "2023-08-10 10:00:00", end = "2023-08-10 18:00:00", increment = 60, lat = 48.0001, lon = 7.8) {
  # create a sequence of times
  
  times <- seq(as.POSIXct(start), as.POSIXct(end), by = increment)
  # create a sequence of sun vectors
  sun_vecs <- t(sapply(times, function(x) {
    oce_pos <- oce::sunAngle(as.POSIXct(x, tz = "CET"), lon, lat)
    return(c(oce_pos$azimuth, 90 - oce_pos$altitude))
  }))
  return(data.frame(time = times, azimuth = sun_vecs[,1], altitude = sun_vecs[,2]))
}

# function to read the ENVI output off less
#
# @param f filepath to the ENVI output
read_less <- function(f){
  caTools::read.ENVI(f)
}

# function to read the less npy output (requires reticulate and a python environment with numpy installed)
#
# @param f filepath to the npy output
read_less_npy <- function(f){
  np <- reticulate::import("numpy")
  np$load(f)
}

# function to calculate the area covered by a sensor
# 
# @param FOV field of view of the sensor in degrees (single value if FOV_X == FOV_Y or vector of length 2)
# @param D distance distance to the object in m
footprint_area <- function(FOV, D){
  if(length(FOV) == 1){
    FOV <- c(FOV,FOV)
  } 
  return(prod(2*D*tan((FOV*(pi/180))/2)))

}

# function to convert dart to less spectra
#
# @param f filepath to the dart spectra
dart_to_less_spectra <- function(f){
  dart <- read.table(f, sep = "\t", comment.char = "*", header = T)
  less_df <- data.frame(wavelength = dart$wavelength*1000, front_reflectance = dart$reflectance/100, back_reflectance = dart$reflectance/100, transmittance = dart$diffuse_transmittance/100)
  # write results to same folder
  write.table(less_df, paste0(tools::file_path_sans_ext(f), "_less.txt"), row.names = F, col.names = F, sep = "\t")
}

# Function to convert PAR surface irradiance to photon flux
#
# @param irradiance surface irradiance
# @param bw wavelength of the sensor in nm (default = 550)
convert_irradiance_to_photon_flux <- function(irradiance, bw = 550) {
  h <- 6.626e-34  # Planck's constant (J·s)
  c <- 3e8        # Speed of light (m/s)
  avogadro_number <- 6.022e23  # Avogadro's number (photons/mol)
  lambda_avg <- bw * 1e-9  # Average wavelength of PAR in meters (550 nm)

  energy_per_photon <- h * c / lambda_avg  # J/photon
  photon_flux <- irradiance / energy_per_photon  # photons/m^2/s
  photon_flux_umol <- (photon_flux / avogadro_number) * 1e6  # ?mol/m^2/s

  return(photon_flux_umol)
}

# helper function to extract the model iteration integer from batch processing in less
#
# @param filepath path to the file
# @param pattern regex pattern to extract the integer
extract_integer <- function(filepath, pattern = "time_([0-9]+)\\.hdr") {
  match <- stringr::str_match(filepath, pattern)
  return(as.integer(match[2]))
}

# ###############
# ## read and preprocess data 
# ###############

# get start time 
t1 <- Sys.time()
las <- lidR::readTLSLAS(filename)

# apply offset for to avoid numerical problems
means <- apply(las@data,2,mean)
maxs  <- apply(las@data,2,max)
mins <- apply(las@data,2,min)
las@data$X <- las@data$X - mins[1]
las@data$Y <- las@data$Y - mins[2]
las@data$Z <- las@data$Z - mins[3]

x_off <- max(las@data$X)

veg_bark <- tree_mesh(las, target_dir = NA, segments = 5)

write_obj_df(veg_bark$leaves, "V:/RayLeaf/less/leaves.obj", mtl_name = "Malus_domestica_leaf", mtl_lib = "Malus_domestica", backfaces = F, mode = "less")
write_obj_df(veg_bark$bark, "V:/RayLeaf/less/bark.obj", mtl_name = "Malus_domestica_bark", mtl_lib = "Malus_domestica", backfaces = F, mode = "less")
# end 
t2 <- Sys.time()
print(t2-t1)


# generate voxel objs with a list of lad, etc. for turbid media in less
library(AMAPVox)

vox_base <- c(-6, -13.0, -1.0)
mins_vox <- c(-4.815,   -13.006,    -1.513)
vox_res <- 1
vxsp_foliage <- AMAPVox::readVoxelSpace("V:/RayLeaf/less/AMAPvox_out/out1/apfelbaum_1_foliage_pad.vox")

# print summary stats
sum(vxsp_foliage@data$pad > 0.0001)
mean(vxsp_foliage@data$pad[vxsp_foliage@data$pad > 0.0001])
sd(vxsp_foliage@data$pad[vxsp_foliage@data$pad > 0.0001])
sum(vxsp_foliage@data$pad)

vxsp_foliage@data$i <- (vxsp_foliage@data$i)*vox_res  +vox_base[1] - mins_vox[1] + vox_res
vxsp_foliage@data$j <- (vxsp_foliage@data$j)*vox_res  +vox_base[2] - mins_vox[2] + vox_res
vxsp_foliage@data$k <- (vxsp_foliage@data$k)*vox_res  +vox_base[3] - mins_vox[3] 

write_obj_vox(vxsp_foliage@data[,1:3],vox_res, "V:/RayLeaf/less/AMAPvox_out/out1/vox_foliage1.obj", mtl_name = "Malus_domestica_leaf", mtl_lib = "materials", backfaces = F, mode = "less")
write.csv(data.frame(vox = paste0("voxel_vxsp_foliage@data", 1:nrow(vxsp_foliage@data)), pad = vxsp_foliage@data$pad), "V:/RayLeaf/less/AMAPvox_out/out1/vox_foliage1.csv")

vox_res <- 0.2
vxsp_foliage <- AMAPVox::readVoxelSpace("V:/RayLeaf/less/AMAPvox_out/out1/apfelbaum_02_foliage_pad.vox")
# print summary stats
sum(vxsp_foliage@data$pad > 0.0001)
mean(vxsp_foliage@data$pad[vxsp_foliage@data$pad > 0.0001])
sd(vxsp_foliage@data$pad[vxsp_foliage@data$pad > 0.0001])
sum(vxsp_foliage@data$pad)*(0.2^3)


vxsp_foliage@data$i <- (vxsp_foliage@data$i)*vox_res  +vox_base[1] - mins_vox[1] + vox_res
vxsp_foliage@data$j <- (vxsp_foliage@data$j)*vox_res  +vox_base[2] - mins_vox[2] + vox_res
vxsp_foliage@data$k <- (vxsp_foliage@data$k)*vox_res  +vox_base[3] - mins_vox[3] 

write_obj_vox(vxsp_foliage@data[,1:3],vox_res, "V:/RayLeaf/less/AMAPvox_out/out1/vox_foliage02.obj", mtl_name = "Malus_domestica_leaf", mtl_lib = "materials", backfaces = F, mode = "less")
write.csv(data.frame(vox = paste0("voxel_vxsp_foliage@data", 1:nrow(vxsp_foliage@data)), pad = vxsp_foliage@data$pad), "V:/RayLeaf/less/AMAPvox_out/out1/vox_foliage02.csv")
