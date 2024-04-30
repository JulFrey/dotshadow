# dotshadow

This R script takes a point cloud of vegetation from terrestrail laser scanning as input and computes polygons represensiting the leaves and woody components. Therefore the point cloud is split into voxel cells and for every cell an optimal fitting plane is calculated based on the points in the voxel. The points itself get projected to the plane and a convex hull is fit to the planar set of points. Th econvex hull is buffered to the theoretical distance between to points on a plane due to th edownsampling distance of the point cloud, to avoid gaps between the polygons. 

How to cite:
