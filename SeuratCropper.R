library(Seurat)
library(data.table)
library(sf)
library(geojsonsf)

PolygonCropFOV <- function(fov, polygon, pyxel_size = 0.2125){
  #' Crop an FOV using a polygon
  #' 
  #' @description This function crops an FOV including centroids, segmentations and molecules.
  #'
  #' The FOV must contain centroids as those are used to determine if the cells are in the polygon or not. 
  #
  #' @param fov fov Object The FOV to be cropped.
  #' @param polygon sf POLYGON object.
  #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
  #' @usage PolygonCropFOV(polygon, pyxel_size = 0.2125)
  #' @return A new cropped FOV object.
  #' @details If no cells are found the unchanged fov is returned with a Warning.
  #' @note TODO: if no centroids are present try to use segmentation or fail with an Error.
    
    # Use the centroids to identify the cells that should be kept
    message("Start processing centroids...")
    polygon <- st_set_crs(polygon, NA)
    centroids <- st_as_sfc(fov$centroids, forceMulti = FALSE)
    centroids <- st_set_crs(centroids, NA)
    centroids <- st_combine(centroids)
    centroids <- centroids / pyxel_size
    centroids <- st_cast(centroids, "POINT")
    contained <- st_contains_properly(polygon, centroids, sparse = FALSE)
    cells_to_keep <- Cells(fov)[contained]
    message("Done processing centroids...")
    
    # If no cells are found within the polygon, return the FOV unchanged with a warning
    if (!any(contained)){
        warning("No cells, found inside the polygon. Returnin g original object")
        return(fov)
    }
    
    # Build a new molecule object with only the molecules inside the polygon
    if (!is.null(fov$molecule)){
        message("Start processing molecules...")
        molecules <- rbindlist(lapply(fov$molecule, function(x){
            as.data.table(x@coords)
        }), use.names = T, idcol = "gene")
        mols <- st_multipoint(as.matrix(molecules[, .(x, y)]))
        mols <- mols / pyxel_size
        mols <- st_cast(st_sfc(mols), "POINT")
        mols <- st_set_crs(mols, NA)
        contained <- st_contains_properly(polygon, mols, sparse = FALSE)
        molecules <- molecules[as.logical(contained)]
        molecules <- CreateMolecules(as.data.frame(molecules))
        message("Done processing molecules...")
    }
    
    # Subset the input FOV to create a new FOV
    message("Preparing new fov...")
    new_fov <- subset(x = fov, cells = cells_to_keep) # This subsets the centroids and segmentation
    if (!is.null(fov$molecule)){
         new_fov[["molecules"]] <- molecules # Replace the molecules with the cropped ones
         new_fov[["molecule"]] <- molecules
    }
    message("Done cropping fov.")
    
    return(new_fov)
}

PolygonCropSeurat <- function(object, polygons, pyxel_size = 0.2125){
  #' Crop a Seurat Object using a list of Polygons
  #' 
  #' @description This function crops all FOVs included in the list, other FOVs included unchanged
  #' 
  #' The FOVs to be cropped must contain centroids as those are used to determine if the cells are in the polygon or not. 
  #
  #' @param object Seurat object to be cropped
  #' @param polygons a named list of sf POLYGON object, the names must match those of the FOVs to be cropped.
  #' @param pixel_size ratio between the FOV coordinates and the sf Objects coordinates.
  #' @usage PolygonCropSeurat(object polygons, pyxel_size = 0.2125)
  #' @return A new Seurat object, including the retained cells for the cropped FOVs and all cells for the other FOVs.
  #' @details It relies on `PolygonCropFOV` so if a cropped FOV would contain no cells, the uncropped FOV is included instead
    
    # Crop all FOVs for which a polygon has been provided
    message("Start cropping fovs...")
    cropped_fovs <- lapply(names(polygons), function(name){
        fov <- object@images[[name]]
        fov <- PolygonCropFOV(fov, polygons[[name]], pyxel_size)
        message("Done: ", name)
        return(fov)
    })
    names(cropped_fovs) <- names(polygons)
    message("Done cropping fovs...")

    # Make a list of all FOVs, cropped and non-cropped
    message("Gathering fovs...")
    untouched_fovs <- object@images[!(names(object@images) %in% names(cropped_fovs))]
    all_fovs <- c(untouched_fovs, cropped_fovs)
    message("Found fovs: ", names(all_fovs))

    # Get the cells to keep and make a subset of the object
    message("Making new object fovs...")
    cells_to_keep <- unlist(sapply(all_fovs, Cells))
    new_object <- subset(object, cells = cells_to_keep)
    new_object@images <- all_fovs # Replace the old FOVs
    message("Done cropping Seurat object")

    return(new_object)
}
