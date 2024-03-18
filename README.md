# SeuratCropper
Crops a Seurat Object using a sf polygon. Only cells and molecules inside the polygon are retained. \
If no cells remain after the cropping the FOV is returned unchanged with a warning.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Depends on:
- `Seurat`
- `sf`
- `data.table`

Recommended:
- `geojsonsf`

## Usage

Setup:
```
library(Seurat)
library(data.table)
library(sf)
library(geojsonsf)
```

Read in a polygon from a geojson file:
```
Polygon_A <- geojson_sf("cropping/B.geojson")
Polygon_B <- geojson_sf("cropping/A.geojson")
```

Subset the Seurat Objecet:
```
subsetted_object <- PolygonCropSeurat(object, list(FOV_A = Polygon_A, FOV_B = Polygon_B))
```

## TODO
To do:
- Setup so cell polygons can be used instead of centroid.
- Fail with an Error if no centroids are available


