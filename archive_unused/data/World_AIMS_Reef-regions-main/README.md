# Global reef regions
This dataset represents broad boundaries around all global reefs with a buffer of 
approximately 30 km around reefs. This boundary map is intended to capture all
shallow regions where tropical coral reefs might occur and is intended to assist with
determining areas where satellite imagery needs to be processes to map global reefs.
The features are split into regions with the purpose of organising global reef mapping 
into managable areas. 

This dataset contains two shapefiles:
1. Global reef regions (`World_AIMS_Reef-regions_V1.shp`) - broad boundary around reefs from -30 degrees to +30 degree latitude.
2. Sentinel 2 tiles that intersect with each region (`World_AIMS_Reef-regions_Sentinel2_V1.shp`)

This dataset was derived from on the regionalisation and identification of reef locations
performed by the Allen Coral Atlas. 

![Preview map of reef regions dataset](media/Reef-regions-preview-map.png)
Preview of the Global reef regions dataset (`World_AIMS_Reef-regions_V1.shp`).

![Preview map of sentinel 2 tiles](media/Reef-regions-sentinel2-preview-map.png)
Preview of the intersecting Senintel 2 tiles (`World_AIMS_Reef-regions_Sentinel2_V1.shp`).

![Detailed map of sentinel 2 tiles](media/Reef-regions-sentinel-2-detailed.jpg)
Map showing the reef regions and labelled Sentinel 2 tiles. Note this image is large.

## Change log
## 2022-02-10 - V1 - Initial release of the dataset. 
This dataset was developed sufficiently to allow the organisation of files for the 
[Coral Sea Mapping project](https://eatlas.org.au/projects-other/coral-sea-reef-mapping). As such
there has been a greater level of adjustment and optimisation of the boundaries in the Australia
regions.
The regional boundaries were heavily based on the Allen Atlas boundaries, but with greater division 
around australia. The original Timor and Arafura sea region was split in two, and the Timor Sea region extended
to cover the banks and shoals along the edge of the continental shelf of North West Australia. To minimise
the number of Sentinel 2 tiles that need to be processed the boundaries were manually reviewed to look for
areas where the boundary just touched a Sentinel 2 tile. In these cases these tiles would simply be open
water. The boundaries were then adjusted in tighter around the reefs. The main focus of this work was
in the Coral Sea region, reducing the number of tiles by 8. Little review has been undertaken for 
most of the non-Australian regions and so these should be considered as draft quality.


## Method 

The boundaries were created by roughly digitising the outer boundary (plus a buffer) of
the Allen Coral Atlas satellite imagery. These boundaries were then sliced and clustered
to align with the existing named "Mapped Areas" of the Allen Coral Atlas. The regions were
remapped from the satellite imagery as the "Mapped Areas" did cover all the satellite imagery.

### 1. Adding Allen Atlas Satellite imagery to QGIS
The Allen Atlas satellite imagery was added to QGIS by creating two XYZ Tiles layers with the 
following settings:
Layer 1:
- *Title:* Allen Coral Atlas 0 - 12
- *URL:* https://storage.googleapis.com/production-coral-tiles/visual/coral_reefs_2020_visual_v1_mosaic/{z}/{x}/{y}.png
- *Min Zoom Level:* 0
- *Max Zoom Level:* 12

Layer 2:
- *Title:* Allen Coral Atlas 13 - 15
- *URL:* https://allencoralatlas.org/tiles/planet/visual/2020/{z}/{x}/{y}
- *Min Zoom Level:* 13
- *Max Zoom Level:* 15

Note that this method of connecting to the imagery is not an official method for accessing the 
Allen Coral Atlas imagery. As a result these URL are likely to break. These URL paths were
obtained using the web browser developer tools to inspect the URLs requests on the 
[mapping page](https://www.allencoralatlas.org/atlas/).

### 2. Exporting the satellite imagery as a single raster
To get the extent of the satellite imagery we first export from the tiles to a single raster,
then convert it to a polygon. 
1. Turn on the Allen COral Atlas 0 - 12
2. Create a `New Print Layout`. Then adjust the page size to be 400 x 400 mm to match the shape of the
world in WGS84 / Pseudo-Mercator. 
3. Add a map to our page, ensuring that it covers the entire page. 
4. Set the Extents to be -20026376.39 -20048966.10, 20026376.39, 20048966.10. This is from 
[EPSG:3857 Projected Bounds](https://epsg.io/3857)
4. Make sure that the map shows the entire world. The satellite imagery should cover a strip in the middle.
5. Export as an image (`02-World-Allen-Coral-Atlas-Satellite-2020.png`), ensuring that 'Generate world file' is turned on. 
Set the export resolution to 8096 x 8096 pixels. 
6. Add the newly created image back into QGIS as a layer, then export it as a GeoTiff. This will embed the
projection information into the image, allowing us to perform reprojections on the data.

### 3. Convert the satellite imagery to a polygon boundary.
Before we convert to a polygon we need to turn the satellite imagery into a binary mask.
1. Use the `Raster\Raster Calculator`, with expression 
 ```
 "02-World-Allen-Coral-Atlas-Satellite-2020@1" + 
 "02-World-Allen-Coral-Atlas-Satellite-2020@2" + 
 "02-World-Allen-Coral-Atlas-Satellite-2020@3"<(255+255+255)
 ```

 This is intended to mask all areas other than white.
 Set the Output CRS to EPSG:3857 and save it as a GeoTiff, `03-World-satellite-mask`.
 Note: This generates a 256 MB image and so is not included in the repository.

2. Export the `03-World-satellite-mask.tif` as `03-World-satellite-mask-LZW.tif`, setting
the 'Create Options' to COMPRESS, LZW. This shrinks the file 100x.

3. Use the `Raster\Conversion\Polygonize`, with default values, to create a matching polygon. 
Save to a temporary file.

4. On the new vector layer set the Query Builder to `"DN" = 1` to remove the unwanted areas.

5. Expand the boundaries by applying a 20 km buffer, followed by a dissolve, followed by -5 km buffer.

6. Apply a Simplify with 'Distance (Douglas-Peucker)' with a tolerance of 5 km. Save as `03-World-reef-mask.shp`.

### 4. Adjust the mask to include new reefal areas not covered by existing boundaries
The Allen Coral Atlas doesn't include all areas where reefs exists. As we process these areas we need to expand
the mask. To do this we used the Australia Bathymetry and Topography (GA, 2009) to estimate where additional
shoals around Australia might be found. We also adjusted the smoothness of the regions.

### 5. Split and name regions
Use the Multipart to Singleparts tool to separate all the polygons into separate features, saving to 
`05-World-reef-regions.shp`. 
1. Where multiple regions were joined together, such as along the coast of Australia, the features were
split. The boundaries used closed followed those used by the Allen Coral Atlas and used the boundaries
of the seas and oceans as additional guides (IHO 23-3rd: Limits of Oceans and Seas).

 Allen Coral Atlas (2022). Imagery, maps and monitoring of the world's tropical coral reefs. doi.org/10.5281/zenodo.3833242
IHO 23-3rd: Limits of Oceans and Seas, Special Publication 23, 3rd Edition 1953, published by the 
International Hydrographic Organization. (https://epic.awi.de/id/eprint/29772/1/IHO1953a.pdf)

2. Add `Region` and `RegionName` attributes to the shapefile. Copy the region names from the 
Coral Allen Atlas for where there are matches. Additional names were allocated for regions not
covered by the exiting Allen Atlas regions. 

3. Export the final edited shapefile as `finaldata\World_AIMS_Reef-regions_V1.shp`. Set the CRS to EPSG:4326.

### 6. Merge the Reef Regions with the Sentinel 2 regions
1. Load `src-data\Sentinel-2-Shapefile-Index-master\derived\sentinel_2_index_shapefile_33S-33N-extract.shp`.
2. Use `Vector\Data Management Tools\Join Attributes by Location ...` with 
 *Base Layer:* `sentinel_2_index_shapefile_33S-33N-extract.shp`
 *Join Layer:* `World_AIMS_Reef-regions_V1.shp`
 *Geometric predicate:* intersects
 *Join type:* Create separate feature for each marching feature (on-to-many)
Create temporary layer
3. For the new layer set the Query builder to remove all the non matches (`"Region"`).
4. Export to `finaldata\World_AIMS_Reef-regions_Sentinel2_V1.shp`


## Data Dictionary
`Region`: Full name of the region.

`RegionName`: Folder name to store files associated with this region. This name is of limited length
to limit problems with file path length.

