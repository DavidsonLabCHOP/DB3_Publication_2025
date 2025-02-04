// Set Image Type
setImageType('FLUORESCENCE');

// Part for changing the color
setChannelColors(
getColorRGB(0, 0, 255), // Blue
getColorRGB(0, 255, 0), // Green
getColorRGB(255, 0, 0)) // Red



// Clear existing objects
selectDetections()
clearDetections()
selectAnnotations()
clearSelectedObjects()
// Get all cell objects from the current image
def cells = getCellObjects()

// Remove all cell objects
removeObjects(cells, true)

// Update the hierarchy to reflect changes
fireHierarchyUpdate()

// Part for renaming channels
setChannelNames(
'Blue',
'Green',
'Red')


//Select the all annotations
createFullImageAnnotation(true)
selectAnnotations()

// Detect cells using the Red channel
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"Blue","backgroundByReconstruction":true,"backgroundRadius":15.0,"medianRadius":0.0,"sigma":10.0,"minArea":50.0,"maxArea":10000.0,"threshold":5.0,"watershedPostProcess":true,"cellExpansion":5.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')

// Run subcellular puncta detection
selectCells();
runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]":-1.0,"detection[Channel 2]":25.0,"detection[Channel 3]":25.0,"detection[Channel 4]":25.0,"doSmoothing":false,"splitByIntensity":true,"splitByShape":true,"spotSizeMicrons":1.0,"minSpotSizeMicrons":0.05,"maxSpotSizeMicrons":20.0,"includeClusters":true}')

// Compute Intensity Features
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"downsample":1.0,"region":"ROI","tileSizePixels":200.0,"channel1":false,"channel2":true,"channel3":true,"doMean":true,"doStdDev":true,"doMinMax":true,"doMedian":true,"doHaralick":false,"haralickMin":NaN,"haralickMax":NaN,"haralickDistance":1,"haralickBins":32}')