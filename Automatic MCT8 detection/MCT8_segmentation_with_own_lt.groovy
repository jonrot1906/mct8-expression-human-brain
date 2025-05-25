//MIT License

//Copyright (c) 2025 Jonas Rotter

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.



import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.gui.commands.Commands;
import qupath.lib.roi.RoiTools.CombineOp;
import qupath.lib.objects.PathAnnotationObject
import qupath.lib.roi.RectangleROI
import qupath.lib.scripting.QP

def project = QP.getProject()

// Function to create an annotation around a cell
def createAnnotationAroundCell(cellObject) {
    print cellObject
    def ml = cellObject.getMeasurementList()
    def roi = cellObject.getROI()
    def halfWidth = 100.0 // Modify this value to adjust the annotation size

    // Calculate the coordinates of the bounding box
    def x = roi.getCentroidX() - halfWidth
    def y = roi.getCentroidY() - halfWidth
    def width = halfWidth * 2
    def height = halfWidth * 2
    def rectangles = []
    // Create a rectangle ROI representing the bounding box
    def rectangleROI = new RectangleROI(x, y, width, height)
    rectangles << rectangleROI
    // Create a new annotation object
    def mct8Backgroundclass = getPathClass('MCT8-Background')    
    def backgroundAnnotations = rectangles.collect {PathObjects.createAnnotationObject(it)}
    for (backgroundAnnotation in backgroundAnnotations) {
        backgroundAnnotation.setPathClass(mct8Backgroundclass)
    }
    cellObject.addPathObjects(backgroundAnnotations)
    fireHierarchyUpdate()
}

// Function to measure intensity within each annotation
def measureIntensityInAnnotations() {
    selectObjectsByClassification("MCT8-Background");
    runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0, "channel1": false, "channel2": false, "channel3": true,  "doMean": true,  "doStdDev": true,  "doMinMax": true,  "doMedian": false,  "doHaralick": false,  "haralickDistance": 1,  "haralickBins": 32}');
}

for (entry in project.getImageList()) {
    def name = entry.getImageName()
    print name

setImageType('FLUORESCENCE');
setChannelNames('DAPI', 'NeuN', 'MCT8');

//tresholding the tissue
int channel = 2 // 0-based index for the channel to threshold
double threshold = 1000 // Threshold value
int level = 2 // 0-based resolution level for the image pyramid (choosing 0 may be slow)
def belowClass = getPathClass('Ignore*') // Class for pixels below the threshold
def aboveClass = getPathClass('MCT8') // Class for pixels above the threshold

// Create a single-resolution server at the desired level, if required
def server = getCurrentServer()
if (level != 0) {
  server = qupath.lib.images.servers.ImageServers.pyramidalize(server, server.getDownsampleForResolution(level))
}

// Create a thresholded image
def thresholdServer = PixelClassifierTools.createThresholdServer(server, channel, threshold, belowClass, aboveClass)

// Create annotations and add to the current object hierarchy
def hierarchy = getCurrentHierarchy()
PixelClassifierTools.createAnnotationsFromPixelClassifier(hierarchy, thresholdServer, 5000, 500)

println("Tissue detection done!")

//Get single object of each class
ignore_annotations = getAnnotationObjects().findAll {
  p -> p.getPathClass() == getPathClass("Ignore*")
}
mct8_annotation = getAnnotationObjects().findAll {
  p -> p.getPathClass() == getPathClass("MCT8")
} [0]

for (ignore_annotation in ignore_annotations) {
  getCurrentHierarchy().getSelectionModel().setSelectedObject(ignore_annotation, true);
  getCurrentHierarchy().getSelectionModel().setSelectedObject(mct8_annotation, true);
  Commands.combineSelectedAnnotations(getCurrentImageData(), CombineOp.SUBTRACT);
  fireHierarchyUpdate()
}

removeObjects(ignore_annotations, true)
selectAnnotations();

runPlugin('qupath.lib.algorithms.TilerPlugin', '{"tileSizeMicrons": 1500, "trimToROI": true, "makeAnnotations": true, "removeParentAnnotation": true}');

println("Tissue tiling done!")


//define StarDist model

def modelPath = "dsb2018_heavy_augment.pb"
def dapi = getPathClass('DAPI')
def neun = getPathClass('NeuN')


  hierarchy = getCurrentHierarchy()
def annotations = getAnnotationObjects()

//loop through tiles and do cell classification
for (annotation in annotations) {
    print(annotation)
    hierarchy.getSelectionModel().clearSelection();
    selectObjects{p -> p == annotation}
  
  def stardistDAPI = StarDist2D
    .builder(modelPath)
    .channels('DAPI') // Extract channel called 'DAPI'
    .normalizePercentiles(1, 99) // Percentile normalization
    .threshold(0.6) // Probability (detection) threshold
    .pixelSize(0.5) // Resolution for detection
    .cellExpansion(5) // Expand nuclei to approximate cell boundaries
    .measureShape() // Add shape measurements
    .measureIntensity() // Add cell measurements (in all compartments)
    .build()

  def imageData = QP.getCurrentImageData()
  def selectedObjects = QP.getSelectedObjects()
  // Detect objects using StarDist for DAPI channel
  stardistDAPI.detectObjects(imageData, selectedObjects)

  var dapiCells = annotation.getChildObjects().findAll {
    d -> d.isDetection()
  }

  if (dapiCells.size() > 0) {

    def average = dapiCells.stream().mapToDouble(cell -> cell.getMeasurementList().getMeasurementValue("NeuN: Nucleus: Mean")).average().get()
    //def channelNames = getCurrentServer().getMetadata().getChannels().collect{ c -> c.name }
    //println channelNames

    dapiCells.each {
      cell ->
        def dapiSignal = cell.getMeasurementList().getMeasurementValue("DAPI: Nucleus: Mean")
      def neuNSignal = cell.getMeasurementList().getMeasurementValue("NeuN: Nucleus: Mean")

      if (neuNSignal > average) {
        cell.setPathClass(neun)
      } else {
        cell.setPathClass(dapi)
      }
    }

    print("DAPI/NeuN classification done")
    resetSelection()
  }
  stardistDAPI.close()

def tile_detections = annotation.getChildObjects().findAll{d-> d.isDetection()}

// Create an annotation around each cell
for (tile_detection in tile_detections) {
    createAnnotationAroundCell(tile_detection)
}
//cellObjects.each { cellObject ->
//    createAnnotationAroundCell(cellObject)
//}

// Measure intensity within each annotation
measureIntensityInAnnotations()

String mct8BackgroundMean = "ROI: 2.00 µm per pixel: MCT8: Mean"
String mct8BackgroundDev = "ROI: 2.00 µm per pixel: MCT8: Std.dev."

//loop through all cells and perform subcellular detection
for (tile_detection in tile_detections) {
    hierarchy = getCurrentHierarchy()
    def cell_annotations = getAnnotationObjects()
    def mct8Background = "MCT8-Background"
    def mct8RoiAverage = cell_annotations.findAll{it -> it.getPathClass() == getPathClass(mct8Background)}.stream().mapToDouble(mct8RoiAverageAnnotation -> mct8RoiAverageAnnotation.getMeasurementList().getMeasurementValue(mct8BackgroundMean) + (mct8RoiAverageAnnotation.getMeasurementList().getMeasurementValue(mct8BackgroundDev)*4)).average().get()
    selectObjects{p -> p == tile_detection}
    runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1.0,  "detection[Channel 2]": -1,  "detection[Channel 3]": ' + mct8RoiAverage + ',  "doSmoothing": true,  "splitByIntensity": true,  "splitByShape": true,  "spotSizeMicrons": 5,  "minSpotSizeMicrons": 1,  "maxSpotSizeMicrons": 10.0,  "includeClusters": true}');
    resetSelection();
    fireHierarchyUpdate()
}

//end with MCT8 classification
}

//proceed with subcellular detection

// Loop through all available cells
for (cell in getCellObjects()) {
    // Determine the current classification, with any intensity component removed
    def pathClass = cell.getPathClass()
    // Determine how many single spots and clusters we have
    double numSingleSpots = measurement(cell, "Subcellular: Channel 3: Num single spots")
    print numSingleSpots
    double numClusters = measurement(cell, "Subcellular: Channel 3: Num clusters")
    print numClusters
    // Apply the new classification
    if (numClusters >= 2 || numSingleSpots >= 5)
        cell.setPathClass(getDerivedPathClass(pathClass, "MCT8-Positive"))
    else
        cell.setPathClass(getDerivedPathClass(pathClass, "MCT8-Negative"))
}
// Ensure the GUI is updated
fireHierarchyUpdate()

//export data

Set annotationMeasurements = []

getDetectionObjects().each{it.getMeasurementList().getMeasurementNames().each{annotationMeasurements << it}}
//println(annotationMeasurements)

annotationMeasurements.each{ removeMeasurements(qupath.lib.objects.PathCellObject, it);}

// write to file
boolean prettyPrint = false // false results in smaller file sizes and thus faster loading times, at the cost of nice formating
def gson = GsonTools.getInstance(prettyPrint)
def all_annotations = getDetectionObjects()
//println gson.toJson(annotations) // you can check here but this will be HUGE and take a long time to parse


// automatic output filename, otherwise set explicitly
String imageLocation = getCurrentImageData().getServer().getPath()
outfname = imageLocation.split("file:/")[1]+".json"



File file = new File(outfname)
file.withWriter('UTF-8') {
    gson.toJson(all_annotations,it)
}

}