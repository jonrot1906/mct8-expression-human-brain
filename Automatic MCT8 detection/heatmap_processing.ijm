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

SourceDir = getDirectory("Set a Directory");

TargetDir = getDirectory("Choose Destination Directory ");

masterFileList = getFileList(SourceDir);

// Takes a directory containing sub-directories, loads them in sorted order and does an operation on them
for(j = 0; j < masterFileList.length; j++) {
     showProgress(j+1, masterFileList.length);
     subdir = SourceDir + masterFileList[j];
     createOverlayHeatmap(subdir);
}

function createOverlayHeatmap(subdir) {
    folder_name = File.getName(subdir);
    // Get the name of the selected folder

    // Get list of files in the selected folder
    list = getFileList(subdir);

    open(subdir + folder_name + "_TissueDetection.tif");
    tissueDetectionImage = getTitle();

    // Loop through each file in the list
    for (i = 0; i < list.length; i++) {
        // Get the full path of the current file
        file_path = subdir + list[i];

        // Check if the file contains "dmap" in its name
        if (indexOf(list[i], "dmap") >= 0) {
            // Open the current image
            open(subdir + list[i]);
            heatmapImage = getTitle();

            // Apply the "Fire" LUT
            run("Fire");

            // Multiply ImageCalculator with "TissueDetection"
            imageCalculator("Multiply create", heatmapImage, tissueDetectionImage);

            // Optionally, you can save the processed image
            // Uncomment the line below to save each processed image
            // saveAs("Tiff", output_folder + "processed_" + File.nameWithoutExtension);

            // Optionally, you can close the original image
            // Uncomment the line below to close each original image after processing
            // close();
        }
    }

    open(subdir + folder_name + ".tif");
    originalImage = getTitle();
    imageCalculator("Multiply create stack", originalImage, tissueDetectionImage);
    selectImage("Result of " + originalImage);
    tissueImage = getTitle();
    run("Stack to RGB");
    selectImage(tissueImage + " (RGB)");
    RGBimage = getTitle();
    run("Enhance Contrast", "saturated=0.25");
    //create good looking originalImage with tissueBorders
    close(tissueDetectionImage);
    outputFolder=TargetDir + File.separator + folder_name;
    File.makeDirectory(outputFolder);

    //loop through Heatmap images to create heatmap-original overlays
    for (i = 0; i < list.length; i++) {
        if (indexOf(list[i], "dmap") >= 0) {
            // Open the current "dmap" image
            iteratingImage = "Result of " + list[i];
            selectImage(tissueImage + " (RGB)");
            run("Duplicate...", " ");
            selectImage(getTitle());
            // Overlay the "dmap" image on the directory-named image with 50% transparency
            run("Add Image...", "image=&iteratingImage x=0 y=0 opacity=50");
            saveAs("PNG", outputFolder + "/" + list[i]);
        }
    }
    selectImage(RGBimage);
    saveAs("PNG", outputFolder + "/" + folder_name);
    
    close("*");

}