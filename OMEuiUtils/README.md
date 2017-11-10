OME ui utils

Copied from the minimal Maven project for connecting to OMERO.

This repository has been created to hold Graphical utilities intended for use by OMERO clients and by tools that access files using bio-formats. 

Currently it holds OMEROImageChooser which graphically allows the selection of one or more object from an OMERO server with a 'look and feel' modelled on the OMERO.insight client.

This tool can be used from either Matlab or Java.

To select a single object, the constructor has the form: OMEROImageChooser(omero.client omeroclient, int selectedType )  
selectedType determines which type of object can be selected, 0 = image,1 = dataset, 2 = plate.

Other constructors allow the selection of multiple objects (currently only enabled for when selectedType == image)

To get the selected object/s call the appropriate routine from the following:

  Image[] getSelectedImages()  
 
  Dataset getSelectedDataset()  
  
  Plate getSelectedPlate()  
      
  
