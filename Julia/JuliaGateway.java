/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author yalexand
 */




import loci.formats.Modulo;

//import ome.xml.model.enums.handlers.DimensionOrderEnumHandler;
//import ome.xml.model.enums.handlers.PixelTypeEnumHandler;
//
//import ome.xml.model.enums.Enumeration;
//import ome.xml.model.enums.PixelType;
//import ome.xml.model.enums.DimensionOrder;
import ome.xml.model.primitives.PositiveInteger;
//
import loci.formats.ImageWriter;
import loci.formats.IFormatWriter;
import loci.formats.out.TiffWriter;
//
//import loci.formats.ChannelFiller;
//import loci.formats.ChannelSeparator;
//import loci.formats.services.OMEXMLServiceImpl;
import loci.formats.ome.OMEXMLMetadata;
//import loci.formats.services.OMEXMLService;
//import loci.formats.meta.MetadataStore;
//import loci.formats.FormatReader;
//import loci.formats.IFormatReader;
import loci.common.DebugTools;
//import loci.formats.Modulo;
import loci.common.DataTools;
//import loci.formats.FormatTools;

public class JuliaGateway 
{    
    public static String getParentDimension(Modulo modlo)  {    return modlo.parentDimension; }
    public static double getStart(Modulo modlo)            {    return modlo.start; }
    public static double getStep(Modulo modlo)             {    return modlo.step; }
    public static double getEnd(Modulo modlo)              {    return modlo.end; }
    public static String getParentType(Modulo modlo)       {    return modlo.parentType; }
    public static String getType(Modulo modlo)             {    return modlo.type; }
    public static String getTypeDescription(Modulo modlo)  {    return modlo.typeDescription; }
    public static String getUnit(Modulo modlo)             {    return modlo.unit; }
    public static String[] getLabels(Modulo modlo)         {    return modlo.labels; }
    
    public static void  setParentDimension(Modulo modlo,String rhs)  {    modlo.parentDimension = rhs; }
    public static void  setStart(Modulo modlo,double rhs)            {    modlo.start = rhs;  }
    public static void  setStep(Modulo modlo,double rhs)             {    modlo.step = rhs; }
    public static void  setEnd(Modulo modlo,double rhs)              {    modlo.end = rhs; }
    public static void  setParentType(Modulo modlo,String rhs)       {    modlo.parentType = rhs; }
    public static void  setType(Modulo modlo,String rhs)             {    modlo.type = rhs; }
    public static void  setTypeDescription(Modulo modlo,String rhs)  {    modlo.typeDescription = rhs; }
    public static void  setUnit(Modulo modlo,String rhs)             {    modlo.unit = rhs; }
    public static void  setLabels(Modulo modlo,String[] rhs)         {    modlo.labels = rhs; }    
    
    public static void  setPixelsSizeX(OMEXMLMetadata metadata, int rhs) { metadata.setPixelsSizeX(new PositiveInteger(rhs),0); }
    public static void  setPixelsSizeY(OMEXMLMetadata metadata, int rhs) { metadata.setPixelsSizeY(new PositiveInteger(rhs),0); }
    public static void  setPixelsSizeZ(OMEXMLMetadata metadata, int rhs) { metadata.setPixelsSizeZ(new PositiveInteger(rhs),0); }
    public static void  setPixelsSizeC(OMEXMLMetadata metadata, int rhs) { metadata.setPixelsSizeC(new PositiveInteger(rhs),0); }
    public static void  setPixelsSizeT(OMEXMLMetadata metadata, int rhs) { metadata.setPixelsSizeT(new PositiveInteger(rhs),0); }
    
    public static void  setChannelSamplesPerPixel(OMEXMLMetadata metadata, int nchannels, int index)
    { metadata.setChannelSamplesPerPixel(new PositiveInteger(nchannels), 0, index); }
    
    public static void  setPixelsBinDataBigEndian(OMEXMLMetadata metadata, boolean flag) { metadata.setPixelsBinDataBigEndian(flag, 0, 0); }
    
    public static void  setMetadataRetrieve(ImageWriter writer, OMEXMLMetadata metadata)    { writer.setMetadataRetrieve(metadata); }
    
    public static int getPixelsSizeX(OMEXMLMetadata metadata)   { return metadata.getPixelsSizeX(0).getValue(); }
    public static int getPixelsSizeY(OMEXMLMetadata metadata)   { return metadata.getPixelsSizeY(0).getValue(); }
    public static int getPixelsSizeZ(OMEXMLMetadata metadata)   { return metadata.getPixelsSizeZ(0).getValue(); }
    public static int getPixelsSizeC(OMEXMLMetadata metadata)   { return metadata.getPixelsSizeC(0).getValue(); }
    public static int getPixelsSizeT(OMEXMLMetadata metadata)   { return metadata.getPixelsSizeT(0).getValue(); }
    
    public static byte[] shortsToBytes(short[] arr, boolean little) { return loci.common.DataTools.shortsToBytes(arr,little); }        
    public static byte[] intsToBytes(int[] arr, boolean little) { return loci.common.DataTools.intsToBytes(arr,little); }
    public static byte[] floatsToBytes(float[] arr, boolean little) { return loci.common.DataTools.floatsToBytes(arr,little); }
    public static byte[] doublesToBytes(double[] arr, boolean little) { return loci.common.DataTools.doublesToBytes(arr,little); }
    
    public static void enableLogging(String spec) { loci.common.DebugTools.enableLogging(spec); }
    
}












