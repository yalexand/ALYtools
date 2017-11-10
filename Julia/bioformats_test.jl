include("bioformats.jl")
include("tomography.jl")

# this code opens tomographic projections file,
# reads pixels and angle metadata,
# makes OPT reconstruction and saves 3D volume to file

id = "..\\TestData\\fluor.OME.tiff"
r = bfGetReader(id);

I = bfGetVolume(r)

angles = bfGetModulo(r,"Z")

jcall(r, "close", Void, ()) # reader not needed anymore

sizeX,sizeY,sizeZ,sizeC,sizeT = size(I)

z_slice = iradon(I[1,:,:,1,1],angles)
w,h = size(z_slice)

reconstruction = zeros(w,h,sizeX,1,1)

p = Progress(sizeX, "OPT reconstruction..")
for x = 1:sizeX
  reconstruction[:,:,x,1,1]  = iradon(I[x,:,:,1,1],angles)
  next!(p)
end

bfsaveAsOMEtiff(reconstruction,
                "..\\TestData\\reconstruction_fluor.OME.tiff",
                [],
                "Float64", # output pixel type
                "LZW", # compession
                true) # last argument - bigTiff
