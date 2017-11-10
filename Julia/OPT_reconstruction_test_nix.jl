include("bioformats.jl")
include("tomography.jl")


# e.g. "/home/yalexand/OPT_JULIA/TestData/fluor.OME.tiff"
src_id = ARGS[1]
dst_id = ARGS[2]


r = bfGetReader(src_id);
I = bfGetVolume(r)
angles = bfGetModulo(r,"Z")

jcall(r, "close", Void, ()) # reader not needed anymore

sizeX,sizeY,sizeZ,sizeC,sizeT = size(I)

z_slice = iradon(I[1,:,:,1,1],angles)
w,h = size(z_slice)

reconstruction = zeros(w,h,sizeX,1,1)

p = Progress(sizeX, "OPT reconstruction..")
#for x = 1:4 
for x = 1:sizeX
  reconstruction[:,:,x,1,1]  = iradon(I[x,:,:,1,1],angles)
  next!(p)
end

bfsaveAsOMEtiff(reconstruction,
                dst_id,
                [],
                "Float64",
                "LZW",
                true) # last argument - bigTiff

println(ARGS[1])
println(ARGS[2])
