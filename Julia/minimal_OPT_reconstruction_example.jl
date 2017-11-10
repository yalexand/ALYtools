include("bioformats.jl")
include("tomography.jl")

src_img_filename = "..\\TestData\\fluor.OME.tiff"
dst_img_filename = "..\\TestData\\reconstruction_fluor.tiff"

r = bfGetReader(src_img_filename)
angles = bfGetModulo(r,"Z")
jcall(r, "close", Void, ())

proj = convert(Array{Float64,3},load(src_img_filename))
sizeX,sizeY,sizeZ = size(proj)

reconstruction = []

p = Progress(sizeY, "OPT reconstruction..")
for y = 1:sizeY
  slice = iradon(proj[:,y,:],angles)
    if isempty(reconstruction)
      reconstruction = zeros(w,h,sizeY)
    end
  reconstruction[:,:,y] = slice
  next!(p)
end

save(dst_img_filename,
colorview(Gray,AxisArray(normedview(convert(Array{Normed{UInt16,16},3},reconstruction)),:x,:y,:z)))
