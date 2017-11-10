include("tomography.jl")
include("utils.jl")

z = shepp_logan(1000; highContrast=true)
ImageView.imshow(z)

angles = (0:1:359) # [degree]

tic();
sinogram = radon(z,angles)
display(toc())
display_image(sinogram)

tic();
reconstruction = iradon(sinogram,angles)
display(toc())
display_image(reconstruction)
