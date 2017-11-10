using Images, ImageView, CoordinateTransformations, ImageTransformations

#################################
function display_image(I)
  v1 = minimum(I)
  v2 = maximum(I)
  todisplay = (I-v1)/(v2-v1)
  ImageView.imshow(todisplay)
end

###################################
# rotate image "I" around its centre by "angle" [degrees]
# returns or "full" or "central" part of the rotated image
function imrotate(I,angle,spec::String="central")

  tform = recenter(RotMatrix(angle/180*pi), center(I))
  Irot = parent(warp(I,tform))

  if (spec=="full")
    Irot
  elseif (spec=="central")
    Ifin = Irot[UnitRange.(indices(I))...]
  else
    I
  end

end

#############################
function polyfit(x, y, n)
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  A \ y
end
