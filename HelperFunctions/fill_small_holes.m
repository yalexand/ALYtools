function filled = fill_small_holes(original,max_size)

    filled = imfill(original, 'holes');
    holes = filled & ~original;
    bigholes = bwareaopen(holes, max_size);
    smallholes = holes & ~bigholes;
    filled = original | smallholes;

end

