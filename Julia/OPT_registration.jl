include("tomography.jl")
include("utils.jl")

using StatsBase, ProgressMeter

#############################
                            # YA March 13 - presume I is floating point
function M1_do_registration!(I,angles,OPT_AXIS_INDEX)

      sizeX,sizeY,n_planes = size(I)

      M1_brightness_quantile_threshold = 0.96 # 0.8
      M1_max_shift = 20
      M1_window = 50

  #                 %
  #                 % select the slices of the image that are bright enough to work with.
  #                 % "brightEnough" calculation - different from Sam's

      if 2 == OPT_AXIS_INDEX
        sample = zeros(sizeY)
        for y=1:sizeY
          sino = I[:,y,:]
          sample[y] = mean(sino)
        end
      else # 2 == OPT_AXIS_INDEX
        sample = zeros(sizeX)
        for x=1:sizeX
          sino = I[x,:,:]
          sample[x] = mean(sino)
        end
      end

      T = quantile(sample, M1_brightness_quantile_threshold);

      brightEnough = sample .> T

  #                 % The shift correction and spearman correlation for individual slices are then
  #                 % found
  #
      shift = fill(NaN,length(sample))
      r = fill(NaN,length(sample))

      p = Progress(length(brightEnough),1,"gathering the data to calculate corrections.. please wait..")
      for n = 1:length(brightEnough)
          if 1 == brightEnough[n]
            sino = ( 2 == OPT_AXIS_INDEX ) ? I[:,n,:] : I[n,:,:]
            (shift[n],r[n]) = M1_quickMidindex(sino,M1_max_shift,angles)
          end
          next!(p)
      end

      # filter out shifts which imply large image rotation
      for n = 1:2
        # delt = abs(shift-nanmean(shift));
        delt = abs(shift-mean(shift[!isnan(shift)]))
        # shift(delt>9) = NaN;
        shift[delt.>9].=NaN
      end

      # slice numbers relative to centre of image, filtered, and then cropped
      ns = (1:length(brightEnough))'-length(brightEnough)/2;
      ns = ns[!isnan(shift)]

      shift = shift[!isnan(shift)]

      # fit shift and rotation
      p = polyfit(ns,shift,1)
      hshift = round(p[2])
      rotation = atan(p[1])

      NUMBER = (2 == OPT_AXIS_INDEX) ? sizeY : sizeX
      if rotation < NUMBER
        rotation = 0
      end

        vshift = hshift
        hshift = 0;
        # ?

        vshift = Int64(round(vshift))
        hshift = Int64(round(hshift))

        p = Progress(n_planes, "introducing corrections..")
        for k = 1:n_planes
          u = I[:,:,k]
          Ishift = M1_imshift(u,hshift,vshift,rotation);
          I[:,:,k] = Ishift;
          next!(p)
        end

end

#############################
function M1_quickMidindex(sino,maxshift,angles)

    spectrum_peak2 = zeros(maxshift+1)

    for i = 1:(maxshift+1)
        shiftsino = M1_sinoshift(sino',i,maxshift+1,0,0)
        slice = iradon(shiftsino',angles)
        # %figure
        # %imshow(slice,[]);
        # %icy_imshow(slice);
        sz_slice_X, sz_slice_Y = size(slice)
        slice = reshape(slice,1,sz_slice_X*sz_slice_Y)
        spectrum_peak2[i] = sqrt(sum(abs(slice.^2)));
    end

    sorted2 = sort(spectrum_peak2,rev=true)
    I2 = sortperm(spectrum_peak2,rev=true)

    diff = abs(I2-I2[1])

    r = corspearman(diff, sorted2)

    hshift = maxshift - 2*I2[1] + 2

    return (hshift, r)

end

#############################
function M1_sinoshift(sino,n,steps,shift,split)

        numOfAngularProjections, numOfParallelProjections = size(sino)

        if split == 1
            sino1 = sino[shift:round(numOfAngularProjections/2+shift-1),n:(numOfParallelProjections+n-steps)]
            sino2 = reverse(sino1)
            shiftedSino = [sino1 sino2]
        else
            shiftedSino = sino[:,n:(numOfParallelProjections+n-steps)]
        end
end

#############################
function M1_imshift(img,hShift,vShift,rotation)

        vsize,hsize = size(img)

        if hShift < 0
            img = img[:,(1-hShift):hsize]
        else
            img = img[:,1:(hsize-hShift)]
        end

        if vShift < 0
            img = img[(1-vShift):vsize,:]
        else
            img = img[1:(vsize-vShift),:]
        end

        if !isempty(rotation)
             if tan(rotation) > 1/min(hsize,vsize)
               # rotate around centre and get same-size central part of the rotated image
               img = imrotate(img,rotation*180/pi)
             end
        end

        shiftedImage = img;
end

# Matlab

# %-------------------------------------------------------------------------%
#         function M1_do_registration(obj,~) % by Samuel Davis
#             %
#              [sizeX,sizeY,n_planes] = size(obj.proj);
#              wait_handle = waitbar(0,'Ini proj memmap...');
#              [mapfile_name_proj,memmap_PROJ] = initialize_memmap([sizeX,sizeY,n_planes],1,'pixels',class(obj.proj),'ini_data',obj.proj);
#              close(wait_handle);
#
#              obj.proj = [];
#
#              PROJ = memmap_PROJ.Data.pixels; % reference
#
#                 %
#                 M1_brightness_quantile_threshold = 0.8;
#                 M1_max_shift = 20;
#                 %M1_window = 50;
#                 %
#                 % select the slices of the image that are bright enough to work with.
#                 % "brightEnough" calculation - different from Sam's
#                 sample = zeros(1,sizeY);
#                 for y=1:sizeY
#                     sino = squeeze(cast(PROJ(:,y,:),'single'));
#                     sample(y) = mean(sino(:));
#                 end
#                 T = quantile(sample, M1_brightness_quantile_threshold);
#                 brightEnough = sample > T;
#
#                 % main loop - starts
#
#                 % The shift correction and spearman correlation for individual slices are then
#                 % found
#
#                 if obj.isGPU
#                     shift = gpuArray(NaN(length(brightEnough),1));
#                     r = gpuArray(NaN(length(brightEnough),1));
#                 else
#                     shift = NaN(length(brightEnough),1);
#                     r = NaN(length(brightEnough),1);
#                 end
#
#                 waitmsg = 'gathering the data to calculate corrections.. please wait..';
#                 hw = waitbar(0,waitmsg);
#                 for n = 1:length(brightEnough)
#                     if brightEnough(n)
#                         sino = squeeze(cast(PROJ(:,n,:),'single'));
#                         [shift(n),r(n)] = obj.M1_quickMidindex(sino,M1_max_shift);
#                         %
#                         if ~isempty(hw), waitbar(n/length(brightEnough),hw); drawnow, end;
#                     end
#                 end
#                 if ~isempty(hw), delete(hw), drawnow, end;
#
#                 if obj.isGPU
#                     shift = gather(shift);
#                 end
#                 %
#                 % filter out shifts which imply large image rotation
#                 for n = 1:2
#                     delt = abs(shift-nanmean(shift));
#                     shift(delt>9) = NaN;
#                 end
#
#                 % slice numbers relative to centre of image, filtered, and then cropped
#                 ns = (1:length(brightEnough))'-length(brightEnough)/2;
#                 ns = ns(~isnan(shift));
#                 %ns = ns(M1_window/2:end-M1_window/2);
#                 shift = shift(~isnan(shift));
#                 %shift = shift(M1_window/2:end-M1_window/2);
#
#                 % fit shift and rotation
#                 p = polyfit(ns,shift,1);
#                 hshift = round(p(2));
#                 rotation = atan(p(1));
#                 if rotation < sizeY
#                     rotation = 0;
#                 end
#             %
#             vshift = hshift;
#             hshift = 0;
#             % ?
#             %
#                 obj.M1_hshift = hshift; % keep it for the case of re-usage for FLIM
#                 obj.M1_vshift = vshift;
#                 obj.M1_rotation = rotation;
#             %
#             waitmsg = 'introducing corrections..';
#             hw = waitbar(0,waitmsg);
#             for k = 1:n_planes
#                 I = PROJ(:,:,k);
#                 Ishift = obj.M1_imshift(I,hshift,vshift,rotation);
#                     if isempty(obj.proj)
#                         [szx,szy] = size(Ishift);
#                         obj.proj = zeros(szx,szy,n_planes,class(Ishift));
#                     end
#                 obj.proj(:,:,k) = Ishift;
#                 if ~isempty(hw), waitbar(k/n_planes,hw); drawnow, end;
#             end
#             if ~isempty(hw), delete(hw), drawnow, end;
#             %
#             clear('memmap_PROJ');
#             delete(mapfile_name_proj);
#
#         end
# %-------------------------------------------------------------------------%
# function [hshift, r] = M1_quickMidindex(obj,sino,maxshift)
#
#     if obj.isGPU
#         sino = gpuArray(sino);
#     end
#
#     for i = 1:(maxshift+1)
#         shiftsino = obj.M1_sinoshift(sino',i,maxshift+1,0,0);
#         slice = iradon(shiftsino',obj.angles,'linear','Hann');
#         %figure
#         %imshow(slice,[])
#         %icy_imshow(slice);
#         spectrum_peak2(i) = sqrt(sum(abs(slice(:).^2)));
#     end
#
#     [sorted2, I2] = sort(spectrum_peak2,'descend');
#
#     diff = abs(I2-I2(1));
#
#     if obj.isGPU
#         r = corr(gather(diff(:)),gather(sorted2(:)),'type','Spearman');
#     else
#         r = corr(diff(:),sorted2(:),'type','Spearman');
#     end
#
#     %figure
#     %plot(-maxshift:2:maxshift,spectrum_peak2);
#     %ylabel('Peak Intensity/a.u.')
#     %xlabel('rotation axis shift/pixels')
#     %box off
#     %drawnow
#
#     %if r < -0.6
#         hshift = maxshift - 2*I2(1) + 2;
#     %else
#     %    hshift = NaN;
#     %end
# end
# %-------------------------------------------------------------------------%
# function shiftedSino = M1_sinoshift(obj,sino,n,steps,shift,split)
#
#         numOfParallelProjections = size(sino,2);
#         numOfAngularProjections = size(sino,1);
#
#         if split == 1
#             sino1 = sino(shift:(numOfAngularProjections/2+shift-1),n:(numOfParallelProjections+n-steps));
#             sino2 = fliplr(sino1);
#             shiftedSino = cat(1,sino1,sino2);
#         else
#             shiftedSino = sino(:,n:(numOfParallelProjections+n-steps));
#         end
# end
# %-------------------------------------------------------------------------%
# function shiftedImage = M1_imshift(obj,img,hShift,vShift,rotation)
#
#         hsize = size(img,2);
#         vsize = size(img,1);
#
#         if hShift < 0
#             img = img(:,(1-hShift):hsize);
#         else
#             img = img(:,1:(hsize-hShift));
#         end
#
#         if vShift < 0
#             img = img((1-vShift):vsize,:);
#         else
#             img = img(1:(vsize-vShift),:);
#         end
#         if nargin > 3
#             if tan(rotation) > 1/min(hsize,vsize)
#                 img = imrotate(img,rotation*180/pi,'bicubic');
#             end
#         end
#         shiftedImage = img;
# end
# %-------------------------------------------------------------------------%
