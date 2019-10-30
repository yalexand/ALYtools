function Vnew = resizeVolume(Volume,newSize)
% Function resizeVolume scales a volume to the new size.
% Volume: A 3D volume which has to be scaled
% newsize: A vector containing the new x y and z dimensions [x y z]

curSize = size(Volume);

% First resize X,Y
tempResult = zeros(newSize(1),newSize(2),curSize(3));
for ii = 1 : curSize(3)
    tempResult(:,:,ii) = imresize(Volume(:,:,ii),[newSize(1) newSize(2)],'lanczos3');
end

% Now resize X,Z
Vnew = zeros(newSize);
for ii = 1 : newSize(2)
    Vnew(:,ii,:) = imresize(squeeze(tempResult(:,ii,:)),[newSize(1) newSize(3)],'lanczos3');
end

end


