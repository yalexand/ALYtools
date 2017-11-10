function [microns_per_pixel_x,microns_per_pixel_z,I] = bfopen_v(id, varargin)
% Open microscopy images using Bio-Formats.
%
% SYNOPSIS r = bfopen(id)
%          r = bfopen(id, x, y, w, h)
%
% Input
%    r - the reader object (e.g. the output bfGetReader)
%
%    x - (Optional) A scalar giving the x-origin of the tile.
%    Default: 1
%
%    y - (Optional) A scalar giving the y-origin of the tile.
%    Default: 1
%
%    w - (Optional) A scalar giving the width of the tile. 
%    Set to the width of the plane by default.
%
%    h - (Optional) A scalar giving the height of the tile.
%    Set to the height of the plane by default.
%
% Output
%
%    result - a cell array of cell arrays of (matrix, label) pairs, 
%    with each matrix representing a single image plane, and each inner 
%    list of matrices representing an image series.
%
% Portions of this code were adapted from:
% http://www.mathworks.com/support/solutions/en/data/1-2WPAYR/
%
% This method is ~1.5x-2.5x slower than Bio-Formats's command line
% showinf tool (MATLAB 7.0.4.365 R14 SP2 vs. java 1.6.0_20),
% due to overhead from copying arrays.
%
% Thanks to all who offered suggestions and improvements:
%     * Ville Rantanen
%     * Brett Shoelson
%     * Martin Offterdinger
%     * Tony Collins
%     * Cris Luengo
%     * Arnon Lieber
%     * Jimmy Fong
%
% NB: Internet Explorer sometimes erroneously renames the Bio-Formats library
%     to bioformats_package.zip. If this happens, rename it back to
%     bioformats_package.jar.
%
% For many examples of how to use the bfopen function, please see:
%     http://trac.openmicroscopy.org.uk/ome/wiki/BioFormats-Matlab

% OME Bio-Formats package for reading and converting biological file formats.
%
% Copyright (C) 2007 - 2014 Open Microscopy Environment:
%   - Board of Regents of the University of Wisconsin-Madison
%   - Glencoe Software, Inc.
%   - University of Dundee
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% -- Configuration - customize this section to your liking --

% Toggle the autoloadBioFormats flag to control automatic loading
% of the Bio-Formats library using the javaaddpath command.
%
% For static loading, you can add the library to MATLAB's class path:
%     1. Type "edit classpath.txt" at the MATLAB prompt.
%     2. Go to the end of the file, and add the path to your JAR file
%        (e.g., C:/Program Files/MATLAB/work/bioformats_package.jar).
%     3. Save the file and restart MATLAB.
%
% There are advantages to using the static approach over javaaddpath:
%     1. If you use bfopen within a loop, it saves on overhead
%        to avoid calling the javaaddpath command repeatedly.
%     2. Calling 'javaaddpath' may erase certain global parameters.
autoloadBioFormats = 1;

% Toggle the stitchFiles flag to control grouping of similarly
% named files into a single dataset based on file numbering.
stitchFiles = 0;

% To work with compressed Evotec Flex, fill in your LuraWave license code.
%lurawaveLicense = 'xxxxxx-xxxxxxx';

% -- Main function - no need to edit anything past this point --

% load the Bio-Formats library into the MATLAB environment
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

% initialize logging
loci.common.DebugTools.enableLogging('INFO');

% Get the channel filler
r = bfGetReader(id, stitchFiles);

sizeZ = r.getSizeZ();
sizeC = r.getSizeC();
sizeT = r.getSizeT();
                
arr = bfGetPlane(r, 1, varargin{:});
[sizeX,sizeY] = size(arr);
I = zeros(sizeX,sizeY,sizeZ,sizeC,sizeT,class(arr));

    r.setSeries(0);
    numImages = r.getImageCount();
    
    for i = 1:numImages
        if mod(i, 72) == 1
            fprintf('\n    ');
        end
        fprintf('.');
        
        zct = r.getZCTCoords(i - 1);
        z = zct(1)+1;
        c = zct(2)+1;
        t = zct(3)+1;
        
        arr = bfGetPlane_comp_5_10_2017(r, i, varargin{:});
        I(:,:,z,c,t) = double(arr);
    end
            
            microns_per_pixel_x = [];
            microns_per_pixel_z = [];            
            try
                omeMeta = r.getMetadataStore();  
                phszx = omeMeta.getPixelsPhysicalSizeX(0).getValue;
                phszy = omeMeta.getPixelsPhysicalSizeY(0).getValue;
                phszz = omeMeta.getPixelsPhysicalSizeZ(0).getValue;
                if ~isempty(phszx) && ~isempty(phszy) && phszx == phszy && ~isempty(phszz)
                    microns_per_pixel_x = phszx;
                    microns_per_pixel_z = phszz;
                end
            catch
            end
                
    fprintf('\n');
    r.close();
end
