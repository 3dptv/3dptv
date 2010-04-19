% DETECTION_PROC_MATLAB
% Detection procedire in Matlab
% a different type of the particle centroid detection
% using the same subroutine that allows to detect breakage
% the particles are identified as bright spots within a larger, grey spot,
% see COLLOID_IMAGE_SEGMENTATION_INPOLY.M for details
%
% Author: Alex Liberzon

% Copyright (c) 2010, Alex Liberzon
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

close all
clear all

first = 10000;
last  = 10004;
n_img = 4;

name_root=['C:\PTV\Software\origo-ptv\test_for_trunk\'];

% change verbose to 1 if you want to display the process
verbose = 0;
for img = 1:n_img
    for filenum = first:last
        % read the file
        f1 = imread([name_root,'img\cam',int2str(img),'.',int2str(filenum)]);

        % Process the file - single liner at the moment:
        stats = regionprops(bwlabel(bwareaopen(imfill(imclose(imextendedmax(adapthisteq(f1),60),...
            ones(3,3)), 'holes'), 10)),{'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
        for i = 1:length(stats)
            stats(i).sumg = sum(f1(stats(i).PixelIdxList));
        end
        % Write the target file
        fname = [name_root,'img\cam',int2str(img),'.',int2str(filenum),'_targets'];
        write_targets(fname,stats);
    end
end
