function detection_proc_matlab_db(directory,n_img,first,last)
% DETECTION_PROC_MATLAB(DIRECTORY,N_IMG,FIRST,LAST)
% 
% this routine is written to binary filter round objects from an image
% the objects should have a radius that is always slightly larger than
% what is defined in 'min_pixel_radius'


%   Example:
%   >> directory = 'D:\technical\codes\IfU-IPG-PTV\working_folder_Dumbbell_b';
%   >> n_img = 4;
%   >> first = 1;
%   >> last = 115;
%   >> detection_proc_matlab_db('D:\technical\codes\IfU-IPG-PTV\working_folder_Dumbbell_b',4,1,115,1)

% Authors: Alex Liberzon and modified by Beat Lüthi
%
% Copyright (c) 2010, Alex Liberzon
% Copyright (c) 2010, Beat Lüthi
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

min_pixel_radius=30;

for img = 1:n_img
    for filenum = first:last
        [img filenum]
        % read the file
        f1 = imread([directory,'\img\db',int2str(img),'.',int2str(filenum)]);

        % put on the disk glasses and then only see disks with our radius
        se = strel('disk',min_pixel_radius=30);
        f2 = imopen(f1,se);
        f3=imextendedmax(f2,10);
        %imshow(f3);
        stats = regionprops(bwlabel(f3),{'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
        for i = 1:length(stats)
            stats(i).sumg = sum(f1(stats(i).PixelIdxList));
        end
        % Write the target file
        fname = [directory,'\img\db',int2str(img),'.',int2str(filenum),'_targets'];
        write_targets(fname,stats);
    end
end
