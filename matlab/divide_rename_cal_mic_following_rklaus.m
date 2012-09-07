function [] = split_image
% Reads all the image between the first and last, provided
% by the user, subtracts from each image the mean image
% Inputs:         none
% Outputs:        none
%
% Usage:
%        >> subtract_mean_image
% point to the desired first and last files
%
% See also: HELP UIGETFILE, IMREAD, IMWRITE

% Author: Alex Liberzon
% Copyright (c) 2004, Alex Liberzon, IHW, ETHZ
% Last modified at: June 17, 2004
%  Version 1.01, on VideoPC
% - takes average of 100 files, than convert and divide all

%first image
[filename1,pathname] = uigetfile('*.tif','First file');
wd = cd;
cd(pathname);

[pathstr,name1,ext1,versn] = fileparts(filename1);

first = eval(name1(end-4:end));

%second image
[filename2,pathname2] = uigetfile('*.tif','Last file');
[pathstr2,name2,ext2,versn2] = fileparts(filename2);
last = eval(name2(end-4:end));


%cicle images were renamed such run210002.tif where run2 is the experiment
%and 10002 is the image number
for i =  first:last
    [tmp,map] = imread([name1(1:end-5),sprintf('%5d',i),'.tif']);
    % tmp = imsubtract(tmp, uint8(meana));
    imwrite(tmp(1:512,1:512),['Cam1.',sprintf('%5d',i)],'tiff','compression','none');
    imwrite(tmp(1:512,513:1024),['Cam2.',sprintf('%5d',i)],'tiff','compression','none');
    imwrite(tmp(513:1024,513:1024),['Cam3.',sprintf('%5d',i)],'tiff','compression','none');
    imwrite(tmp(513:1024,1:512),['Cam4.',sprintf('%5d',i)],'tiff','compression','none');
    i
end

cd(wd);
disp('Done ...')

stop
%calibration image


[filename1,pathname] = uigetfile('*.tif','First file');
wd = cd;
cd(pathname);
[pathstr,name1,ext1,versn] = fileparts(filename1);


[tmp,map] = imread([name1,'.tif']);
    % tmp = imsubtract(tmp, uint8(meana));
    imwrite(tmp(1:512,1:512),[name1,'Cam1.tif'],'tiff','compression','none');
    imwrite(tmp(1:512,513:1024),[name1,'Cam2.tif'],'tiff','compression','none');
    imwrite(tmp(513:1024,513:1024),[name1,'Cam3.tif'],'tiff','compression','none');
    imwrite(tmp(513:1024,1:512),[name1,'Cam4.tif'],'tiff','compression','none');

    