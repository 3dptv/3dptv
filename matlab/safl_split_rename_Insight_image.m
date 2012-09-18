% safl_split_rename_Insight_image
%
% Reads all the image between the first and last, provided
% by the user, subtracts from each image the mean image
% Inputs:         none
% Outputs:        none
%
%
% See also: HELP UIGETFILE, IMREAD, IMWRITE

% Author: Alex Liberzon
% Copyright (c) 2012, alex.liberzon@gmail.com
% Last modified at: Sep 13, 2012, at SAFL
% - reads now a list of A,B images captured by Insight 4G at 180 Hz
% - renames and splits the iamges into 4 sub-images
% - useful only for the present SAFL installation


% First image in the sequence
% NOTE: user might choose A or B, ignores the choice

[filename1,pathname] = uigetfile({'*.tif'; '*.TIF'},'Pick the FIRST file');
wd = cd;
cd(pathname);

% Last image in the sequence
[filename2,pathname2] = uigetfile({'*.tif'; '*.TIF'},'LAST file');



% Split the image into few useful parts:
% base_name is what is the experimental name in Insight
% digits will be used later just for proper sequence
[first_name,last_part] = strtok(filename1,'.'); % the actual name and the number before the first dot
digits_in_the_name = regexpi(first_name,'[0-9]');
base_name = first_name(1:digits_in_the_name(1)-1);
first_image = str2double(first_name(digits_in_the_name));
num_digits = length(digits_in_the_name);


% Split the last image to get the number
last_name = strtok(filename2,'.');
last_image = str2double(last_name(digits_in_the_name));



% how to reconstruct the name:
% all, i.e A and B files with this name and number

format_string  = ['%s%0',num2str(num_digits),'d*'];

% running counter of images
counter = 10^num_digits;


% Safety - if this directory already has cam1.1***** images,
% then ask the user if to overwrite or append
d = dir('cam1.*');
if ~isempty(d)
    [~,last_existing] = strtok(d(end).name,'.');
    last_existing = str2num(last_existing(2:end)); % . is the first character
    btn = questdlg('Overwrite ?', ...
        'Cam1.*** files exist ',...
        'Yes','No','No');
    switch btn,
        case 'No',
            counter = last_existing+1;
        case 'Yes'
            continue
    end % switch
end




% In order to be sure that the files were selected in ordered mode:

if last_image < first_image
    tmp = first_image;
    first_image = last_image;
    last_image = tmp;
end

% Loop through the images, remember there are A and B ones in PIV mode
for i =  first_image:last_image
    
    % Reconstruct the name, get the names of the A, B ones
    images_to_process = dir(sprintf(format_string,base_name,i));
    
    for j = 1:length(images_to_process) % there are maybe LA and LB
        imname = images_to_process(j).name;
        
        
        [tmp,map] = imread(imname);
        
        dim = size(tmp)/2;
        
        % Split and rename, using the running counter, add the basename to the
        % TIFF file as a description
        
        imwrite(tmp(1:dim(1),1:dim(2)),['cam1.',sprintf('%05d',counter)],'tiff','compression','none','Description',imname);
        imwrite(tmp(1:dim(1),dim(2)+1:2*dim(2)),['cam2.',sprintf('%05d',counter)],'tiff','compression','none','Description',imname);
        imwrite(tmp(dim(1)+1:2*dim(1) ,dim(2)+1:2*dim(2)),['cam3.',sprintf('%05d',counter)],'tiff','compression','none','Description',imname);
        imwrite(tmp(dim(1)+1:2*dim(1),1:dim(2)),['cam4.',sprintf('%05d',counter)],'tiff','compression','none','Description',imname);
        counter = counter + 1;
    end
end


cd(wd);
disp('Done ...')


% %% calibration image
%
%
% [filename1,pathname] = uigetfile('*.tif','First file');
% wd = cd;
% cd(pathname);
% [pathstr,name1,ext1,versn] = fileparts(filename1);
%
%
% [tmp,map] = imread([name1,'.tif']);
% dim=size(tmp)/2;
% % tmp = imsubtract(tmp, uint8(meana));
% imwrite(tmp(1:dim(1),1:dim(2)),['Cam1.tif'],'tiff','compression','none');
% imwrite(tmp(1:dim(1),dim(2)+1:2*dim(2)),['Cam2.tif'],'tiff','compression','none');
% imwrite(tmp(dim(1)+1:2*dim(1) ,dim(2)+1:2*dim(2)),['Cam3.tif'],'tiff','compression','none');
% imwrite(tmp(dim(1)+1:2*dim(1),1:dim(2)),['Cam4.tif'],'tiff','compression','none');
%
% %check
% [c1,map] =imread('Cam1.tif');
% size(c1)
% [c2,map] =imread('Cam2.tif');
% size(c2)
% [c3,map] =imread('Cam3.tif');
% size(c3)
% [c4,map] =imread('Cam4.tif');
% size(c4)
% %ok
%
%
% % One must test the stuff. I guess we somehow messed up the
% % definitions?
%
% figure, subplot(221),imshow(Cam1),title('Cam 1'),subplot(222),imshow(Cam2),title('Cam 2'),subplot(223),imshow(Cam4),title('Cam 4'),subplot(224),imshow(Cam3),title('Cam 3')
% figure, imshow(tmp)


