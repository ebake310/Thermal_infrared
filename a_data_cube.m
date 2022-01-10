%% Organizing image data into a data cube
% By Caroline Aubry-Wake

% all images have to be in a folder as asci files
% Calls all the images in a folder and organizeds them in a data cube with them

cd 'C:\Users\Emily\Google Drive\Peru_2016\Peru_2016_asc' ;% directory of all files

filesStruct = dir('*.asc'); % the parameters for each file that is imported
files = {filesStruct.name}; % the name of each file in the cube

dataCell = cell(numel(files), 1);

 delimiterIn = '\t';     % tab delimiter
 headerlinesIn = 1;      % headerline from column title

for i=1:numel(files)
    a = importdata(files{i},delimiterIn,headerlinesIn);
    dataCell{i} =a.data(:, 2:1025); % because column 1 is just the number of the pixel
end

data = cat(3, dataCell{:});

save IRdata_Aug2016 data files filesStruct -v7.3
% save data as the dataset to used
