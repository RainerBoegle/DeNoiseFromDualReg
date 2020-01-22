function [newFileList,stdNoiseFileList] = subtractDeSignalDataFromOrgData(pathDeSignalFileList,pathOrgFileList,outputPathNewFileList,pathMaskFile)
% This function will subtract the DeSignaled/Noise data from the Original data as 
% a means of noise removal.
% 
% The FileList of the DeSignaled/Noise Data must be provided, as well as the
% the FileList of the Original Data, and these must have the same length, of course.
%
% Also, an output path for a new file list should be provided, if not provied 
% or left empty, then no new filelist will be created and the user 
% has to do this manually later.
% A mask can be provided or created automatically, as well.
%
% The results will be saved in the path of the Original Data and the filenames
% will get the prefix "DeNoised_".
% An additional file with the prefix "NoiseStd_" will be stored in the 
% Original Data folder which contains the standard deviation over time of
% the DeSignaled/Noise Data. 
% --> This can later be used on independent components that were not normalized 
%     by the error variance to normalize them for the following mixture model inference.
%
%Steps:
% For both data the mean over time will be removed before subtraction and
% later, the mean over time from the Original data will be added in again.
%
%
%
%Usage:
%      [newFileList,stdNoiseFileList,subtractionReport] = subtractDeSignalDataFromOrgData(pathDeSignalFileList,pathOrgFileList,outputPathNewFileList,pathMaskFile);
%      newFileList = subtractDeSignalDataFromOrgData(pathDeSignalFileList,pathOrgFileList,outputPathNewFileList); %auto-mask 
%      newFileList = subtractDeSignalDataFromOrgData(pathDeSignalFileList,pathOrgFileList); %auto-mask & do not output new file list (this is not suggested!) 
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: (07.09.2018): initial implementation.

%% check inputs
assert(exist(pathDeSignalFileList,'file')~=0,['Error(pathDeSignalFileList): Could not find "',getFName(pathDeSignalFileList),'" in folder "',fileparts(pathDeSignalFileList),'"! Check file path of DeSignal Data FileList!']);
assert(exist(pathOrgFileList,'file')~=0,     ['Error(pathOrgFileList): Could not find "',     getFName(pathOrgFileList),     '" in folder "',fileparts(pathOrgFileList),     '"! Check file path of Original Data FileList!']);
%check that both FileLists have the same length
fileListDeSignalData = importdata(pathDeSignalFileList);
fileListOriginalData = importdata(pathOrgFileList);
if(length(fileListDeSignalData)~=length(fileListOriginalData))
    error(['FileList of DeSignal Data has length==',num2str(length(fileListDeSignalData)),' while FileList Original Data has length==',num2str(length(fileListOriginalData)),', but should be equal! Check inputs!']);
end

%outputPathNewFileList
if(~exist('outputPathNewFileList','var'))
    outputPathNewFileList = [];
end
if(isempty(outputPathNewFileList))
    disp('Will not create a new file list of the DeNoised data, you will have to do this manually then.');
end

%pathMaskFile
if(~exist('pathMaskFile','var'))
    pathMaskFile = [];
end
if(isempty(pathMaskFile))
    disp('No mask provided, will create mask automatically from the data.');
elseif(~exist(deSPMpath(pathMaskFile),'file'))%is it SPM-style path with ,# added at the end? --> make sure this doesn't interfere here
    error(['Error(pathMaskFile): Could not find "',getFName(deSPMpath(pathMaskFile)),'" in folder "',fileparts(pathMaskFile),'"! Check file path of Mask!']);
else
    disp(['Will use mask "',getFName(pathMaskFile),'" in folder "',fileparts(pathMaskFile),'".']);
end

%% do all subtractions
nFiles = length(fileListOriginalData);
newFileList      = cell(nFiles,1);
stdNoiseFileList = cell(nFiles,1);
mask = prepMask(pathMaskFile,fileListOriginalData{1});
for indFile = 1:nFiles
    disp(['File ',num2str(indFile),'of',num2str(nFiles),': Subtracting "',getFName(fileListDeSignalData{indFile}),'" from "',getFName(fileListOriginalData{indFile}),'"...']);
    [newFileList{indFile},stdNoiseFileList{indFile}] = doSubtraction(fileListDeSignalData{indFile},fileListOriginalData{indFile},mask);
end

%% write out FileList if possible
if(~isempty(outputPathNewFileList))
    writeNewFileList(newFileList,outputPathNewFileList);
    subtractionReport = struct('fileListDeSignalData',{fileListDeSignalData},'fileListOriginalData',{fileListOriginalData},'newFileList',{newFileList},'pathMaskFile',pathMaskFile,'outputPathNewFileList',outputPathNewFileList);
    save([fileparts(subtractionReport.outputPathNewFileList),filesep,'subtractionReport.mat'],'subtractionReport')
end

%% Done.
disp(' ');
disp('All Done.');
disp(' ');

end

%% subfunction
%% writeNewFileList
function [] = writeNewFileList(newFileList,outputPath)
% Write the new filelist to path

%% write new file list
[baseDirNewFileList,fNameNewFileList,extNewFileList] = fileparts(outputPath); 
disp(['Writing out new file list "',fNameNewFileList,extNewFileList,'" to directory "',baseDirNewFileList,'".']);
nFiles = length(newFileList);
fileID = fopen(outputPath,'w'); %open file for new filelist for writing
try
    for indFile = 1:nFiles
        %%write line to new denoise_output.filelist file in OutputDir
        fprintf(fileID,'%s\n',newFileList{indFile}); %write line to file denoise_output.filelist in OutputDir
    end
catch CATCH_exception
    fclose(fileID);
    assignin('base','CATCH_exception_filelistWrite',CATCH_exception);
    error(['Error: could not write out new file list "',fNameNewFileList,extNewFileList,'" to directory "',baseDirNewFileList,'". (file was successfully closed though, and the exception was saved in the base workspace.)']);
end
%% close file 
fclose(fileID);

end

%% doSubtraction
function [outputFilePath,outputFilePathStdev] = doSubtraction(pathDeSignalData,pathOriginalData,mask)
% Do the subtraction of a original file using the equivalent DeSignal/Noise file.

%% read data & check dims
vols_DeSignal = spm_vol(pathDeSignalData);
vols_OrgData = spm_vol(pathOriginalData);
assert(length(vols_DeSignal)==length(vols_OrgData),['Error: DeSignal Data has ',num2str(length(vols_DeSignal)),' volumes, while Original Data has ',num2str(length(vols_OrgData)),' volumes, but should be of equal length! Check inputs!']); 
assert(all(vols_DeSignal(1).dim(:)==vols_OrgData(1).dim(:)),['Error: DeSignal Data Dimensions are [',num2str(vols_DeSignal(1).dim(1)),',',num2str(vols_DeSignal(1).dim(2)),',',num2str(vols_DeSignal(1).dim(3)),'], while Original Data Dimensions are [',num2str(vols_OrgData(1).dim(1)),',',num2str(vols_OrgData(1).dim(2)),',',num2str(vols_OrgData(1).dim(3)),'], but should be equal! Check inputs!']); 

nVols = length(vols_OrgData);
data4D_DeSignal = spm_read_vols(vols_DeSignal);
data4D_Original = spm_read_vols(vols_OrgData);

%% reshape and apply mask
data2D_DeSignal = reshape(data4D_DeSignal,[],nVols);
data2D_Original = reshape(data4D_Original,[],nVols);

data2D_DeSignal = data2D_DeSignal(mask~=0,:);
data2D_Original = data2D_Original(mask~=0,:);

%% clear memory
% clear data4D_DeSignal data4D_Original

%% remove each mean but keep original data mean & stdev of the DeSignaled/Noise Data
mean_Original = mean(data2D_Original,2);
deMeanData2D_DeSignal = data2D_DeSignal - mean(data2D_DeSignal,2);
deMeanData2D_Original = data2D_Original - mean_Original;
std_DeSignal = std(deMeanData2D_DeSignal,0,2);

%% do subtraction 
deNoisedData2D = deMeanData2D_Original - deMeanData2D_DeSignal;

%% add mean of original
deNoisedData2DwMean = deNoisedData2D + mean_Original;

%% write DeNoised Data
vOut = rmfield(vols_OrgData(1),'private');
if(vOut.dt(1)<16)
    vOut.dt(1) = 16;
end

[baseDir,fName,ext] = fileparts(pathOriginalData);
outputFilePath = [baseDir,filesep,'DeNoised_',fName,ext];
vOut.fname = outputFilePath;

Order      = ceil(log10(nVols)); %number of zeros to add for output text
reverseStr = ''; %init reverse string empty and later set it to the length of character IN BACKSPACES ("\b"), that has been printed in the last step, in order to remove them.
for indVol = 1:nVols
    vOut.n(1) = indVol;
    currVol = zeros(vOut.dim(:)');
    currVol(mask~=0) = deNoisedData2DwMean(:,indVol);
    
    %info
    msg = sprintf(['Writting volume %0',num2str(Order),'d of %0',num2str(Order),'d Volumes of file "DeNoised_',fName,ext,'".'], indVol, nVols); %Don't forget this semicolon
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %prep the next reverse string, i.e. a string of BACKSPACES in the length of the message, i.e. delete the message on screen.
    
    %write
    spm_write_vol(vOut,currVol);
end
disp(' ');
disp(['Writting "NoiseStd_',fName,ext,'" to directory "',baseDir,'".']);
vOut = rmfield(vols_OrgData(1),'private');
if(vOut.dt(1)<16)
    vOut.dt(1) = 16;
end
vOut.n(1) = 1;
outputFilePathStdev = [baseDir,filesep,'NoiseStd_',fName,ext];
vOut.fname = outputFilePathStdev;
currVol = zeros(vOut.dim(:)');
currVol(mask~=0) = std_DeSignal;
spm_write_vol(vOut,currVol);

end

%% prepMask
function mask = prepMask(pathMaskFile,pathOriginalData)
%prepare the mask

vols_OrgData = spm_vol(pathOriginalData);

if(isempty(pathMaskFile))
    mask = ones(vols_OrgData(1).dim(:)');
    mask = mask(:);
else
    vol_mask = spm_vol(pathMaskFile);
    if(length(vol_mask)>1)
        disp('mask is a 4D NIFTI--> will take the first volume.');
    end
    assert(all(vols_OrgData(1).dim(:)==vol_mask.dim(:)),['Error: Original Data Dimensions are [',num2str(vols_OrgData(1).dim(1)),',',num2str(vols_OrgData(1).dim(2)),',',num2str(vols_OrgData(1).dim(3)),'], while the Mask has dimensions [',num2str(vol_mask.dim(1)),',',num2str(vol_mask.dim(2)),',',num2str(vol_mask.dim(3)),'], but should be equal! Check inputs!']); 
    mask = spm_read_vols(vol_mask);
    mask = mask(:);
end

end

%% deSPMpath
function path2FileDeSPM = deSPMpath(pathMaybeSPMstyle)
% remove any comma and numbers in the extension of a file

[baseDir,fName,extOrg] = fileparts(pathMaybeSPMstyle);
idx = regexp(extOrg,',');
if(~isempty(idx))
    path2FileDeSPM = [baseDir,filesep,fName,extOrg(1:(idx-1))];
else
    path2FileDeSPM = [baseDir,filesep,fName,extOrg];
end

end

%% getFName
function fName = getFName(path2File)
% get the filename

[~,fName,ext] = fileparts(path2File);

if(~isempty(ext))
    fName = [fName,ext];
end

end
