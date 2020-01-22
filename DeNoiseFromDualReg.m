function [NewFileListPath,NewFileList,DeNoiseInfoStruct] = DeNoiseFromDualReg(DualRegDir,FileListPath,ICsRemove,RemoveRes,OutputDir,pathAltMask)
% This function can be used to remove selected ICs (and residuals) from fMRI data, based on the results
% of dual regression analysis of MELODIC analysis of group data.
% NB: 
%    use dual_regression_residuals to create the dual regression results that also include the residuals
%    after dual regression, i.e. the variability that was not explained by dual regression or
%    in an alternative view the residuals can be seen as the variability not included in the initial ICA.
%
%Inputs:
%       DualRegDir    (string)   The path to the directory containing the dual regression results.
%       FileListPath  (string)   The path to the .filelist file that contains the paths to the original
%                                input data to the ICA, i.e. the data to be denoised.
%       ICsRemove     (vector)   The indices of the ICs to remove from the data, i.e. perform group ICA
%                                and identify the artifact ICs for removing from all data.
%                                NB: the ICs identified on the group level are "adjusted" for each
%                                input data file/session/subject via dual regression, which also
%                                account for the non-artifact parts of the data --> regression will
%                                not damage the data as normal nuisance regression might do.
%       RemoveRes    (bool 1or0) Indicates if the residuals after dual regression should be removed as well. 
%                                (NB: these residuals can also be seen as the variability that was excluded
%                                     from the ICA by the initial PCA step.)
%       OutputDir     (string)   Path to the directory in which to save the output data/denoised data. 
%                                The data will be output as denoise_output_00XXX.nii XXX being the input number 
%                                with maximum trail of 5 zeros.
%                                NB: in this path there will be a "denoise_output.filelist" file that points
%                                    to the denoised data, in case you want to directoy continue with MELODIC on this data.
%       pathAltMask   (string)   Path to alternative mask, if not given or empty just use DualReg mask.
%
%
%Usage:
%      [NewFileListPath,NewFileList,DeNoiseInfoStruct] = DeNoiseFromDualReg(DualRegDir,FileListPath,ICsRemove,RemoveRes,OutputDir,pathAltMask);
%
%
%
%V1.6
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.6(01.06.2018): allow alternative mask V1.5: 05.03.2017): complete reworking of the original function! Memory efficient and (comparatively) faster than before. This should handle high-res fMRI without a problem. NB: V1.4 was a hack that allowed high-res fMRI data but very inefficient, however this inspired the current version. V1.3 (30.12.2016): More outputs that show progress. V1.2: (29.09.2016): From this version on we will add the mean of the timeseries, i.e. the mean image, back onto the denoised data, because pure residuals cause a problem in SPM12 first-level analysis. (NB: This should not matter as it just adds a constant offset to all data, but somehow this really screws with SPM12's estimation of the data mask, ending in "no significant voxels in analysis" error!) (NB2: I had used a different function to add the mean to the data to fix this error, but I choose to avoid it altogether by including this here.) V1.1: (28.04.2016): fixed major bug last version only removed ONE IC now it removes all selected ICs!!! Damn, I should be ashamed. V1.0: (11.04.2016): initial implementation based on earlier version for GVSaging Project

NewMode = 1;

%% check inputs
%DualRegDir
if(~exist('DualRegDir','var'))
    DualRegDir = spm_select(1,'dir','Select Dual Regression Directory...');
    if(isempty(DualRegDir))
        NewFileListPath = [];
        NewFileList     = [];
        disp('Quit');
        return;
    else
        disp(['Will load Dual Regression results from directory "',DualRegDir,'"...']);
    end
elseif(isempty(DualRegDir))
    DualRegDir = spm_select(1,'dir','Select Dual Regression Directory...');
    if(isempty(DualRegDir))
        NewFileListPath = [];
        NewFileList     = [];
        disp('Quit');
        return;
    else
        disp(['Will load Dual Regression results from directory "',DualRegDir,'"...']);
    end
elseif(~exist(DualRegDir,'dir'))
    error(['Could not find Dual Regression Directory "',DualRegDir,'"!']);
else
    disp(['Will load Dual Regression results from directory "',DualRegDir,'"...']);
end

%FileListPath
if(~exist('FileListPath','var'))
    FileListPath = spm_select(1,'any','Select .filelist file pointing to the original data paths...');
    if(isempty(FileListPath))
        NewFileListPath = [];
        NewFileList     = [];
        disp('Quit');
        return;
    else
        disp(['Will use file list "',FileListPath,'"...']);
    end
elseif(isempty(FileListPath))
    FileListPath = spm_select(1,'any','Select .filelist file pointing to the original data paths...');
    if(isempty(FileListPath))
        NewFileListPath = [];
        NewFileList     = [];
        disp('Quit');
        return;
    else
        disp(['Will use file list "',FileListPath,'"...']);
    end
elseif(~exist(FileListPath,'file'))
    error(['Could not find .filelist file at "',FileListPath,'"!']);
else
    disp(['Will use file list "',FileListPath,'"...']);
end
    
%ICsRemove
if(~exist('ICsRemove','var'))
    AllICs = 1:length(cellstr(spm_select('List',DualRegDir,'dr_stage2_ic\d*.nii')));
    [ICsRemove,ok] = listdlg('ListString',cellstr(num2str(AllICs')),'PromptString','Select ICs that should be removed:','CancelString','Quit','Name','ICs for removal');
    if(~ok)
        NewFileList = [];
        return;
    end
end
NICsTotal = length(cellstr(spm_select('List',DualRegDir,'dr_stage2_ic\d*.nii')));
if(~isempty(ICsRemove))
    if(any(ICsRemove>NICsTotal))
        error(['ICsRemove contains a IC number (',num2str(max(ICsRemove)),') higher than the maximum number of ICs available (',num2str(NICsTotal),')!']);
    else
        if(any(ICsRemove<1))
            error(['ICsRemove must be numbers between 1 and the maximum number of ICs available (',num2str(NICsTotal),')!']);
        end
    end
end

%RemoveRes
if(~exist('RemoveRes','var'))
    if(strcmp('Remove Residuals',questdlg('Remove residuals of dual regression (i.e. equivalent to the variability that was excluded from ICA in the initial PCA step) as well?','Remove residuals?','Remove Residuals','Do NOT Remove Residuals','Remove Residuals')))
        RemoveRes = 1;
    else
        RemoveRes = 0;
    end
elseif(isempty(RemoveRes))
    RemoveRes = 0;
else
    if(RemoveRes<1)
        RemoveRes = 0;
    else
        RemoveRes = 1;
    end
end
%inform user
if(~isempty(ICsRemove))
    if(size(ICsRemove,1)==length(ICsRemove))
        ICsRemove = ICsRemove';
    end
    disp(['Will remove ICs: ',num2str(ICsRemove)]);
    if(RemoveRes==0)
        disp('Will NOT remove residuals after dual regression...');
    else
        disp('Will remove residuals after dual regression...');
    end
else
    if(RemoveRes==1)
        disp('Will not remove any ICs, but only residuals after dual regression, i.e. variability that was not included in ICA.');
    else
        error('You have choosen to NOT remove ICs AND residuals after dual regression! Why even call this function if you do not want to do anything with it?!');
    end
end

%OutputDir
if(~exist('OutputDir','var'))
    OutputDir = spm_select(1,'dir','Select output directory...');
elseif(isempty(OutputDir))
    OutputDir = spm_select(1,'dir','Select output directory...');
else
    if(~exist(OutputDir,'dir'))
        disp(['Output directory "',OutputDir,'" does not exist, will create it.']);
        mkdir(OutputDir);
    end
end
disp(['Will save outputs to directory "',OutputDir,'".']);

%pathAltMask
if(~exist('pathAltMask','var')||isempty(pathAltMask))
    disp('Will use DualReg mask.nii file.');
    MaskFile = spm_select('FPList',DualRegDir,'mask.nii');
else
    disp(['Will use alternative mask file "',pathAltMask,'".']);
    MaskFile = pathAltMask; %spm_select('FPList',DualRegDir,'mask.nii');
end

%% get the needed files
FileList = importdata(FileListPath); %get filelist.
[~,~,TempTest] = fileparts(FileList{1}); %check extension of first entry and assume all others are the same.
if(isempty(TempTest))
    FileList            = cellfun(@(x) [x,'.nii'],importdata(FileListPath),'UniformOutput',0); %add extension ".nii"
elseif(ischar(TempTest))
    if(strcmpi(TempTest,'.nii')||strcmpi(TempTest,'.hdr'))
        disp('NIFTI-inputs recognized.');
    else
        disp(['WARNING: Input file extension "',TempTest,'" is unknown. Will try to continue, but this could lead to an error later.']);
    end
else
    disp(['WARNING: Input file extension "',TempTest,'" is unknown. Will try to continue, but this could lead to an error later.']);
end

DualRegBetasPerInput    = cellstr(spm_select('FPList',DualRegDir,'dr_stage2_subject\d*.nii'));
DualRegRegsPerInput     = cellstr(spm_select('FPList',DualRegDir,'dr_stage1_subject\d*.txt'));
if(RemoveRes==1)
    if(~isempty(spm_select('List',DualRegDir,'dr_residuals_stage2_subject\d*.nii')))
        DualRegResidualsPerInput= cellstr(spm_select('FPList',DualRegDir,'dr_residuals_stage2_subject\d*.nii'));
    else
        error(['Could not find Residuals of dual regression ("dr_residuals_stage2_subject\d*.nii" files) in folder "',DualRegDir,'". Please check. (Did you really use dual_regression_residuals or the original dual_regression?)']);
    end
else
    if(~isempty(spm_select('List',DualRegDir,'dr_residuals_stage2_subject\d*.nii')))
        DualRegResidualsPerInput= cellstr(spm_select('List',DualRegDir,'dr_residuals_stage2_subject\d*.nii'));
    else
        DualRegResidualsPerInput = 'Not available and not used, therefore no harm done.';
    end
end

%% prepare mask
VMask   = spm_vol(MaskFile); %get mask volume
Mask3D  = VMask.private.dat(:,:,:);
mask    = Mask3D(:)~=0;
NVoxMask= length(find(mask~=0));
disp([num2str(NVoxMask),' are in the mask.']);

%% basic checks of needed files
if(length(FileList)~=length(DualRegBetasPerInput))
    error('"FileList" does not match the number of files listed in "DualRegBetasPerInput"!');
end
if(length(FileList)~=length(DualRegRegsPerInput))
    error('"FileList" does not match the number of files listed in "DualRegRegsPerInput"!');
end
if(length(DualRegBetasPerInput)~=length(DualRegRegsPerInput))
    error('"DualRegBetasPerInput" does not match the number of files listed in "DualRegRegsPerInput"!');
end
if(RemoveRes==1)
    if(length(FileList)~=length(DualRegResidualsPerInput))
        error('"FileList" does not match the number of files listed in "DualRegResidualsPerInput"!');
    end
    if(length(DualRegBetasPerInput)~=length(DualRegResidualsPerInput))
        error('"DualRegBetasPerInput" does not match the number of files listed in "DualRegResidualsPerInput"!');
    end
    if(length(DualRegResidualsPerInput)~=length(DualRegRegsPerInput))
        error('"DualRegResidualsPerInput" does not match the number of files listed in "DualRegRegsPerInput"!');
    end
end
NInputs = length(FileList); %number of inputs

VRawData   = spm_vol(FileList{1});
if((VRawData(1).private.dat.dim(1)~=size(Mask3D,1))||(VRawData(1).private.dat.dim(2)~=size(Mask3D,2))||(VRawData(1).private.dat.dim(3)~=size(Mask3D,3)))
    error('Mask does not fit fMRI data!');
else
    clear VRawData %clear this again to avoid any complications later...
end
  
%% do calculations
NewFileList  = cell(length(FileList),1); %this will be filled with the new file list and then written to a txt file.
StartTimeStr = datestr(now,'yyyymmmdd_HHMM'); %starting time of cleanup 
for IndInput = 1:NInputs
    disp('---------------------------');
    disp(['Treating input ',num2str(IndInput,['%0',num2str(max([2; ceil(log10(NInputs))])),'d']),'of',num2str(NInputs,['%0',num2str(max([2; ceil(log10(NInputs))])),'d'])]);
    
    %% get spm volume structure for the input.
    disp('loading data...');
    VRawData   = spm_vol(FileList{IndInput});
    RawDataDims= VRawData(1).private.dat.dim;
       
    %% remove selected ICs for this input
    RemoveData2D = zeros(NVoxMask,RawDataDims(4)); %voxels in mask -x- time points, i.e., TRs
    if(~isempty(ICsRemove))
        disp('removing selected ICs...');
        %% get betas for this input
        VBetas = spm_vol(DualRegBetasPerInput{IndInput});
        if(VBetas(1).private.dat.dim(4)~=NICsTotal)
            error(['Dual Regression Beta for input ',num2str(IndInput),' seems to have a different number of ICs (',num2str(VBetas(1).private.dat.dim(4)),') than previously determined (',num2str(NICsTotal),') from the dual regression stage 2 files in the dual regression folder?!']);
        end
        Betas2D= reshape(VBetas(1).private.dat(:,:,:,:),[],VBetas(1).private.dat.dim(4));
        Betas2D= Betas2D(mask,:); %only use voxels in mask
        
        %% get regs for this inputs (apply design normalization --> zscore)
        RegressorsDesNorm = zscore(load(DualRegRegsPerInput{IndInput}));
        if(size(RegressorsDesNorm,1)~=RawDataDims(4))
            error(['Dual Regression Regressors for input ',num2str(IndInput),' seems to have a different number of TRs (',num2str(size(RegressorsDesNorm,1)),') than previously determined (',num2str(VRawData(1).private.dat.dim(4)),') from the 4D Raw Data input?!']);
        else
            if(size(RegressorsDesNorm,2)~=NICsTotal)
                error(['Dual Regression Regressors for input ',num2str(IndInput),' seems to have a different number of ICs (',num2str(size(RegressorsDesNorm,2)),') than previously determined (',num2str(NICsTotal),') from the dual regression stage 2 files in the dual regression folder?!']);
            end
        end
        
        %% create remove data from this
        if(NewMode)
            RemoveData2D = Betas2D(:,ICsRemove)*RegressorsDesNorm(:,ICsRemove)';
        else
            reverseStr = ''; %init reverse string empty and later set it to the length of character IN BACKSPACES ("\b"), that has been printed in the last step, in order to remove them.
            for ICrem = 1:length(ICsRemove)
                %info
                msg = sprintf('%03.1f percent done...', 100*(ICrem/length(ICsRemove))); %percentage rounded to one decimal place
                fprintf([reverseStr, msg]); %delete the last message and print current message.
                reverseStr = repmat(sprintf('\b'), 1, length(msg)); %prep the next reverse string, i.e. a string of BACKSPACES in the length of the message, i.e. delete the message on screen.
                
                RemoveData2D = RemoveData2D + repmat(Betas2D(:,ICsRemove(ICrem)),1,size(RegressorsDesNorm,1)).*repmat(RegressorsDesNorm(:,ICsRemove(ICrem)),1,NVoxMask)';
            end
            disp(' ');
        end
    end
    
    %% clear data from memory
    clear Betas2D RegressorsDesNorm
    
    %% remove residuals for this input
    if(RemoveRes==1)
        disp('removing residuals...');
        %% get residuals data
        VRes = spm_vol(DualRegResidualsPerInput{IndInput});
        if((VRes(1).private.dat.dim(1)~=VRawData(1).private.dat.dim(1))||(VRes(1).private.dat.dim(2)~=VRawData(1).private.dat.dim(2))||(VRes(1).private.dat.dim(3)~=VRawData(1).private.dat.dim(3))||(VRes(1).private.dat.dim(4)~=VRawData(1).private.dat.dim(4)))
            error(['Dimensions of residuals (',num2str(VRes(1).private.dat.dim(1)),',',num2str(VRes(1).private.dat.dim(2)),',',num2str(VRes(1).private.dat.dim(3)),',',num2str(VRes(1).private.dat.dim(4)),') does not match the dimensions of the input data (',num2str(VRawData(1).private.dat.dim(1)),',',num2str(VRawData(1).private.dat.dim(2)),',',num2str(VRawData(1).private.dat.dim(3)),',',num2str(VRawData(1).private.dat.dim(4)),')!']);
        end
        
        %% include this in remove data
        Res2D = reshape(VRes(1).private.dat(:,:,:,:),[],RawDataDims(4));
        Res2D = Res2D(mask,:);
        RemoveData2D = RemoveData2D + Res2D;
        clear Res2D %free memory
    end
    
    %% remove selected ICs & residuals == remove data from raw data
    disp('Adjusting input data...');
    RawData2D = reshape(VRawData(1).private.dat(:,:,:,:),[],RawDataDims(4));
    RawData2D = RawData2D(mask,:);
    NewData2D = RawData2D - RemoveData2D;
    
    %% clear data from memory
    clear RawData2D RemoveData2D 
    
    %% apply mask i.e. bring it back to 4D from 2D
    disp('applying mask...');
    NewData = zeros(prod(RawDataDims(1:3)),RawDataDims(4));
    NewData(mask,:) = NewData2D;
    clear NewData2D %free memory
    NewData = reshape(NewData,RawDataDims);
    
    %% make new volume structure for output in OutputDir (file name "denoise_output00XXX.nii")
    disp('writing out result...');
    Vout = rmfield(VRawData,'private'); %init and remove reference
    Order      = ceil(log10(length(Vout))); %number of zeros to add for output text
    reverseStr = ''; %init reverse string empty and later set it to the length of character IN BACKSPACES ("\b"), that has been printed in the last step, in order to remove them.
    for IndFrame = 1:length(Vout)
        Vout(IndFrame).fname = [OutputDir,filesep,'denoise_output',num2str(IndInput,'%05d'),'.nii']; %all get the same name as it is 4D
        if(Vout(IndFrame).dt(1)<16) %check to make sure data is encoded high enough.
            Vout(IndFrame).dt(1) = 16;
        end
        if(Vout(IndFrame).n(1)~=IndFrame)
            error(['4D-Volume index mismatch at Frame/TR ',num2str(IndFrame),' for Input ',num2str(IndInput)]);
        end
        Vout(IndFrame).descrip = ['denoise_output ',StartTimeStr];
        
        %info
        msg = sprintf(['Writting volume %0',num2str(Order),'d of %0',num2str(Order),'d Volumes.'], IndFrame, length(Vout)); %Don't forget this semicolon
        fprintf([reverseStr, msg]); %delete the last message and print current message.
        reverseStr = repmat(sprintf('\b'), 1, length(msg)); %prep the next reverse string, i.e. a string of BACKSPACES in the length of the message, i.e. delete the message on screen.
        
        %% write out each volume for this input (4D) to file "denoise_output00XXX.nii" in OutputDir
        spm_write_vol(Vout(IndFrame),NewData(:,:,:,IndFrame));
    end
    disp(' ');
    
    %% assign file path to NewFileList
    NewFileList{IndInput} = [OutputDir,filesep,'denoise_output',num2str(IndInput,'%05d'),'.nii'];
    
    %% next
    disp(' ');
end
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++');

%% write new file list
[NewFileListPath,FileListName] = CreateNewFileListPath(OutputDir); %try if there is already a filelist, if there is then extend the filename.
disp(['Writing out new file list "',FileListName,'" to directory "',OutputDir,'".']);
fileID = fopen(NewFileListPath,'w'); %open file for new filelist for writing
for IndInput = 1:NInputs
    %%write line to new denoise_output.filelist file in OutputDir
    fprintf(fileID,'%s\n',NewFileList{IndInput}); %write line to file denoise_output.filelist in OutputDir
end
%% close file "denoise_output.filelist"
fclose(fileID);

%% save denoise info structure
disp('Done with denoising each input, will write out "DeNoiseInfoStruct" as well...');
DeNoiseInfoStruct.StartTimeStr             = StartTimeStr;
DeNoiseInfoStruct.DualRegDir               = DualRegDir;
DeNoiseInfoStruct.FileListPath             = FileListPath;
DeNoiseInfoStruct.FileList                 = FileList;
DeNoiseInfoStruct.NInputs                  = NInputs;
DeNoiseInfoStruct.DualRegBetasPerInput     = DualRegBetasPerInput;
DeNoiseInfoStruct.DualRegRegsPerInput      = DualRegRegsPerInput;
DeNoiseInfoStruct.DualRegResidualsPerInput = DualRegResidualsPerInput;
DeNoiseInfoStruct.ICsRemove                = ICsRemove;
DeNoiseInfoStruct.NICsTotal                = NICsTotal;
DeNoiseInfoStruct.RemoveRes                = RemoveRes;
DeNoiseInfoStruct.MaskFile                 = MaskFile;
DeNoiseInfoStruct.OutputDir                = OutputDir;
DeNoiseInfoStruct.NewFileList              = NewFileList;
DeNoiseInfoStruct.NewFileListPath          = NewFileListPath;

save([OutputDir,filesep,'DeNoiseInfoStruct.mat'],'DeNoiseInfoStruct');

%% Done.
disp(' ');
disp('ALL DONE.');
disp(' ');

end

%% subfunctions
%% 
function [NewFileListPath,FileListName] = CreateNewFileListPath(OutputDir)
FileListName= 'denoise_output.filelist';

NewFileListPath = [OutputDir,filesep,FileListName];
if(exist(NewFileListPath,'file'))
    Count = 2;
    while(exist(NewFileListPath,'file')~=0)
        disp(['FileList Name "',FileListName,'", does already exist, will raise output number suffix...']);
        FileListName    = ['denoise_output.filelist',num2str(Count)];
        NewFileListPath = [OutputDir,filesep,FileListName];
        Count = Count + 1;
    end
end


end