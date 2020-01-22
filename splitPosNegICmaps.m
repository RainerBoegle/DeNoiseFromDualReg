function [pathMelodic_IC_splitPosNeg,varargout] = splitPosNegICmaps(pathMelodic_IC,useAbs)
% This function will split the spatial ICs in a melodic_IC.nii file into
% positive and negative components.
% An input melodic_IC.nii with N Components will be split into 2N
% Components with the first N being the positive amplitudes of the N
% components and zeros for the negative parts and the next N are the
% negative amplitudes with zeros for the positive parts.
%
% As an extra option, each of the final 2N components (N positive, N
% negative) can be transformed to be strictly positive, i.e., applying the
% absolute function.
%
%Usage:
%       pathMelodic_IC_splitPosNeg = splitPosNegICmaps(pathMelodic_IC,useAbs);
%       pathMelodic_IC_splitPosNeg = splitPosNegICmaps(pathMelodic_IC);
%       pathMelodic_IC_splitPosNeg = splitPosNegICmaps(pathMelodic_IC,0); %same as above.
%       pathMelodic_IC_splitPosNeg = splitPosNegICmaps(pathMelodic_IC,1); %after splitting positive and negative parts apply abs function.
%
%NB: If there are no positive or negative values for a certain component,
%    then the resulting split will create a component with only zeros, this
%    has to be be avoided in dual regression, therefore this function will
%    also output a copy of the split with the zero components removed, -automatically.
%    
%    [pathMelodic_IC_splitPosNeg,pathMelodic_IC_splitPosNegWithOutZeroComps,nonZeroInds] = splitPosNegICmaps(pathMelodic_IC,useAbs);
%
%
%Author: Rainer.Boegle@googlemail.com
%V1.0
%Comment: V1.0(03.08.2018): initial implementation

%% check inputs
assert(~isempty(pathMelodic_IC),'Error: no input file provided!');
assert(ischar(pathMelodic_IC),'Error: input path must be a char-vector indicating the file location!');
if(~exist(pathMelodic_IC,'file'))
    [baseDir,fName,ext] = fileparts(pathMelodic_IC);
    if(~isempty(regexp(ext,',\d'))) %spm input
        pathMelodic_IC = [baseDir,filesep,fName,ext(1:(regexp(ext,',')-1))];
    end
end 
assert(exist(pathMelodic_IC,'file')~=0,['Error: could not find "',getFName(pathMelodic_IC),'" in directory "',fileparts(pathMelodic_IC),'".']);
[baseDir,fName,ext] = fileparts(pathMelodic_IC);

if(~exist('useAbs','var'))
    useAbs = 0;
end

%% get 4D data
disp(['Loading "',fName,ext,'" from "',baseDir,'"...']);
vols = spm_vol(pathMelodic_IC);
data4D = spm_read_vols(vols);
nICs = size(data4D,4);

%% generate new volume data
disp(['Splitting ',num2str(nICs),' ICs into ',num2str(2*nICs),' positive and negative components...']);
newData4D = zeros(size(data4D,1),size(data4D,2),size(data4D,3),2*nICs);
newData4D(:,:,:,1:nICs) = data4D.*double(sign(data4D)>0);
newData4D(:,:,:,nICs+(1:nICs)) = data4D.*double(sign(data4D)<0);

%% abs?
if(useAbs)
    disp('Applying ABS...');
    newData4D = abs(newData4D);
    absStr = 'ABS';
else
    absStr = '';
end

%% check data
zeroInds = [];
for indIC = 1:(2*nICs)  
    currData = newData4D(:,:,:,indIC);
    if(all(currData(:)))
        zeroInds = [zeroInds; indIC];
        if(indIC<=nICs)
            disp(['Warning(Index ',num2str(indIC),'): Component ',num2str(indIC),' (positive) has no non-zero Voxels!!!']);
        else
            disp(['Warning(Index ',num2str(indIC),'): Component ',num2str(indIC-nICs),' (negative) has no non-zero Voxels!!!']);
        end
    end
end

%% write out new 4D data
disp(['Writting "',fName,'_splitPosNeg',absStr,ext,'" to "',baseDir,'"...']);
pathMelodic_IC_splitPosNeg = [baseDir,filesep,fName,'_splitPosNeg',absStr,ext];
vOut = rmfield(vols(1),'private');
vOut.fname = pathMelodic_IC_splitPosNeg;
if(vOut.dt(1)<16)
    vOut.dt(1) = 16;
end
Order      = ceil(log10(2*nICs)); %number of zeros to add for output text
reverseStr = ''; %init reverse string empty and later set it to the length of character IN BACKSPACES ("\b"), that has been printed in the last step, in order to remove them.
for indIC = 1:(2*nICs)  
    %info
    msg = sprintf(['Writting volume %0',num2str(Order),'d of %0',num2str(Order),'d Volumes. (1:',num2str(nICs),'==Pos && Neg==',num2str(nICs+1),':',num2str(2*nICs),')'], indIC, 2*nICs); %Don't forget this semicolon
    fprintf([reverseStr, msg]); %delete the last message and print current message.
    reverseStr = repmat(sprintf('\b'), 1, length(msg)); %prep the next reverse string, i.e. a string of BACKSPACES in the length of the message, i.e. delete the message on screen.
    
    %write
    vOut.n(1) = indIC;
    spm_write_vol(vOut,newData4D(:,:,:,indIC));
end

%% any empty components?
if(~isempty(zeroInds))
    nonZeroInds = (1:2*nICs)';
    nonZeroInds(zeroInds) = [];
    newData4D(:,:,:,zeroInds) = [];
    
    disp(['Writting "',fName,'_splitPosNeg',absStr,'WithoutZeroComps',ext,'" to "',baseDir,'"...']);
    varargout{1} = [baseDir,filesep,fName,'_splitPosNeg',absStr,'WithoutZeroComps',ext];
    vOut = rmfield(vols(1),'private');
    vOut.fname = varargout{1};
    if(vOut.dt(1)<16)
        vOut.dt(1) = 16;
    end
    Order      = ceil(log10(size(newData4D,4))); %number of zeros to add for output text
    reverseStr = ''; %init reverse string empty and later set it to the length of character IN BACKSPACES ("\b"), that has been printed in the last step, in order to remove them.
    for indIC = 1:size(newData4D,4)
        %info
        msg = sprintf(['Writting volume %0',num2str(Order),'d of %0',num2str(Order),'d Volumes.'], indIC, size(newData4D,4)); %Don't forget this semicolon
        fprintf([reverseStr, msg]); %delete the last message and print current message.
        reverseStr = repmat(sprintf('\b'), 1, length(msg)); %prep the next reverse string, i.e. a string of BACKSPACES in the length of the message, i.e. delete the message on screen.
        
        %write
        vOut.n(1) = indIC;
        spm_write_vol(vOut,newData4D(:,:,:,indIC));
    end
    
    varargout{2} = nonZeroInds;
end

%% done.
disp(' ');
disp('Done.');
disp(' ');
end

%% subfunction
function fName = getFName(pathStr)
% get only fileName

[~,fName,ext] = fileparts(pathStr);
if(~isempty(ext))
    fName = [fName,ext];
end
end