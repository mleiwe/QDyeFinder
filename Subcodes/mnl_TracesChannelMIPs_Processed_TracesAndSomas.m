function []=mnl_TracesChannelMIPs_Processed_TracesAndSomas(dim,efPxTrace,Somas)
%function that creates MIPs for each channel based on the traces from
%neurolucida
% Inputs
%  dim - dimensions of the original image [x y c z]
%  efPxTrace - structure containg the relevant inforation
%  data - the image matrix
nTrace=size(efPxTrace,2);
nSoma=size(Somas,2);
hMeanRaw=figure('Name','Mean Raw Values');
hBkgSub=figure('Name','BkgSub Values');
hMxNorm=figure('Name','MaxNorm Raw Values');
hVecNorm=figure('Name','VecNorm Mean Raw Values');
%Per Channel
for i=1:dim(3)
    BW_RawMean=zeros(dim(2),dim(1)); %Create a BW mask
    BW_BkgSub=zeros(dim(2),dim(1)); %Create a BW mask
    BW_MxNorm=zeros(dim(2),dim(1)); %Create a BW mask
    BW_VecNorm=zeros(dim(2),dim(1)); %Create a BW mask
    AllVx=[]; %The list of the key xy positions
    %Per Trace - get the important voxels
    for j=1:nTrace
        tAllVx=efPxTrace(j).AllVoxels(:,1:2);
        fAllVx=unique(tAllVx,'rows');
        numVx=size(fAllVx,1);
        RawMeanVal=efPxTrace(j).RawMean(i);
        BkgSubVal=efPxTrace(j).BkgRemovedMean(i);
        MaxNormVal=efPxTrace(j).NormMean(i);
        VecNormVal=efPxTrace(j).VecNormMean(i);
        for k=1:numVx
            %Raw Mean
            tRawVal=BW_RawMean(fAllVx(k,2),fAllVx(k,1));
            if tRawVal<RawMeanVal
                BW_RawMean(fAllVx(k,2),fAllVx(k,1))=RawMeanVal;
            end
            %BkgSub
            tBkgSubVal=BW_BkgSub(fAllVx(k,2),fAllVx(k,1));
            if tBkgSubVal<BkgSubVal
                BW_BkgSub(fAllVx(k,2),fAllVx(k,1))=BkgSubVal;
            end
            %Max Norm
            tMaxNormVal=BW_MxNorm(fAllVx(k,2),fAllVx(k,1));
            if tMaxNormVal<MaxNormVal
                BW_MxNorm(fAllVx(k,2),fAllVx(k,1))=MaxNormVal;
            end
            %VecNorm
            tVecNorm=BW_VecNorm(fAllVx(k,2),fAllVx(k,1));
            if tVecNorm<VecNormVal
                BW_VecNorm(fAllVx(k,2),fAllVx(k,1))=VecNormVal;
            end
        end
        clear tAllVx fAllVx
    end
    %Per Soma - get the important voxels
    for j=1:nSoma
        tAllVx=Somas(j).VxList(:,1:2);
        fAllVx=unique(tAllVx,'rows');
        numVx=size(fAllVx,1);
        RawMeanVal=Somas(j).RawMean(i);
        BkgSubVal=Somas(j).BkgSubMean(i);
        MaxNormVal=Somas(j).NormMean(i);
        VecNormVal=Somas(j).VecNormMean(i);
        for k=1:numVx
            %Raw Mean
            tRawVal=BW_RawMean(fAllVx(k,2),fAllVx(k,1));
            if tRawVal<RawMeanVal
                BW_RawMean(fAllVx(k,2),fAllVx(k,1))=RawMeanVal;
            end
            %BkgSub
            tBkgSubVal=BW_BkgSub(fAllVx(k,2),fAllVx(k,1));
            if tBkgSubVal<BkgSubVal
                BW_BkgSub(fAllVx(k,2),fAllVx(k,1))=BkgSubVal;
            end
            %Max Norm
            tMaxNormVal=BW_MxNorm(fAllVx(k,2),fAllVx(k,1));
            if tMaxNormVal<MaxNormVal
                BW_MxNorm(fAllVx(k,2),fAllVx(k,1))=MaxNormVal;
            end
            %VecNorm
            tVecNorm=BW_VecNorm(fAllVx(k,2),fAllVx(k,1));
            if tVecNorm<VecNormVal
                BW_VecNorm(fAllVx(k,2),fAllVx(k,1))=VecNormVal;
            end
        end
        clear tAllVx fAllVx
    end
    %Make Tiffs
    fn=sprintf('%s%d%s','TracesMIP_RawMean_Channel',i,'.tiff');
    imwrite(uint16(BW_RawMean*2^5),fn);
    
    fn=sprintf('%s%d%s','TracesMIP_BkgSub_Channel',i,'.tiff');
    imwrite(uint16(BW_BkgSub*(2^5)),fn);
    
    fn=sprintf('%s%d%s','TracesMIP_MxNorm_Channel',i,'.tiff');
    imwrite(BW_MxNorm,fn);
    
    fn=sprintf('%s%d%s','TracesMIP_VecNorm_Channel',i,'.tiff');
    imwrite(BW_VecNorm,fn);
    
    %plot in MATLAB
    figure(hMeanRaw)
    subplot(ceil(sqrt(dim(3))),ceil(sqrt(dim(3))),i)
    imagesc(BW_RawMean)
    axis equal
    axis off
    tn=sprintf('%s%d','Channel ',i);
    title(tn)
    
    figure(hBkgSub)
    subplot(ceil(sqrt(dim(3))),ceil(sqrt(dim(3))),i)
    imagesc(BW_BkgSub)
    axis equal
    axis off
    tn=sprintf('%s%d','Channel ',i);
    title(tn)
    
    figure(hMxNorm)
    subplot(ceil(sqrt(dim(3))),ceil(sqrt(dim(3))),i)
    imagesc(BW_MxNorm)
    axis equal
    axis off
    tn=sprintf('%s%d','Channel ',i);
    title(tn)
    
    figure(hVecNorm)
    subplot(ceil(sqrt(dim(3))),ceil(sqrt(dim(3))),i)
    imagesc(BW_VecNorm)
    axis equal
    axis off
    tn=sprintf('%s%d','Channel ',i);
    title(tn)
end
end
        
        
        
        
        


