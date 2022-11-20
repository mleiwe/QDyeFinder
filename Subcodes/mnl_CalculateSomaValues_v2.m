function [Cell]=mnl_CalculateSomaValues_v2(Cell,Scale,cData,BkgMean)
%% Step 1 - Get the voxels involved in each soma
dim=size(cData);
Scount=size(Cell,2);
for i=1:Scount
    NumPoints=size(Cell(i).dxyz,1);
    for j=1:NumPoints
        %Convert the surface points to voxels
        dxyz(j,1)=round(Cell(i).dxyz(j,1));
        dxyz(j,2)=(round(Cell(i).dxyz(j,2)/Scale(1)))+1; %Add one to compensate
        dxyz(j,3)=(round(Cell(i).dxyz(j,3)/Scale(2))*-1)+1;%Flip the Y axis, add one to compensate
        dxyz(j,4)=(round(Cell(i).dxyz(j,4)/Scale(3))*-1)+1;%Flip the Z axis, add one to compensate
        if dxyz(j,1)<1
            dxyz(j,1)=1;
        end
        if dxyz(j,2)>dim(1)
            dxyz(j,2)=dim(1);
        elseif dxyz(j,2)<1
            dxyz(j,2)=1;
        end
        if dxyz(j,3)>dim(2)
            dxyz(j,3)=dim(2);
        elseif dxyz(j,3)<1
            dxyz(j,3)=1;
        end
        if dxyz(j,4)>dim(4)
            dxyz(j,4)=dim(4);
        elseif dxyz(j,4)<1
            dxyz(j,4)=1;
        end
    end
    Zlist=unique(dxyz(:,4));
    numZ=size(Zlist,1);
    VxN=1;
    for j=1:numZ
        idx=find(dxyz(:,4)==Zlist(j));
        %szI=size(idx,1);
        xPos=dxyz(idx,2);
        yPos=dxyz(idx,3);
        %Create the Mask
        maxlims=max(dxyz);
        BW=poly2mask(yPos,xPos,maxlims(2)+10,maxlims(3)+10); %Add a little padding to be safe
        %Find the Key Voxels
        [r,c]=find(BW==1);
        NumVx=size(r,1);
        TempVxList=ones(NumVx,3);
        TempVxList(:,3)=TempVxList(:,3).*dxyz(1,4);
        TempVxList(:,1)=r;
        TempVxList(:,2)=c;
        %Now Assign them to the full Vx List
        VxN2=VxN-1+NumVx;
        VxList(VxN:VxN2,:)=TempVxList;
        VxN=VxN2+1;        
        clear BW maxlims r c TempVxList x y    
    end
    %Put the Vx into the structure
    Cell(i).VxList=VxList;
    Cell(i).Scale=Scale;
    Cell(i).dim=size(cData);
    clear VxList dxyz
end
%% Step 2 - Now calculate the Mean Colour Values of the Soma, and also the BkgSub
ImDim=size(cData);
for i=1:Scount
    VxList=Cell(i).VxList;
    nVx=size(VxList,1);
    ValList=nan(nVx,ImDim(3));
    for j=1:nVx
        ValList(j,:)=cData(VxList(j,2),VxList(j,1),:,VxList(j,3));
    end
    Cell(i).RawMean=mean(ValList,'omitnan');
    %Quick switching to 0 if negative
    smV=mean(ValList,'omitnan')-BkgMean;
    for j=1:dim(3)
        if smV(j)<0
            smV(j)=0;
        end
    end
    Cell(i).BkgSubMean=smV;
    CellMatrixRaw(i,:)=mean(ValList,'omitnan');
    CellMatrixSub(i,:)=smV;
end
%% Step 3 - Max Normalise and then Vector Normalise
MxVal=max(CellMatrixSub);
CellMatrixMxNorm=CellMatrixSub./MxVal;
[VecNormMatrix]=mnl_NormaliseVectors(CellMatrixMxNorm);
for i=1:Scount
    Cell(i).NormMean=CellMatrixMxNorm(i,:);
    Cell(i).VecNormMean=VecNormMatrix(i,:);
end
figure('Name','Cell Soma Colour Values')
subplot(1,4,1)
imagesc(CellMatrixRaw)
title('Raw Mean')
subplot(1,4,2)
imagesc(CellMatrixSub)
title('Bkg Sub')
subplot(1,4,3)
imagesc(CellMatrixMxNorm,[0 1])
title('Max Normalised')
subplot(1,4,4)
imagesc(VecNormMatrix,[0 1])
title('Vector Normalised')
colormap(magma)
if ImDim(3)==7
    figure('Name','Trial RGB representation')
    [testMatrix]=mnl_RepresentInRGB(VecNormMatrix);
    image(testMatrix)    
end
end
function [Matrix]=mnl_RepresentInRGB(InputMatrix)
%Establish the chosen scoring system
ChannelColours(1,:)=[125 0 255];
ChannelColours(2,:)=[0 0 255];
ChannelColours(3,:)=[0 255 255];
ChannelColours(4,:)=[0 255 0];
ChannelColours(5,:)=[255 255 0];
ChannelColours(6,:)=[255 125 0];
ChannelColours(7,:)=[255 0 0];
%
szIm=size(InputMatrix);
Matrix=zeros(szIm(1),szIm(2),3); %Make an RGB matrix
for i=1:szIm(1)
    for j=1:szIm(2)
        Value=InputMatrix(i,j);
        Matrix(i,j,:)=Value*ChannelColours(j,:);
    end
end
Matrix=Matrix/255;
end