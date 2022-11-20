function [Cells]=mnl_ReadXML_Somas(fname)
%Function to quickly read all the somas only, based on the functio mnl_ReadXML_DendritesAndSomas4
filename=sprintf('%s%s',fname,'.xml');
[s]=xml2struct(filename);
%% Step 1 - Store the somas
if isfield(s.mbf,'contour')==1 %if somas exist
    nContours=size(s.mbf.contour,2); %All the contours are stored individually not by soma so we need to change the grouping
    %Find where the soma splits are
    Scount=1;
    for i=1:nContours
        if i==1
            SomaNum=s.mbf.contour{1,i}.Attributes.name;
            CellSplits(i,1)=i;
            CellNames{Scount}=SomaNum;
        end
        chk=strcmp(s.mbf.contour{1,i}.Attributes.name,SomaNum);
        if chk==0 %If it is different
            %% Double Check to see if manually traced somas are included
            if strcmp(s.mbf.contour{1,i}.Attributes.name,'CellBody')==1
                UserDefId=s.mbf.contour{1,i}.property{1,2}.s.Text;
                if strcmp(SomaNum,UserDefId)==0 %If it is different
                    CellSplits(Scount,2)=i-1;
                    %Start the next one
                    Scount=Scount+1;
                    SomaNum=UserDefId;
                    CellSplits(Scount,1)=i;
                    CellNames{Scount}=SomaNum;
                end
            else
                %% Update the CellSplits matrix
                CellSplits(Scount,2)=i-1;
                %Now Start the next one
                Scount=Scount+1;
                SomaNum=s.mbf.contour{1,i}.Attributes.name;
                CellSplits(Scount,1)=i;
                CellNames{Scount}=SomaNum;
            end
        end
        if i==nContours
            CellSplits(Scount,2)=nContours;
        end
    end
    [Cells]=mnl_ExtractSomaVoxels(s,Scount,CellSplits); %Extracts the key voxels and also the set name if specified
else
    Cells=[]; %If there are no somas then make it empty
end

%% Step Three make a 2D plot
prompt='Do you want a figure? (y/n)';
UserInput=input(prompt,'s');
FigYN=strcmp(UserInput,'y');
if FigYN==1
    %Find out how many colours to use in the plot
    nSoma=size(Cells,2);
    SetIdList={};
    SomaSetIdList={};
    nColours=0;
    if nSoma>0 %For somas
        SetIdList={};
        for i=1:nSoma
            SetId=Cells(i).SetId;
            if ~isempty(SetId) && sum(strcmp(SetId,SetIdList))==0
               nColours=nColours+1;
               SetIdList{nColours}=SetId;
               SomaSetIdList{nColours}=SetId;
            end
        end
        if ~isempty(SomaSetIdList)
            SomaSetIdList=unique(SomaSetIdList);
        end
    end
    if nColours==0
        nColours=nSoma;
    end
    [cmap]=mnl_GenerateShuffledColourmap(nColours);
    h=figure('Name','The correct traces according to Neurolucida');
    
    if isempty(Cells)==0
        for i=1:nSoma
            if isempty(Cells(i).dxyz)==0
                SetId=Cells(i).SetId;
                %Plot the soma
                Points=Cells(i).dxyz(:,2:3);
                [k,~] = convhull(Points);
                %scatter(Points(:,1),Points(:,2),'.','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:))
                 patch(Points(k,1),Points(k,2),cmap(i,:),'EdgeColor','none')
                hold on
                
            end
        end
    end
    axis equal
    xlimval=xlim;
    ylimval=ylim;
    xlim([0 xlimval(2)])
    ylim([ylimval(1) 0])
    grid on
    mnl_ExportEPSdense(h,'NeurolucidaSomas')
    savefig('NeurolucidaSomas.fig')
end
end
%% Sub-functions
function [Cell]=mnl_ExtractSomaVoxels(s,Scount,CellSplits)
%Function to group the perimeters and also asign the set name
for i=1:Scount
    st=CellSplits(i,1);
    ed=CellSplits(i,2);
    n=1;
    for j=st:ed %J is the Z plane of each Soma
        NumPoints=size(s.mbf.contour{1,j}.point,2);
        dxyz=nan(NumPoints,4);
        for k=1:NumPoints
            if NumPoints~=1
                %Extract the Surface Points of the soma
                dxyz(k,1)=str2double(s.mbf.contour{1,j}.point{1,k}.Attributes.d);
                dxyz(k,2)=str2double(s.mbf.contour{1,j}.point{1,k}.Attributes.x);
                dxyz(k,3)=str2double(s.mbf.contour{1,j}.point{1,k}.Attributes.y);
                dxyz(k,4)=str2double(s.mbf.contour{1,j}.point{1,k}.Attributes.z);
            else
                %Extract the Surface Points of the soma
                dxyz(k,1)=str2double(s.mbf.contour{1,j}.point.Attributes.d);
                dxyz(k,2)=str2double(s.mbf.contour{1,j}.point.Attributes.x);
                dxyz(k,3)=str2double(s.mbf.contour{1,j}.point.Attributes.y);
                dxyz(k,4)=str2double(s.mbf.contour{1,j}.point.Attributes.z);
            end
        end
        n2=n-1+NumPoints;
        Cell(i).dxyz(n:n2,:)=dxyz;
        clear dxyz
        n=n2+1;
        %Now assign the manual set (if present)
        if j==st
            if size(s.mbf.contour{1,j}.property,2)>=3
                if isfield(s.mbf.contour{1,j}.property{1,3}.s,'Text')==1
                    SetName=s.mbf.contour{1,j}.property{1,3}.s.Text;
                    Cell(i).SetId=SetName;
                else
                    Cell(i).SetId=[];
                end
            else
                if size(s.mbf.contour{1,j}.property,2)==2
                    if isfield(s.mbf.contour{1,j}.property{1,2},'s')==1
                        if isfield(s.mbf.contour{1,j}.property{1,2}.s,'Text')==1
                            SetName=s.mbf.contour{1,j}.property{1,2}.s.Text;
                            Cell(i).SetId=SetName;
                        else
                            Cell(i).SetId=[];
                        end
                    else
                        Cell(i).SetId=[];
                    end
                end
            end
        end
    end  
end
end
function [points,diameter]=mnl_GetPoints(TempTrace)
n=1; %Points per branch counter
NumPoints=size(TempTrace.point,2);
if NumPoints>1
    for j=1:NumPoints
        x=str2double(TempTrace.point{1,j}.Attributes.x);
        y=str2double(TempTrace.point{1,j}.Attributes.y);
        z=str2double(TempTrace.point{1,j}.Attributes.z);
        d=str2double(TempTrace.point{1,j}.Attributes.d);
        points(n,:)=[x y z];
        diameter(n)=d;
        n=n+1;
    end
elseif NumPoints==1
    x=str2double(TempTrace.point.Attributes.x);
    y=str2double(TempTrace.point.Attributes.y);
    z=str2double(TempTrace.point.Attributes.z);
    d=str2double(TempTrace.point.Attributes.d);
    points(n,:)=[x y z];
    diameter(n)=d;
end
end
function [FinalBranches,Branches]=mnl_GetBranches(TempTrace)
nB=size(TempTrace.branch,2);
c=1;
for i=1:nB
    InitialBranchCounter=c;
    %Allocate branches
    if i==1
        Branches{1,c}=TempTrace.branch;
    else
        Branches{1,c}=TempTrace.branch{1,i};
    end
    %Are there branches?
    if isfield(Branches{1,c},'branch')==1
        [FinalTwigs,Twigs]=mnl_GetBranches(Branches{1,c});
        sz1=size(Twigs,2);
        for j=1:sz1
            Branches{1,c+1}=Twigs{1,j};
            c=c+1;
        end
    else
        %Check for forks -i.e. substructures with no points or branches
        chkP=isfield(Branches{1,c},'point');
        chkB=isfield(Branches{1,c},'branch');
        if chkP==0 && chkB==0 %If it is still empty
            tc=1;
            szF=size(Branches{1,c}); %A fork should be more than one
            for j=1:szF(2)
                OriginalBranchCounter=tc;
                tBranches{1,tc}=Branches{1,c}{1,j};
                %Does it have branches
                if isfield(tBranches{1,tc},'branch')==1
                    t=tBranches{1,tc};
                    [FinalTwigs,Twigs]=mnl_GetBranches(t);
                    clear t
                    sz1=size(FinalTwigs,2);
                    for k=1:sz1
                        tBranches{1,tc+1}=FinalTwigs{1,k};
                        tc=tc+1;
                    end
                end
                %Check to see if there are points
                if isfield(tBranches{1,OriginalBranchCounter},'point')==1
                    [points,diameter]=mnl_GetPoints(tBranches{1,OriginalBranchCounter});
                    tBranches{1,OriginalBranchCounter}.point=tBranches{1,OriginalBranchCounter}.point;
                    tBranches{1,OriginalBranchCounter}.Points=points;
                    tBranches{1,OriginalBranchCounter}.Diameter=diameter;
                    tBranches{1,OriginalBranchCounter}.SetId=[];
                    clear points diameter
                end
                if OriginalBranchCounter==tc %If there were no branches we need to shift tc on by 1
                    tc=tc+1;
                end
            end
            szTemp=size(tBranches,2);
            FinalBranches{1,c}=[];
            for j=1:szTemp
                FinalBranches{1,c-1+j}.Points=tBranches{1,j}.Points;
                if isfield(tBranches{1,j},'point')==1
                    FinalBranches{1,c-1+j}.point=tBranches{1,j}.point;
                end
                FinalBranches{1,c-1+j}.Diameter=tBranches{1,j}.Diameter;
                FinalBranches{1,c-1+j}.SetId=tBranches{1,j}.SetId;
            end
            c=c+szTemp;
        end
        clear chkP chkB
    end       
    %Check for Points
    if isfield(Branches{1,InitialBranchCounter},'point')==1 %Have Points?
        [points,diameter]=mnl_GetPoints(TempTrace);
        FinalBranches{1,InitialBranchCounter}.Points=points;
        FinalBranches{1,InitialBranchCounter}.Diameter=diameter;
        FinalBranches{1,InitialBranchCounter}.points=TempTrace.point;
        if isfield(TempTrace,'property')==1
            szB=size(TempTrace.property,2);
            if szB>1
                FinalBranches{1,InitialBranchCounter}.SetId=TempTrace.property{1,1}.s.Text;
            else
                FinalBranches{1,InitialBranchCounter}.SetId=TempTrace.property.s.Text;
            end
        else
            FinalBranches{1,InitialBranchCounter}.SetId=[];
        end
    end
end
end
function mnl_PlotSingleTrace2D(TracePoints,Diameters,Colour)
%Inputs
% TracePoints - The position of the points
% Diameters - The diameter of the trace
% Colour - the RGB values that you want to plot the trace in
%% Now make the volume for each trace
PointList=TracePoints(:,1:2);
szD=size(Diameters);
[nD,p]=max(szD);
%Work if the diameters are encoded in rows (p=1) or columns (p=2)
if p==1
    for j=1:nD
        Diameter(j,1)=Diameters(j,1);
    end
elseif p==2
    for j=1:nD
        Diameter(j,1)=Diameters(1,j);
    end
end
%Left Side
Xleft=PointList(:,1)-(Diameter/2);
%Right Side
Xright=PointList(:,1)+(Diameter/2);
X=[Xleft;flipud(Xright)];
%Sort Out Y
Y=[PointList(:,2);flipud(PointList(:,2))];
%Merge
patch(X,Y,Colour,'EdgeColor','none')
end