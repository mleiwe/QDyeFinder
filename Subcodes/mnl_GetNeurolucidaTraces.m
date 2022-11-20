function [Trace]=mnl_GetNeurolucidaTraces(s)
NumTraces=size(s.mbf.tree,2);
Tree=s.mbf.tree;
nT=1; %Final Number of Fragments/Traces
counter=1;
while counter<=NumTraces
    TempTrace=Tree{1,counter};
    %% Are there points?
    if isfield(TempTrace,'point')==1
        [points,diameter]=mnl_GetPoints(TempTrace);
        Trace(nT).Points=points;
        Trace(nT).Diameter=diameter;
        Trace(nT).OriginalTrace=counter;
        TypeId=TempTrace.Attributes.type;
        Trace(nT).TypeId=TypeId;
        a=isfield(TempTrace,'property');
        if a==1
            sz=size(TempTrace.property,2);
            if sz>1
                SetId=TempTrace.property{1,1}.s.Text;
            else
                SetId=TempTrace.property.s.Text;
            end
            Trace(nT).SetId=SetId;
        else
            SetId=[];
            Trace(nT).SetId=[];
        end
        clear points diameter
        nT=nT+1;
    end
    %Allocate Checker
    CheckBranch=isfield(TempTrace,'branch');
    %% Are there branches
    if CheckBranch==1
        [FinalBranches]=mnl_GetBranches(TempTrace);%Get Branches
        %NB Here FinalBranches is my desired structure, while Branches are a cell array
        szFB=size(FinalBranches,2);
        for i=1:szFB
            Trace(nT).Points=FinalBranches(i).Points;
            Trace(nT).Diameter=FinalBranches(i).Diameter;
            Trace(nT).SetId=SetId;
            Trace(nT).TypeId=TypeId;
            Trace(nT).OriginalTrace=counter;
            nT=nT+1;
        end
    end
    counter=counter+1;
end
end
%% Subfunctions
%% Get the Points
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
%% Get the Branches
function [FinalBranches]=mnl_GetBranches(TempTrace)
nB=size(TempTrace.branch,2); %The number of branches
%% For a standard branch set up
FinalCounter=1;
for i=1:nB
    InitialBranchCounter=FinalCounter;
    %Allocate branches
    if nB==1
        ChosenBranch=TempTrace.branch;
    else
        ChosenBranch=TempTrace.branch{1,i};
    end
    CheckPoints(i)=isfield(ChosenBranch,'point');
    CheckBranch(i)=isfield(ChosenBranch,'branch');
    %% Check for Points
    if CheckPoints(i)==1 
        %Have Points?
        [points,diameter]=mnl_GetPoints(ChosenBranch);
        FinalBranches(FinalCounter).Points=points;
        FinalBranches(FinalCounter).Diameter=diameter;
        %Search for the Id if present
        if isfield(ChosenBranch,'property')==1
            szB=size(ChosenBranch.property,2);
            if szB>1
                FinalBranches(FinalCounter).SetId=ChosenBranch.property{1,1}.s.Text;
            else
                FinalBranches(FinalCounter).SetId=ChosenBranch.property.s.Text;
            end
        else
            FinalBranches(FinalCounter).SetId=[];
        end
        FinalCounter=FinalCounter+1;
    end   
    %% Are there branches?
    if CheckBranch(i)==1
        [FinalTwigs]=mnl_GetBranches(ChosenBranch);
        sz1=size(FinalTwigs,2);
        for j=1:sz1
            FinalBranches(FinalCounter).SetId=FinalTwigs(j).SetId;
            FinalBranches(FinalCounter).Points=FinalTwigs(j).Points;
            FinalBranches(FinalCounter).Diameter=FinalTwigs(j).Diameter;
            FinalCounter=FinalCounter+1;
        end
    end
end
end