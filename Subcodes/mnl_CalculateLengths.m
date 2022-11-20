function [Trace]=mnl_CalculateLengths(Trace,Scale)
%Calculate how long each Trace is
%
% Inputs
% Trace - Trace structure
% Scale - The scale
%
% Outputs
% Trace - Trace Structure with distances added

szT=size(Trace);
for i=1:szT(2)
    clear PointsList
    clear VoxList
    PointsList=Trace(i).Points;
    if isempty(PointsList)==0
        %Calculate in microns
        [Length_real]=mnl_MeasureDistaceOfPoints(PointsList);
        Trace(i).LengthReal=Length_real;
        %Calculate in voxel distances
        VxList=PointsList./Scale;
        [Length_vx]=mnl_MeasureDistaceOfPoints(VxList);
        Trace(i).LengthVx=Length_vx;
    else
        Trace(i).LengthReal=[];
        Trace(i).RealPoints=[];
        Trace(i).LengthVx=[];
    end
end
end
%% Nested Functions
function [dist]=mnl_MeasureDistaceOfPoints(PointsList)
nPoints=size(PointsList,1);
dist=0;
for i=1:nPoints-1
    Pos1=PointsList(i,:);
    Pos2=PointsList(i+1,:);
    [D]=mnl_EuclideanDistance(Pos1,Pos2);
    distPoints(i)=D;
    dist=dist+D;   
end
end
function [EuD]=mnl_EuclideanDistance(X,Y)
EuD=sqrt(sum((X-Y).^2));
end