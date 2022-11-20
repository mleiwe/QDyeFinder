function mnl_PlotClustersFromMultipleClusters(efPxTrace,Cmap,FinalClusterIDs)
figure('Name','Clusters Per Section')
szT=length(efPxTrace);
ImNum=nan(szT,1);
for i=1:szT
    ImNum(i)=efPxTrace(i).ImageNum;
end
nCol=max(ImNum);
for i=1:szT
    subplot(1,nCol,ImNum(i))
    Scale=efPxTrace(i).Scale;
    TracePoints=efPxTrace(i).Points(:,1:2);
    Diameters=efPxTrace(i).Diameter./Scale(1);
    ClustN=FinalClusterIDs(i);
    mnl_PlotSingleTrace2D(TracePoints,Diameters,Cmap(ClustN,:))
    hold on
    axis equal
    axis ij
    dim=efPxTrace(i).dim;
    xlim([0 dim(1)])
    ylim([0 dim(2)])
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