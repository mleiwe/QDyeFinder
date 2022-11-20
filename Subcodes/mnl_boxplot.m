function h=mnl_boxplot(data,ColumnNames,YvalueLabel)
% Marcus' version of boxplot where outliers are not recorded, and raw data
% points are plotted too
%
% Red line is the mean, red whiskers are the standard deviation, blue line
% is the median, blue box is the interquartile range, and the blue whiskers
% are the maximum and minimum values.
% 
% Inputs
% data=data points*groups set of data
% ColumnNames=Xticklabels for each group. NB Input as a cell {'Group
% 1','Group2',....}
% YvalueLabel=String containing Y axis label
%
% Outputs
% h=gcf
% In addition figure is produced

%% Base Information
sz=size(data);
%Pre-allocate
Means=zeros(sz(2),1);
Medians=zeros(sz(2),1);
StdDevs=zeros(sz(2),1);
LowQuarts=zeros(sz(2),1);
HighQuarts=zeros(sz(2),1);
Mins=zeros(sz(2),1);
Maxs=zeros(sz(2),1);
xtickval=1:sz(2);
%Calculate Values
for i=1:sz(2)
    Means(i)=nanmean(data(:,i));
    Medians(i)=nanmedian(data(:,i));
    StdDevs(i)=nanstd(data(:,i));
    LowQuarts(i)=prctile(data(:,i),25);
    HighQuarts(i)=prctile(data(:,i),75);
    Mins(i)=nanmin(data(:,i));
    Maxs(i)=nanmax(data(:,i));
end
% Create Figure
xticks(xtickval)
xticklabels(ColumnNames)
ylabel(YvalueLabel)

%% Now draw boxplots
prompt='Do you want to include the standard deviation (in red).... y/n';
result=input(prompt,'s');
for i=1:sz(2)
    if isnan(Means(i))==0
        %Draw IQRs
        drawIQRbox(i,LowQuarts(i),(HighQuarts(i)-LowQuarts(i)));
        hold on
        %Draw Median Line
        line([i-0.33 i+0.33],[Medians(i) Medians(i)],'Color','blue')
        %Draw Mean Line
        line([i-0.33 i+0.33],[Means(i) Means(i)],'Color','red')
        if result=='y'
            %Draw StdDev lines
            line([i i],[Means(i)-StdDevs(i) Means(i)+StdDevs(i)],'Color','red')%Draw main bar
            line([i-0.2 i+0.2],[Means(i)-StdDevs(i) Means(i)-StdDevs(i)],'Color','red')%Draw bottom whisker
            line([i-0.2 i+0.2],[Means(i)+StdDevs(i) Means(i)+StdDevs(i)],'Color','red')%Draw top whisker
        end
    end
end
%% Plot individual data points
prompt='Plot Individual Points? y/n';
result=input(prompt,'s');
for i=1:sz(2)
    if result=='y'
        %Beeswarm the plots
        [x,y]=mnl_BeeswarmPlots(data(:,i));
        plot(i-1+x,y,'ok','MarkerSize',2)
        hold on
    elseif result=='n'
        %Draw min
        line([i i],[LowQuarts(i) Mins(i)],'Color','blue')%Draw main bar
        line([i-0.2 i+0.2],[Mins(i) Mins(i)],'Color','blue')%Draw whisker
        %Draw max
        line([i i],[HighQuarts(i) Maxs(i)],'Color','blue')%Draw main bar
        line([i-0.2 i+0.2],[Maxs(i) Maxs(i)],'Color','blue')%Draw whisker
    else
        disp('Unrecognised Input....skipping step')
    end
end
xticks(xtickval)
xticklabels(ColumnNames)
ylabel(YvalueLabel)
h=gcf;
end

function drawIQRbox(x,y,h)
x1=x-0.33;
w=0.66;
rectangle('Position',[x1 y w h],'EdgeColor','b')
end