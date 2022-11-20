function h=mnl_boxplot2(data,ColumnNames,YvalueLabel,IncSD,IncIP)
% Marcus' version of boxplot where outliers are not recorded, and raw data
% points are plotted too
%
% Red line is the mean, red whiskers are the standard deviation, blue line
% is the median, blue box is the interquartile range, and the blue whiskers
% are the maximum and minimum values.
% 
% Inputs
% data - data points*groups set of data
% ColumnNames - Xticklabels for each group. NB Input as a cell {'Group 1','Group2',....}
% YvalueLabel - String containing Y axis label
% IncSD - include SD in the plots, y/n
% IncIP - include individual points in the plots, y/n
%
% Outputs
% h=gcf
% Marcus Leiwe, Kyushu University, Dec 2021

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
    Means(i)=mean(data(:,i),'omitnan');
    Medians(i)=median(data(:,i),'omitnan');
    StdDevs(i)=std(data(:,i),'omitnan');
    LowQuarts(i)=prctile(data(:,i),25);
    HighQuarts(i)=prctile(data(:,i),75);
    Mins(i)=min(data(:,i),[],'omitnan');
    Maxs(i)=max(data(:,i),[],'omitnan');
end
% Create Figure
xticks(xtickval)
xticklabels(ColumnNames)
ylabel(YvalueLabel)
%% Plot individual data points
for i=1:sz(2)
    if strcmp(IncIP,'y')==1
        %Beeswarm the plots
        [x,y]=mnl_BeeswarmPlots(data(:,i));
        plot(i-1+x,y,'.k','MarkerSize',1)
        hold on
    elseif strcmp(IncIP,'n')==1
        %Draw min
        line([i i],[LowQuarts(i) Mins(i)],'Color','blue','LineWidth',2)%Draw main bar
        line([i-0.2 i+0.2],[Mins(i) Mins(i)],'Color','blue','LineWidth',2)%Draw whisker
        %Draw max
        line([i i],[HighQuarts(i) Maxs(i)],'Color','blue','LineWidth',2)%Draw main bar
        line([i-0.2 i+0.2],[Maxs(i) Maxs(i)],'Color','blue','LineWidth',2)%Draw whisker
    else
        disp('Unrecognised Input....skipping step')
    end
end

%% Now draw boxplots
for i=1:sz(2)
    if ~isnan(Means(i))
        %Draw IQRs
        drawIQRbox(i,LowQuarts(i),(HighQuarts(i)-LowQuarts(i)));
        hold on
        %Draw Median Line
        line([i-0.33 i+0.33],[Medians(i) Medians(i)],'Color','blue','LineWidth',2)
        %Draw Mean Line
        line([i-0.33 i+0.33],[Means(i) Means(i)],'Color','red','LineWidth',2)
        if strcmp(IncSD,'y')==1
            %Draw StdDev lines
            line([i i],[Means(i)-StdDevs(i) Means(i)+StdDevs(i)],'Color','red','LineWidth',2)%Draw main bar
            line([i-0.2 i+0.2],[Means(i)-StdDevs(i) Means(i)-StdDevs(i)],'Color','red','LineWidth',2)%Draw bottom whisker
            line([i-0.2 i+0.2],[Means(i)+StdDevs(i) Means(i)+StdDevs(i)],'Color','red','LineWidth',2)%Draw top whisker
        end
    end
end

xticks(xtickval)
xticklabels(ColumnNames)
ylabel(YvalueLabel)
h=gcf;
xlim([0 sz(2)+1])
end

function drawIQRbox(x,y,h)
x1=x-0.33;
w=0.66;
rectangle('Position',[x1 y w h],'EdgeColor','b','LineWidth',2)
end