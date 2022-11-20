function []=mnl_CumulativePlot3(varargin)
%mnl_CumulativePlot is designed to plot cumulative distribution plots of
%multiple groups to compare the differences in the spread of the data.

%Inputs
%The data needs to be arranged into different columns - the sorting is
%done automatically (This is Name1 and Name2)

cmap=colormap(jet(nargin));
for i=1:nargin
    mxVals(i)=max(varargin{i});
    mnVals(i)=min(varargin{i});
end
mxVal=max(mxVals);
mnVal=min(mnVals);
for i=1:nargin
    %% Determine the size of the Data Setc
    tempCol=[];
    Y=[];
    tempCol=varargin{i};
    s1=length(tempCol(~isnan(tempCol)));    
    if s1~=0
        
        %% Sort Data
        StempCol=sort(tempCol);
        StempCol=StempCol(1:s1,1);
        %% Get Y Values
        %ID and remove NaNs
        index=[];
        index=isnan(StempCol)~=1; %Index of real values
        TotVals=sum(index);
        n=1;
        for i2=1:s1(1)
            if index(i2)==1
                Y(n)=n*(1/TotVals)*100;
                n=n+1;
            end
        end
        if isempty(Y)==1 %Basically if all are NaNs
            Y(1:2)=0;
            StempCol=[mnVal;mxVal];
        end
        %% Plot Data
        plot(StempCol',Y,'MarkerSize',5,'LineWidth',2,'Color',cmap(i,:,:))
        hold on
    end
end
hold off
ylabel('Cumulative Percent(%)')
end