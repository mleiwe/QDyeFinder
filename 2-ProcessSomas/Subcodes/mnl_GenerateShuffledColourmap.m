function [cmap]=mnl_GenerateShuffledColourmap(n)
%Generates a shuffled hsv map so that each point would be a unique colour
% n=number of colours wanted
%
% cmap - the colourmap outputted
tcmap=colormap(jet(n));
% h=gcf;
% close(h)
order=randperm(n);
cmap=nan(n,3);
for i=1:n
    cmap(i,:)=tcmap(order(i),:);
end
end
