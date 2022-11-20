function mnl_ExportEPSdense(h,fn)
% Code to switch back the renderer to a vector format from a bmp (image
% format)
% Input
% h = figure number (h=gcf)
% fn = filename
%
% Output
% .eps file

set(h,'renderer','Painters')
print(fn,'-depsc','-tiff','-r300','-painters')
end