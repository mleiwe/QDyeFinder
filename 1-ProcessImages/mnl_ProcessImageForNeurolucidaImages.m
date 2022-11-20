function [cData,Scale,dim]=mnl_ProcessImageForNeurolucidaImages
%The first step of the Neurolucida Imaging

% Inputs
% User prompts to define images and other properties
%
% Outputs
% -cData - The 4D image that has been corrected for chromatic abberation
% -Scale - The XYZ scale 
% -dim - The dimensions of the image.
%
% What to do next? First save all the outputs. Then auto trace the image in
% Neurolucida save the .xml file then run the code
% mnl_ProcessTracesFromNeurolucida
%
%Marcus Leiwe, Kyushu University - 21st Nov 2019

%% Load the unmixed data - NB make sure the metadata is still attached
%[Data,Scale,dim,metadata]=mnl_Load4Dimage;
[Data,Scale,dim,~]=mnl_Load4Dimage;
%mnl_MakeMaxNormalisedMIPs(Data,100);%Makes MIPs normalised to the Max value
%NB pick which one you need
%% Chromatic Abberation
prompt='Do you want to correct for chromatic aberations? y/n';
Cacorr=input(prompt,'s');
if strcmp(Cacorr,'y')==1
    %Tell the laser for each channel
    for i=1:dim(3)
        txt=sprintf('%s%d%s','Please input the Excitation Laser Wavelength for Channel ',i,' in nm');
        prompt=txt;
        LaserEx(i)=input(prompt);
    end
    %Which Channel is our base
    prompt='Which Channel should we correct to?';
    n=input(prompt);
    BaseLaser=LaserEx(n);
    %Now find the required equations from the chromatic corrections
    disp('Please Load the Chromatic Abberation Corrections File')
    [Wkspaces]=uipickfiles; %Load the Chromatic corrections file - NB This relies on the variable still being named "ChromaticCorrections"
    load(Wkspaces{1,1});
    sz1=size(ChromaticCorrections);
    for i=1:sz1(2)
        %Does the chosen 'Base Laser' match the 'ChromaticCorrections(i)'
        %Laser?
        if BaseLaser==ChromaticCorrections(i).ToWhichLaser
            WhichLasers=ChromaticCorrections(i).ForWhichLaser;
            Laser_mValues=ChromaticCorrections(i).mValues;
            Laser_cValues=ChromaticCorrections(1).cValues;
        end
    end
    if exist('WhichLasers','var')==0
        error('You do not have information for Chromatic Abberatiorns with this laser')
    end
    %Create the equations per image channel
    sz=size(Data);
    for i=1:sz(3)
        for j=1:sz1(2)
            if LaserEx(i)==WhichLasers(j)
                mValues(i)=Laser_mValues(j);
                cValues(i)=Laser_cValues(j);
            end
        end
    end
    %Now Apply the Chromatic Abberations
    [cData]=mnl_ChromaticAbberationCorrectionZ4(Data,Scale,mValues,cValues);
else
    cData=Data;
end
mnl_MakeMaxNormalisedMIPs(cData,100)%Makes MIPs normalised to the Max value
%% Create Grayscale
gsData=mnl_4DtoGrayscale('cData',cData,dim);
%Now make a gray MIP
gsMIP=max(gsData,[],3);
imwrite(gsMIP,'Grayscale_MIP.tiff');
%% Save Workspace
save('ProcessedImage.mat','cData','Scale','dim');
end