function [A,B,C] = SimulateGrating(Material,Thickness_range,GratingSize,RectangleLocation, RectangleArea,z_elevation,NA,threshold)
%Outputs two figures. A is comparison of image contrast vs specified range
%of thickness. B is an example of the amplitude of light transmissed at the
%midrange of thickness.
%Material- give numerical data or choose from list.
%thickness range is a 1-D array
%Grating size is two numbers length and width in nm
%Rectangle Location is an array with one element per rectangle. Each
%element contains the x and y coordinate for the center of the rectangle.
%Rectangle area is a 2 column by n row matrix containing the x and y size of the
%rectangles
%zElevation means if it is not sitting on the bottom layer of material, how
%much material in nm came along before it? it is one number
%NA is apurature in nm

%Author Mack Sowers
%date: June 2021
%%
%EUV only
   lambda=13.5;
   
%Define the Space
width=GratingSize(1);
height=GratingSize(2);
x = linspace(0, width, width*4);
y = linspace(0, height, height*4);

[X, Y] = meshgrid(x, y);
idealMask = ones(size(X));
wave = ones(size(X));
size_rec=RectangleArea;
%Image SIM
    F0x = 1/width; % Smallest spatial frequency
    F0y = 1/height; % Smallest spatial frequency
        
    fCutoff = 1/(lambda/NA); % Largest frequency supported by lens
    
    filtMask = sqrt((Y*F0y).^2 + (X*F0x).^2) < fCutoff;
    
    filtMask = filtMask | flipud(filtMask);
    filtMask = filtMask | fliplr(filtMask);
%%
%Populate for each thickness value
%combine rectangle area and thicknesss for use with Ryan's function

for k=1:length(Thickness_range)
    size_rec(:,3)=Thickness_range(k);
    Size_Cell{k}=size_rec;
    
    %Populate wave functions
    for n=1:4
        [wave_new, idealMask] = createRectangleTransmission(X, Y, RectangleLocation(n,:), z_elevation, Size_Cell{k}(n,:), Material);
        wave=wave.*wave_new;
    end
    %Image SIM
    objSpectrum = fft2(wave);
    filtered{k} = ifft2(objSpectrum.*filtMask);
    %Compute the contrast
    %(max - min)/(max + min)
    %For these absorbers, max is 1. in vaccum region it is 1, but since edge
    %effects will show a wave vale grater than one due to ringing, I wont even
    %ask what the max value is.
    min_sim = min(abs(filtered{k}),[],'all'); %returns the smallest element of the simulation.
    Contrast(k)=(1-min_sim)/(1);
    if Contrast(k)>=0.95
        k;
    end
end
%%
%plot thickness vs contrast
A=figure(1);
subplot(1,3,1)
plot(Thickness_range,Contrast)
axis auto
title('Image Contrast vs Material Thickness')


%plot the image simulation
% if default_thick==[]
% %find middle thickness, 
%     thick= Thickness_range(round(length(Thickness_range)/2));
% else
%     thick=default_thick;
% end
thick= round(length(Thickness_range)/2);

B=figure(2);
subplot(1,3,2)
z=abs(filtered{thick});
surf(x,y,z), shading interp
%put in overhead view
view(2);
axis image
title(sprintf('Simulated Image at thickness %dnm\n and numerical aperature %d.',Thickness_range(thick),NA))

C=figure(3);
subplot(1,3,3)
imgIntens = abs(filtered{thick}).^2;

imagesc(imageToDeprotection(imgIntens, [.210, .200], 50, 10, 40, threshold))
axis image
axis xy
title('Latent Image')
end


