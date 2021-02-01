%January 31,2021 J. Schwartzman


% Quantitative measurement of cells colonizing hydrogel particle perimeter.
% Prior to this anaysis: processed raw data (2-channel .czi files) by splitting them, and writing separate
% TIF files for the autofluorescence _AF or brightfield _BF channels.

%For each image read in, the script creates an image showing what was segmented as cells on the
%perimeter of the bead, and a file containing the following columns:

%Areaperimeter
%Areacells_perimeter
%Fraction_perimeter_area

clear
close all 
 
myFolder ='/Users/j/Desktop/Pro_segmentation/FINAL_ANALYSIS/images_from_dec';
filePatternBF = fullfile(myFolder,'*_BF.tif');
filePatternAF = fullfile(myFolder,'*_AF.tif');
theFiles = dir(filePatternBF);

%define threshold above which cell signal is visible
%for data collected in december, used a value of 10
%for data collected in january, used a value of 50
%for data collected in october, used a value of 15
threshold=10;

%Iterate through the files
for a = 1 : length(theFiles)
    FileName = theFiles(a).name;
    baseFileName = FileName(1:end-6)
    fullFileName_BF = fullfile(myFolder, FileName);
    AF=('AF.tif');
    AFname=strcat(baseFileName,AF);
    fullFileName_AF = fullfile(myFolder,AFname);
    fprintf(1, 'Now reading %s\n', baseFileName);

    %define mask for beads
    beads = double(imread(fullFileName_BF));
    beads = double(imgaussfilt(beads,2));
    clear maskb
    [Gmag,Gdir] = imgradient(beads);
    Gmag =imfill(Gmag,8,'holes');
    se = strel('disk',5);
    Gmag = imerode(Gmag,se);
    maskb =double(Gmag>mean2(Gmag));
    maskb =imfill(maskb);
    
    % remove small segmented objects (not particles)
    se = strel('disk',5);
    maskb = imerode(maskb,se);
    maskb = imdilate(maskb,se);
    clear labels
    labels=bwlabel(maskb>0);

    %Sometimes things that aren't beads are found above. This step removes
    %small objects from labels, by assigning them the value 0  
    clear stats
    stats= regionprops(labels,'Area');
        for b = 1:length(stats);
            junk=stats(b);
        if junk.Area <3000;
            labels(labels==b)=0;
        end
        end
        
    %labels now contains the area occupied by beads. 
    labels=bwlabel(labels>0);
    %this part defines the perimeter of the bead (zone) and the middle
    %(center) Signal coming from the middle can be from cells getting 
    %stuck in pockets in the beads, so we want to measure these separately.
    clear centers
    centers=labels;
    er = strel('disk',15);
    centers = imerode(centers,er);
    clear zone
    zone=labels-centers;

%export an image of the cells in the zone- this is useful to check the
%segmentation
    clear signal
    signal = double(imread(fullFileName_AF));
    zonemask=imgradient(signal);
    zonemask=double(zonemask>threshold*mean2(zonemask));
    zonemask=imfill(zonemask);
    zonemask=zone.*zonemask;
    cells_on_beads=imfuse(zonemask,beads);
    zonelabel='_perimeter.png';
    zonelabel=strcat(baseFileName,zonelabel);
    imwrite(cells_on_beads,zonelabel);
    
%%%%This section of code gets the area of cells in the perimeter zone and
%%%%the center

%Measure the area of the zone
    stats2= regionprops(zone,'Area');
    Areaperimeter=[];
    for c = 1:length(stats2);
        junk3=stats2(c);
        Areaperimeter(c) = junk3.Area;
    end   
%Quantify the area occupied by autofluorescent cells in that zone.
Areacells_perimeter=[];
    for d = 1:length(stats2)
        temp=double(zone==d);
        clear cellsignal
        cellsignal=imgradient(signal);
        cellmask=double(cellsignal>threshold*mean2(cellsignal));
        cellmask=imfill(cellmask);
        cellmask=temp.*cellmask;
        clear labels2
        labels2 = bwlabel(cellmask);    
        labels2=bwlabel(labels2>0);
    %Measure the area of each patch of cells, and the total area of the
    %bead covered by cells.
        clear stats4
        stats4= regionprops(labels2,'Area');
    %The find the area of the patches
    if length(stats4)>0
        for e = 1:length(stats4);
            junk2=stats4(e);
            if junk2.Area>1
                Areacells(e) = junk2.Area;
            end
            if junk2.Area==1
                Areacells(e)= 0;
            end
        end
    else
        Areacells=0;
    end
        Areacells_perimeter=[Areacells_perimeter sum(Areacells)];
    end
    end
    
%%%%%% This section writes the data to its own file. Two files will be written- one for the perimeters and one for the centers. 
%Individual data files are processed with a separate script to combine the
%measurements.
    clear Fraction_perimeter_area
    clear Fraction_center_area
Fraction_perimeter_area=(Areacells_perimeter./Areaperimeter);
data_end='_perimeter.csv';
name = strcat(baseFileName,data_end);
%The output is in columns, variables listed below:
titles = {'perimeter_area','cells_in_perimeter_area','fraction_of_cells_in_perimeter area'};
data = table(Areaperimeter',Areacells_perimeter',Fraction_perimeter_area','VariableNames',titles);
writetable(data,name)  
    end






