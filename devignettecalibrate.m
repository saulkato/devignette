%function devignettecalibrate(masterFolder)

%create calibration functions for devignette.m
%
%Saul Kato
%2013/04/17
%Update 2014/02/14
%
%how to use:
%put calibration images in subfolders, one folder for each brightness
%level, within masterFolder

useBorderSplinesFlag=true;  %set this to false if there isn't much good edge data, so a spline along each border wont be created
gridSize=17; %number of splines in x and y direction.

%if nargin<1
    masterFolder=pwd;
    %imagefolder='/Users/skato/Desktop/Dropbox/WormData3/Yifan130411-lowmagcalibration/20130401_01_fluorescine_5xDemag_35msExp/';
%end

%get all sub-folders within the folder
files=dir(masterFolder);
j=0;
for i=1:length(files)
    if files(i).isdir && ~strcmp(files(i).name,'.') && ~strcmp(files(i).name,'..') 
        j=j+1;        
        imageDir{j}=[filesep files(i).name];
    end
end

if j==0 
    disp('no subfolders found.  Looking for mult-image tifs');
    
    files=dir([masterFolder filesep '*.tif']);

    numLevels=length(files);
       
    for i=1:numLevels
          
        imageDir{i}='';     
        imageCount(i)=length(imfinfo(files(i).name));
        for j=1:length(imfinfo(files(i).name))
              imagefiles(i).f(j).name=files(i).name;  
        end
        
    end
    
    
else  %directories of .tif images

    numLevels=length(imageDir); 
    
    for i=1:numLevels
        imagefiles(i).f=dir([masterFolder filesep imageDir{i} filesep '*.tif']);
        imageCount(i)=length(imagefiles(i).f);
    end
    
end



%load tiffs and make mean images
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
for i=1:numLevels   %number of calibration levels
      
    for j=1:imageCount(i)  %all images for a single calibration level
        X=tiffread2([masterFolder imageDir{i} filesep imagefiles(i).f(j).name],j);
        if j==1
            if i==1
                width=X.width;
                height=X.height;
            end
            accumdata=double(X.data);
        else
            accumdata=accumdata+double(X.data);
        end
     
    end
    M(i).meanImage=uint16(accumdata/j);

    
end
%%load thresholds
k=1;
for i=1:numLevels
    if ~exist([masterFolder imageDir{i} filesep 'thresholds.txt'],'file')
        th(k)=600; th(k+1)=4000;
    else
        th=load([masterFolder imageDir{i} filesep 'thresholds.txt']);
    end
    
    lowerthresh(i)=th(k);
    upperthresh(i)=th(k+1);
    k=k+2;

end


%%load geometries
k=1;
for i=1:numLevels
    if ~exist([masterFolder imageDir{i} filesep 'geometries.txt'],'file')
        ge(k)=669;
        ge(k+1)=160;
        ge(k+2)=27;
        ge(k+3)=74;
        ge(k+4)=1;

    else
        ge=load([masterFolder imageDir{i} filesep 'geometries.txt']);
    end
    

    centerx(i)=ge(k);
    centery(i)=ge(k+1);
    radius(i)=ge(k+2);
    spacing(i)=ge(k+3);
    rotation(i)=ge(k+4);
    
    
    
    k=k+5;

end

%%make masks from thresholds
for i=1:numLevels

    M(i).mask = (M(i).meanImage > lowerthresh(i)) & (M(i).meanImage < upperthresh(i));
    
    %add hex posts here.
    
    
    %1 where there is no post [good value]
    %0 where there is a post  
    
    %add masks from geometry
    M(i).mask = (M(i).mask + 2*(1-makeHexPostMask(width,height,centerx(i),centery(i),radius(i),spacing(i),rotation(i)))) >1 ;
    
    %grow mask
    %M(i).mask = ~bwmorph(~M(i).mask,'thicken',4);  %grow mask

end 




nr=ceil(sqrt(numLevels)); %number of rows for plotting
figure('Position',[0 0 1024 1024]);

%plot mean images
for i=1:numLevels
    subtightplot(nr,nr,i);
    imagesc(M(i).meanImage);
    axis image;
    axis off;  
    drawnow;
end
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');

mtit('averaged calibration images');

figure('Position',[0 0 1024 1024]);
for i=1:numLevels
    subtightplot(nr,nr,i);
    imagesc(uint16(M(i).mask));
    axis image;
    axis off;
    
    drawnow;
end

mtit('masks');
   
    
%compute overall mean levels
for level=1:numLevels
     meanlevel(level)=sum(sum(M(level).meanImage))/sum(sum(M(level).mask));
end


%sort images by mean levels
[meanlevel_sorted sortedlevel_index] =sort(meanlevel);

% figure;plot(meanlevel_sorted,'.-');drawnow;
% title('mean brightness levels, sorted');


%%
%fill AL (All Levels) matrix with sorted mean images
AL=zeros(numLevels,width,height,'uint16');
for level=1:numLevels
    AL(level,:,:)=M(sortedlevel_index(level)).meanImage;
end

%add bad pixels manually to mask here


% for level=1:numLevels
%     j=1;
%     for i=fs(level):fe(level)
%         M(level).diag1(j,:)=diag(X.data.*uint16(XmaskGrown));
%         M(level).diag2(j,:)=diag(X.data(end:-1:1,:).*uint16(XmaskGrown(end:-1:1,:)));
%         j=j+1;
%     end
% end

%%fit smooth surfaces to each brightness level
%made from a grid of evenly spaced 1D splines
%figure('Position',[1024 0 1024 1024]);

plots=1;
%if (plots) figure('Position',[0 300 1024 1024]); end

polyOrder=2;

pixstep=floor(width/(gridSize-1));
for level=1:numLevels
    
    if (plots) figure('Position',[0 300 1024 1024]); end
    
    if useBorderSplinesFlag
        splineIndex=1:gridSize;
    else
        splineIndex=2:gridSize-1;
    end
    
    nr2=ceil(sqrt(gridSize));
    for i=splineIndex
        if (plots) subplot(nr2,nr2,i); end
        
        pixval=1+pixstep*(i-1);
        if i==gridSize
            pixval=pixval-1;
        end
        xl=AL(level,:,pixval);
        yl=AL(level,pixval,:);
        
        ALgrid(i)=pixval;
        
        splX=(xl(:)).*uint16(M(sortedlevel_index(level)).mask(:,pixval));
        splY=(yl(:)').*uint16(M(sortedlevel_index(level)).mask(pixval,:));
        [splHx{i},~,splHy{i}]=find(splX);
        [~,splVx{i},splVy{i}]=find(splY);
        
        if (plots) plot(yl(:),'b'); end
        if (plots) hold on; end        
        if (plots) plot(splVx{i},splVy{i},'.','Color',color(i,gridSize,'jet')); end

        
        pV=polyfit(double(splVx{i}'),double(splVy{i}'),polyOrder);
        
        if (plots) plot(1:width,polyval(pV,1:width),'Color',color(i,gridSize,'jet')); end
        
        ALreducedV(level,:,i)=polyval(pV,1:pixstep:height+1);
        
        %plot(splHx{i},splHy{i},'.','Color',color(i,gridSize,'jet'));

        pH=polyfit(double(splHx{i}'),double(splHy{i}'),polyOrder);
        
        %plot(1:height,polyval(pH,1:height),'Color',color(i,gridSize,'jet'));
        
        ALreducedH(level,i,:)=polyval(pH,1:pixstep:width+1);
        
        xlim([1 height]);
        
    end
% size(ALreducedH)
% size(ALreducedV)
    ALreducedM(level,:,:)=(ALreducedH(level,:,:)+ALreducedV(level,:,:))/2;
    
    %Vq = interp2(V,Xq,Yq) assumes X=1:N and Y=1:M where [M,N]=SIZE(V).

    [MGpx,MGpy]=meshgrid(1:pixstep:width+1,1:pixstep:height+1);
    [MGx,MGy] = meshgrid(1:width,1:height);
    ALsm=interp2(MGpx,MGpy,squeeze(ALreducedM(level,:,:)),MGx,MGy,'cubic');

    ALsmooth(level,:,:)=ALsm(1:width,1:height);

end

figure('Position',[500 0 1024 1024]);
for i=1:numLevels
   subtightplot(nr,nr,i);
   imagesc(squeeze(ALsmooth(i,:,:))');
   axis image;
   axis off;
end
mtit('smoothed/patched calibration images');




%% test code for visualization of smoothing
%{
for level=1:numLevels
    figure;
    subplot(3,2,1);
    imagesc(squeeze(ALreducedV(level,:,:)));

    subplot(3,2,2);
    imagesc(squeeze(ALreducedH(level,:,:)));

    subplot(3,2,3);
    imagesc(squeeze(ALreducedH(level,:,:))-squeeze(ALreducedV(level,:,:)));

    subplot(3,2,4);
    imagesc(squeeze(ALreducedM(level,:,:)));

    subplot(3,2,5);
    imagesc(squeeze(ALsmooth(level,:,:)));
end
%}


%this would be how to start a 2D surface fit
% points=find(~Xmask);
% xind=repmat(1:512,1,512);
% yind=repmat(1:512,512,1);
% yindu=yind(:)';
% sf = griddata( [xind(points),yindu(points)], M(1).meandata(points),xind,yindu,'cubic');


%% create intensity curves for each pixel 

disp('computing per-pixel interpolation curves.');
for i=1:width
    for j=1:height
         %v = double(AL(:,i,j)) \ [ones(size(meanlevel')) log(meanlevel')];
         %p1(i,j)=v(1);
         %p2(i,j)=v(2);
         
         p(:,i,j)=polyfit(meanlevel_sorted',double(ALsmooth(:,i,j)),1);
    end
    
    if mod(i,10)==1
        disp([num2str(i) '/' num2str(width) ' rows computed.']);
    end
end


%save calibration matrix
save([masterFolder filesep 'devignettecal.mat'],'p');

%% plot test spline
figure;

numPoints=20;
pixInc=floor(width/numPoints);

for j=1:numPoints
    
    i=j*pixInc;
    
    subplot(5,5,j);
    plot(meanlevel_sorted,AL(:,i,i),'r.');
    xlim([0 4500]);
    %ylim([0 4500]);

    hold on;
    %ft=log(p1(i,i))*(meanlevel).^p2(i,i);
    mv=0:meanlevel_sorted(numLevels);
    ft=polyval(p(:,i,i),mv);
    plot(mv,ft,'b--');
    
end

mtit('test correction curves');

%%plot of diagonal cross sections of mean levels
%{
figure;
for level=1:numLevels
    meanlevel(level)=sum(sum(M(level).meandata))/sum(sum(Xmask));   
    hold on;
    plot(mean(M(level).diag1,1)/(meanlevel(level)),'.','MarkerSize',1,'Color',color(level,numLevels,'jet'));
    %plot(mean(M(level).diag2,1),'r');  %other diagonal  
    k=mean(M(level).diag1,1)/meanlevel(level); 
end
%}

    