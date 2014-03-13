%vignettetests
%create calibration functions for devignette.m
%

%Saul Kato
%04/17/13

imagedir='/Users/skato/Desktop/Dropbox/WormData3/Yifan130411-lowmagcalibration/20130401_01_fluorescine_5xDemag_35msExp/';
cd(imagedir);

reftifname='tmpStream8027.tif';

%test
%load .stk.
%load mandrill;
%X = tiffread2([imagedir '021011_kyex2865_1Mgly1mMtet_1secflicker_W13T1.stk'],190,300);
%load 021011-pointpath.xls;
%cps=X021011_pointpath;

%X= tiffread2([imagedir stackname],firstimageindex,lastimageindex);
Xref=tiffread2([imagedir reftifname]);
width=Xref.width;
height=Xref.height;


imagesc(Xref.data);
%% create hexpost mask
%1 where there is no post
%0 where there is a post
Xmask=Xref.data>1500;
figure;
imagesc(Xmask);

%grow mask by a few pixels
XmaskGrown=~bwmorph(~Xmask,'thicken',4);
figure;
imagesc(XmaskGrown);

%add bad pixels to mask


%% average frames for each brightness level
fs=[643 1863 3060 4277 5473 6680 7861];
fe=[900 2100 3290 4489 5680 6870 8060];
AL=zeros(7,512,512,'uint32');
for level=1:7
    accumdata=zeros(width,height,'uint32');
    j=1;
    for i=fs(level):fe(level)

        tifname=['tmpStream' num2str(i,'%04i') '.tif'];
        X=tiffread2([imagedir tifname]);
        accumdata=accumdata+uint32(X.data);
        M(level).diag1(j,:)=diag(X.data.*uint16(XmaskGrown));
        M(level).diag2(j,:)=diag(X.data(end:-1:1,:).*uint16(XmaskGrown(end:-1:1,:)));
        j=j+1;
    end
    M(level).meandata=accumdata/(fe(level)-fs(level)+1);
    AL(level,:,:)=M(level).meandata;
    disp(['level ' num2str(level) ' done.']);
end

%%compute mean levels

for level=1:7
    meanlevel(level)=sum(sum(M(level).meandata))/sum(sum(XmaskGrown));
end

%% fit smooth surfaces to each brightness level
%lets do a 17x17 grid by first fitting evenly spaced 1D splines
figure;
pixstep=32;
for level=1:7
    figure;
    for i=1:17
        pixval=1+pixstep*(i-1);
        if i==17
            pixval=pixval-1;
        end
        xl=AL(level,:,pixval);
        yl=AL(level,pixval,:);
        
        splX=(xl(:)).*uint32(XmaskGrown(:,pixval));
        splY=(yl(:)').*uint32(XmaskGrown(pixval,:));
        [splHx{i},~,splHy{i}]=find(splX);
        [~,splVx{i},splVy{i}]=find(splY);
        
        
        %plot(splVx{i},splVy{i},'.','Color',color(i,17,'jet'));
        hold on;
        pV=polyfit(double(splVx{i}'),double(splVy{i}'),4);
        plot(1:512,polyval(pV,1:512),'Color',color(i,17,'jet'));
        
        ALreducedV(level,:,i)=polyval(pV,1:pixstep:513);
        
        plot(splHx{i},splHy{i},'.','Color',color(i,17,'jet'));

        pH=polyfit(double(splHx{i}'),double(splHy{i}'),4);
        plot(1:512,polyval(pH,1:512),'Color',color(i,17,'jet'));
        
        ALreducedH(level,i,:)=polyval(pH,1:pixstep:513);
        
        
    end
    i=17;

    ALreducedM(level,:,:)=(ALreducedH(level,:,:)+ALreducedV(level,:,:))/2;
    
    ALsm=interp2(squeeze(ALreducedM(level,:,:)),5,'cubic');
    ALsmooth(level,:,:)=ALsm(1:512,1:512);
end


%% test code for visualization of smoothing
for level=1:7
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


%this is the start of a 2D surface fit
% points=find(~Xmask);
% xind=repmat(1:512,1,512);
% yind=repmat(1:512,512,1);
% yindu=yind(:)';
% sf = griddata( [xind(points),yindu(points)], M(1).meandata(points),xind,yindu,'cubic');


%% create intensity curves for each pixel 
for i=1:512
    for j=1:512
         %v = double(AL(:,i,j)) \ [ones(size(meanlevel')) log(meanlevel')];
         %p1(i,j)=v(1);
         %p2(i,j)=v(2);
         
         p(:,i,j)=polyfit(meanlevel',double(ALsmooth(:,i,j)),1);
    end
    if mod(i,10)==1
        disp(i)
    end
end

%%
save([imagedir '/devignettecal.mat'],'p');

%% plot test spline
figure;
j=1;

for i=1:20:500
    
subplot(5,5,j);
plot(meanlevel,AL(:,i,i),'r.');
xlim([0 4500]);
%ylim([0 4500]);

hold on;
%ft=log(p1(i,i))*(meanlevel).^p2(i,i);
mv=0:meanlevel(7);
ft=polyval(p(:,i,i),mv);
plot(mv,ft,'b');
j=j+1;
end


%% plot of diagonal cross sections of mean levels

figure;
for level=1:7 
    meanlevel(level)=sum(sum(M(level).meandata))/sum(sum(Xmask));   
    hold on;
    plot(mean(M(level).diag1,1)/(meanlevel(level)),'.','MarkerSize',1,'Color',color(level,7,'jet'));
    %plot(mean(M(level).diag2,1),'r');  %other diagonal  
    k=mean(M(level).diag1,1)/meanlevel(level); 
end

    