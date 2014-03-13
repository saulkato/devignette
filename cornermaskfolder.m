function cornermaskfolder(folder,cornerLengths)
%run through a folder of tifs and zero out the corners on all of them.
%creates a new file "foo-cornermasked.tif" for every file "foo.tif"
%
%you need to specify the size of the corners.
%
%cornerLengths can be:
%- a single number between 0 and 1:  a fraction of the whole image edge length
%- four numbers between 0 and 1: one fraction for each corner NW,NE,SE,SW
%- a single number > 1:  absolute pixel corner length for all corners
%- four numbers between 0 and 1:  absolute pixel corner length NW,NE,SE,SW
%
%
%Saul Kato

if nargin<2
    cornerLengths=2*[.1 .15 .2 .25]; % 5 percent
end

%turn cornerLengths into a 4 vector
if length(cornerLengths)==1
    cornerLengths=[cornerLengths cornerLengths cornerLengths cornerLengths];
end


if nargin<1
    folder=pwd;
end

%look for tifs in directory, both absolute and relative path
fnames=dir([folder filesep '*.tif']);
if isempty(fnames)  
    fnames=dir([pwd filesep folder filesep '*.tif']);
    if isempty(fnames)
        disp('no .tifs found.  Quitting.');
        return;
    end
end

X=tiffread2([folder filesep fnames(1).name]);

%convert fractional cornerLengths into absolute pixel cornerlengths
for i=1:4   
    if cornerLengths(i)<1
        cornerLengths(i)=floor(cornerLengths(i)*X.width);
    end
end


cornerLengths

warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir([folder filesep 'cornermasked']);
warning('on','MATLAB:MKDIR:DirectoryExists');

for i=1:length(fnames)
    tifname=fnames(i).name;   
    X=tiffread2([folder filesep tifname]);
    imout=cornermask(X,cornerLengths);
    outputfilename=[folder filesep 'cornermasked/' tifname(1:end-4) '-cornermasked.tif'];
    imwrite(imout,outputfilename,'tif','Compression','none')
end

end  %main function

function imgout=cornermask(X,cornerLengths)

    imgout=X.data;
    
    %corner 1
    for j=1:cornerLengths(1)   
            imgout(j,1:cornerLengths(1)-j+1)=0;
    end

    %corner 2
    for j=1:cornerLengths(2)    
            imgout(j,X.width-cornerLengths(2)+j-1:X.width)=0;
    end
    
    %corner 3
    for j=X.height-cornerLengths(3) : X.height
            imgout(j,X.width-  (j-X.height+cornerLengths(3) ) :X.width)=0;
    end
    
    %corner 4
    for j=X.height-cornerLengths(4) : X.height
            imgout(j,1:  1+j+(cornerLengths(4)-X.height)     )=0;
    end

end