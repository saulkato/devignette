function maskImage=makeHexPostMask(width,height,centerx,centery,radius,spacing,rotation)

    %currently only support vertical orientation
    if (nargin<7) rotation=1; end %degrees counterclockwise
    if (nargin<6) spacing=74; end
    if (nargin<5) radius=27; end
    if (nargin<4) centery=160; end
    if (nargin<3) centerx=669; end
    if (nargin<2) height=1024; end
    if (nargin<1) width=1024; end
    

    bufferImage=zeros(width+4*radius,height+4*radius);
    
    
    centerList=MakeCenterList(2*radius+centerx,2*radius+centery,1+radius,1+radius,width+3*radius,height+3*radius,spacing);
    

    [centerList.x, centerlist.y]=applyRotation((centerList.x)',(centerList.y)',rotation);

    for i=1:length(centerList.x)
        DrawCircle(round(centerList.x(i)),round(centerList.y(i)));
    end
    
    %crop bufferImage
    maskImage=bufferImage(1+2*radius:width+2*radius,1+2*radius:height+2*radius)';
    
    %plot for debugging
    %figure; imagesc(maskImage); axis image;
    
%end main
    


    
    %nested functions
  
    function list=MakeCenterList(x,y,xMin,yMin,xMax,yMax,spacing)
    
        xspacing=spacing*sqrt(3)/2;
        
        
        list.x=x;
        list.y=y;
        
        listAdd=MakeStripList(x,y,yMin,yMax,spacing);
        list.x=[list.x round(listAdd.x)];
        list.y=[list.y listAdd.y];        
        
        
        parity=1;
        xStep=x;
        while (xStep>xMin+xspacing)
            
            
            xStep=xStep-xspacing;
                   
            if parity
                yStep=y-spacing/2;
            else
                yStep=y;
            end
            
            listAdd=MakeStripList(xStep,yStep,yMin,yMax,spacing);
            list.x=[list.x round(listAdd.x)];
            list.y=[list.y listAdd.y];

            
            parity=1-parity;
        end
        
        parity=1;
        xStep=x;
        while (xStep<xMax-xspacing)
 
            xStep=xStep+xspacing;
                        
            if parity
                yStep=y-spacing/2;
            else
                yStep=y;
            end
            
            listAdd=MakeStripList(xStep,yStep,yMin,yMax,spacing);
            list.x=[list.x round(listAdd.x)];
            list.y=[list.y listAdd.y];

            parity=1-parity;
        end
        
  
        
    end


    function list=MakeStripList(x,y,yMin,yMax,spacing)
        
        list.y=y;
        
        %make above list
        yStep=y;
        while(yStep>=yMin)
            list.y=[list.y yStep];
            yStep=yStep-spacing;
        end
           
        %make below list
        yStep=y;        
        while (yStep<=yMax)
            list.y=[list.y yStep];
            yStep=yStep+spacing;
            
        end
          
        list.x=x*ones(size(list.y));
        
    end



    function DrawCircle(x,y)
        
         circImage=circularmask(radius,2*radius+1);   

         bufferImage(x-radius:x+radius,y-radius:y+radius)=circImage;
    end


    function [xcoordsrot, ycoordsrot]=applyRotation(xcoords,ycoords,degrees)
        
        rotMat=computeRotationMatrix(degrees/180*pi);
        res=(rotMat*[xcoords ycoords ]')';
        xcoordsrot=res(:,1);
        ycoordsrot=res(:,2);

    
    end

    function mat=computeRotationMatrix(theta)

        %this is a counterclockwise rotation
        mat=[cos(theta) -sin(theta); sin(theta) cos(theta) ]';

    end
    

end %main function