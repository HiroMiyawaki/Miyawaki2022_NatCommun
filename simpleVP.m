function simpleVP(x,vp,edgeColor,faceColor,lineColor,width,type,rawRad,jitter)
    
    if ~exist('width','var') || isempty(width)
        width = 0.8;
    end
    
    if ~exist('type','var')
        type='';
    end
    
    if ~exist('rawRad','var')
        rawRad=[];
    end
    if ~exist('jitter','var')||isempty(jitter)
        jitter=0.2;
    end


    fill(x+vp.x*width,vp.y,faceColor,'EdgeColor',edgeColor)
%     plot(x+[vp.x,vp.x(1)]*width,[vp.y,vp.y(1)],'Color',faceColor)

    if ~isempty(rawRad)
        for ii =1:length(vp.raw)
            jit=jitter*(rand()-0.5);
            rectangle('Position',[x+jit-rawRad(1),vp.raw(ii)-rawRad(2),rawRad(1)*2,rawRad(2)*2],...
                'Curvature',[1,1],'LineStyle','-','FaceColor','none','EdgeColor','k','LineWidth',0.25)
        end    
    end
    
    if strcmpi(type,'line') || strcmpi(type,'l')
        plot(x+vp.medDen*[-1,1]*width,vp.median+[0,0],'-','Color',lineColor,'LineWidth',0.5)        
    elseif strcmpi(type,'dot') || strcmpi(type,'d')
        plot(x+[0,0],vp.quartile,'-','Color',lineColor,'LineWidth',0.25)
        plot(x,vp.median,'.','color',lineColor,'MarkerSize',3)
    end

    