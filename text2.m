% function h=text2(Xpos,Ypos,String,Axis,varargin)
% wrapper of text to specify text as in relative position
%
% Hiro Miyawaki at UWM, 2011-
%
function h=text2(Xpos,Ypos,String,Axis,varargin)

    if nargin<4
        Axis=axis();
    end
    
    % just for backward compatibility
    if length(varargin)==1 && iscell(varargin)
        Property=varargin{1};
    else
        Property=varargin;
    end

    if strcmpi(get(gca,'XScale'),'log')
        convPosX=exp(log(Axis(1:2))*[1-Xpos;Xpos]);
    else
        convPosX=Axis(1:2)*[1-Xpos;Xpos];
    end
    if strcmpi(get(gca,'yScale'),'log')
        convPosY=exp(log(Axis(3:4))*[1-Ypos;Ypos]);
    else
        convPosY=Axis(3:4)*[1-Ypos;Ypos];      
    end
    H=text(convPosX,convPosY,String,Property{:});
    if nargout>0
        h=H;
    end
    
    
    
    
    
    