function y=joinNumVec(x)
    temp=join(arrayfun(@(x) num2str(x),x,'UniformOutput',false),',');
    y=temp{1};