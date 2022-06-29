function [dout ] = dbt_deriv(din, order, axis, scale)
% First derivative with a maximally linear FIR differentiator
% REF: http://ahay.org/blog/2012/05/01/program-of-the-month-sfderiv/
% By Yangkang Chen
% Aug 26, 2019
%  
if nargin==1
   order=6;
   axis=1;
   scale=1;
end

if nargin==2
   axis=1;
   scale=1;
end

if nargin==3
    scale=1;
end
[n1,n2,n3]=size(din);
dout=zeros(n1,n2,n3);
c1=0;

switch axis
    case 1
        for i3=1:n3
            for i2=1:n2
                tmp=din(:,i2,i3);
                tmp2 = der(tmp,n1,order,c1)*scale;
                dout(:,i2,i3)=tmp2(:);
            end
        end
    case 2
         for i3=1:n3
            for i1=1:n1
                tmp=din(i1,:,i3);
                tmp2 = der(tmp,n2,order,c1)*scale;
                dout(i1,:,i3)=tmp2(:);
            end
        end       
        
    case 3
          for i2=1:n2
            for i1=1:n1
                tmp=din(i1,i2,:);
                tmp2 = der(tmp,n3,order,c1)*scale;
                dout(i1,i2,:)=tmp2(:);
            end
        end        
    otherwise 
    error("Invalid parameter");
end


return



function trace2 = der(trace,nt,n,c1)
% derivative operator
% nt: transform length
% n: order
% c1: filter parameter

c=1./(2*sqrt(c1));
c2=c*c;c2=c2*1.0;

h=trace;%length of trace is nt
trace2=zeros(size(trace));

for ii=n:-1:1
    for it=1:nt-2
    trace2(it+1)=h(it+1) - 0.5*(h(it+2)+h(it));
    end
    
    trace2(1)=trace2(2);
    trace2(nt)=trace2(nt-1);
    
    for it=0:nt-1
        h(it+1)=trace(it+1) + trace2(it+1)*ii/(2*ii+1);
    end
    
end

trace2(1)=h(2)-h(1);
for it=1:nt-2
    trace2(it+1) = 0.5*(h(it+2)-h(it));
end


trace2(nt)=h(nt)-h(nt-1);

return