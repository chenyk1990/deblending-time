function [dout]=dbt_mutter(din,x0,t0,t1)
%
%
% INPUT
% din
% x0: (mutter point 1: (x0,t0);
% t0: defaut: 0
% t1: (mutter point 2: (1,t1);
% [x0,t0,t1] are all indices
% 
% threshold curve: 
% t = (x0-x)/(x0-1)*(t1-t0)+t0 (left)
% t = (x-x0)/(x0-1)*(t1-t0)+t0 (right)
% 
% OUTPUT
% dout
%
% DEMO
% test_yc_mutter.m

[n1,n2]=size(din);
dout=din;
for i2=1:n2
    for i1=1:n1
        
        if i2<=x0
            
            if i1< ((x0-i2)/(x0-1)*(t1-t0)+t0)
                dout(i1,i2)=0;
            end
        else
            if i1< ((i2-x0)/(x0-1)*(t1-t0)+t0)
                dout(i1,i2)=0;
            end

        end
        
        
    end
end





return