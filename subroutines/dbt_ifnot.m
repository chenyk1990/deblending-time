function v = yc_ifnot(yes, v1, v2)
% yc_ifnot: equivalent to C grammar (v=yes?v1:v2)
% Yangkang Chen
% July, 22, 2020

if yes
    v=v1;
else
    v=v2;
end

return