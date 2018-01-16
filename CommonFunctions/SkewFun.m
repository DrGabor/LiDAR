function [A] = SkewFun(a)
if length(a) == 3
    A = [0 -a(3) a(2)
        a(3)  0 -a(1)
        -a(2) a(1) 0];
end
if length(a) == 2
    A = [a(2); -a(1)];
end
end

