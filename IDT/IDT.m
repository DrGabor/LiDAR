% I = imread('C:\Users\Ë∫\Desktop\±œ“µ…Ëº∆\Automatic online calibration\right1.jpg');
function  D = IDT(I,r,a)
if nargin == 1
    r = 0.98;
    a = 1/3;
end
if nargin == 0
    I = imread('G:\20171023morning\Calib_Source\edge_5.jpg');
    r = 0.98;
    a = 1/3;
end
[m,n] = size(I);
E = double(I);
I = double(I);
for i = 2:m
    for j = 2:n-1
        E(i,j)  = max(max(max(E(i-1:i,j-1:j+1)))*r,E(i,j));
    end
end

for i = m-1:-1:1
    for j = n-1:-1:2
        E(i,j) = max(max(max(E(i:i+1,j-1:j+1)))*r,E(i,j));
    end
end
D = (1-a)*E+a*I;
D = uint8(D);