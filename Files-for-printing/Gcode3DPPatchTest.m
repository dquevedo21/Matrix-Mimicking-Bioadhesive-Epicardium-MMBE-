%Description: This script generates a txt file that can be used as gcode
%to generate the route of the nozzle to print a sinusoidal fiber patch.
%Written by Diego Quevedo-Moreno (07/08/2020)
clc, clear all, close all
yv=linspace(0,10*pi,500);
Ns=1;%Number of Samples
A=1.5;%Amplitude(mm)
P=.25;%Density 
WL=(2*pi)/10;%Wave Length(2*pi/WL=wavelength)(mm) 
W=.150;%Fiber Width(mm)
SL=50;%Sample Length(mm)
SW=5;%Sample Width(mm)
nL=10;%#Layers
Lt=.025;%Layer thickness
iL=.1;%Initial Layer Thickness
dp=.05;%Distance Between points
F= 150;%Printing Speed(mm/min)
nf=(SW/W)*P;%#Fibers
% df=(SW-(2*A)-W)/(nf-1);%Distance between fibers
df=(SW-W)/(nf-1);%Distance between fibers
Offset=[-20 0];%%Offset if needed [X Y] (mm)
x=0:dp:SL;
Vy=0:0.01:(SW-(2*A));
for j=1:nf
    if rem(j,2)==0
        for i=1:length(x)
            y(i,j)=A*sin(x(i)*WL)+(j-1)*df;
        end
    else
        for i=1:length(x)
            y(i,j)=A*sin(x(i)*WL)+(j-1)*df;
        end
    end
end
%Sinusoidal Pattern Generator
for j=0:nf-1
    if rem(j,2)==0 %Odd direction 
        for i=1:length(x)
            xt(i+(length(x)*j))=x(i);
            yt(i+(length(x)*j))=y(i,(j+1));
        end
    else %Even direction
        for i=1:length(x)
            xt(i+(length(x)*j))= x((SL/dp+1)-i+1);
            yt(i+(length(x)*j))=y(((SL/dp+1)-i+1),(j+1));
        end
    end
end

%%Skirt 
svx=[-4 -4 -3 -3 -2 -2 -1 -1];
svy=[0 SW SW 0 0 SW SW 0];
svz=[Lt Lt Lt Lt Lt Lt Lt Lt];

xtf=flip(xt);
ytf=flip(yt);
for k=1:nL
    if rem(k,2)==1
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xt(i);
            gy(i+length(xt)*(k-1))=yt(i);
            gz(i+length(xt)*(k-1))= k*Lt;
        end
    else
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xtf(i);
            gy(i+length(xt)*(k-1))=ytf(i);
            gz(i+length(xt)*(k-1))=k*Lt;
        end
    end
end

%Grip Zone
gzvx=0:W/2:10;
for k=1:2:length(gzvx)
gzx(k)=gzvx(k);
gzy(k)=SW+A;
gzz(k)=Lt*nL;
if k>1
gzx(k-1)=gzvx(k);
gzy(k-1)=-A;
gzz(k-1)=Lt*nL;
end
end

gx=[svx gx gzx gzx+SL-10];
gy=[svy gy gzy gzy];
gz=[svz gz gzz gzz];

% gx=[svx gx];
% gy=[svy gy];
% gz=[svz gz];

 
%%Offset (if needed)
gx=gx + Offset(1);
gy=gy + Offset(2);



for j=1:Ns
    for i=1:length(gx)
        gtx(i+(length(gx)*(j-1)))=gx(i);
        gty(i+(length(gx)*(j-1)))=(gy(i)+(SW+5)*(j-1));
        gtz(i+(length(gx)*(j-1)))=gz(i);
    end 
end
figure
axis equal 
grid on
hold on
plot3(gtx,gty,gtz,'.-');
xlabel('mm');
ylabel('mm');


%% File saving
% Nomenclature
% p(density)_a(amplitude)_size_F(printing speed)_dp(distance between points)_L(# of layers)_#samples_fw(fiber width)_Lt(Layer Thickness)
currentFolder = pwd;
path=strcat(currentFolder,'\SampleTest_p100_a0_7x7mm_F275_dp.5mm_L5_1_fw150_Lt.02mm.gcode');
fid = fopen(path,'w');
fprintf(fid,'F%8.6f\r\n',F);

for j = 1:length(gtx)

    fprintf(fid,'LINEAR X%8.6f Y%8.6f Z%8.6f\r\n',gtx(j),gty(j),gtz(j));

end

fclose('all');
