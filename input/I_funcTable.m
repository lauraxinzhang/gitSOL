% This script computes the integral function I(x;R) on a 2D X-R grid and
% outputs to a text file for future interpolation purposes

% ====== Define the grid ======

xMin = -5;
xMax = 10;
numX = 50;

RMin = 1;
RMax = 5;
numR = 20;

xList = linspace(xMin, xMax, numX);
RList = linspace(RMin, RMax, numR);

I = zeros(numX, numR);

% Integration settings
tMax = 100;

% ====== Compute I ======
for i=1:numX
    
    x = xList(i);
    
    for j=1:numR
        
        R = RList(j);
        
        if x == 0
            
            % Compute the no-mirror `Debye' Contribution
            I(i,j) = 0.25;
            
            % Finite mirror contributions
            if R > 1
                
                I(i,j) = I(i,j) - 0.25 * (R-1) * acsch(sqrt(R-1)) / sqrt(R);
                
            end
            
        elseif x > 0
            
            % Compute the no-mirror 'Debye' Contribution
            I(i,j) = 0.125*x*exp(0.5*x)*( besselk(1, 0.5*x)+besselk(0, 0.5*x) );
            
            % Finite mirror contributions
            if R > 1
                
                a = sqrt(x/(R-1));
                f = integral(@(t) cosh(0.5*t).*dawson(a*cosh(0.5*t)).*exp(-0.5*x*cosh(t)),...
                             0,tMax);
                I(i,j) = I(i,j) - 0.25*sqrt( x*(R-1) )*exp(0.5*x)*f;
                 
            end
            
        else
            
            xAbs = abs(x);
            
            % Compute the no-mirror 'Debye' Contribution
            I(i,j) = 0.125*xAbs*exp(-0.5*xAbs)*( besselk(1, 0.5*xAbs)-besselk(0, 0.5*xAbs) );
            
            % Finite mirror contributions
            if R > 1
                
                a = sqrt(xAbs/(R-1));
                f = integral(@(t) sinh(0.5*t).*dawson(a*sinh(0.5*t)).*exp(-0.5*xAbs*cosh(t)),...
                             0,tMax);
                I(i,j) = I(i,j) - 0.25*sqrt( xAbs*(R-1) )*exp(-0.5*xAbs)*f;
                
            end
            
        end
        

    end
    
end

% ====== Plot table ======

xGrid = repmat(xList',[1,numR]);
RGrid = repmat(RList,[numX,1]);

surf(xGrid,RGrid,I)

% % ====== Output table ======
% 
% % Path to data storage location
% dataPath = '/Users/Nick/Documents/MATLAB/Research_Data/SOL';
% 
% % data file name
% fileName = 'I_func_TableData.txt';
% 
% file = fullfile(dataPath,fileName);
% fid = fopen(file,'w');
% 
% % Flatten all arrays into 1-D vectors to store as 3 column file
% xOut = repmat(xList,[1,numR]);
% rOut = reshape( repmat(RList,[numX,1]), [1,numR*numX] );
% IOut = reshape(I, [1,numR*numX]);
% 
% % file key: xvals Rvals Ivals
% fprintf(fid,'%f %f %f \n',[xOut;rOut;IOut]);
% 
% fclose(fid);
