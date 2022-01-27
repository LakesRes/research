%% Bending Round Cosserat Analysis of Experiment

clc,clear, close all
%% Finding the best fit parameters

%This data is for the 111 structures
diam = [80.53 44.3 22.21]';         %Enter diameters from largest to smallest in mm
E_measured = [7.19 9.91 24.6]';     %Enter moduli from BVS bending in MPa from big specimen to small

%Enter in the following values
v = 0.3;    %Poissons ratio

count = 1;  %Setting counting index for final solutions matrix.
for Ec = 6.0:0.1:6.0                %Searching a range of asymptotic bending moduli (MPa) in 0.1 MPa increments.
                                    %This command reads, search from 6.0 to\6.0 in 0.1 increments.  Make start and
                                    %stop values the same to only search one value.
    Output = E_measured./Ec;        %Calculating Omega values for each Ec
    N_index = 0;
    for N = 0.99950:0.00001:0.99958 %Range of possible N values searched at specified increment.  Here, N is searched
                                    %from 0.99950 to 0.99958 in increments of 0.00001
        N_index = N_index+1;
        
        %Fitting function for Cosserat elastic solids with circular cross
        %sections in bending.  Fit parameters can be limited by adjusting
        %options.  First value in square brackets is for del, second value
        %is for p.  These are defined in parameter section just below.
        %Origin of equation is G.V. Krishna Reddy's "On the flexural
        %rigidity of a micropolar elastic circular cylinder"
        
        %p = beta/gamma
        %del = 1/characteristic length of bending
        fitFunc = @(del,p,x) (1+(8.*N^2./(v+1)).*((1-p.^2)./(del.*x./2).^2+((p+v).^2./((8.*N.^2.*(1-v)) ...
            +((del.*x./2).^2.*(del.*(x./2).*besseli(0,(del.*x./2))-besseli(1,(del.*x./2)))./ ...
            (del.*x./2.*besseli(0,(del.*x./2))-2.*besseli(1,(del.*x./2))))))));
        options = fitoptions('Method','NonLinearLeastSquares',...
            'StartPoint',[0.001,0.1],'Lower',[0.001,-0.5],'Upper',[10,0.5],'TolX',1e-10);


        %Here I am fitting the data and trying to get both parameters
        %parameter1 is assigned as beta/gamma
        %parameter2 is assigned as 1/characteristic length
        [fitobject,gof,output] = fit(diam,Output,fitFunc,options);
        parameter1 = fitobject.p;
        parameter2 = fitobject.del;

        %% Creating the plots
        
        %Fitting function for creating plots
        fitFunc1 = @(x) (1+(8.*N^2./(v+1)).*((1-parameter1.^2)./(parameter2.*x./2).^2+((parameter1+v).^2./((8.*N.^2.*(1-v)) ...
            +((parameter2.*x./2).^2.*(parameter2.*(x./2).*besseli(0,(parameter2.*x./2))-besseli(1,(parameter2.*x./2)))./ ...
            (parameter2.*x./2.*besseli(0,(parameter2.*x./2))-2.*besseli(1,(parameter2.*x./2))))))));
        xx = (1:0.1:100)';          %Plot x values
        yy = fitFunc1(xx);          %Plot y values for each trial

        
        %% Calculating R^2 in this section
        
        num = numel(Output);
        sum_obsv = sum(Output);
        yavg = (1/num)*sum_obsv;

        %Sum of squares of residuals
        Predicted_E = fitFunc1(diam);
        diff = @(d) (Output(d)-Predicted_E(d)).^2;
        summation = diff(1:num);
        SSres = sum(summation);

        %Total sum of squares
        diff1 = @(d) (Output(d)-yavg).^2;
        summation1 = diff1(1:num);
        SStot = sum(summation1);

        %Calculating R^2
        R2(N_index,1) = 1-SSres/SStot;
        R2(N_index,2) = N;
        R2(N_index,3) = parameter1;
        R2(N_index,4) = parameter2;
    end
    
    %Recording best-fit best values from N trials at each Ec value. From R2
    %matrix
    [maxi idx] = max(R2(:,1));
    Final(count,1) = maxi;
    Final(count,2) = R2(idx,2);
    Final(count,3) = R2(idx,3);
    Final(count,4) = R2(idx,4);
    count = count + 1;
end

%% Fitting function for best fit parameters from all trials
fitFunc2 = @(x) (1+(8.*(Final(1,2)).^2./(v+1)).*((1-(Final(1,3)).^2)./((Final(1,4)).*x./2).^2+(((Final(1,3))+v).^2./((8.*(Final(1,2)).^2.*(1-v)) ...
            +(((Final(1,4)).*x./2).^2.*((Final(1,4)).*(x./2).*besseli(0,((Final(1,4)).*x./2))-besseli(1,((Final(1,4)).*x./2)))./ ...
            ((Final(1,4)).*x./2.*besseli(0,((Final(1,4)).*x./2))-2.*besseli(1,((Final(1,4)).*x./2))))))));
yy2 = fitFunc2(xx);
        
%% Plotting final graph
hold on
plot(diam,Output,'ko');
plot(xx,yy2,'k-');
plot(xx,ones(size(xx,1)),'r--');
xlim([0 100]);
ylim([0 7]);
title(['Best fit when E = ',num2str(Ec),' MPa, l_b = ',num2str(Final(1,4)^(-1)),' mm, N = ',num2str(Final(1,2)),' R^2 = ',num2str(Final(1,1)),', \beta/\gamma = ',num2str(Final(1,3)),'.'])
xlabel('Diameter (mm)');
ylabel('Relative Stiffness \Omega')
hold off