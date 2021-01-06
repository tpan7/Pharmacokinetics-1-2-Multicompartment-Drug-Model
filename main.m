%% Instructions
%Final Code: You need to input the type of model in the first section and
%fill everything present in the section. Afterward, just hit run. 
%Some functions will ask for more parameters which you can type into the 
%command window.

%% Input variables here

%How many compartments do you want?
%Use '1C' for 1 compartment, '2C' for Two, 'MC' for Multicompartment
numComp = '1C';

%Which route of administration do you want? Use the number coding below
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
route = 1;

%Do you have an active metabolite? 
%0 for No, 1 for Yes 
metabolite = 0;

%What kind of drug formualtion do you have?
%Use 1 for Instant (Injection, inhalation, normal oral pill)
%Use 2 for extended release (Oral extended Tablets, patch, infusion). For
%extended release we assume 0th order release
%Use 3 for timed release (oral timed tablet)
%Use 4 for multiple doses
formulation = 3;

%% Producing plot
%Nothing needs to be changed here

%Checking for valid route
route_vec = zeros(1, 9);
if route <= 0 || route > 9
    error('Invalid route of administration.')
else
    route_vec(route) = 1;
end

%Checking valid metabolite input
if metabolite ~= 0 && metabolite ~= 1
    error('Invalid metabolite option.')   
end

%Checking valid formulation input
if formulation ~= 1 && formulation ~= 2 && formulation ~= 3 && formulation ~= 4
    error('Invalid drug formulation option.')    
end

%Checking that the drug formulation exists
if formulation == 1 && route == 9
    error('There is no instant drug formulation for the route you chose.')   
elseif formulation == 2 && route ~= 1 && route ~= 2 && route ~= 9
    error('There is no extended formulation for the route you chose.')   
elseif formulation == 3 && route ~= 2
    error('There is no timed release formulation for the route you chose.')    
elseif formulation == 4 && route == 9
    error('There is no multiple dose formulation for the route you chose.')
end

%Getting delayed start
if formulation == 3
    prompt_time = 'How long does it take for the drug to release (hr)?: ';
    time_start = input(prompt_time);
end

%Getting multiple doses
if formulation == 4
    prompt_doses = 'How many doses are you giving?: ';
    num_doses = input(prompt_doses);
    if num_doses == 1
        error('There is only one dose. Use the instant model instead.')
    end
end

%Going to the right compartment model
if strcmp(numComp, '1C') %1C model 
    parameters = get1C_para(metabolite);
    if metabolite %1C model with metabolite
        switch formulation
            case 1
                instant1C_Met(parameters, route_vec, 0);
            case 2
                extended1C_Met(parameters, route_vec);
            case 3
                instant1C_Met(parameters, route_vec, time_start);
            case 4
                multiple1C_Met(parameters, route_vec, num_doses);
        end
    else %1C model no metabolite
        switch formulation
            case 1
                instant1C_noMet(parameters, route_vec, 0);
            case 2
                extended1C_noMet(parameters, route_vec);
            case 3
                instant1C_noMet(parameters, route_vec, time_start);
            case 4
                multiple1C_noMet(parameters, route_vec, num_doses);
        end
    end
%2C model
elseif strcmp(numComp, '2C')  
    parameters = get2C_para(metabolite);
    if metabolite %2C model with metabolite
        switch formulation
            case 1
                instant2C_Met(parameters, route_vec, 0);
            case 2
                extended2C_Met(parameters, route_vec);
            case 3
                instant2C_Met(parameters, route_vec, time_start);
            case 4
                multiple2C_Met(parameters, route_vec, num_doses);
        end
    else %2C model no metabolite
        switch formulation
            case 1
                instant2C_noMet(parameters, route_vec, 0);
            case 2
                extended2C_noMet(parameters, route_vec);
            case 3
                instant2C_noMet(parameters, route_vec, time_start);
            case 4
                multiple2C_noMet(parameters, route_vec, num_doses);
        end
    end
%MC model    
elseif strcmp(numComp, 'MC')
    parameters = getMC_para(metabolite);
    if metabolite %MC model with metabolite
        error('Our code cannot handle an active metabolite for the multicompartment model.');
    else %MC model no metabolite
        switch formulation
            case 1
                instantMC_noMet(parameters, route_vec, 0);
            case 2
                extendedMC_noMet(parameters, route_vec);
            case 3
                instantMC_noMet(parameters, route_vec, time_start);
            case 4
                multipleMC_noMet(parameters, route_vec, num_doses);
        end
    end
else
    error('Invalid choice for comparments.');
end