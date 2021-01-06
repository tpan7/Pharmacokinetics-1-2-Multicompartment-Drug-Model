function multiple2C_Met(parameters, route, num_doses)
%Which route of administration
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
[temp, indx] = max(route);

%Blood flow rate
Q = 24;
%The Volume of distribution of first compartment
Vd1 = parameters(1);
%The Volume of distribution of 2nd compartment
Vd2 = parameters(2);
%The bioavailability factor
F = parameters(4);
if indx == 1
    F = 1;
end
%The fraction unbound of drug
fu = parameters(6);
%The end of the plotting
range = parameters(7);
%The Volume of distribution of first compartment
Vd1_met = parameters(8);
%The Volume of distribution of 2nd compartment
Vd2_met = parameters(9);
%the fraction of drug is converted to metabolite
F_met = parameters(11);
%The fraction unbound of drug
fu_met = parameters(12);

dosing_times = [0];
Dose = [];
for i=1:num_doses-1
    ti = input(['At what time did you give dose #' num2str(i+1) ' (hr)? (Assuming dose #1 was at t=0): ']);
    if ti > range
        error('This dose is outside the time you are integrating over.')
    elseif ti <= dosing_times(end)
        error('This dose time is before or at the same time as the previous one.'); 
    end
    dosing_times = [dosing_times ti];
end
for i=1:num_doses
    Dosei = input(['What was dose #' num2str(i) ' (mg)?: ']);
    Dose = [Dose Dosei];
end

%Define IV experiment to solve for ke
ke_ode = @(ke, t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2))];
ke = ke_solve2C(parameters, ke_ode, [Dose(1)/Vd1 0]);

%To solve for ka
if indx ~= 1
    ka_ode = @(ka, t, x) [-ka*x(1);...
            ka*F*x(1)/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3))];
    ka = ka_solve2C(parameters, ka_ode, [Dose(1) 0 0]);
end

%Solve for ke_met
ke_met_ode = @(ke_met, t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2));...
        ke*F_met*fu*x(1)*(Vd1/Vd1_met) + (Q/Vd1_met)*(x(4)-fu_met*x(3)) - ke_met*fu_met*x(3);...
        (Q/Vd2_met)*(fu_met*x(3) - x(4))];
ke_met = ke_met_solve2C(parameters, ke_met_ode, [Dose(1)/Vd1 0 0 0], ke);

if indx == 1
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2));...
        ke*F_met*fu*x(1)*(Vd1/Vd1_met) + (Q/Vd1_met)*(x(4)-fu_met*x(3)) - ke_met*fu_met*x(3);...
        (Q/Vd2_met)*(fu_met*x(3) - x(4))];
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [Dose(1)/Vd1 0 0 0];
        else
            y = [C(end,1)+Dose(i)/Vd1, C(end,2) C(end,3) C(end,4)];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [(Q/Vd1)*(x(2)-fu*x(1)) - ke*fu*x(1);...
        (Q/Vd2)*(fu*x(1) - x(2));...
        ke*F_met*fu*x(1)*(Vd1/Vd1_met) + (Q/Vd1_met)*(x(4)-fu_met*x(3)) - ke_met*fu_met*x(3);...
        (Q/Vd2_met)*(fu_met*x(3) - x(4))];
        tspan = [dosing_times(i) dosing_times(i+1)];
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1)+Dose(num_doses)/Vd1, C(end,2) C(end,3) C(end,4)];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
else 
    t = [];
    C = [[]];
    for i=1:num_doses-1
        sys = @(t, x) [-ka*x(1); ...
            ka*x(1)*F/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3));...
            ke*F_met*fu*x(2)*(Vd1/Vd1_met) + (Q/Vd1_met)*(x(5)-fu_met*x(4)) - ke_met*fu_met*x(4);...
            (Q/Vd2_met)*(fu_met*x(4) - x(5))];
        tspan = [dosing_times(i) dosing_times(i+1)];
        if i == 1
            y = [Dose(i) 0 0 0 0];
        else
            y = [(C(end,1)+Dose(i)) C(end,2) C(end,3) C(end,4) C(end,5)];
        end
        [temp_t, temp_C] = ode45(sys, tspan, y);
        t = [t; temp_t];
        C = [C; temp_C];
    end
    sys = @(t, x) [-ka*x(1); ...
            ka*x(1)*F/Vd1 + (Q/Vd1)*(x(3)-fu*x(2)) - ke*fu*x(2);...
            (Q/Vd2)*(fu*x(2) - x(3));...
            ke*F_met*fu*x(2)*(Vd1/Vd1_met) + (Q/Vd1_met)*(x(5)-fu_met*x(4)) - ke_met*fu_met*x(4);...
            (Q/Vd2_met)*(fu_met*x(4) - x(5))];
        tspan = [dosing_times(i) dosing_times(i+1)];
    tspan = [dosing_times(num_doses) range];
    y = [C(end,1)+Dose(num_doses) C(end,2) C(end,3) C(end,4) C(end,5)];
    [temp_t, temp_C] = ode45(sys, tspan, y);
    t = [t; temp_t];
    C = [C; temp_C];
end

figure(1)
plot(t, C(:,end-3), 'b');
hold on; grid on;
plot(t, C(:,end-2), 'r');
plot(t, C(:,end-1), 'b--');
plot(t, C(:,end), 'r--');
xlabel('Time (hr)')
ylabel('Conc of drug(mg/L)');
switch indx
    case 1
        title('Drug Conc vs Time from IV Administration')
    case 2
        title('Drug Conc vs Time from Oral Administration')
    case 3
        title('Drug Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Conc vs Time from Buccal Administration')
    case 6
        title('Drug Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Conc vs Time from Lung Administration')
    case 8
        title('Drug Conc vs Time from Rectal Administration')
end
legend('Drug 1st Compartment', 'Drug 2nd Compartment', 'Metabolite 1st Compartment', 'Metabolite 2nd Compartment', 'Location', 'best')
xlim([0 range])
hold off

figure(2)
plot(t, log(C(:,end-3)), 'b')
hold on; grid on;
plot(t, log(C(:,end-2)), 'r')
plot(t, log(C(:,end-1)), 'b--')
plot(t, log(C(:,end)), 'r--')
xlabel('Time (hr)')
ylabel('Log conc (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Administration')
    case 2
        title('Drug Log Conc vs Time from Oral Administration')
    case 3
        title('Drug Log Conc vs Time from Intramusclar Administration')
    case 4
        title('Drug Log Conc vs Time from Subcutaneous Administration')
    case 5
        title('Drug Log Conc vs Time from Buccal Administration')
    case 6
        title('Drug Log Conc vs Time from Sublingual Administration')
    case 7
        title('Drug Log Conc vs Time from Lung Administration')
    case 8
        title('Drug Log Conc vs Time from Rectal Administration')
end
legend('Drug 1st Compartment', 'Drug 2nd Compartment', 'Metabolite 1st Compartment', 'Metabolite 2nd Compartment', 'Location', 'best')
hold off
end