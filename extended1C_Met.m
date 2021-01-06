function extended1C_Met(parameters, route)
%Which route of administration
%1 for Intravenous, 2 for Oral, 3 for Intramusclar, 4 for Subcutaneous
%5 for Buccal, 6 for Sublingual, 7 for Lung, 8 for Rectal, 9 for
%transdermal
[temp, indx] = max(route);

%The Volume of distribution
Vd = parameters(1);
%The half life
t_half = parameters(2);
%The bioavailability factor
F = parameters(3);
if indx == 1
    F = 1;
end
%The time where max concentration occurs
t_max = parameters(4);
%The end of the plotting
range = parameters(5);
%the volume of distribution of the metabolite (L/kg)
Vd_met = parameters(6);
%the half-life of the metabolite (hr)
t_half_met = parameters(7);
%the fraction of drug is converted to metabolite
F_met = parameters(8);

%rate of elimination
ke = log(2)/t_half;
%rate of elimination of metabolite
ke_met = log(2)/t_half_met;

prompt_dose = 'What is the total dose of the drug (mg)?: ';
Dose = input(prompt_dose);

prompt_time = 'Over what length of time is the drug being administered (hr)? In other words, what is the extinction time of administration?: ';
t_extinct = input(prompt_time);

%Ka is the same for all oral systems
low_guess_ka = 0; %define some high and low bounds for guessing
high_guess_ka = 100;
ka = (low_guess_ka + high_guess_ka)/2;
if indx ~= 1
    tol = .1;
    T_max = t_max/2;
    while T_max < t_max - tol * t_max || T_max > t_max + tol * t_max
        if ka == 0 || ka == 100
            error('Ka solver hit the bounds. Go to extended1C_noMet.m to change low and high guess.')
        end
        ka = (low_guess_ka + high_guess_ka)/2;
        ode = @(t, y) [-ka*y(1); ...
            ka*F*y(1)/Vd - ke*y(2);...
            ke*y(2)*F_met*(Vd/Vd_met) - ke_met*y(3)];
        [T, Y] = ode45(ode, [0 1000], [Dose, 0, 0]);
        [argmax, maxVal] = max(Y(:,2));
        T_max = T(maxVal);
        if T_max < t_max - tol * t_max
            high_guess_ka = ka;
        elseif T_max > t_max + tol * t_max
            low_guess_ka = ka;
        end
    end
end

%IV infusion or patch
if indx == 1 || indx == 9
    %During infusion or patch
    sys1 = @(t, x) [(Dose/Vd)*F/t_extinct - ke*x(1);...
        ke*x(1)*(Vd/Vd_met)*F_met - ke_met*x(2)];
    tspan1 = [0 t_extinct];
    y1 = [0 0];
    [t1, C1] = ode45(sys1, tspan1, y1);
    
    if range > t_extinct
        sys2 = @(t, x) [-ke*x(1);...
            ke*x(1)*(Vd/Vd_met)*F_met - ke_met*x(2)];
        tspan2 = [t_extinct range];
        y2 = [C1(end,1) C1(end,2)];
        [t2, C2] = ode45(sys2, tspan2, y2);
        
        t1 = [t1; t2];
        C1 = [C1; C2];
    end
else %Oral extended tablet
    sys1 = @(t, x) [Dose/t_extinct-ka*x(1); ...
        ka*F*x(1)/Vd - ke*x(2);...
        ke*x(2)*(Vd/Vd_met)*F_met - ke_met*x(3)];
    tspan1 = [0 t_extinct];
    y1 = [0 0 0];
    [t1, C1] = ode45(sys1, tspan1, y1);
    
    if range > t_extinct
        sys2 = @(t, x) [-ka*x(1); ...
            ka*F*x(1)/Vd - ke*x(2);...
            ke*x(2)*(Vd/Vd_met)*F_met - ke_met*x(3)];
        tspan2 = [t_extinct range];
        y2 = [C1(end, 1) C1(end, 2) C1(end, 3)];
        [t2, C2] = ode45(sys2, tspan2, y2);
        
        t1 = [t1; t2];
        C1 = [C1; C2];
    end
end

figure(1)
plot(t1, C1(:,end-1), 'b');
hold on; grid on;
plot(t1, C1(:,end), 'r');
xlabel('Time (hr)')
ylabel('Conc in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Conc vs Time from IV Infusion')
    case 2
        title('Drug Conc vs Time from Oral Extended Tablet')
    case 9
        title('Drug Conc vs Time from Transdermal Patch')
end
xlim([0 range])
legend('Main Drug', 'Metabolite', 'Location', 'best')
hold off

figure(2)
plot(t1, log(C1(:,end-1)), 'b');
hold on; grid on;
plot(t1, log(C1(:,end)), 'r');
xlabel('Time (hr)')
ylabel('Log conc of drug in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Infusion')
    case 2
        title('Drug Log Conc vs Time from Oral Extended Tablet')
    case 9
        title('Drug Log Conc vs Time from Transdermal Patch')
end
xlim([0 range])
hold off
ylabel('Log conc of drug in central compartment (mg/L)');
switch indx
    case 1
        title('Drug Log Conc vs Time from IV Infusion')
    case 2
        title('Drug Log Conc vs Time from Oral Extended Tablet')
    case 9
        title('Drug Log Conc vs Time from Transdermal Patch')
end
xlim([0 range])
legend('Main Drug', 'Metabolite', 'Location', 'best')
hold off

end