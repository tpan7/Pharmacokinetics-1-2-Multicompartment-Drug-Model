function [parameters] = get2C_para(metabolite)
    prompt1 = 'What is the volume of distribution of the drug in the central compartment (L/kg)?: ';
    parameters(1) = input(prompt1)*70;
    
    prompt2 = 'What is the volume of distribution of the drug in the 2nd compartment (L/kg)?: ';
    parameters(2) = input(prompt2)*70;
    
    prompt3 = 'What is the half-life of the drug from IV experiment (hr)?: ';
    parameters(3) = input(prompt3);
    
    prompt4 = 'What is the bioavailability of the drug?: ';
    parameters(4) = input(prompt4);
    if parameters(4) > 1 || parameters(4) < 0
        error('Invalid bioavailability');
    end
    
    prompt5 = 'What is the tmax of the drug (hr)?: ';
    parameters(5) = input(prompt5);
    
    prompt6 = 'What is the fraction unbound of the drug in plasma?: ';
    parameters(6) = input(prompt6);
    if parameters(6) > 1 || parameters(6) < 0
            error('Invalid fraction unbound.');
    end
    
    prompt7 = 'How long do you want to model concentration (hr)?: ';
    parameters(7) = input(prompt7);
    
    if metabolite
        prompt8 = 'What is the volume of distribution of the metabolite in the central compartment (L/kg)?: ';
        parameters(8) = input(prompt8)*70;

        prompt9 = 'What is the volume of distribution of the metabolite in the 2nd compartment (L/kg)?: ';
        parameters(9) = input(prompt9)*70;
        
        prompt10 = 'What is the half-life of the metabolite from IV experiment (hr)?: ';
        parameters(10) = input(prompt10);

        prompt11 = 'What is the fraction of drug is converted to metabolite?: ';
        parameters(11) = input(prompt11);
        if parameters(11) > 1 || parameters(11) < 0
            error('Invalid fraction converted');
        end
        
        prompt12 = 'What is the fraction unbound of the metabolite in plasma?: ';
        parameters(12) = input(prompt12);
        if parameters(12) > 1 || parameters(12) < 0
            error('Invalid fraction unbound.');
        end
    end
end