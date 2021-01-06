function [parameters] = getMC_para(metabolite)
    prompt1 = 'What is the half-life of the drug from IV experiment (hr)?: ';
    parameters(1) = input(prompt1);
    
    prompt2 = 'What is the bioavailability of the drug?: ';
    parameters(2) = input(prompt2);
    if parameters(2) > 1 || parameters(2) < 0
        error('Invalid bioavailability');
    end
    
    prompt3 = 'What is the tmax of the drug (hr)?: ';
    parameters(3) = input(prompt3);
    
    prompt4 = 'What is the fraction unbound of the drug in plasma?: ';
    parameters(4) = input(prompt4);
    if parameters(4) > 1 || parameters(4) < 0
            error('Invalid fraction unbound.');
    end
    
    prompt5 = 'What is the water octanol partition coefficient for the drug?: ';
    parameters(5) = input(prompt5);
    
    prompt6 = 'What fraction of drug is eliminated in urine unchanged?: ';
    parameters(6) = input(prompt6);
    if parameters(6) > 1 || parameters(6) < 0
            error('Invalid fraction eliminated.');
    end
    
    prompt7 = 'How long do you want to model concentration (hr)?: ';
    parameters(7) = input(prompt7);
    
    if metabolite
        prompt8 = 'What is the half-life of the metabolite from IV experiment (hr)?: ';
        parameters(8) = input(prompt8);
        
        prompt9 = 'What is the water octanol partition coefficient for the drug?: ';
        parameters(9) = input(prompt9);
        
        prompt10 = 'What is the fraction unbound of the metabolite in plasma?: ';
        parameters(10) = input(prompt10);
        if parameters(10) > 1 || parameters(10) < 0
            error('Invalid fraction unbound.');
        end
    end
end