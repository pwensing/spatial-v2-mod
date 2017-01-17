function B = factorFunctions(I, v, number)
    if nargin == 2
        number = 3;
    end
    if number == 1
       B = 1/2 * crf(v)*I; 
    elseif number == 2
       B =  crfinv( I * v) - I * crm(v); 
    else
       B = 1/2 * (crf(v)*I + crfinv( I * v) - I * crm(v)); 
    end
end