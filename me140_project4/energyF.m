% energyF.m
% 4/15/16 - Created Jon Renslo
% 4/22/16 - Adapted from hMix.m
% Returns enthalpy and Gibbs function of a specified number of moles of a selected species at 
% a selected temperature.
% Also returns the mole-specific heat at that temperature

function out = energyF(T,P,species,moles)
    % H in Joules, or Joules per mol if no moles number specified
    % P in Pa
    % cpbar in J/mol-K
    % cp in J/kg-K
    P0 = 101.3e3;
    R = 8.314; %J/mol-K
    
    
    if(~exist('moles','var')) moles = 1; end
    T0 = 273 + 25; %standard conditions 25C
    [co2, h2ovap, h2o, n2, o2, air, airConst,h2] = deal(1,2,3,4,5,6,7,8);
    
    fit{co2} =    [22.26, 5.981*10^-2,   -3.501 *10^-5, 7.469*10^-9  ];
    fit{h2ovap} = [32.24, 0.1923*10^-2,  1.055*10^-5,   -3.595 *10^-9];
    fit{h2o} =    [75.42271 0 0 0]; %calculated from 4.1855 at 15C
    fit{n2} =     [28.90, -0.1517*10^-2, 0.8081*10^-5,  -2.873*10^-9 ];
    fit{o2} =     [25.48, 1.520*10^-2,   -0.7155*10^-5, 1.312*10^-9  ];
    fit{air} =    [28.11, 0.1967*10^-2,  0.4802*10^-5,  -1.966*10^-9];
    fit{airConst} = [27.8715 0 0 0];
    fit{h2} =     [29.11, -.1916e-2,     0.4003e-5     -0.8704e-9];

    
    hf{co2} = -393520; %J/mol
    hf{h2ovap} = -241820; 
    hf{h2o} = -285830; % for liquid water
    hf{n2} = 0;
    hf{o2} = 0;
    hf{air} = 0;
    hf{airConst} = 0;
    hf{h2} = 0;
    
    sf{co2} = 213.8; %J/mol
    sf{h2ovap} = 188.83; 
    sf{h2o} = 69.92; % for liquid water
    sf{n2} = 191.61;
    sf{o2} = 205.04;
    sf{air} = sf{n2}*3.76/4.76 + sf{o2}/4.76; %** cannot do
    sf{airConst} = sf{n2}*3.76/4.76 + sf{o2}/4.76; %** not in table
    sf{h2} = 130.68;
    
    gf{co2} = -394360; %J/mol
    gf{h2ovap} = -228590; 
    gf{h2o} = -237180; % for liquid water
    gf{n2} = 0;
    gf{o2} = 0;
    gf{air} = 0;
    gf{airConst} = 0;
    gf{h2} = 0;
    
    m{co2} = 44; %g/mol
    m{h2ovap} = 18.02;
    m{h2o} = 18.02;
    m{n2} = 28;
    m{o2} = 32;
    m{air} = 28.98;
    m{airConst} = 28.98;
    m{h2} = 2.02;
    
    switch lower(species)
        case 'co2'     
            i = co2;
        case 'h2o'
            i = h2o;
        case'h2ovap'
            i = h2ovap;
        case 'n2'
            i = n2;
        case 'o2'
            i = o2;
        case 'air'
            i = air;
        case 'airconst'
            i = airConst;
        case 'h2'
            i = h2;
        otherwise
            disp 'input a supported species';
    end
    fits = num2cell(fit{i});
    [a, b, c, d] = deal(fits{:});    
%     I1 = R*(a*(T0 - T) + b/2*(T0^2 - T^2) + c/3*(T0^3 - T^3) + d/4*(T0^4 - T^4) + e/5*(T0^5 - T^5));
%     cp_ave = I1/(T0 - T);
%     I2 = a*log(T0/T) + b*(T0 - T) + c/2*(T0^2 - T^2) + d/3*(T0^3 - T^3) + e/4*(T0^4 - T^4);

    gPerKg = 1000;
    out.cpbar = (a + b*T + c*T.^2 + d*T.^3); % J/mol-k
    out.cp = out.cpbar/m{i}*gPerKg; %J/g-k
    
    % integrate cp from t0 to t in j/mol kelvin
    delH = (a*(T - T0) + b/2*(T.^2 - T0.^2) + c/3*(T.^3 - T0.^3) + d/4*(T.^4 - T0.^4)); 
    intCpbarOnT = a*log(T./T0) + b*(T - T0) + c/2*(T.^2 - T0.^2) + d/3*(T.^3 - T0.^3);
    delS = intCpbarOnT - R *log(P/P0); 
    out.S = (sf{i} + delS)*moles;
    out.H = (hf{i} + delH)*moles;
    
    out.G = (delH - T.*delS+ gf{i})*moles; 
    %**double check this calculation
    % delH*moles -T.*out.S + gf{i}*moles seems to produce best results, why?
end
