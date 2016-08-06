function [Tely,eta_tot,V_H2,V_O2,Ustack] =  alkalineElectrolyzer(ely,Iely,Tini,on_off,Tambient,Tcool_in,Vcool_in)
% This is an implementation of advanced alkaline electrolyzer
% and it calculate the Ely operating status at each time interval
% reference: Modeling of advanced alkaline electrolyzers: a system
% simulation approach

Timestep = 1; % This is a hourly calculation
% general information of Ely
Area =ely.Area;
Ncells = ely.Ncells; % number of cells in a stack
Nstack = ely.Nstack; % number of stacks in a electrolyzer unit.

% I-V curve cof

IVcurve_par(1) = ely.r1;
IVcurve_par(2) = ely.r2;
IVcurve_par(3)= ely.s;
IVcurve_par(4)= ely.t1;
IVcurve_par(5) = ely.t2;
IVcurve_par(6) =ely.t3;

% Faraday efficiency cof

Farad_par(1) = ely.a1;
Farad_par(2) =ely.a2;

% thermal model ocf
therm_par(1) =ely.h1;
therm_par(2) = ely.h2;
therm_par(3) = ely.Rt;

%

Pely =1;  % pressure, unit:bar
Ta  = Tambient; %ambient temperature(Ta is assumed to the the ambient temperature here.)
Tcool_in = Tcool_in; % cooling water input temperature
Vcool_in = Vcool_in; % cooling water flow rate
Tely_avg = Tini; % This is the average electrolyzer temperatre over this time interval
Tely = Tini; % This is the first value of Tely



% operating limit
Tmax   = 'to_define';   % maximum allowable operating temperature
Idmax   = 'to_define';  % maximum allowable current density
Ucmin = 'to_define';      % minimum allowable operating cell voltage


% output initialize
eta_e =[]; % energy efficiency
eta_f = []; % Faraday efficieny
v_H2  = []; %Hydrogen production rate
v_O2 =[]; % Oxygen production rate

%% This is the main part of this function
% if Electrolyzer is switch off
if on_off == 0 || Iely<1e-6
    [Tely,V_H2,V_O2,eta_e,eta_f,Qloss] = ElyOff(therm_par,Ta,Ncells,Tini,Timestep);
    eta_tot = eta_e*eta_f;
else
    % if Electrolyzer is witch off
    beta = 0.6; % This is a factor to estimate the average Ely temperature over a time interval, cause curve(T=f(t))  is convex.
    diff = 10; % a criterion to check if the iteration loop converges
    ii = 0;    % a counter to count the number of  iteration
    
    
    while (ii<100) &&(diff>0.001)
        % calculate Gibbs free energy
        [Utn,Urev] = ElyGibbs(Tely_avg, Pely);
        
        %
        [Ucell,Idensity] = ElyElec(IVcurve_par,Area,Iely,Tely_avg,Urev);
        
        %
        eta_f = ElyFarad(Farad_par,Idensity);
        
        %
        [Ustack,Ptot,V_H2,V_O2,eta_tot,eta_e] = ElyStack(Ncells,Iely,eta_f,Ucell,Utn);
        
        %Thermal model
        [Telyfin,Tcool_out,Qstore,Qgen,Qloss,Qcool] =  ElyTherm(therm_par,Iely,Ta,Tcool_in,Vcool_in,Ucell,Utn,Ncells,Tini,Tely_avg,Timestep);
        diff = abs(Tely-Telyfin);  %Tely is the value of last iteration
        Tely = Telyfin;
        Tely_avg = beta*Tely+(1-beta)*Tini;
        ii = ii+1;
        
    end
end

    function [Utn, Urev] = ElyGibbs(Tely,Pely)
        %  Standard conditions:25 centigrade, 1 bar
        Tref = 25; % centigrade
        Pref = 1;
        FaradConst = 96485; %As/mol
        Rgas = 8.3145; %J/(k，mol��
        Nelec = 2;
        
        Href_H2O = -286E3; %L/mol
        Href_H2 = 0;
        Href_O2 = 0;
        Sref_H2O = 70;  %J/(k，mol)
        Sref_H2 = 131;
        Sref_O2 = 205;
        Cp_H2O = 75; %J/(k，mol)
        Cp_H2 = 29;
        Cp_O2 = 29;
        
        % temperature
        [T_H2,T_O2,T_H2O] = deal(Tely);
        % Enthalpy
        H_H2 = Cp_H2*(T_H2-Tref)+Href_H2;
        H_O2 = Cp_O2*(T_O2-Tref)+Href_O2;
        H_H2O = Cp_H2O*(T_H2O -Tref)+Href_H2O;
        dH = H_H2 + 0.5*H_O2 - H_H2O;
        
        %entropy :check logarithm arguments brefore trying them
        if (((T_H2+273.15)/(Tref+273.15))<=0)
            error('The electrolyzer temperature should not be lower than absolute zero');
        end
        S_H2 = Cp_H2 * log((T_H2+273.15)/(Tref+273.15)) - Rgas * log(Pely/Pref) + Sref_H2;
        
        if (((T_O2+273.15)/(Tref+273.15))<=0)
            error('The electrolyzer temperature should not be lower than absolute zero');
        end
        S_O2 = Cp_O2 * log((T_O2+273.15)/(Tref+273.15)) - Rgas * log(Pely/Pref) + Sref_O2;
        
        if (((T_H2O+273.15)/(Tref+273.15))<=0)
            error('The electrolyzer temperature should not be lower than absolute zero');
        end
        S_H2O = Cp_H2O * log((T_H2O+273.15)/(Tref+273.15)) +  Sref_H2O;
        dS = S_H2+ 0.5*S_O2 - S_H2O;
        
        % Gibbs free energy
        dG = dH - (Tely+273.15)*dS;
        
        %Thermoneutral voltage(per cell)
        Utn = dH / (Nelec*FaradConst);
        % Reversible voltage (per cell)
        Urev = dG / (Nelec*FaradConst);
    end

    function[Ucell, Idensity] =   ElyElec(Par, Area, Iely,Tely,Urev)
        r1 = Par(1);
        r2 = Par(2);
        s1 = Par(3);
        t1 = Par(4);
        t2 = Par(5);
        t3 = Par(6);
        
        if (Tely==0)
            error('The electrolyzer temperature equalled zero centigrade')
        end
        
        Uohmic = (r1+r2*Tely)*Iely/Area;
        Uoverpot = s1*log10((t1+t2/Tely+t3/(Tely^2))*Iely/Area+1);
        Ucell = Urev+Uohmic + Uoverpot;
        % current density
        Idensity = Iely/Area/10; % 1 A/squareMeter = 10  mA/squareCentiMeter
    end

    function eta_f =  ElyFarad(par,Idensity)
        % This function calculateds the Farady efficiency
        a1  = par(1);
        a2 = par(2);
        
        eta_f = Idensity^2/(a1+Idensity^2)*a2;
    end

    function [Uely,Ptot,V_H2,V_O2,eta_tot,eta_e] = ElyStack( Ncells,Iely,eta_f,Ucell,Utn)
        FaradConst = 96485; %As/mol
        Rgas = 8.3145; %J/(k，mol��
        Nelec = 2;
        Tstd = 0; %centigrade
        Pstd =1; %bar
        
        % energy efficiency
        eta_e = Utn/Ucell;
        % Overall efficiency
        eta_tot = eta_e*eta_f;
        
        % Hydrogen & Oxygen production
        n_H2 = eta_f*Ncells*Iely/(Nelec*FaradConst); %mol
        n_O2 = 0.5* n_H2;
        % production flowrate
        V_H2 = n_H2 *(Rgas*(Tstd+273.15))/Pstd;  %Nm3/hr, Nm3--standard cubic meters
        V_O2 = n_O2 *(Rgas*(Tstd+273.15))/Pstd;  %Nm3/hr
        
        %stack voltage
        Ustack = Ncells*Ucell;
        % stack power
        Ptot = Uely*Iely;
    end

    function [Telyfin,Tcool_out,Qstore,Qgen,Qloss,Qcool] =  ElyTherm(par,Iely,Ta,Tcool_in,Vcool_in,Ucell,Utn,Ncells,Tini,Tely_avg,Timestep)
        %
        tmode =2;
        
        m_H2O = 18.016; %g/mol
        rho_H2O = 1000; % kg/m3;
        Cpo_H2O = 75; %J/(K，mol)
        Ct = 625000;   %J/centigrade
        
        hx1 = par(1);
        hx2 = par(2);
        Rt = par(3);
        % Generated thermal Energy
        Qgen = Ncells*Iely*(Ucell-Utn);
        % heat loss to ambient
        Qloss = 1/Rt*(Tely_avg-Ta);
        % Auxiliary cooling demand...
        Cp_H2O = Vcool_in*rho_H2O*1000 /m_H2O*Cpo_H2O/3600; %J/K
        
        UA = hx1 + hx2*Iely;
        Tcool_out = Tcool_in + (Tely_avg - Tcool_in)*(1- exp(-UA/Cp_H2O));
        Qcool = Cp_H2O *(Tcool_out-Tcool_in);
        Tau_t = Ct*Rt/3600; %hr
        
        if tmode == 1  % If the timestep is small enough
            Telyfin =Tini+(Timestep*3600) /Ct*(Qgen-Qcool-Qloss);
        elseif tmode ==2 % A complex dumped thermal compacitance model
            a = 1/Ct*(1/Rt+Cp_H2O*(1-exp(-UA/Cp_H2O)));
            b = 1/Ct*(Qgen+Ta/Rt+Cp_H2O*(1-exp(-UA/Cp_H2O))*Tcool_in);
            Telyfin = (Tini-b/a)*exp(-a*Timestep*3600)+b/a;
        end
        
        Qstore = Qgen-Qcool-Qloss;
    end

    function [Tely,V_H2,V_O2,eta_e,eta_f,Qloss] = ElyOff(par,Ta,Ncells,Tini,Timestep)
        
        tmode = 2;
        Rt = par(3);
        Ct = 625000;   %J/centigrade
        Tau_t = Ct*Rt/3600; %hr, time constant
        
        
        Qgen =0;
        Qcool = 0;
        
        beta =0.6;
        diff = 10;
        ii = 0;
        [Tely_avg,Tely] = deal(Tini);
        while (diff>0.001) &&(ii<100)
            if tmode==1
                Qloss = 1/Rt*(Tely_avg-Ta);
                Telyfin = Tini+(Timestep*3600)/Ct*(Qgen-Qcool-Qloss); % 1Hr = 3600s
                diff = Tely - Telyfin;
                Tely = Telyfin;
                Tely_avg = beta*Tely+(1-beta)*Tini;
                ii = ii+1;
            else % A complex dumped capacitance thermal model
                a = 1/Ct*(1/Rt);
                b = 1/Ct*(Qgen +Ta/Rt);
                Telyfin = (Tini - b/a)*exp(-a*Timestep*3600)+b/a;
                diff = Tely - Telyfin;
                Tely = Telyfin;
                Tely_avg = beta*Tely+(1-beta)*Tini;
                Qloss = 1/Rt*(Tely_avg-Ta);
            end
            
            Qstore = Qgen -Qcool-Qloss;
        end
        
        %Output
        V_H2 = 0;
        V_O2 = 0;
        eta_e = 0;
        eta_f = 0;
        
    end
end



