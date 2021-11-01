clear,clc
close all
format compact
format long


% INPUT %
% =========================================================================

Flue_CorrectionFactor = 1.5;
Steam_CorrectionFactor = 2.5;

BankType   = 2;
SolverType = 1;
NuFType    = 3;
NuSType    = 3;
EHEXType   = 1;
ZpType     = 2;

BoilerCapacity = 3140; %(kcal/hr)
SteamFlowRate = 12000.0; %(kg/hr)
ExcessAir = 17.0;

FuelHeatValue = 10.0; %(kcal/m3)

X_METHANE = 96.00;
X_PROPANE = 4.00;

T_gas = 25.0+273.15;
P_gas = 101.325;
T_air = 25.0+273.15;
P_air = 101.325;

Tube_Nr = 18;
Tube_Nc = 16;
Tube_Kw = 50.0;
Tube_Rn = 0.000001;
Tube_Ru = 100.0*0.001;
Tube_Do = 42.2*0.001;
Tube_Wt = 3.0*0.001;
Tube_Ls = 2100.0*0.001;
Tube_St = 84.4*0.001;
Tube_Sl = 100.0*0.001;
Tube_cl = 25.0*0.001;
Tube_Rfi = 0.0000;
Tube_Rfo = 0.0000;

Duct_corbel = 100*0.001;

T_flue_i = 1200.0;
T_steam_i = 225.0;

P_steam_i = 25.0;
P_flue_i = 0.9;

T_flue_o = 300.0;
T_steam_o = 400.0;

FC = 0.95;

MaxIter = 100;
MaxError = 0.0001;
Ws = 0.25;

% INITIALIZE %
% =========================================================================

X_CARBON = (X_METHANE+3.0*X_PROPANE)*0.01;
X_HYDROGEN = (4.0*X_METHANE+8.0*X_PROPANE)*0.01;

M_H = 1.0;
M_O = 16.0;
M_C = 12.0;
M_N = 14.0;

M_O2 = 2.0*M_O;
M_N2 = 2.0*M_N;
M_CO2 = M_C+2.0*M_O;
M_H2O = 2.0*M_H+M_O;

alpha = (1.0+0.01*ExcessAir)*(X_CARBON+0.25*X_HYDROGEN);

n_O2 = alpha-X_CARBON-0.25*X_HYDROGEN;
n_N2 = 3.76*alpha;
n_CO2 = X_CARBON;
n_H2O = 0.5*X_HYDROGEN;

n_gas = n_O2+n_N2+n_CO2+n_H2O;

x_O2 = n_O2/n_gas;
x_N2 = n_N2/n_gas;
x_CO2 = n_CO2/n_gas;
x_H2O = n_H2O/n_gas;

M_gas = x_O2*M_O2+x_N2*M_N2+x_CO2*M_CO2+x_H2O*M_H2O;

y_O2 = x_O2*(M_O2/M_gas);
y_N2 = x_N2*(M_N2/M_gas);
y_CO2 = x_CO2*(M_CO2/M_gas);
y_H2O = x_H2O*(M_H2O/M_gas);

AFRm = (alpha*4.76*M_gas)/(M_C*X_CARBON+M_H*X_HYDROGEN);
AFRv = alpha*4.76;

Q_gas = BoilerCapacity/FuelHeatValue;
Q_air = AFRv*Q_gas;

Rho_gas = refpropm('D','T',T_gas,'P',P_gas,'METHANE','PROPANE',[0.01*X_METHANE 0.01*X_PROPANE]);
Rho_air = refpropm('D','T',T_air,'P',P_air,'AIR.PPF');

Mdot_flue = ((Q_gas*Rho_gas)+(Q_air*Rho_air))/3600.0;
Mdot_steam = (SteamFlowRate/3600.0)/Tube_Nc;

Tube_Di = Tube_Do-2*Tube_Wt;
Tube_Lt1 = Tube_Nc*(Tube_Nr*Tube_Ls+(Tube_Nr-1)*(0.5*pi*Tube_Ru));
Tube_Lt2 = Tube_Nc*Tube_Nr*Tube_Ls;
Tube_Ao = 0.25*pi*Tube_Do*Tube_Do;
Tube_Ai = 0.25*pi*Tube_Di*Tube_Di;
Tube_Sd = sqrt(0.25*Tube_St*Tube_St+Tube_Sl*Tube_Sl);
Tube_to = Tube_St/Tube_Do;
Tube_tl = Tube_St/Tube_Sl;

Duct_Wl = ((Tube_Nc-1)*Tube_St)/BankType+2.0*Tube_cl+Tube_Do;
Duct_Bl = 2.0*Tube_Ru+Tube_Ls+2.0*Tube_cl;
Duct_De = Tube_Do*(((4.0/pi)*(Tube_St/Tube_Do)^2)-1.0);
Duct_Ae = 0.25*pi*Duct_De*Duct_De;
Duct_Ar = (Duct_Wl*(Duct_Bl-Duct_corbel))-((Tube_Do*Tube_Ls*Tube_Nc)/BankType);
Duct_Am = min([(Tube_St-Tube_Do),(2.0*(Tube_Sd-Tube_Do))]);
Duct_Ht = (Tube_Nr+BankType-1)*2.0*Tube_Ru;

if ( BankType == 1 )
    if ( Tube_Nr == 1 )
        F = 0.70;
    elseif ( Tube_Nr == 2 )
        F = 0.80;
    elseif ( Tube_Nr == 3 )
        F = 0.86;
    elseif ( Tube_Nr == 4 )
        F = 0.90;
    elseif ( Tube_Nr == 5 || Tube_Nr == 6 )
        F = 0.93;
    elseif ( Tube_Nr == 7 || Tube_Nr == 8 || Tube_Nr == 9 )
        F = 0.96;
    elseif ( Tube_Nr == 10 || Tube_Nr == 11 || Tube_Nr == 12 )
        F = 0.98;
    elseif ( Tube_Nr == 13 || Tube_Nr == 14 || Tube_Nr == 15 )
        F = 0.99;
    elseif ( Tube_Nr > 15 )
        F = 1.0;
    end
else
    if ( Tube_Nr == 1 )
        F = 0.64;
    elseif ( Tube_Nr == 2 )
        F = 0.76;
    elseif ( Tube_Nr == 3 )
        F = 0.84;
    elseif ( Tube_Nr == 4 )
        F = 089;
    elseif ( Tube_Nr == 5 || Tube_Nr == 6 )
        F = 0.93;
    elseif ( Tube_Nr == 7 || Tube_Nr == 8 || Tube_Nr == 9 )
        F = 0.96;
    elseif ( Tube_Nr == 10 || Tube_Nr == 11 || Tube_Nr == 12 )
        F = 0.98;
    elseif ( Tube_Nr == 13 || Tube_Nr == 14 || Tube_Nr == 15 )
        F = 0.99;
    elseif ( Tube_Nr > 15 )
        F = 1.0;
    end
end

% SOLVER %
% =========================================================================
I = 1;
Error = 1.0;
fprintf('% 10d : %16.12f\n',I,Error);
while( Error > MaxError )
    
    T_flue_b = 0.5*(T_flue_i+T_flue_o);
    T_steam_b = 0.5*(T_steam_i+T_steam_o);
    
    Cp_steam = refpropm('C','T',T_steam_b+273.15,'P',100.0*P_steam_i,'WATER');
    Rho_steam = refpropm('D','T',T_steam_b+273.15,'P',100.0*P_steam_i,'WATER');
    Mu_steam = refpropm('V','T',T_steam_b+273.15,'P',100.0*P_steam_i,'WATER');
    K_steam = refpropm('L','T',T_steam_b+273.15,'P',100.0*P_steam_i,'WATER');
    Pr_steam = refpropm('^','T',T_steam_b+273.15,'P',100.0*P_steam_i,'WATER');
    
    [ data_flue ] = FlueProp( T_flue_b,P_flue_i,y_O2,y_N2,y_CO2,y_H2O );
    Cp_flue = data_flue.cp;
    Rho_flue = data_flue.rho;
    Mu_flue = data_flue.mu;
    K_flue = data_flue.k;
    Pr_flue = data_flue.pr;
    
    C_hot = Mdot_flue*Cp_flue;
    C_cold = Mdot_steam*Cp_steam*Tube_Nc;
    
    C_min = C_hot;
    C_max = C_cold;
    if ( C_cold < C_hot )
        C_min = C_cold;
        C_max = C_hot;
    end
    R = C_min/C_max;
    
    V_flue = Mdot_flue/(Duct_Ar*Rho_flue);
    if ( BankType == 1 )
        V_flue_max = (Tube_St/(Tube_St-Tube_Do))*V_flue;
    elseif ( BankType == 2 )
        cond = (2.0*Tube_Sd)/(Tube_St+Tube_Do);
        if ( cond >= 1.0 )
            V_flue_max = (Tube_St/(Tube_St-Tube_Do))*V_flue;
        else
            V_flue_max = 0.5*(Tube_St/(Tube_Sd-Tube_Do))*V_flue;
        end
    end
    V_steam = Mdot_steam/(Tube_Ai*Rho_steam);
    
    Re_flue = (V_flue_max*Tube_Do)/(Mu_flue/Rho_flue);
    Re_steam = (V_steam*Tube_Di)/(Mu_steam/Rho_steam);
    
    if ( BankType == 1 )
        if ( Re_flue <= 100 )
            C = 0.90;
            m = 0.40;
            n = 0.36;
        elseif ( 100 < Re_flue && Re_flue <= 1000 )
            C = 0.52;
            m = 0.50;
            n = 0.36;
        elseif ( 1000 < Re_flue && Re_flue <= 200000 )
            C = 0.27;
            m = 0.63;
            n = 0.36;
        elseif ( 200000 < Re_flue && Re_flue < 2000000 )
            C = 0.033;
            m = 0.80;
            n = 0.40;
        end
    else
        if ( Re_flue <= 500 )
            C = 1.04;
            m = 0.40;
            n = 0.36;
        elseif ( 500 < Re_flue && Re_flue <= 1000 )
            C = 0.71;
            m = 0.50;
            n = 0.36;
        elseif ( 1000 < Re_flue && Re_flue <= 200000 )
            C = 0.35;
            m = 0.60;
            n = 0.36;
        elseif ( 200000 < Re_flue && Re_flue < 2000000 )
            C = 0.031;
            m = 0.80;
            n = 0.36;
        end
    end
    
    w = 0.0;
    if ( BankType == 1 && Re_flue > 1000 )
        w = 0.20;
    end
    
    if ( NuFType == 1 )
        Nusselt_flue = F*C*((Tube_St/Tube_Sl)^w)*(Re_flue^m)*(Pr_flue^n);
    elseif ( NuFType == 2 )
        Nusselt_flue = 2.0+(0.4*(Re_flue^(1/2))+0.06*(Re_flue^(2/3)));
    elseif ( NuFType == 3 )
        Nusselt_flue = 0.3+((0.62*(Re_flue^(1/2))*(Pr_flue^(1/3)))/((1.0+((0.4/Pr_flue)^(2/3)))^(1/4)))*((1.0+(Re_flue/282000)^(5/8))^(4/5));
    end
    
    if ( NuSType == 1 )
        a = 1.0/(1.0+(Re_steam/2712.0)^8.4);
        b = 1.0/(1.0+(Re_steam/(150.0/(Tube_Rn/Tube_Di)))^1.8);
        fe = ((64.0/Re_steam)^a)*((0.75*log(Re_steam/5.37))^(2.0*(a-1.0)*b))*((0.88*log(3.41/(Tube_Rn/Tube_Di)))^(2.0*(a-1.0)*(1.0-b)));
        Nusselt_steam = (0.125*fe*(Re_steam-1000.0)*Pr_steam)/(1.0+12.7*sqrt(0.125*fe)*((Pr_steam^(2.0/3.0))-1.0));
    elseif ( NuSType == 2 )
        syms f
        fe = double(solve((1/sqrt(f)) + 2.0*log10(((Tube_Rn/Tube_Di)/3.7)+(2.51/(Re_steam*sqrt(f))))));
        Nusselt_steam = 0.023*(Re_steam^0.80)*(Pr_steam^0.40);
    elseif ( NuSType == 3 )
        fe = (0.79*log(Re_steam)-1.64)^(-2.0);
        Nusselt_steam = 0.125*fe*Re_steam*(Pr_steam^(1/3));
    end
    
    h_flue = (Nusselt_flue*K_flue)/(Tube_Do);
    h_steam = (Nusselt_steam*K_steam)/(Tube_Di);
    
    R_O = 1.0/(h_flue*pi*Tube_Do*Tube_Lt2);
    R_I = 1.0/(h_steam*pi*Tube_Di*Tube_Lt2);
    R_W = (log(Tube_Do/Tube_Di))/(2.0*pi*Tube_Lt2*Tube_Kw);
    UA = 1.0/(R_O+R_I+R_W+Tube_Rfi+Tube_Rfo);
    
    if ( SolverType == 1 )
        NTU = UA/C_min;
        if ( EHEXType == 1)
            E = (1.0-exp(-NTU*(1.0-R)))/(1.0-R*exp(-NTU*(1.0-R)));
        elseif ( EHEXType == 2)
            E = 1.0-exp((NTU^0.22)*(exp(-R*(NTU^0.78))-1.0)*(1.0/R));
        end
        Q = E*C_min*(T_flue_i-T_steam_i);
    elseif ( SolverType == 2 )
        DTM = 0.5*(T_flue_o+T_flue_i)-0.5*(T_steam_o+T_steam_i);
        LMTD = ((T_flue_o-T_steam_i)-(T_flue_i-T_steam_o))/log((T_flue_o-T_steam_i)/(T_flue_i-T_steam_o));
        Q = FC*UA*LMTD;
    end
    
    T_flue = T_flue_i-(Q/C_hot);
    T_steam = T_steam_i+(Q/C_cold);
    
    Th = T_flue_o;
    Tc = T_steam_o;
    T_flue_o = 0.5*(T_flue_o+T_flue*(1.00+0.01*Flue_CorrectionFactor));
    T_steam_o = 0.5*(T_steam_o+T_steam*(1.00+0.01*Steam_CorrectionFactor));
    
    Error = max(abs([(Th-T_flue_o),(Tc-T_steam_o)]));
    I = I+1;
    if ( I > MaxIter )
        break;
    end
    
    fprintf('% 10d : %16.12f\n',I,Error);
    
end

% PRINT %
% =========================================================================

fprintf('=============================================================\n');
fprintf('Flue Temperature  (IN - OUT)  :  %10.2f  - %10.2f (0C) \n',T_flue_i,T_flue_o);
fprintf('Steam Temperature (IN - OUT)  :  %10.2f  - %10.2f (0C) \n',T_steam_i,T_steam_o);
fprintf('Process Duty                  :  %10.2f  (kW) \n',Q/1000.0);
fprintf('=============================================================\n\n');






