function [STATES, CONSTANTS] = initConsts(protocol)
    
if (protocol==0 | protocol==1)
    flag = 0;
else
    flag = 1;
end

VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];

switch(flag)
    case 0
    switch(protocol)
        case 0
        CONSTANTS(:,1) = 1;          % Specific membrane capacitance - Cm [uF/cm2]
        case 1
        CONSTANTS(:,1) = 1.22;       % Specific membrane capacitance - Cm [uF/cm2]
    end
    CONSTANTS(:,2) = 2.2e-5;         % Myoplasmic volume - Vmyo [uL]
    CONSTANTS(:,3) = 1.022e-7;       % Junctional SR volume - VJSR [uL]
    CONSTANTS(:,4) = 1.786e-6;       % Network SR volume VNSR [uL]
    CONSTANTS(:,5) = 1.264e-9;       % Subspace volume - Vss [uL]
    CONSTANTS(:,6) = 0.00013866;     % Capacitive membrane area - Acap [cm2]
    CONSTANTS(:,7) = 5400;           % Extracellular K+ concentration - Ko [uM]
    CONSTANTS(:,8) = 134000;         % Extracellular Na+ concentration - Nao [uM]
    CONSTANTS(:,9) = 1400;           % Extracellular Ca+ concentration - Cao [uM]
    CONSTANTS(:,10) = 8.314;         % Ideal gas constant - R [J/molK]
    CONSTANTS(:,11) = 308.00;        % Absolute temperature - T [K]
    CONSTANTS(:,12) = 96.5;          % Faraday constant - F [C/mmol]
    CONSTANTS(:,13) = 0;             % Stimulation offset - stim_offset [ms]
    CONSTANTS(:,14) = 1000;          % Stimulation period - stim_period [ms]
    CONSTANTS(:,15) = 5;             % Stimulation duration - stim_duration [ms]
    CONSTANTS(:,16) = -15;           % Stimulation amplitude - stim_amplitude [pA/pF]
    CONSTANTS(:,17) = 109;           % Total myoplasmic calmodulin concentration - CMDNtot [uM]
    CONSTANTS(:,18) = 15000;         % Total junctional SR calsequestrin concentration - CSQNtot [uM]
    CONSTANTS(:,19) = 0.6;           % Ca2+ half-saturation constant for calmodulin - K_m_CMDN [uM]
    CONSTANTS(:,20) = 800;           % Ca2+ half-saturation constant for calsequestrin - K_m_CSQN [uM]
    CONSTANTS(:,21) = 4.5;           % Maximum RyR channel Ca2+ permeability - v1 [1/ms] 
    CONSTANTS(:,22) = 20;            % Time constant for transfer from NSR to JSR - tau_tr [ms] 
    CONSTANTS(:,23) = 1e-5;          % Ca2+ leak rate constant from the NSR - v2 [1/ms]
    CONSTANTS(:,24) = 8;             % Time constant for tansfer from subspace to myoplasm - tau_xfer [ms]
    CONSTANTS(:,25) = 0.05;          % on_rate [1/ms]
    CONSTANTS(:,26) = 0.00068;       % off_rate [1/ms]
    CONSTANTS(:,27) = 0.412;         % Half-saturation constant for SR Ca2+-ATPase pump - Km_up [uM]
    CONSTANTS(:,28) = 4;             % Normalization constant for L-type Ca2+ current - ICaL_max [pA/pF]
    CONSTANTS(:,29) = 0.3158;        % vmup_init [uM/ms]
    CONSTANTS(:,30) = -0.09;         % P_RyR_const1 [1/ms]
    CONSTANTS(:,31) = -0.225;        % P_RyR_const2 [1/ms]
    CONSTANTS(:,32) = 0.006075;      % RyR: Pc1-Po1 rate constant - K+_a [uM4/ms]
    CONSTANTS(:,33) = 0.07125;       % RyR: Po1-Pc1 rate constant - K-_a [1/ms]
    CONSTANTS(:,34) = 0.00405;       % RyR: Po1-Po2 rate constant - K+_b [uM3/ms]
    CONSTANTS(:,35) = 0.965;         % RyR: Po2-Po1 rate constant - K-_b [1/ms]
    CONSTANTS(:,36) = 0.009;         % RyR: Po1-Pc2 rate constant - K+_c [1/ms]
    CONSTANTS(:,37) = 0.0008;        % RyR: Pc2-Po1 rate constant - K-_c [1/ms]
    CONSTANTS(:,38) = 3;             % RyR: Ca2+ cooperative parameter Po1-Po2 - m [-]
    CONSTANTS(:,39) = 4;             % RyR: Ca2+ cooperative parameter Pc1-Po1 - n [-]
    CONSTANTS(:,40) = 19.1078;       % P_CaL [1/ms]    
    CONSTANTS(:,41) = -2;            % VL [mV]
    CONSTANTS(:,42) = 7.0671;        % delta_VL [mV]
    CONSTANTS(:,43) = 1.1683;        % Time switching between C and O states - tL [ms]
    CONSTANTS(:,44) = 1.6411;        % Proportion of time closed in open mode - phiL [-]
    CONSTANTS(:,45) = 0.07;          % Biasing to make inactivation function of V - a [-]
    CONSTANTS(:,46) = 14;            % Biasing to make inactivation function of V - b [-]
    CONSTANTS(:,47) = 972.9715;      % Inactivation time - tauL [ms]
    CONSTANTS(:,48) = 0.0964;        % Concentration at inactivation - KL [uM]
    CONSTANTS(:,49) = 6.6755;        % C_2 [mV]
    CONSTANTS(:,50) = 0.9;           % Maximum Ca2+ pump current - IpCa_max [pA/pF]
    CONSTANTS(:,51) = 0.4;           % Ca2+ half-saturation  constant for Ca2+ pump current - K_mpCa [uM]
    CONSTANTS(:,52) = 772.8991;      % Scaling factor of Na+/Ca2+ exchange - k_NaCa [pA/pF]
    CONSTANTS(:,53) = 86500;         % Na+ half-saturation constant for Na+/ca2+ exchange - K_mNa [uM]
    CONSTANTS(:,54) = 1380;          % Ca2+ half-saturation constant for Na+/ca2+ exchange - K_mCa [uM]
    CONSTANTS(:,55) = 0.1;           % Na+/Ca2+ exchange saturation factor - k_sat [-]
    CONSTANTS(:,56) = 0.35;          % Controls voltage dependence of Na+/Ca2+ exchange - eta [-]
    CONSTANTS(:,57) = 0.00088;       % ICab: Max conductance - GCab [mS/uF]
    CONSTANTS(:,58) = 13;            % INa: Max conductance - GNa [mS/uF]
    CONSTANTS(:,59) = 0.0026;        % INab: Max conductance - GNab [mS/uF]
    CONSTANTS(:,60) = 0.4067;        % Iktof: Max conductance - Gktof [mS/uF]
    CONSTANTS(:,61) = 0;             % Max conductance - GKtos [mS/uF]
    CONSTANTS(:,62) = 0.2938;        % IK1: Max conductance - GK1 [mS/uF]
    CONSTANTS(:,63) = 0.00575;       % IKs: Max conductance - GKs [mS/uF]
    CONSTANTS(:,64) = 0.16;          % IKur: Max conductance - GKur [mS/uF]
    CONSTANTS(:,65) = 0.05;          % IKss: Max conductance - GKss [mS/uF]
    CONSTANTS(:,66) = 0.078;         % IKr: Max conductance - GKr [mS/uF]
    CONSTANTS(:,67) = 0.036778;      % Rate constant for rapid delayed-rectifier K+ current - kb [1/ms]
    CONSTANTS(:,68) = 0.023761;      % Rate constant for rapid delayed-rectifier K+ current - kf [1/ms]
    CONSTANTS(:,69) = 1.66;          % Max exchange current - INaK_max [pA/pF]
    CONSTANTS(:,70) = 21000;         % Na+ half-saturation constant for Na+/K+ exchange current - K_mNai [uM]
    CONSTANTS(:,71) = 1500;          % K+ half-saturation constant for Na+/K+ exchange current - K_mKo [uM]
    CONSTANTS(:,72) = 10;            % ICaCl: Max conductance - G_ClCa [mS/mF]
    CONSTANTS(:,73) = -40;           % ICaCl: Reversal potential E_Cl [mV]
    CONSTANTS(:,74) = 10;            % ICaCl: Half saturation constant K_mCl [uM]
    CONSTANTS(:,75) = CONSTANTS(:,44)./CONSTANTS(:,43);                % alpha_m [1/ms]
    CONSTANTS(:,76) = (1.0./7.0).*(exp(CONSTANTS(:,8)./67300.0)-1.0);  % sigma [-]
    CONSTANTS(:,76) = 0.0;

    STATES(:,1) = -85.64004;         % Membrane potential - V [mV]
    STATES(:,2) = 0.1040595;         % Myoplasmic Ca2+ concentration - Cai [uM]
    STATES(:,3) = 0.1043777;         % Subspace SR Ca2+ concentration - Cass [uM]
    STATES(:,4) = 730.0589;          % JSR Ca2+ concentration - CaJSR [uM]
    STATES(:,5) = 841.106;           % NSR Ca2+ concentration - CaNSR [uM]
    STATES(:,6) = 2.290355e-9;       % RyR modulation factor - P_RyR [-]
    STATES(:,7) = 0.08989079;        % CaMKt [-]
    STATES(:,8) = 0.003825599;       % Fraction of RyR channels in state Po1 - P_O1 [-]
    STATES(:,9) = 1.835831e-8;       % Fraction of RyR channels in state Po2 - P_O2 [-]
    STATES(:,10) = 0.3797679;        % Fraction of RyR channels in state Pc2 - P_C2 [-]
    STATES(:,11) = 4.373318e-6;      % L-type Ca2+ channel conducting state - O [-]
    STATES(:,12) = 0.009171979;      % I [-]
    STATES(:,13) = 0.8876797;        % y_gate [-]
    STATES(:,14) = 16522.45;         % Myoplasmic Na+ concentration - Nai [uM]
    STATES(:,15) = 2.639399e-7;      % Open state of fast Na+ channel - O_Na [-]
    STATES(:,16) = 0.0001581035;     % Closed state of fast Na+ channel - C_Na1 [-]
    STATES(:,17) = 0.01702105;       % Closed state of fast Na+ channel - C_Na2 [-]
    STATES(:,18) = 0.00001799179;    % Slow inactivated state 1 of fast Na+ channel - I1_Na [-]
    STATES(:,19) = 0.000005460299;   % Slow inactivated state 2 of fast Na+ channel - I2_Na [-]
    STATES(:,20) = 0.0000556206;     % Fast inactivated state of fast Na+ channel - IF_Na [-]
    STATES(:,21) = 0.005985434;      % Closed-inactivated state of fast Na+ channel - IC_Na2 [-]
    STATES(:,22) = 0.2543133;        % Closed-inactivated state of fast Na+ channel - IC_Na3 [-]
    STATES(:,23) = 141474;           % Myoplasmic K+ concentration - Ki [uM]
    STATES(:,24) = 0.001937245;      % Gating variable for transient outward K+ current - ato_f [-]
    STATES(:,25) = 0.9999985;        % Gating variable for transient outward K+ current - ito_f [-]
    STATES(:,26) = 0.02000568;       % Gating variable for transient outward K+ current - ato_s [-]
    STATES(:,27) = 0.9308568;        % Gating variable for transient outward K+ current - ito_s [-]
    STATES(:,28) = 0.002206261;      % Gating variable for slow delayed-rectifier K+ current - nKs [-]
    STATES(:,29) = 0.02000568;       % Gating variable for ultrarapidly activating  delayed-rectifier K+ current - aur [-]
    STATES(:,30) = 0.9822006;        % Gating variable for ultrarapidly activating  delayed-rectifier K+ current - iur [-]
    STATES(:,31) = 0.8883113;        % Gating variable for noninactivating steady-state K+ current - aKss [-]
    STATES(:,32) = 1;                % Gating variable for noninactivating steady-state K+ current - iKss [-]
    STATES(:,33) = 0.0004858865;     % mERG channel open state - O_K [-]
    STATES(:,34) = 0.0007799137;     % mERG channel closed state - C_K1 [-]
    STATES(:,35) = 0.0005301217;     % mERG channel closed state - C_K2 [-]
    STATES(:,36) = 0.00007519518;    % mERG channel inactivated state - I_K [-]

    case 1
    CONSTANTS(:,1) = 1.22;           % Specific membrane capacitance - Cm0 [uF/cm2]
    CONSTANTS(:,2) = 2.2e-5;         % Myoplasmic volume - Vmyo [uL]
    CONSTANTS(:,3) = 1.022e-7;       % Junctional SR volume - VJSR [uL]
    CONSTANTS(:,4) = 1.786e-6;       % Network SR volume VNSR [uL]
    CONSTANTS(:,5) = 1.264e-9;       % Subspace volume - Vss [uL]
    CONSTANTS(:,6) = 0.00013866;     % Capacitive membrane area - Acap [cm2]
    CONSTANTS(:,7) = 5400;           % Extracellular K+ concentration - Ko [uM]
    CONSTANTS(:,8) = 134000;         % Extracellular Na+ concentration - Nao [uM]
    CONSTANTS(:,9) = 2000;           % Extracellular Ca+ concentration - Cao [uM]
    CONSTANTS(:,10) = 8.314;         % Ideal gas constant - R [J/molK]
    CONSTANTS(:,11) = 308.00;        % Absolute temperature - T [K]
    CONSTANTS(:,12) = 96.5;          % Faraday constant - F [C/mmol]
    CONSTANTS(:,13) = 0;             % Stimulation offset - stim_offset [ms]
    CONSTANTS(:,14) = 1000;          % Stimulation period - stim_period [ms]
    CONSTANTS(:,15) = 0;             % Stimulation duration - stim_duration [ms]
    CONSTANTS(:,16) = 0;             % Stimulation amplitude - stim_amplitude [pA/pF]
    CONSTANTS(:,17) = 109;           % Total myoplasmic calmodulin concentration - CMDNtot [uM]
    CONSTANTS(:,18) = 15000;         % Total junctional SR calsequestrin concentration - CSQNtot [uM]
    CONSTANTS(:,19) = 0.6;           % Ca2+ half-saturation constant for calmodulin - K_m_CMDN [uM]
    CONSTANTS(:,20) = 800;           % Ca2+ half-saturation constant for calsequestrin - K_m_CSQN [uM]
    CONSTANTS(:,21) = 4.5;           % Maximum RyR channel Ca2+ permeability - v1 [1/ms]
    CONSTANTS(:,22) = 20;            % Time constant for transfer from NSR to JSR - tau_tr [ms] 
    CONSTANTS(:,23) = 1e-5;          % Ca2+ leak rate constant from the NSR - v2 [1/ms]
    CONSTANTS(:,24) = 8;             % Time constant for tansfer from subspace to myoplasm - tau_xfer [ms]
    CONSTANTS(:,25) = 0.05;          % on_rate [1/ms]
    CONSTANTS(:,26) = 0.00068;       % off_rate [1/ms]
    CONSTANTS(:,27) = 0.412;         % Half-saturation constant for SR Ca2+-ATPase pump - Km_up [uM]
    CONSTANTS(:,28) = 4;             % Normalization constant for L-type Ca2+ current - ICaL_max [pA/pF]
    CONSTANTS(:,29) = 0.3158;        % vmup_init [uM/ms]
    CONSTANTS(:,30) = -0.09;         % P_RyR_const1 [1/ms]
    CONSTANTS(:,31) = -0.225;        % P_RyR_const2 [1/ms]
    CONSTANTS(:,32) = 0.006075;      % RyR: Pc1-Po1 rate constant - K+_a [uM4/ms]
    CONSTANTS(:,33) = 0.07125;       % RyR: Po1-Pc1 rate constant - K-_a [1/ms]
    CONSTANTS(:,34) = 0.00405;       % RyR: Po1-Po2 rate constant - K+_b [uM3/ms]
    CONSTANTS(:,35) = 0.965;         % RyR: Po2-Po1 rate constant - K-_b [1/ms]
    CONSTANTS(:,36) = 0.009;         % RyR: Po1-Pc2 rate constant - K+_c [1/ms]
    CONSTANTS(:,37) = 0.0008;        % RyR: Pc2-Po1 rate constant - K-_c [1/ms]
    CONSTANTS(:,38) = 3;             % RyR: Ca2+ cooperative parameter Po1-Po2 - m [-]
    CONSTANTS(:,39) = 4;             % RyR: Ca2+ cooperative parameter Pc1-Po1 - n [-]
    CONSTANTS(:,40) = 19.1078;       % P_CaL [1/ms]
    CONSTANTS(:,41) = -2;            % VL [mV]
    CONSTANTS(:,42) = 7.0671;        % delta_VL [mV]
    CONSTANTS(:,43) = 1.1683;        % Time switching between C and O states - tL [ms]
    CONSTANTS(:,44) = 1.6411;        % Proportion of time closed in open mode - phiL [-]
    CONSTANTS(:,45) = 0.07;          % Biasing to make inactivation function of V - a [-]
    CONSTANTS(:,46) = 14;            % Biasing to make inactivation function of V - b [-]
    CONSTANTS(:,47) = 972.9715;      % Inactivation time - tauL [ms]
    CONSTANTS(:,48) = 0.0964;        % Concentration at inactivation - KL [uM]
    CONSTANTS(:,49) = 6.6755;        % C_2 [mV]
    CONSTANTS(:,50) = 0.9;           % Maximum Ca2+ pump current - IpCa_max [pA/pF]
    CONSTANTS(:,51) = 0.4;           % Ca2+ half-saturation  constant for Ca2+ pump current - K_mpCa [uM]
    CONSTANTS(:,52) = 772.8991;      % Scaling factor of Na+/Ca2+ exchange - k_NaCa [pA/pF]
    CONSTANTS(:,53) = 86500;         % Na+ half-saturation constant for Na+/ca2+ exchange - K_mNa [uM]
    CONSTANTS(:,54) = 1380;          % Ca2+ half-saturation constant for Na+/ca2+ exchange - K_mCa [uM]
    CONSTANTS(:,55) = 0.1;           % Na+/Ca2+ exchange saturation factor - k_sat [-]
    CONSTANTS(:,56) = 0.35;          % Controls voltage dependence of Na+/Ca2+ exchange - eta [-]
    CONSTANTS(:,57) = 0.00088;       % ICab: Max conductance - GCab [mS/uF]
    CONSTANTS(:,58) = 13;            % INa: Max conductance - GNa [mS/uF]
    CONSTANTS(:,59) = 0.0026;        % INab: Max conductance - GNab [mS/uF]
    CONSTANTS(:,60) = 0.4067;        % Iktof: Max conductance - Gktof [mS/uF]
    CONSTANTS(:,61) = 0;             % Max conductance - GKtos [mS/uF]
    CONSTANTS(:,62) = 0.2938;        % IK1: Max conductance - GK1 [mS/uF]
    CONSTANTS(:,63) = 0.00575;       % IKs: Max conductance - GKs [mS/uF]
    CONSTANTS(:,64) = 0.16;          % IKur: Max conductance - GKur [mS/uF]
    CONSTANTS(:,65) = 0.05;          % IKss: Max conductance - GKss [mS/uF]
    CONSTANTS(:,66) = 0.078;         % IKr: Max conductance - GKr [mS/uF]
    CONSTANTS(:,67) = 0.036778;      % Rate constant for rapid delayed-rectifier K+ current - kb [1/ms]
    CONSTANTS(:,68) = 0.023761;      % Rate constant for rapid delayed-rectifier K+ current - kf [1/ms]
    CONSTANTS(:,69) = 1.66;          % Max exchange current - INaK_max [pA/pF]
    CONSTANTS(:,70) = 21000;         % Na+ half-saturation constant for Na+/K+ exchange current - K_mNai [uM]
    CONSTANTS(:,71) = 1500;          % K+ half-saturation constant for Na+/K+ exchange current - K_mKo [uM]
    CONSTANTS(:,72) = 10;            % ICaCl: Max conductance - G_ClCa [mS/mF]
    CONSTANTS(:,73) = -40;           % ICaCl: Reversal potential E_Cl [mV]
    CONSTANTS(:,74) = 10;            % ICaCl: Half saturation constant K_mCl [uM]
    CONSTANTS(:,75) = CONSTANTS(:,44)./CONSTANTS(:,43);                % alpha_m [1/ms]
    CONSTANTS(:,76) = (1.0./7.0).*(exp(CONSTANTS(:,8)./67300.0)-1.0);  % sigma [-]
    CONSTANTS(:,76) = 0.0;

    switch(protocol)
        case 2
        CONSTANTS(:,77) = CONSTANTS(:,13);                                 % t0 [ms]
        CONSTANTS(:,78) = CONSTANTS(:,13)+10;                              % t1 [ms]
        CONSTANTS(:,79) = CONSTANTS(:,13)+11;                              % t2 [ms]
        CONSTANTS(:,80) = CONSTANTS(:,13)+20;                              % t3 [ms]
        CONSTANTS(:,81) = CONSTANTS(:,13)+120;                             % t4 [ms]
        CONSTANTS(:,82) = 0.93;                                            % Amax [-]
        CONSTANTS(:,83) = 10;                                              % a
        CONSTANTS(:,84) = 15;                                              % b
        CONSTANTS(:,85) = 6;                                               % c
    
        % Evaluation of coefficients
        matrix = [CONSTANTS(:,78)^3    CONSTANTS(:,78)^2  CONSTANTS(:,78)  1  0                   0                 0;
                  3*CONSTANTS(:,78)^2  2*CONSTANTS(:,78)  1                0  0                   0                 0;
                  CONSTANTS(:,79)^3    CONSTANTS(:,79)^2  CONSTANTS(:,79)  1  -CONSTANTS(:,79)^2  -CONSTANTS(:,79)  -1;
                  3*CONSTANTS(:,79)^2  2*CONSTANTS(:,79)  1                0  -2*CONSTANTS(:,79)  -1                0;
                  6*CONSTANTS(:,79)    2                  0                0  -2                  0                 0;
                  0                    0                  0                0  CONSTANTS(:,80)^2   CONSTANTS(:,80)   1;
                  0                    0                  0                0  2*CONSTANTS(:,80)   1                 0];
        rhs = [1.105; -0.0065625; 0; 0; 0; 1.07; 0];
        sol = matrix\rhs;
    
        CONSTANTS(:,86) = sol(1);                                          % d
        CONSTANTS(:,87) = sol(2);                                          % e
        CONSTANTS(:,88) = sol(3);                                          % f
        CONSTANTS(:,89) = sol(4);                                          % g
        CONSTANTS(:,90) = sol(5);                                          % h
        CONSTANTS(:,91) = sol(6);                                          % i
        CONSTANTS(:,92) = sol(7);                                          % l
    
        % SACs
        CONSTANTS(:,93) = 1.0015;                                               % Lambda1
        CONSTANTS(:,94) = 1.01;                                                 % Lambda2
        CONSTANTS(:,95) = 0.2203;                                               % Beta_ns [-]
        CONSTANTS(:,96) = 1.411;                                                % r
        CONSTANTS(:,97) = 3.6121e-04;                                           % G_ns [mS]
        CONSTANTS(:,98) = 0.8742;                                               % Beta_Ko [-]
        CONSTANTS(:,99) = 4.4873e-04;                                           % G_Ko
        CONSTANTS(:,100) = 4.0386e-04;                                          % G_CaP
        CONSTANTS(:,101) = 74;                                                  % Diel_const (water 37 °C)
        CONSTANTS(:,102) = 1.82*10^6*(CONSTANTS(:,101)*CONSTANTS(:,11))^(-1.5); % Const_A
    
        CONSTANTS(:,103) = CONSTANTS(:,13)+10;                                  % t_opening
        CONSTANTS(:,104) = CONSTANTS(:,13)+65;                                  % t_start_stretch
        CONSTANTS(:,105) = CONSTANTS(:,13)+100;                                 % t_peak_stretch
        CONSTANTS(:,106) = CONSTANTS(:,13)+300;                                 % t_end_relax
        
        case 3
        CONSTANTS(:,77) = CONSTANTS(:,13);                                 % t0 [ms]
        CONSTANTS(:,78) = CONSTANTS(:,13)+10;                              % t1 [ms]
        CONSTANTS(:,79) = CONSTANTS(:,13)+11;                              % t2 [ms]
        CONSTANTS(:,80) = CONSTANTS(:,13)+200;                             % t3 [ms]
        CONSTANTS(:,81) = CONSTANTS(:,13)+300;                             % t4 [ms]
        CONSTANTS(:,82) = 0.93;                                            % Amax [-]
        CONSTANTS(:,83) = 10;                                              % a
        CONSTANTS(:,84) = 15;                                              % b
        CONSTANTS(:,85) = 6;                                               % c
    
        % Evaluation of coefficients
        matrix = [CONSTANTS(:,78)^3    CONSTANTS(:,78)^2  CONSTANTS(:,78)  1  0                   0                 0;
                  3*CONSTANTS(:,78)^2  2*CONSTANTS(:,78)  1                0  0                   0                 0;
                  CONSTANTS(:,79)^3    CONSTANTS(:,79)^2  CONSTANTS(:,79)  1  -CONSTANTS(:,79)^2  -CONSTANTS(:,79)  -1;
                  3*CONSTANTS(:,79)^2  2*CONSTANTS(:,79)  1                0  -2*CONSTANTS(:,79)  -1                0;
                  6*CONSTANTS(:,79)    2                  0                0  -2                  0                 0;
                  0                    0                  0                0  CONSTANTS(:,80)^2   CONSTANTS(:,80)   1;
                  0                    0                  0                0  2*CONSTANTS(:,80)   1                 0];
        rhs = [1.105; -0.0065625; 0; 0; 0; 1.07; 0];
        sol = matrix\rhs;
    
        CONSTANTS(:,86) = sol(1);                                          % d
        CONSTANTS(:,87) = sol(2);                                          % e
        CONSTANTS(:,88) = sol(3);                                          % f
        CONSTANTS(:,89) = sol(4);                                          % g
        CONSTANTS(:,90) = sol(5);                                          % h
        CONSTANTS(:,91) = sol(6);                                          % i
        CONSTANTS(:,92) = sol(7);                                          % l
    
        % SACs
        CONSTANTS(:,93) = 1.0015;                                                % Lambda1
        CONSTANTS(:,94) = 1.01;                                                % Lambda2
        CONSTANTS(:,95) = 0.2203;                                               % Beta_ns [-]
        CONSTANTS(:,96) = 1.411;                                                % r
        CONSTANTS(:,97) = 3.6121e-04;                                           % G_ns [mS]
        CONSTANTS(:,98) = 0.8742;                                               % Beta_Ko [-]
        CONSTANTS(:,99) = 4.4873e-04;                                           % G_Ko
        CONSTANTS(:,100) = 4.0386e-04;                                          % G_CaP
        CONSTANTS(:,101) = 74;                                                  % Diel_const (water 37 °C)
        CONSTANTS(:,102) = 1.82*10^6*(CONSTANTS(:,101)*CONSTANTS(:,11))^(-1.5); % Const_A
    
        CONSTANTS(:,103) = CONSTANTS(:,13)+10;                                  % t_opening
        CONSTANTS(:,104) = CONSTANTS(:,13)+210;                                 % t_start_stretch
        CONSTANTS(:,105) = CONSTANTS(:,13)+220;                                 % t_peak_stretch
        CONSTANTS(:,106) = CONSTANTS(:,13)+450;                                 % t_end_relax
    end

    STATES(:,1) = -85.64004;         % Membrane potential - V [mV]
    STATES(:,2) = 0.1040595;         % Myoplasmic Ca2+ concentration - Cai [uM]
    STATES(:,3) = 0.1043777;         % Subspace SR Ca2+ concentration - Cass [uM]
    STATES(:,4) = 730.0589;          % JSR Ca2+ concentration - CaJSR [uM]
    STATES(:,5) = 841.106;           % NSR Ca2+ concentration - CaNSR [uM]
    STATES(:,6) = 2.290355e-9;       % RyR modulation factor - P_RyR [-]
    STATES(:,7) = 0.08989079;        % CaMKt [-]
    STATES(:,8) = 0.003825599;       % Fraction of RyR channels in state Po1 - P_O1 [-]
    STATES(:,9) = 1.835831e-8;       % Fraction of RyR channels in state Po2 - P_O2 [-]
    STATES(:,10) = 0.3797679;        % Fraction of RyR channels in state Pc2 - P_C2 [-]
    STATES(:,11) = 4.373318e-6;      % L-type Ca2+ channel conducting state - O [-]
    STATES(:,12) = 0.009171979;      % I [-]
    STATES(:,13) = 0.8876797;        % y_gate [-]
    STATES(:,14) = 16522.45;         % Myoplasmic Na+ concentration - Nai [uM]
    STATES(:,15) = 2.639399e-7;      % Open state of fast Na+ channel - O_Na [-]
    STATES(:,16) = 0.0001581035;     % Closed state of fast Na+ channel - C_Na1 [-]
    STATES(:,17) = 0.01702105;       % Closed state of fast Na+ channel - C_Na2 [-]
    STATES(:,18) = 0.00001799179;    % Slow inactivated state 1 of fast Na+ channel - I1_Na [-]
    STATES(:,19) = 0.000005460299;   % Slow inactivated state 2 of fast Na+ channel - I2_Na [-]
    STATES(:,20) = 0.0000556206;     % Fast inactivated state of fast Na+ channel - IF_Na [-]
    STATES(:,21) = 0.005985434;      % Closed-inactivated state of fast Na+ channel - IC_Na2 [-]
    STATES(:,22) = 0.2543133;        % Closed-inactivated state of fast Na+ channel - IC_Na3 [-]
    STATES(:,23) = 141474;           % Myoplasmic K+ concentration - Ki [uM]
    STATES(:,24) = 0.001937245;      % Gating variable for transient outward K+ current - ato_f [-]
    STATES(:,25) = 0.9999985;        % Gating variable for transient outward K+ current - ito_f [-]
    STATES(:,26) = 0.02000568;       % Gating variable for transient outward K+ current - ato_s [-]
    STATES(:,27) = 0.9308568;        % Gating variable for transient outward K+ current - ito_s [-]
    STATES(:,28) = 0.002206261;      % Gating variable for slow delayed-rectifier K+ current - nKs [-]
    STATES(:,29) = 0.02000568;       % Gating variable for ultrarapidly activating  delayed-rectifier K+ current - aur [-]
    STATES(:,30) = 0.9822006;        % Gating variable for ultrarapidly activating  delayed-rectifier K+ current - iur [-]
    STATES(:,31) = 0.8883113;        % Gating variable for noninactivating steady-state K+ current - aKss [-]
    STATES(:,32) = 1;                % Gating variable for noninactivating steady-state K+ current - iKss [-]
    STATES(:,33) = 0.0004858865;     % mERG channel open state - O_K [-]
    STATES(:,34) = 0.0007799137;     % mERG channel closed state - C_K1 [-]
    STATES(:,35) = 0.0005301217;     % mERG channel closed state - C_K2 [-]
    STATES(:,36) = 0.00007519518;    % mERG channel inactivated state - I_K [-]

end

if (isempty(STATES))
    warning('Initial values for states not set');
end

end