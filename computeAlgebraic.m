function ALGEBRAIC = computeAlgebraic(ALGEBRAIC,CONSTANTS,STATES,VOI,protocol)
   
if (protocol==0 | protocol==1)
    flag = 0;
else
    flag = 1;
end

statesSize = size(STATES);
statesColumnCount = statesSize(2);
if (statesColumnCount == 1)
    STATES = STATES';
    utilOnes = 1;
else
    statesRowCount = statesSize(1);
    utilOnes = ones(statesRowCount,1);
end

switch(flag)
    case 0

    ALGEBRAIC(:,1) = floor(VOI./CONSTANTS(:,14)).*CONSTANTS(:,14);                                                 % past [ms]
    if (VOI-ALGEBRAIC(:,1)>=CONSTANTS(:,13) & VOI-ALGEBRAIC(:,1)<=CONSTANTS(:,13)+CONSTANTS(:,15))
        ALGEBRAIC(:,13) = CONSTANTS(:,16);                                                                         % Istim [pA/pF]
    else
        ALGEBRAIC(:,13) = 0.0;                                                                                     % Istim [pA/pF]
    end
    
    % ------------------------
    %   CALCIUM DYNAMICS
    % ------------------------
    
    % Calcium concentration
    ALGEBRAIC(:,27) = power(1.0+(CONSTANTS(:,17).*CONSTANTS(:,19))./power(CONSTANTS(:,19)+STATES(:,2),2.0),-1.0);  % Bi [-]
    ALGEBRAIC(:,32) = power(1.0+(CONSTANTS(:,17).*CONSTANTS(:,19))./power(CONSTANTS(:,19)+STATES(:,3),2.0),-1.0);  % Bss [-]
    ALGEBRAIC(:,36) = power(1.0+(CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,4),2.0),-1.0);  % BJSR [-]
    
    % Calcium fluxes
    ALGEBRAIC(:,39) = CONSTANTS(:,21).*(STATES(:,8)+STATES(:,9)).*(STATES(:,4)-STATES(:,3)).*STATES(:,6);          % Jrel [uM/ms]
    ALGEBRAIC(:,41) = (STATES(:,5)-STATES(:,4))./CONSTANTS(:,22);                                                  % Jtr [uM/ms]
    ALGEBRAIC(:,43) = (STATES(:,3)-STATES(:,2))./CONSTANTS(:,24);                                                  % Jxfer [uM/ms]
    ALGEBRAIC(:,45) = CONSTANTS(:,23).*(STATES(:,5)-STATES(:,2));                                                  % Jleak [uM/ms]
    ALGEBRAIC(:,47) = (0.05.*(1.0-STATES(:,7)).*1.0)./(1.0+0.7./STATES(:,3));                                      % CaMKb [-]
    ALGEBRAIC(:,49) = ALGEBRAIC(:,47)+STATES(:,7);                                                                 % CaMKa [-]
    ALGEBRAIC(:,51) = ((3.1512.*power(ALGEBRAIC(:,49),2.2062))./(power(0.1588,2.2062)+...
        power(ALGEBRAIC(:,49),2.2062))+1.0).*CONSTANTS(:,29);                                                      % vmup [uM/ms]
    ALGEBRAIC(:,53) = (ALGEBRAIC(:,51).*power(STATES(:,2), 2.0))./(power(CONSTANTS(:,27),2.0)+...
        power(STATES(:,2),2.0));                                                                                   % Jserca [uM/ms]
    
    % Ryanodine receptors
    ALGEBRAIC(:,2) = 1.0-(STATES(:,10)+STATES(:,8)+STATES(:,9));                                                   % Function of RyR channels in state Pc1 - PC1 [-]
        
    % Calcium L-type current
    ALGEBRAIC(:,55) = (CONSTANTS(:,12).*STATES(:,1))./( CONSTANTS(:,10).*CONSTANTS(:,11));                         % FVRT [-]
    ALGEBRAIC(:,57) = 2.0.*ALGEBRAIC(:,55);                                                                        % FVRT_Ca [-]
    ALGEBRAIC(:,4) = exp((STATES(:,1)-CONSTANTS(:,41))./CONSTANTS(:,42));                                          % exp_VL [-]
    ALGEBRAIC(:,15) = ALGEBRAIC(:,4)./(CONSTANTS(:,43).*(ALGEBRAIC(:,4)+1.0));                                     % alpha_p [1/ms]
    ALGEBRAIC(:,28) = (ALGEBRAIC(:,4)+CONSTANTS(:,45))./(CONSTANTS(:,47).*CONSTANTS(:,48).*(ALGEBRAIC(:,4)+1.0));  % epsilon_p [1/uMms]
    ALGEBRAIC(:,33) = (CONSTANTS(:,46).*(ALGEBRAIC(:,4)+CONSTANTS(:,45)))./...
        (CONSTANTS(:,47).*(CONSTANTS(:,46).*ALGEBRAIC(:,4)+CONSTANTS(:,45)));                                      % epsilon_m [1/ms]
    ALGEBRAIC(:,3) = 1.0./(1.0+exp((STATES(:,1)+16.6577)./CONSTANTS(:,49)))+0.1./...
        (1.0+exp((-STATES(:,1)+40.0)./6.0));                                                                       % y_inf [-]
    ALGEBRAIC(:,14) = 20.0+600.0./(1.0+exp((STATES(:,1)+30.0)./9.6));                                              % tau_y [ms]
    ALGEBRAIC(:,37) = (1.0-STATES(:,11))-STATES(:,12);                                                             % C [-]
    if (abs(ALGEBRAIC(:,57))>1.0e-05)
        ALGEBRAIC(:,59) = ((((-CONSTANTS(:,40).*2.0.*CONSTANTS(:,5).*...
            CONSTANTS(:,12))./(CONSTANTS(:,6).*CONSTANTS(:,1))).*STATES(:,11).*STATES(:,13).*ALGEBRAIC(:,57))./...
            (1.0-exp(-ALGEBRAIC(:,57)))).*(CONSTANTS(:,9).*exp(-ALGEBRAIC(:,57))-STATES(:,3));                     % ICaL [pA/pF]
    else
        ALGEBRAIC(:,59) = ((((-CONSTANTS(:,40).*2.0.*CONSTANTS(:,5).*CONSTANTS(:,12))./(CONSTANTS(:,6).*...
            CONSTANTS(:,1))).*STATES(:,11).*STATES(:,13).*1.0e-05)./(1.0-exp(-1.0e-05))).*(CONSTANTS(:,9).*...
            exp(-1.0e-05)-STATES(:,3));                                                                            % ICaL [pA/pF]
    end
    
    % Calcium pump current
    ALGEBRAIC(:,60) = (CONSTANTS(:,50).*power(STATES(:,2),2.0))./(power(CONSTANTS(:,51),2.0)+...
        power(STATES(:,2),2.0));                                                                                   % IpCa [pA/pF]
    ALGEBRAIC(:,62) = (-ALGEBRAIC(:,60).*CONSTANTS(:,6).*CONSTANTS(:,1))./...
        (2.0.*CONSTANTS(:,2).*CONSTANTS(:,12));                                                                    % JpCa [uM/ms]
      
    % Sodium/calcium exchange current
    ALGEBRAIC(:,61) = ((((((CONSTANTS(:,52).*1.0)./(power(CONSTANTS(:,53),3.0)+power(CONSTANTS(:,8),3.0))).*1.0)./...
        (CONSTANTS(:,54)+CONSTANTS(:,9))).*1.0)./(1.0+CONSTANTS(:,55).*exp(((CONSTANTS(:,56)-1.0).*STATES(:,1).*...
        CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*(exp((CONSTANTS(:,56).*STATES(:,1).*CONSTANTS(:,12))./...
        (CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,14),3.0).*CONSTANTS(:,9)-exp(((CONSTANTS(:,56)-1.0).*...
        STATES(:,1).*CONSTANTS(:,12))./(CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8),3.0).*...
        STATES(:,2));                                                                                              % INaCa [pA/pF]
    ALGEBRAIC(:,64) = (ALGEBRAIC(:,61).*CONSTANTS(:,6).*CONSTANTS(:,1))./(CONSTANTS(:,2).*CONSTANTS(:,12));        % Jncx [uM/ms]
    
    % Calcium background current
    ALGEBRAIC(:,63) = ((CONSTANTS(:,10).*CONSTANTS(:,11))./(2.0.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./...
        STATES(:,2));                                                                                              % ECaN [mV]
    ALGEBRAIC(:,65) = CONSTANTS(:,57).*(STATES(:,1)-ALGEBRAIC(:,63));                                              % ICab [pA/pF]
    ALGEBRAIC(:,67) = (-ALGEBRAIC(:,65).*CONSTANTS(:,6).*CONSTANTS(:,1))./(2.0.*CONSTANTS(:,2).*CONSTANTS(:,12));  % JCab [uM/ms]
    
    % ------------------------
    %   SODIUM DYNAMICS
    % ------------------------
    
    % Fast Na+ current
    ALGEBRAIC(:,66) = ((CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log((0.9.*CONSTANTS(:,8)+0.1.*...
        CONSTANTS(:,7))./(0.9.*STATES(:,14)+0.1.*STATES(:,23)));                                                   % ENa [mV]
    ALGEBRAIC(:,68) = CONSTANTS(:,58).*STATES(:,15).*(STATES(:,1)-ALGEBRAIC(:,66));                                % INa [pA/pF]
    ALGEBRAIC(:,5) = 1.0-(STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,20)+STATES(:,18)+STATES(:,19)+...
        STATES(:,21)+STATES(:,22));                                                                                % Closed state of fast Na+ channel - C_Na3 [-] 
    ALGEBRAIC(:,16) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./17.0)+0.2.*exp(-(STATES(:,1)+2.5)./150.0));          % alpha_Na11 [1/ms]
    ALGEBRAIC(:,29) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./15.0)+0.23.*exp(-(STATES(:,1)+2.5)./150.0));         % alpha_Na12 [1/ms]
    ALGEBRAIC(:,34) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./12.0)+0.25.*exp(-(STATES(:,1)+2.5)./150.0));         % alpha_Na13 [1/ms]
    ALGEBRAIC(:,38) = 0.1917.*exp(-(STATES(:,1)+2.5)./20.3);                                                       % beta_Na11 [1/ms]
    ALGEBRAIC(:,40) = 0.2.*exp(-(STATES(:,1)-2.5)./20.3);                                                          % beta_Na12 [1/ms]
    ALGEBRAIC(:,42) = 0.22.*exp(-(STATES(:,1)-7.5)./20.3);                                                         % beta_Na13 [1/ms]
    ALGEBRAIC(:,44) = 7.0e-07.*exp(-(STATES(:,1)+7.0)./7.7);                                                       % alfa_Na3 [1/ms]
    ALGEBRAIC(:,46) = 0.0084+2.0e-05.*(STATES(:,1)+7.0);                                                           % beta_Na3 [1/ms]
    ALGEBRAIC(:,48) = 1.0./(0.188495.*exp(-(STATES(:,1)+7.0)./16.6)+0.393956);                                     % alfa_Na2 [1/ms]
    ALGEBRAIC(:,50) = (ALGEBRAIC(:,34).*ALGEBRAIC(:,48).*ALGEBRAIC(:,44))./( ALGEBRAIC(:,42).*ALGEBRAIC(:,46));    % beta_Na2 [1/ms]
    ALGEBRAIC(:,52) = ALGEBRAIC(:,48)./1000.0;                                                                     % alfa_Na4 [1/ms]
    ALGEBRAIC(:,54) = ALGEBRAIC(:,44);                                                                             % beta_Na4 [1/ms]
    ALGEBRAIC(:,56) = ALGEBRAIC(:,48)./95000.0;                                                                    % alfa_Na5 [1/ms]
    ALGEBRAIC(:,58) = ALGEBRAIC(:,44)./50.0;                                                                       % beta_Na5 [1/ms]
    
    % Background Na+ current
    ALGEBRAIC(:,69) =  CONSTANTS(:,59).*(STATES(:,1)-ALGEBRAIC(:,66));                                             % INab [pA/pF]
    
    % ------------------------
    %   POTASSIUM DYNAMICS
    % ------------------------
    
    % Transient outward K+ current IKto_f
    ALGEBRAIC(:,70) =  ((CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,23));   % EK [mV]
    ALGEBRAIC(:,71) = CONSTANTS(:,60).*power(STATES(:,24),3.0).*STATES(:,25).*(STATES(:,1)-ALGEBRAIC(:,70));       % IKtof [pA/pF]
    ALGEBRAIC(:,6) = 0.18064.*exp(0.03577.*(STATES(:,1)+30.0));                                                    % alpha_a [1/ms]
    ALGEBRAIC(:,17) = 0.3956.*exp(-0.06237.*(STATES(:,1)+30.0));                                                   % beta_a [1/ms]
    ALGEBRAIC(:,7) = (0.000152.*exp(-(STATES(:,1)+13.5)./7.0))./(0.0067083.*exp(-(STATES(:,1)+33.5)./7.0)+1.0);    % alpha_i [1/ms]
    ALGEBRAIC(:,18) = (0.00095.*exp((STATES(:,1)+33.5)./7.0))./(0.051335.*exp((STATES(:,1)+33.5)./7.0)+1.0);       % beta_i [1/ms]
    
    % Transient outward K+ current IKto_s
    ALGEBRAIC(:,72) =  CONSTANTS(:,61).*STATES(:,26).*STATES(:,27).*(STATES(:,1)-ALGEBRAIC(:,70));                 % IKtos [pA/pF]
    ALGEBRAIC(:,8) = 1.0./(1.0+exp(-(STATES(:,1)+22.5)./7.7));                                                     % a_ss [-]
    ALGEBRAIC(:,9) = 1.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                                      % i_ss [-]
    ALGEBRAIC(:,19) = 0.493.*exp(-0.0629.*STATES(:,1))+2.058;                                                      % tau_tas [ms]
    ALGEBRAIC(:,20) = 270.0+1050.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                            % tau_tis [ms]
    
    % Time-indipendent K+ current
    ALGEBRAIC(:,73) = (((CONSTANTS(:,62).*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.0)).*(STATES(:,1)-...
        ALGEBRAIC(:,70)))./(1.0+exp(0.0896.*(STATES(:,1)-ALGEBRAIC(:,70))));                                       % IK1 [pA/pF]
    
    % Slow delayed rectifier K+ current
    ALGEBRAIC(:,74) =  CONSTANTS(:,63).*power(STATES(:,28),2.0).*(STATES(:,1)-ALGEBRAIC(:,70));                    % IKs [pA/pF]
    ALGEBRAIC(:,10) = (4.81333e-06.*(STATES(:,1)+26.5))./(1.0-exp(-0.128.*(STATES(:,1)+26.5)));                    % alpha_n [1/ms]
    ALGEBRAIC(:,21) = 9.53333e-05.*exp(-0.038.*(STATES(:,1)+26.5));                                                % beta_n [1/ms]
    
    % Ultrarapidly activating delayed rectifier K+ current
    ALGEBRAIC(:,75) =  CONSTANTS(:,64).*STATES(:,29).*STATES(:,30).*(STATES(:,1)-ALGEBRAIC(:,70));                 % IKur [pA/pF]
    ALGEBRAIC(:,22) = 0.493.*exp(-0.0629.*STATES(:,1))+2.058;                                                      % tau_aur [ms]
    ALGEBRAIC(:,23) = 1200.0-170.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                            % tau_iur [ms]
    
    % Noninactivating steady-state K+ current
    ALGEBRAIC(:,76) = CONSTANTS(:,65).*STATES(:,31).*STATES(:,32).*(STATES(:,1)-ALGEBRAIC(:,70));                  % IKss [pA/pF]
    ALGEBRAIC(:,24) = 39.3.*exp(-0.0862.*STATES(:,1))+13.17;                                                       % tau_Kss [ms]
    
    % Rapid delayed rectifier K+ current (mERG)  
    ALGEBRAIC(:,77) =  CONSTANTS(:,66).*STATES(:,33).*(STATES(:,1)-((CONSTANTS(:,10).*CONSTANTS(:,11))./...
        CONSTANTS(:,12)).*log((0.98.*CONSTANTS(:,7)+ 0.02.*CONSTANTS(:,8))./(0.98.*STATES(:,23)+...
        0.02.*STATES(:,14))));                                                                                     % IKr [pA/pF]
    ALGEBRAIC(:,11) = 1.0-(STATES(:,34)+STATES(:,35)+STATES(:,33)+STATES(:,36));                                   % mERG channel closed state - C_K0 [-] 
    ALGEBRAIC(:,25) = 0.022348.*exp(0.01176.*STATES(:,1));                                                         % alpha_a0 [1/ms]
    ALGEBRAIC(:,30) = 0.047002.*exp(-0.0631.*STATES(:,1));                                                         % beta_a0 [1/ms]
    ALGEBRAIC(:,12) = 0.013733.*exp(0.038198.*STATES(:,1));                                                        % alpha_a1 [1/ms]
    ALGEBRAIC(:,26) = 6.89e-05.*exp(-0.04178.*STATES(:,1));                                                        % beta_a1 [1/ms]
    ALGEBRAIC(:,31) = 0.090821.*exp(0.023391.*(STATES(:,1)+5.0));                                                  % alpha_i [1/ms]    
    ALGEBRAIC(:,35) = 0.006497.*exp(-0.03268.*(STATES(:,1)+5.0));                                                  % beta_i [1/ms]
        
    % Na+/K+ pump
    ALGEBRAIC(:,78) = 1.0./(1.0+0.124500.*exp((-0.1.*STATES(:,1).*CONSTANTS(:,12))./(CONSTANTS(:,10).*...
        CONSTANTS(:,11)))+0.0365.*CONSTANTS(:,76).*exp((-STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*...
        CONSTANTS(:,11))));                                                                                        % fNaK [-]
    ALGEBRAIC(:,79) = (((CONSTANTS(:,69).*ALGEBRAIC(:,78).*1.0)./(1.0+power(CONSTANTS(:,70)./STATES(:,14),1.5))).*...
        CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,71));                                                         % INaK [pA/pF]
    
    % Ca2+ activated Cl- current
    ALGEBRAIC(:,80) = 0.2./(1.0+exp(-(STATES(:,1)-46.7)./7.8));                                                    % O_ClCa [-]
    ALGEBRAIC(:,81) = ((CONSTANTS(:,72).*ALGEBRAIC(:,80).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,74))).*...
        (STATES(:,1)-CONSTANTS(:,73));                                                                             % IClCa [pA/pF]

    case 1
    
    ALGEBRAIC(:,1) = floor(VOI./CONSTANTS(:,14)).*CONSTANTS(:,14);                                                 % past [ms]
    if (VOI-ALGEBRAIC(:,1)>=CONSTANTS(:,13) & VOI-ALGEBRAIC(:,1)<=CONSTANTS(:,13)+CONSTANTS(:,15))
        ALGEBRAIC(:,13) = CONSTANTS(:,16);                                                                         % Istim [pA/pF]
    else
        ALGEBRAIC(:,13) = 0.0;                                                                                     % Istim [pA/pF]
    end

    ALGEBRAIC(:,82) = ((VOI-ALGEBRAIC(:,1))-CONSTANTS(:,77))./((CONSTANTS(:,78).*2)-CONSTANTS(:,77));              % chi_1 [-]
    ALGEBRAIC(:,83) = 1.22+(CONSTANTS(:,82)-1.0)*(power(ALGEBRAIC(:,82),3.0)).*(CONSTANTS(:,83)-CONSTANTS(:,84)*...
        ALGEBRAIC(:,82)+CONSTANTS(:,85)*power(ALGEBRAIC(:,82),2.0));                                               % Cm_1 [pF]    
    ALGEBRAIC(:,84) = CONSTANTS(:,86).*(power((VOI-ALGEBRAIC(:,1)),3.0))+CONSTANTS(:,87).*...
        (power((VOI-ALGEBRAIC(:,1)),2.0))+CONSTANTS(:,88).*(VOI-ALGEBRAIC(:,1))+CONSTANTS(:,89);                   % Cm_2 [pF]
    ALGEBRAIC(:,85) = CONSTANTS(:,90).*(power((VOI-ALGEBRAIC(:,1)),2.0))+CONSTANTS(:,91).*...
        (VOI-ALGEBRAIC(:,1))+CONSTANTS(:,92);                                                                      % Cm_3 [pF]
    ALGEBRAIC(:,86) = ((VOI-ALGEBRAIC(:,1))-CONSTANTS(:,80))./(CONSTANTS(:,81)-CONSTANTS(:,80));                   % chi_2 [-]
    ALGEBRAIC(:,87) = CONSTANTS(:,82)+0.22+(1-CONSTANTS(:,82))*(power(ALGEBRAIC(:,86),3.0)).*(CONSTANTS(:,83)-...
        CONSTANTS(:,84)*ALGEBRAIC(:,86)+CONSTANTS(:,85)*power(ALGEBRAIC(:,86),2.0));                               % Cm_4 [pF]
        
    ALGEBRAIC((VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,78),88) = ALGEBRAIC((VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,78),83);
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,78) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,79)),88) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,78) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,79)),84);           
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,79) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,80)),88) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,79) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,80)),85);
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,80) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,81)),88) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,80) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,81)),87);
    ALGEBRAIC((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,81),88) = 1.22;                                                     % Cm [pF]                       
    
    ALGEBRAIC(:,89) = (CONSTANTS(:,82)-1.0)*power(ALGEBRAIC(:,82),2.0).*(5.0*CONSTANTS(:,85).*...
        power(ALGEBRAIC(:,82),2.0)-4.0*CONSTANTS(:,84)*ALGEBRAIC(:,82)+3.0*CONSTANTS(:,83))/...
        ((CONSTANTS(:,78).*2.0)-CONSTANTS(:,77));                                                                  % d/dt Cm_1 [-]
    ALGEBRAIC(:,90) = 3.0*CONSTANTS(:,86).*(power((VOI-ALGEBRAIC(:,1)),2.0))+2.0*CONSTANTS(:,87).*...
        (VOI-ALGEBRAIC(:,1))+CONSTANTS(:,88);                                                                      % d/dt Cm_2 [-]
    ALGEBRAIC(:,91) = 2.0*CONSTANTS(:,90).*(VOI-ALGEBRAIC(:,1))+CONSTANTS(:,91);                                   % d/dt Cm_3 [-]
    ALGEBRAIC(:,92) = (1.0-CONSTANTS(:,82))*power(ALGEBRAIC(:,86),2.0).*(5.0*CONSTANTS(:,85).*...
        power(ALGEBRAIC(:,86),2.0)-4.0*CONSTANTS(:,84)*ALGEBRAIC(:,86)+3.0*CONSTANTS(:,83))/...
        (CONSTANTS(:,81)-CONSTANTS(:,80));                                                                         % d/dt Cm_4 [-]

    ALGEBRAIC((VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,78),93) = ALGEBRAIC((VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,78),89);
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,78) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,79)),93) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,78) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,79)),90);
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,79) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,80)),93) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,79) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,80)),91); 
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,80) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,81)),93) = ALGEBRAIC(((VOI-...
        ALGEBRAIC(:,1))>CONSTANTS(:,80) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,81)),92);
    ALGEBRAIC((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,81),93) = 0.0;                                                      % d/dt Cm [-]
    
    % ------------------------
    %   CALCIUM DYNAMICS
    % ------------------------
    
    % Calcium concentration
    ALGEBRAIC(:,27) = power(1.0+(CONSTANTS(:,17).*CONSTANTS(:,19))./power(CONSTANTS(:,19)+STATES(:,2),2.0),-1.0);  % Bi [-]
    ALGEBRAIC(:,32) = power(1.0+(CONSTANTS(:,17).*CONSTANTS(:,19))./power(CONSTANTS(:,19)+STATES(:,3),2.0),-1.0);  % Bss [-]
    ALGEBRAIC(:,36) = power(1.0+(CONSTANTS(:,18).*CONSTANTS(:,20))./power(CONSTANTS(:,20)+STATES(:,4),2.0),-1.0);  % BJSR [-]
    
    % Calcium fluxes
    ALGEBRAIC(:,39) = CONSTANTS(:,21).*(STATES(:,8)+STATES(:,9)).*(STATES(:,4)-STATES(:,3)).*STATES(:,6);          % Jrel [uM/ms]
    ALGEBRAIC(:,41) = (STATES(:,5)-STATES(:,4))./CONSTANTS(:,22);                                                  % Jtr [uM/ms]
    ALGEBRAIC(:,43) = (STATES(:,3)-STATES(:,2))./CONSTANTS(:,24);                                                  % Jxfer [uM/ms]
    ALGEBRAIC(:,45) = CONSTANTS(:,23).*(STATES(:,5)-STATES(:,2));                                                  % Jleak [uM/ms]
    ALGEBRAIC(:,47) = (0.05.*(1.0-STATES(:,7)).*1.0)./(1.0+0.7./STATES(:,3));                                      % CaMKb [-]
    ALGEBRAIC(:,49) = ALGEBRAIC(:,47)+STATES(:,7);                                                                 % CaMKa [-]
    ALGEBRAIC(:,51) = ((3.1512.*power(ALGEBRAIC(:,49),2.2062))./(power(0.1588,2.2062)+...
        power(ALGEBRAIC(:,49),2.2062))+1.0).*CONSTANTS(:,29);                                                      % vmup [uM/ms]
    ALGEBRAIC(:,53) = (ALGEBRAIC(:,51).*power(STATES(:,2), 2.0))./(power(CONSTANTS(:,27),2.0)+...
        power(STATES(:,2),2.0));                                                                                   % Jserca [uM/ms]
    
    % Ryanodine receptors
    ALGEBRAIC(:,2) = 1.0-(STATES(:,10)+STATES(:,8)+STATES(:,9));                                                   % Function of RyR channels in state Pc1 - PC1 [-]
        
    % Calcium L-type current
    ALGEBRAIC(:,55) = (CONSTANTS(:,12).*STATES(:,1))./( CONSTANTS(:,10).*CONSTANTS(:,11));                         % FVRT [-]
    ALGEBRAIC(:,57) = 2.0.*ALGEBRAIC(:,55);                                                                        % FVRT_Ca [-]
    ALGEBRAIC(:,4) = exp((STATES(:,1)-CONSTANTS(:,41))./CONSTANTS(:,42));                                          % exp_VL [-]
    ALGEBRAIC(:,15) = ALGEBRAIC(:,4)./(CONSTANTS(:,43).*(ALGEBRAIC(:,4)+1.0));                                     % alpha_p [1/ms]
    ALGEBRAIC(:,28) = (ALGEBRAIC(:,4)+CONSTANTS(:,45))./(CONSTANTS(:,47).*CONSTANTS(:,48).*(ALGEBRAIC(:,4)+1.0));  % epsilon_p [1/uMms]
    ALGEBRAIC(:,33) = (CONSTANTS(:,46).*(ALGEBRAIC(:,4)+CONSTANTS(:,45)))./...
        (CONSTANTS(:,47).*(CONSTANTS(:,46).*ALGEBRAIC(:,4)+CONSTANTS(:,45)));                                      % epsilon_m [1/ms]
    ALGEBRAIC(:,3) = 1.0./(1.0+exp((STATES(:,1)+16.6577)./CONSTANTS(:,49)))+0.1./...
        (1.0+exp((-STATES(:,1)+40.0)./6.0));                                                                       % y_inf [-]
    ALGEBRAIC(:,14) = 20.0+600.0./(1.0+exp((STATES(:,1)+30.0)./9.6));                                              % tau_y [ms]
    ALGEBRAIC(:,37) = (1.0-STATES(:,11))-STATES(:,12);                                                             % C [-]
    if (abs(ALGEBRAIC(:,57))>1.0e-05)
        ALGEBRAIC(:,59) = ((((-CONSTANTS(:,40).*2.0.*CONSTANTS(:,5).*...
            CONSTANTS(:,12))./(CONSTANTS(:,6).*ALGEBRAIC(:,88))).*STATES(:,11).*STATES(:,13).*ALGEBRAIC(:,57))./...
            (1.0-exp(-ALGEBRAIC(:,57)))).*(CONSTANTS(:,9).*exp(-ALGEBRAIC(:,57))-STATES(:,3));                     % ICaL [pA/pF]
    else
        ALGEBRAIC(:,59) = ((((-CONSTANTS(:,40).*2.0.*CONSTANTS(:,5).*CONSTANTS(:,12))./(CONSTANTS(:,6).*...
                ALGEBRAIC(:,88))).*STATES(:,11).*STATES(:,13).*1.0e-05)./(1.0-exp(-1.0e-05))).*(CONSTANTS(:,9).*...
                exp(-1.0e-05)-STATES(:,3));                                                                        % ICaL [pA/pF]
    end
    
    % Calcium pump current
    ALGEBRAIC(:,60) = (CONSTANTS(:,50).*power(STATES(:,2),2.0))./(power(CONSTANTS(:,51),2.0)+...
        power(STATES(:,2),2.0));                                                                                   % IpCa [pA/pF]
    ALGEBRAIC(:,62) = (-ALGEBRAIC(:,60).*CONSTANTS(:,6).*ALGEBRAIC(:,86))./...
        (2.0.*CONSTANTS(:,2).*CONSTANTS(:,12));                                                                    % JpCa [uM/ms]
      
    % Sodium/calcium exchange current
    ALGEBRAIC(:,61) = ((((((CONSTANTS(:,52).*1.0)./(power(CONSTANTS(:,53),3.0)+power(CONSTANTS(:,8),3.0))).*1.0)./...
        (CONSTANTS(:,54)+CONSTANTS(:,9))).*1.0)./(1.0+CONSTANTS(:,55).*exp(((CONSTANTS(:,56)-1.0).*STATES(:,1).*...
        CONSTANTS(:,12))./( CONSTANTS(:,10).*CONSTANTS(:,11))))).*(exp((CONSTANTS(:,56).*STATES(:,1).*CONSTANTS(:,12))./...
        (CONSTANTS(:,10).*CONSTANTS(:,11))).*power(STATES(:,14),3.0).*CONSTANTS(:,9)-exp(((CONSTANTS(:,56)-1.0).*...
        STATES(:,1).*CONSTANTS(:,12))./(CONSTANTS(:,10).*CONSTANTS(:,11))).*power(CONSTANTS(:,8),3.0).*...
        STATES(:,2));                                                                                              % INaCa [pA/pF]
    ALGEBRAIC(:,64) = (ALGEBRAIC(:,61).*CONSTANTS(:,6).*ALGEBRAIC(:,86))./(CONSTANTS(:,2).*CONSTANTS(:,12));       % Jncx [uM/ms]
    
    % Calcium background current
    ALGEBRAIC(:,63) = ((CONSTANTS(:,10).*CONSTANTS(:,11))./(2.0.*CONSTANTS(:,12))).*log(CONSTANTS(:,9)./...
        STATES(:,2));                                                                                              % ECaN [mV]
    ALGEBRAIC(:,65) = CONSTANTS(:,57).*(STATES(:,1)-ALGEBRAIC(:,63));                                              % ICab [pA/pF]
    ALGEBRAIC(:,67) = (-ALGEBRAIC(:,65).*CONSTANTS(:,6).*ALGEBRAIC(:,86))./(2.0.*CONSTANTS(:,2).*CONSTANTS(:,12)); % JCab [uM/ms]
    
    % ------------------------
    %   SODIUM DYNAMICS
    % ------------------------
    
    % Fast Na+ current
    ALGEBRAIC(:,66) = ((CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log((0.9.*CONSTANTS(:,8)+0.1.*...
        CONSTANTS(:,7))./(0.9.*STATES(:,14)+0.1.*STATES(:,23)));                                                   % ENa [mV]
    ALGEBRAIC(:,68) = CONSTANTS(:,58).*STATES(:,15).*(STATES(:,1)-ALGEBRAIC(:,66));                                % INa [pA/pF]
    ALGEBRAIC(:,5) = 1.0-(STATES(:,15)+STATES(:,16)+STATES(:,17)+STATES(:,20)+STATES(:,18)+STATES(:,19)+...
        STATES(:,21)+STATES(:,22));                                                                                % Closed state of fast Na+ channel - C_Na3 [-] 
    ALGEBRAIC(:,16) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./17.0)+0.2.*exp(-(STATES(:,1)+2.5)./150.0));          % alpha_Na11 [1/ms]
    ALGEBRAIC(:,29) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./15.0)+0.23.*exp(-(STATES(:,1)+2.5)./150.0));         % alpha_Na12 [1/ms]
    ALGEBRAIC(:,34) = 3.802./(0.1027.*exp(-(STATES(:,1)+2.5)./12.0)+0.25.*exp(-(STATES(:,1)+2.5)./150.0));         % alpha_Na13 [1/ms]
    ALGEBRAIC(:,38) = 0.1917.*exp(-(STATES(:,1)+2.5)./20.3);                                                       % beta_Na11 [1/ms]
    ALGEBRAIC(:,40) = 0.2.*exp(-(STATES(:,1)-2.5)./20.3);                                                          % beta_Na12 [1/ms]
    ALGEBRAIC(:,42) = 0.22.*exp(-(STATES(:,1)-7.5)./20.3);                                                         % beta_Na13 [1/ms]
    ALGEBRAIC(:,44) = 7.0e-07.*exp(-(STATES(:,1)+7.0)./7.7);                                                       % alfa_Na3 [1/ms]
    ALGEBRAIC(:,46) = 0.0084+2.0e-05.*(STATES(:,1)+7.0);                                                           % beta_Na3 [1/ms]
    ALGEBRAIC(:,48) = 1.0./(0.188495.*exp(-(STATES(:,1)+7.0)./16.6)+0.393956);                                     % alfa_Na2 [1/ms]
    ALGEBRAIC(:,50) = (ALGEBRAIC(:,34).*ALGEBRAIC(:,48).*ALGEBRAIC(:,44))./( ALGEBRAIC(:,42).*ALGEBRAIC(:,46));    % beta_Na2 [1/ms]
    ALGEBRAIC(:,52) = ALGEBRAIC(:,48)./1000.0;                                                                     % alfa_Na4 [1/ms]
    ALGEBRAIC(:,54) = ALGEBRAIC(:,44);                                                                             % beta_Na4 [1/ms]
    ALGEBRAIC(:,56) = ALGEBRAIC(:,48)./95000.0;                                                                    % alfa_Na5 [1/ms]
    ALGEBRAIC(:,58) = ALGEBRAIC(:,44)./50.0;                                                                       % beta_Na5 [1/ms]
    
    % Background Na+ current
    ALGEBRAIC(:,69) =  CONSTANTS(:,59).*(STATES(:,1)-ALGEBRAIC(:,66));                                             % INab [pA/pF]
    
    % ------------------------
    %   POTASSIUM DYNAMICS
    % ------------------------
    
    % Transient outward K+ current IKto_f
    ALGEBRAIC(:,70) =  ((CONSTANTS(:,10).*CONSTANTS(:,11))./CONSTANTS(:,12)).*log(CONSTANTS(:,7)./STATES(:,23));   % EK [mV]
    ALGEBRAIC(:,71) = CONSTANTS(:,60).*power(STATES(:,24),3.0).*STATES(:,25).*(STATES(:,1)-ALGEBRAIC(:,70));       % IKtof [pA/pF]
    ALGEBRAIC(:,6) = 0.18064.*exp(0.03577.*(STATES(:,1)+30.0));                                                    % alpha_a [1/ms]
    ALGEBRAIC(:,17) = 0.3956.*exp(-0.06237.*(STATES(:,1)+30.0));                                                   % beta_a [1/ms]
    ALGEBRAIC(:,7) = (0.000152.*exp(-(STATES(:,1)+13.5)./7.0))./(0.0067083.*exp(-(STATES(:,1)+33.5)./7.0)+1.0);    % alpha_i [1/ms]
    ALGEBRAIC(:,18) = (0.00095.*exp((STATES(:,1)+33.5)./7.0))./(0.051335.*exp((STATES(:,1)+33.5)./7.0)+1.0);       % beta_i [1/ms]
    
    % Transient outward K+ current IKto_s
    ALGEBRAIC(:,72) =  CONSTANTS(:,61).*STATES(:,26).*STATES(:,27).*(STATES(:,1)-ALGEBRAIC(:,70));                 % IKtos [pA/pF]
    ALGEBRAIC(:,8) = 1.0./(1.0+exp(-(STATES(:,1)+22.5)./7.7));                                                     % a_ss [-]
    ALGEBRAIC(:,9) = 1.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                                      % i_ss [-]
    ALGEBRAIC(:,19) = 0.493.*exp(-0.0629.*STATES(:,1))+2.058;                                                      % tau_tas [ms]
    ALGEBRAIC(:,20) = 270.0+1050.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                            % tau_tis [ms]
    
    % Time-indipendent K+ current
    ALGEBRAIC(:,73) = (((CONSTANTS(:,62).*CONSTANTS(:,7))./(CONSTANTS(:,7)+210.0)).*(STATES(:,1)-...
        ALGEBRAIC(:,70)))./(1.0+exp(0.0896.*(STATES(:,1)-ALGEBRAIC(:,70))));                                       % IK1 [pA/pF]
    
    % Slow delayed rectifier K+ current
    ALGEBRAIC(:,74) = CONSTANTS(:,63).*power(STATES(:,28),2.0).*(STATES(:,1)-ALGEBRAIC(:,70));                     % IKs [pA/pF]
    ALGEBRAIC(:,10) = (4.81333e-06.*(STATES(:,1)+26.5))./(1.0-exp(-0.128.*(STATES(:,1)+26.5)));                    % alpha_n [1/ms]
    ALGEBRAIC(:,21) = 9.53333e-05.*exp(-0.038.*(STATES(:,1)+26.5));                                                % beta_n [1/ms]
    
    % Ultrarapidly activating delayed rectifier K+ current
    ALGEBRAIC(:,75) =  CONSTANTS(:,64).*STATES(:,29).*STATES(:,30).*(STATES(:,1)-ALGEBRAIC(:,70));                 % IKur [pA/pF]
    ALGEBRAIC(:,22) = 0.493.*exp(-0.0629.*STATES(:,1))+2.058;                                                      % tau_aur [ms]
    ALGEBRAIC(:,23) = 1200.0-170.0./(1.0+exp((STATES(:,1)+45.2)./5.7));                                            % tau_iur [ms]
    
    % Noninactivating steady-state K+ current
    ALGEBRAIC(:,76) = CONSTANTS(:,65).*STATES(:,31).*STATES(:,32).*(STATES(:,1)-ALGEBRAIC(:,70));                  % IKss [pA/pF]
    ALGEBRAIC(:,24) = 39.3.*exp(-0.0862.*STATES(:,1))+13.17;                                                       % tau_Kss [ms]
    
    % Rapid delayed rectifier K+ current (mERG)  
    ALGEBRAIC(:,77) =  CONSTANTS(:,66).*STATES(:,33).*(STATES(:,1)-((CONSTANTS(:,10).*CONSTANTS(:,11))./...
        CONSTANTS(:,12)).*log((0.98.*CONSTANTS(:,7)+ 0.02.*CONSTANTS(:,8))./(0.98.*STATES(:,23)+...
        0.02.*STATES(:,14))));                                                                                     % IKr [pA/pF]
    ALGEBRAIC(:,11) = 1.0-(STATES(:,34)+STATES(:,35)+STATES(:,33)+STATES(:,36));                                   % mERG channel closed state - C_K0 [-] 
    ALGEBRAIC(:,25) = 0.022348.*exp(0.01176.*STATES(:,1));                                                         % alpha_a0 [1/ms]
    ALGEBRAIC(:,30) = 0.047002.*exp(-0.0631.*STATES(:,1));                                                         % beta_a0 [1/ms]
    ALGEBRAIC(:,12) = 0.013733.*exp(0.038198.*STATES(:,1));                                                        % alpha_a1 [1/ms]
    ALGEBRAIC(:,26) = 6.89e-05.*exp(-0.04178.*STATES(:,1));                                                        % beta_a1 [1/ms]
    ALGEBRAIC(:,31) = 0.090821.*exp(0.023391.*(STATES(:,1)+5.0));                                                  % alpha_i [1/ms]    
    ALGEBRAIC(:,35) = 0.006497.*exp(-0.03268.*(STATES(:,1)+5.0));                                                  % beta_i [1/ms]
        
    % Na+/K+ pump
    ALGEBRAIC(:,78) = 1.0./(1.0+0.124500.*exp((-0.1.*STATES(:,1).*CONSTANTS(:,12))./(CONSTANTS(:,10).*...
        CONSTANTS(:,11)))+0.0365.*CONSTANTS(:,76).*exp((-STATES(:,1).*CONSTANTS(:,12))./( CONSTANTS(:,10).*...
        CONSTANTS(:,11))));                                                                                        % fNaK [-]
    ALGEBRAIC(:,79) = (((CONSTANTS(:,69).*ALGEBRAIC(:,78).*1.0)./(1.0+power(CONSTANTS(:,70)./STATES(:,14),1.5))).*...
        CONSTANTS(:,7))./(CONSTANTS(:,7)+CONSTANTS(:,71));                                                         % INaK [pA/pF]
    
    % Ca2+ activated Cl- current
    ALGEBRAIC(:,80) = 0.2./(1.0+exp(-(STATES(:,1)-46.7)./7.8));                                                    % O_ClCa [-]
    ALGEBRAIC(:,81) = ((CONSTANTS(:,72).*ALGEBRAIC(:,80).*STATES(:,2))./(STATES(:,2)+CONSTANTS(:,74))).*...
        (STATES(:,1)-CONSTANTS(:,73));                                                                             % IClCa [pA/pF]

    % ------------------------
    %   SACs
    % ------------------------
  
    ALGEBRAIC(:,94) = (CONSTANTS(:,93)-1)*(((VOI-ALGEBRAIC(:,1))-CONSTANTS(:,103))/(CONSTANTS(:,104)-...
        CONSTANTS(:,103)));                                                                                        % x1
    ALGEBRAIC(:,95) = ((CONSTANTS(:,94)-CONSTANTS(:,93))*(((VOI-ALGEBRAIC(:,1))-CONSTANTS(:,104))/...
        (CONSTANTS(:,105)-CONSTANTS(:,104))))+(CONSTANTS(:,93)-1);                                                 % x2
    ALGEBRAIC(:,96) = ((1-CONSTANTS(:,94))*(((VOI-ALGEBRAIC(:,1))-CONSTANTS(:,105))/(CONSTANTS(:,106)-...
        CONSTANTS(:,105))))+(CONSTANTS(:,94)-1);                                                                   % x3

    ALGEBRAIC((VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,103),97) = 0.0;
    ALGEBRAIC((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,103) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,104),97) =...
        ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,103) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,104)),94);  
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,104) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,105)),97) =...
        ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,104) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,105)),95);
    ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,105) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,106)),97) =...
        ALGEBRAIC(((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,105) & (VOI-ALGEBRAIC(:,1))<=CONSTANTS(:,106)),96);
    ALGEBRAIC((VOI-ALGEBRAIC(:,1))>CONSTANTS(:,106),97) = 0.0;                                                     % x

    ALGEBRAIC(:,98) = CONSTANTS(:,95)*ALGEBRAIC(:,97);                                                             % Gamma_SL_ns
    ALGEBRAIC(:,99) = (CONSTANTS(:,96)*CONSTANTS(:,97)*ALGEBRAIC(:,98).*(STATES(:,1)-ALGEBRAIC(:,66)))./...
        ALGEBRAIC(:,88);                                                                                           % I_ns_Na
    ALGEBRAIC(:,100) = (CONSTANTS(:,97)*ALGEBRAIC(:,98).*(STATES(:,1)-ALGEBRAIC(:,70)))./ALGEBRAIC(:,88);          % I_ns_K
    ALGEBRAIC(:,101) = ALGEBRAIC(:,99)+ALGEBRAIC(:,100);                                                           % I_ns

    ALGEBRAIC(:,102) = CONSTANTS(:,98)*ALGEBRAIC(:,97)+0.7;                                                        % Gamma_SL_Ko
    ALGEBRAIC(:,103) = (CONSTANTS(:,99)*(ALGEBRAIC(:,102)./(1+exp(-(10+STATES(:,1))/45.0))).*(STATES(:,1)-...
        ALGEBRAIC(:,70)))./ALGEBRAIC(:,88);                                                                        % I_Ko

    ALGEBRAIC(:,104) = CONSTANTS(:,100)*ALGEBRAIC(:,97);                                                           % P_CaP

    ALGEBRAIC(:,105) = 0.5*(STATES(:,14)+STATES(:,23)+4*STATES(:,2))/1000000;                                      % Ii
    ALGEBRAIC(:,106) = 0.5*(CONSTANTS(:,8)+CONSTANTS(:,7)+4*CONSTANTS(:,9))/1000000;                               % Io
    ALGEBRAIC(:,107) = exp(-CONSTANTS(:,102)*4*(sqrt(ALGEBRAIC(:,105))./(1+sqrt(ALGEBRAIC(:,105)))-...
        0.3*ALGEBRAIC(:,105)));                                                                                    % Gamma_Cai
    ALGEBRAIC(:,108) = exp(-CONSTANTS(:,102)*4*(sqrt(ALGEBRAIC(:,106))./(1+sqrt(ALGEBRAIC(:,106)))-...
        0.3*ALGEBRAIC(:,106)));                                                                                    % Gamma_Cao
    ALGEBRAIC(:,109) = 4.0*ALGEBRAIC(:,104).*((STATES(:,1)*(CONSTANTS(:,12)^2))./(CONSTANTS(:,10)*CONSTANTS(:,11))).*...
        (ALGEBRAIC(:,107).*STATES(:,2).*exp(2.0*((STATES(:,1)*CONSTANTS(:,12))./(CONSTANTS(:,10)*CONSTANTS(:,11))))-...
        ALGEBRAIC(:,108)*CONSTANTS(:,9))./(exp(2.0*((STATES(:,1)*CONSTANTS(:,12))./(CONSTANTS(:,10)*...
        CONSTANTS(:,11))))-1.0);                                                                                   % I_Cap

end
end