%% Coronary flow modelling - Model-2.1 -  Vessel radius constant - Blunt Stenosis present - No branches
% author: Muneeza Azmat
% date : 3/16/2020
% Using navier-stokes equation for incompressible newtonian fluid in a
% circular pipe. We assume the flow is fully developed, steady and laminar with local seperation.
% using semi empirical formulas for pressure drop across a stenosis
% all units are in SI
% Results valid for % stenosis > 75 
% Input : (Lenght of stenosis (m) , diameter of stenosis (m) , Average velocity in m/s)

function [FFR]=LaminarSeperation(L_1,D_1,V_avg,D_0)
    
% Model parameters
rho = 1060;                         % kg/m3
mu = 2.78e-3;                       % Pa.s (kg/m.s)
Pa = 119.7*133.322;                 % Static arterial pressure (Pa)
%D_0 = 4.6e-3;                      % unobstructed dia (m)
R_0 = D_0/2;                        % m
Qt = pi*(R_0)^2*(V_avg);            % m3/s
mdot = 1000*Qt*rho;                 % g/s

% Parametrize vessel-seg data by s
ds = 0.1*1e-3;
s  = [0:ds:140*1e-3];               % starts from 0 

%% -------------------------- Vessel Geometry -----------------------------%
R_start = R_0;
taper_ratio = 1;                     % for constant radius set taper_ratio=1
R_end = taper_ratio*R_start;
Rs = R_start+((R_end-R_start)/(s(end)-s(1)))*(s-s(1));
A_0   = (pi/4)*(D_0)^2;              % unobstructed crosssection area
Re_0  = (rho*V_avg*D_0)/mu;          % unobstructed reynolds number

%% ------------------------- Stenosis Geometry ----------------------------%
st_max_id = round(7e-2/ds);                          % location of max stenosis (7cm)
num_st = round( ((L_1)/(2*ds)) );
st_id = st_max_id + [-num_st : 1 : num_st];           % indices of stenosed elements

R_1 = D_1/2;
A_1 = pi*(R_1)^2;
Re_1 = (rho*(Qt/A_1)*D_1)/mu;                        % stenosed reynolds number

% change radius of stenosed elements
Rs(st_id) = R_1;

%% ------------------------- Calculate pressure ----------------------------%
% % Varying lumen area
% L_a    = L_1;
% A_ratio = 0.75*(A_0/A_1)+0.25;

% Blunt plug
A_ratio = A_0/A_1;
L_a    = 0.83*L_1 + 1.64*D_1;

K_t   = 1.52;
K_nu  = 32*(L_a/D_0)*(A_ratio)^2;   

DP_st = rho*(V_avg^2)*( (K_nu/Re_0) + (K_t/2)*(A_ratio -1 )^2 );

Ps    = zeros(1,length(s));
Ps(1) = Pa;
% before stenosis
for i = 2:st_id(1)-1
    Ps(i)= Pa - norm(s(i)-s(1))*(Qt*8*mu./(pi*Rs(i)^4));
end
% at stenosis
P1 = Ps(st_id(1)-1);
for i = st_id(1):st_id(end)
    x = s(i) - s(st_id(1));                 %horizontal distance from start of stenosis
    DP_st_x = (DP_st/L_1)*x;
    Ps(i) = (P1 - DP_st_x) + (rho*V_avg^2/2)*(1 - ((Rs(1)^4)/(Rs(i)^4)) );
end
% after stenosis
P2 = Ps(st_id(1)-1) - DP_st;
for i= st_id(end)+1:length(s)
    Ps(i)= P2 - norm(s(i)-s(st_id(end)))*(Qt*8*mu./(pi*Rs(i)^4));
end

%% ------------------------- Calculate FFR ----------------------------%
L_D = L_1/D_0;
if L_D < 1.5
    Pa_loc = 6e-2;                     % FFR proximal location
    Pd_loc = 8e-2;                     % FFR distal location
else
    Pa_loc = 4e-2;
    Pd_loc = 10e-2;
end
Pa_idx = round(Pa_loc/ds);
Pd_idx = round(Pd_loc/ds);
FFR = Ps(Pd_idx +2)/Ps(Pa_idx-2);

%% ------------------------- Plots & Figures ----------------------------%
% if FFR <0
%     fig = figure('visible', 'off');
%     subplot(2,1,1)
%     plot(s,Rs,'k', 'linewidth', 1), hold on, plot(s,-Rs,'k', 'linewidth', 1)
%     ylabel('Rs - vessel radius (m)')
%     ylim([-10*R_0  10*R_0])
%     xlim([s(1) s(end)])
%     grid on
% 
%     subplot(2,1,2)
%     title('Pressure (P)')
%     yyaxis right
%     plot(s,(0.00750062)*Ps, 'linewidth', 1.5), hold on , scatter(s(Pa_idx),(0.00750062)*Ps(Pa_idx),'k','filled') , scatter(s(Pd_idx),(0.00750062)*Ps(Pd_idx),'k','filled')
%     ylabel('mmHg')
%     yyaxis left
%     plot(s,Ps,'--')
%     ylabel('Pa')
%     xlabel('s (m)')
%     grid on
%     annotation(fig,'textbox',...
%         [0.734415233415235 0.383597882714843 0.16461916817787 0.0590828932992575],...
%         'String',{strcat('FFR = ',num2str(FFR))},...
%         'FitBoxToText','on');
%     ax.YAxis(1).Exponent = 0;
%     saveas(fig,strcat(num2str(abs(FFR)),'.jpg'))
%     close(fig)
% end

end





