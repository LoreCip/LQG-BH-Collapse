%----------------------------------------------------------
% Matlab Code for effective LQG collapse of dust matter
% V. Husain, J.G. Kelly, R. Santacruz and E. Wilson-Ewing
%----------------------------------------------------------
% Description:
%   This code solves the differential conservation equation
%   stated in the paper by making use of the Godunov Scheme.
%   Initial conditions are set up such we have an infalling
%   gaussian pulse of dust matter. 
%   The evolution is such that the dust field bounces and  
%   and the matter escapes to infinity.
%----------------------------------------------------------

% LQG B equation:
%   dot{B} + 1/(2gammaDelta)( x^3 sin^2(sqrt{Delta}B/x^2) )' = 0

% In dimensioness variables:
%   B-> B/sqrt(Delta), (x,t) -> (x/sqrt(Delta),t/sqrt(Delta))
%   dot{B} + 1/(2 gamma)( x^3 sin^2(B/x^2) )' = 0
%   gamma = 1
 

syms x y m;
 

% Initial parameter for gaussian pulse
sig = 1/2;

%% Gaussian
%Rho0(x,m) = exp(-(x-2.5*m)^2/sig^2);

%% Tanh
radius0=15;
m=5;
%Rho0(x) = 3*m/(4*pi*radius0^3)*(1-tanh(x-radius0))/2;
%Rho0(x) = 3*m/(4*pi*radius0^3)*(pi/2-atan(x-radius0))/2;
Rho0(x) = 3*m/(4*pi*radius0^3)*(1-heaviside(x-radius0));

%% Double peak  
%Rho0(x) = exp(-2*(x-1.5*m)^2)/5 + exp(-4*(x-3.5*m)^2)/10;
   
%% Initial Conditions for Mass, B and Density
%Mfn(x,m) =  m*int(y^2*Rho0(y,m),[0,x])/int(y^2*Rho0(y,m),[0,Inf]); 
%Bfn(x,m) =  -x^2*acos(1-4*Mfn(x,m)/x^3 )/2;

Mfn(x) =  4*pi*int(y^2*Rho0(y),[0,x]); 
Bfn(x) =  -x^2*acos(1-4*Mfn(x)/x^3 )/2;
dx=0.1;

%% Set mass values in array
%  mass = massA(kk);
  mass = 5;  
  % radial axis
  xm = 20.00;
 X = 0.0:dx:xm;
  XN = numel(X);



% Spatial Lattice Spacing


tic;
 
 
%%  Begin mass runs: 
%   We set up the previously stated initial conditions in arrays and define
%   a limit of time for the computation after horizons disappear. Then we
%   start our computation to solve the equation of motion. Once we evolve
%   the solution past the Time limit, the computation stops and we set up a
%   different mass. In this way we get several points of the lifetime vs.
%   mass of the collapse.

%   In what follows, we have set the loop to do just one mass, but it can be changed if
%   needed.

fprintf('\nmass      lifetime\n');

for  kk = 1:1  
    
  %mass = massA(kk);
  %mass = 7;  
  % radial axis
  xm = 50.00;
  X = 0.0:dx:xm;
  XN = numel(X);

  % set scale for the vertical axis in the output plots
  ym = 2.0;

  % set mass run time t_final 
  T_final = 8.5*mass^2 + 40*mass;

  % initial data vectors
  
  B=zeros(1,XN);
  zeroline=zeros(1,XN); % This is the y = 0 line, to check the formation of apparent horizons.
 
  BB = matlabFunction(Bfn(x));

  B=BB(X);
  B(1)=B(2);
   
  %% Video file

  % v = VideoWriter('LTB_Movie.avi');   
  % open (v);
   
  a = [];

  %set(gcf,'WindowState','fullscreen');


  %% Evolution loop for each mass

  counter = 0;

  bh = false;
  t=0;
  ti=0;
  tf=0;
  y_max = -1;

  % This loops keeps repeating until we reached the max computation time
  while t < T_final
    
    % Set dt using velocity of characteristics
    vel = X.*sin(B./X.^2).*cos(B./X.^2);
    v_min = min(vel);
    v_max = max(vel);
    v_abs = max(-v_min, v_max);

    % CFL condition
    dt = 0.45*dx/abs(v_abs);
    if dt > dx / 20
       dt = dx / 20; % largest timestep allowed
    end
    t = t + dt;
    
    Bt = B;
    B(1)=0; % inner boundary condition, since B = x*b, at x=0 we have B=0
    
    % flux from B(1) into B(2) is 0
    B(2) = Bt(2) - dt/dx * flux(Bt(2)/X(2)^2,Bt(3)/X(3)^2,(X(2)+X(3))/2);
    
    for j = 3:XN-1
        B(j) = Bt(j) - dt/dx * (flux(Bt(j)/X(j)^2, Bt(j+1)/X(j+1)^2, (X(j)+X(j+1))/2) - flux(Bt(j-1)/X(j-1)^2,Bt(j)/X(j)^2,(X(j-1)+X(j))/2));         
    end
    
    % outer boundary condition: no dynamics for B(XN) (see limits on the
    % for loop above) implies no incoming matter from beyond XN.
     
       
    theta  = 1 - X.^2/4.*sin(2*B./X.^2).^2 ; % from Eq. (5.15)
    theta(1) = 1; % to avoid NaN
    
    % There is a spike in Theta in the outgoing shock, it
    % can be removed with the following two lines to give a plot
    % that looks smoother, but this is not essential.
    [MT I] = min(theta);
    theta(I) = theta(I+1);

    % Horizon test
    bhformed = any(theta<0) && bh == false;
    bhgone = all(theta>0) && bh==true ;

    Rho1 = (Bt-B) ./ X .^ 2 / (4*pi*dt); % from Eq. (4.15)
    Rho1(1) = Rho1(4);
    Rho1(2) = Rho1(4);
    Rho1(3) = Rho1(4);

    %% Figures
    %  Update figures every 20 timesteps, or when the horizons are gone
    if (mod(counter,20)== 0 || bhgone)
    
    %    tiledlayout(3,1) 
        tiledlayout(1,1) 
        
        set(gcf,'color','w');
     
    
        % p1 = nexttile;
        % plot(p1,X,Rho1,'LineWidth',2);  
        % grid on;
        % xlim([0 xm]);
        % 
        % if max(Rho1*1.1) > y_max 
        %     y_max = max(Rho1)*1.1;
        % end
        % ylim([0 y_max]);
        % ylabel('\rho  ','Rotation',0, 'fontsize',20);
      
     
        p4 = nexttile;
        plot(p4,X,B,'LineWidth',2);  
        grid on;
        xlim([0 xm]);
        ylim([min(B)*1.1,0]);
        ylabel('B  ','Rotation',0, 'fontsize',20);
        % 
        % 
        % Time box in first plot
        delete(a);
        dim = [.445 .875 .145 .05];
        a = annotation('textbox',dim,'String',['t = ', ...
            num2str(t,'%.2f')], 'Horizontalalignment','center',...
            'fontsize',16);
        % 
        % p2 = nexttile; 
        % plot(p2,X,zeroline,X,theta,'LineWidth',2);  
        % grid on;
        % xlim([0 xm]);
        % ylim([-1.5 1.5]);
        % ylabel('\Theta_+  ','Rotation',0,'fontsize',18);
        %xlabel('r','fontsize',14);
        %set(gca,'FontSize',12);

%         p3 = nexttile;
%         plot(p3,X,B ./ X .^ 2,'LineWidth',2);  
%         grid on;
%         xlim([0 xm]);
%         ylim([min(B ./ X .^ 2)*1.1 0]);
%         ylabel('\beta  ','Rotation',0, 'fontsize',14);
%         %xlabel('r','fontsize',14);
%         %set(gca,'FontSize',12);
    
        drawnow;
        
        
        
        % frame = getframe (gcf);
        % writeVideo (v, frame);   
        hold off;

%         if counter == 0
%             pause;
%         end
    
    end

    counter = counter + 1; 
    
    %% Black hole lifetime 
    
    if bhformed
        ti = t;
        bh = true;
    elseif bhgone
        tf = t;
        break
    end
      
  end
  % End of one mass run

  delete(a);


  %% Print data for mass loop  and continue outer loop with next mass
  close (v)
  fprintf('%2.2f      %2.2f\n',mass,tf-ti);

end
% End of mass for loop

toc;

%% Flux function used in the evolution of the solution
function F = flux(uL,uR,x)
 
        FL = (x^3)/2*sin(uL)^2;
        FR = (x^3)/2*sin(uR)^2;
        
        if uL <= uR 
            F = min(FL,FR);
        elseif uL > uR 
            if (uR > -pi/2 || uL < -pi/2)
                F = max(FL,FR); 
            else
                F = (x^3)/2;
            end
        end 

end
