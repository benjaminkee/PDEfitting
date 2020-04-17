function fitCRE2
%%	Initial stuff
%	Load files
% 	load('bzy20-425-450relax.mat','times','resis');
	load('bzy20_2.mat','times','resis');
%	Adjust time to start at zero
	xdata = times'-times(1);
	ydata = resis;

%   Geom Params
	geom.A = 0.45;
	geom.th = 0.0015;
	
%	Temperatures
	T.init = 650;	% C
	T.set  = 675;	% C

%	Indexes
	vars = {'D_oho_0','D_vao_0','Ea_oho','Ea_vao','Rb'};
	nVars = length(vars); 
	id = struct();
	for i=1:nVars
		id.(vars{i}) = i;
	end

%%  Initial guesses
	x0 = zeros(nVars,1);
%   Diffusion coefficients
	x0(id.D_oho_0) = 1.5e-3/100^2; % m^2/s
	x0(id.D_vao_0) = 1.0e-2/100^2; % m^2/s
%   Activation energies
	x0(id.Ea_oho) = 45000000/1e6;  % MJ/kmol
	x0(id.Ea_vao) = 86000000/1e6;  % MJ/kmol
%   Surface concentration offset
	x0(id.Rb) = 0;
        
%   Works better fitting offset for surface concentration, not thermo
% 		DH = -93300000; % J/kmol
% 		DS = -103200;	% J/mol-K


%%  Set bounds
	lb = zeros(nVars,1)  ;	ub				=   lb;
	lb(id.D_oho_0)	=   0;	ub(id.D_oho_0)	=    1; 
	lb(id.D_vao_0)	=   0;	ub(id.D_vao_0)	=    1;
	lb(id.Ea_oho)	=   0;	ub(id.Ea_oho)	= 1000;
	lb(id.Ea_vao)	=   0;	ub(id.Ea_vao)	= 1000;
	lb(id.Rb)		= -20;	ub(id.Rb)		=   20;
		
%%	Least square fitting
%	Define minimization function
	lsqfunc = @(x0,x_data) BZY20CRE(x0,xdata,geom,id,T);
		
%	Solve
	sol = lsqcurvefit(lsqfunc,x0,xdata,ydata,lb,ub);
		
%	Retrieve curve from fit variables
	ySol = BZY20CRE(sol,xdata,geom,id,T);
	close all;
	figure
	plot(xdata,ySol,xdata,ydata)
		
	sol
		
end