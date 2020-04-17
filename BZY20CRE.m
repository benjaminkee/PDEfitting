function [R] = BZY20CRE(x0,tspan,geom,id,T)

%%	Unpackaging inputs
%	Unpackage geometry params
	th		= geom.th;
	A		= geom.A ;
%	Unpackage variables
	D_oho_0 = x0(id.D_oho_0 );
	D_vao_0 = x0(id.D_vao_0 );
	Ea_oho	= x0(id.Ea_oho	)*1e6;
	Ea_vao	= x0(id.Ea_vao	)*1e6;
	Rb		= x0(id.Rb		);

%%	Initial state conditions
%	Temperature
	T_init		= T.init;							% deg C
	T_K_init	= T_init+273.15;					% K
%	Bubbler conditions
	T_bub		= 25;								% deg C
	T_bub_K		= T_bub+273.15;						% K
%	Pressure
	P_h2o_init	= exp(20.386-5132/T_bub_K)/760/1;	% atm
	P_h2		= 0.03*oneatm;						% atm
%	Dopant amount
	Dop			= 0.2;								% Lattice
%	Calculate initial proton concentration
	OHO_init = calcOHO(P_h2o_init,T_K_init,Dop)+Rb;	% Lattice
	
%%	Final state conditions
	T_set = T.set;
	T_K_set = T_set+273.15;
	P_h2o_set = exp(20.386-5132/T_bub_K)/760/1;
	OHO_set = calcOHO(P_h2o_set,T_K_set,Dop)+Rb;

%	Calculate diffusion coefficients at final state temperature	
	D_oho = D_oho_0 * exp(-Ea_oho/gasconstant/T_K_set);
	D_vao = D_vao_0 * exp(-Ea_vao/gasconstant/T_K_set);

%%	Setup PDE solver
	m = 0;
	pde = @(x,t,u,dudx) pdefunc(x,t,u,dudx,Dop,D_oho,D_vao);
	ic = @(x) pdeic(x,OHO_init);
	bc = @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,OHO_set);
	xmesh = linspace(0,th,200);
%	Solve PDE
	sol = pdepe(m,pde,ic,bc,xmesh,tspan);
	
% 	u = sol(:,:);
% 	close all
% 	surf(xmesh,tspan,u)
% 	shading interp
% 	xlabel('Thickness (m)')
% 	ylabel('Time (s)')
% 	zlabel('Proton Conc(x,t)')
% 	fixplot

%%	Calculate conductivity
%	Charges of species
	z_oho = 1;
	z_vao = -2;
%	Lattice to cubic meters
	molV = 0.0000449426990603*1000000*1000/100^3;
	
	
%	Find mesh spacing, in cm
	dx = (xmesh(2)-xmesh(1))*100;

%	Nernst-Einstein relationship for get conductivity
	OHO_sol = sol;
	VAO_sol = (Dop-OHO_sol)/2;
	FRT = Faraday^2/gasconstant/T_K_set;
	sigma_oho = z_oho^2*FRT*OHO_sol*D_oho/molV*10;
	sigma_vao = z_vao^2*FRT*VAO_sol*D_vao/molV*10;
	sigma_total = (sigma_oho+sigma_vao); %mS/cm
	
	sig_avg = (sigma_total(:,2:end)+sigma_total(:,1:end-1))/2;
	Ra = sum(dx./sig_avg,2)/0.45*1000;
	
%	Maximum works better than integral over all cells
	aaa= max(sigma_total,[],2);
	R = th./aaa/A*100000;
	
	
% 	figure
% 	plot(tspan,sigma_total,tspan,sigma_oho,tspan,sigma_vao)
% 	xlabel('Time (s)')
% 	ylabel('\sigma (mS/cm)')
% 	fixplot
% 	
% 	figure
% 	plot(tspan,OHO_sol,tspan,VAO_sol)
% 	xlabel('Time (s)')
% 	ylabel('Concentration')
% 	fixplot

end

%%	Define the pde function
function [c,f,s] = pdefunc(x,t,u,dudx,Dop,D_oho,D_vao)
%	Degree of hydration
	S = u/Dop;
%	Water diffusion coefficient (Kreuer 1999)
	D_h2o = (2-S)*D_oho*D_vao/(S*D_oho+2*(1-S)*D_vao);
%	Coefficients for PDE
%	In the form c*du/dt = x^-m*d/dx(x^m*f)+s
%	Where f is a flux term and s is a source term
%	u is the concentration
	c = 1;
	f = D_h2o*dudx;
	s = 0;
%	Simplified, this gives Fick's law
end

%%	PDE initial condition
function u0 = pdeic(x,OHO_init)
	u0 = OHO_init;
end

%%	PDE boundary conditions
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t,OHO)
%	Set the BC to proton concentrations, same on both sides
	pl = ul-OHO;
	ql = 0;
	pr = ur-OHO;
	qr = 0;
end

%%	Calculate proton concentration
function OHO = calcOHO(P_h2o,T_K,Dop)
%	Kreuer 1994 proton concentration based on p_h2o
%	Use Delta H and Delta S from BZY20 data (zhu 2018)
	DH = -93300000; % J/kmol
	DS = -103200;	% J/mol-K
	K_hyd = exp(-DH/(gasconstant*T_K)+DS/(gasconstant));
	Kpp = K_hyd*P_h2o;
	OHO = (3*Kpp-sqrt(Kpp*((Kpp-4)*(Dop-3)^2+36)))/(Kpp-4);
end






