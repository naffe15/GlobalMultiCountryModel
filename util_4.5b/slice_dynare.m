% ---------------------------------------------------------------------------------
% SINGLE-VARIABLE SLICE SAMPLING as in Neal (2003),
% "Slice Sampling", Annals of Statistics 31, 705-67 
% Project DMM - G.Fiorentini, C.Planas, and A.Rossi (Aug 2011)
% INPUT: 
%        1. it:  position of theta(i) into the theta vector
%        2. theta (nt x 1) all parameters 
%        4. thetaprior: (4 x 1) 1,2 hyperparametrs; 3,4 Lower and upper
%                        bounds for theta(it)
% OTPUT: 
%        1. NEVAL: # of log-posterior evaluations
%        2. XSIM:  simulated value for theta(it) 
%
% NOTE: 
% log_posterior = -dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
% ---------------------------------------------------------------------------------
function [NEVAL,XSIM] = slice_dynare(it,theta,thetaprior,dataset_,dataset_info, options_, M_, estim_params_,bayestopt_,bounds, oo_)
NEVAL = 0;
XOLD  = theta(it);
XLB   = thetaprior(3);
XUB   = thetaprior(4);

% -------------------------------------------------------
% 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
%    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
% -------------------------------------------------------
theta(it) = XOLD;
FXOLD     = -dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
NEVAL     = NEVAL + 1;
Z         = FXOLD + log(rand(1,1));

% -------------------------------------------------------------
% 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE 
%    STEPPING-OUT PROCEDURE
%    W = an estimate of the scale of SC    
%    M = Limit on steps (-1 = +INF)
% -------------------------------------------------------------
M = -1;
W = max((XUB-XLB)/10,1); 
U = rand(1,1);
L = XOLD - W*U;
R = XOLD + W - W*U;  % L + W
if (M == -1)
    while (L > XLB)
        theta(it) = L;
        FXL=-dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
        NEVAL = NEVAL + 1;
        if (FXL <= Z) 
            break;
        end
        L = L - W;
    end
    while (R < XUB)
        theta(it) = R;
        FXR=-dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
        NEVAL = NEVAL + 1;
        if (FXR <= Z) 
            break;
        end
        R = R + W;
    end
else
% never used into DMMM but can be coded if needed
end
if (L < XLB) 
    L = XLB;
end
if (R > XUB) 
    R = XUB;
end
	
% ------------------------------------------------------
% 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
% ------------------------------------------------------
OK = 0;
while (OK == 0)
    U = rand(1,1);
    XSIM = L + U*(R - L);
    theta(it) = XSIM;
    FXSIM=-dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
    NEVAL = NEVAL + 1;  
    if (FXSIM >= Z) 
        OK = 1;
    end
    if (XSIM > XOLD) 
        R = XSIM;
    else
        L = XSIM;
    end
end
theta(it) = XOLD;     