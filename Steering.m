%%  STEERING    Calculates the maximal steering quantity over all assemblages
%   This function has one required input argument:
%     M: a cellarray consisting of a collection of ensembles. 
%
%   STEERING_OPT = STEERING(M) gives the optimal / maximal value achievable
%   via SDP in a steering scenario. 
%
%   In a standard steering scenario, Alice's measurements and state are
%   unknown, while Bob's measurements and state are known.  
%
%   URL: http://www.qetlab.com/Steering

%   requires: CVX (http://cvxr.com/cvx/)
%   author: Vincent Russo (vincentrusso1@gmail.com)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 26, 2015

function steering_opt = Steering(M)

    % check to ensure that M is non-empty
    if isempty(M)
        error('Steering:EmptyMeasurementSet','M must be a non-empty cellarray.');
    end
    
    % check to ensure all elements of M have the same dimension and also 
    % check that elements of M are square
    [~,sz] = size(M);
    [~,dim] = size(M{1});
    for i = 1:sz
        if size(M{i}) ~= dim
            error('Steering:DimensionError','The elements of M must all have the same dimension.');
        end
        if ndims(M{i}) ~= 2
            error('Steering:DimensionError', 'The elements of M must be square.');
        end
    end
    
    tol = eps^(3/4); % numerical tolerance used

    % Steering SDP
    cvx_begin sdp quiet
        %#ok<*VUNUS> % suppress MATLAB warnings for equality in CVX.
        %#ok<*EQEFF> % suppress MATLAB warnings for inequality in CVX.
        cvx_precision(tol);
        variable rho(dim,dim,dim,dim) semidefinite 
        variable rho_hat(dim,dim) semidefinite     

        % construct objective function
        obj_fun = 0;
        m_count = 1;
        for x = 1:dim
           for a = 1:dim
               obj_fun = obj_fun + ip( M{m_count}, rho(:,:,x,a) );
               m_count = m_count + 1;
           end
        end

        maximize obj_fun

        subject to

            rho_x_sum = sum(rho,3);
            for a = 1:dim
                rho_x_sum(:,:,:,a) == rho_hat;
            end

            trace(rho_hat) == 1;

    cvx_end

    steering_opt = cvx_optval;

end