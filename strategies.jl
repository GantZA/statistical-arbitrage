"""
Online ZBCRP and Anti-ZBCRP Expert Generating Algorithms
Inputs:
       X: Matrix of price relatives
       k: Window size parameter
       L: Free-Parameter
       Anti: Boolean indicating which algorithm to use. TRUE = Anti
       absExp: Boolean indicating whether the expert is absolute
       (TRUE) or active(FALSE)
       risk_aver_vec: Not in use
Outputs:
       b: Vector of an individual expert's controls
Workflow:
The online ZBCRP agent generating algorithms choose stocks by
searching for trends (ZBCRP) or mean-reversion (Anti-ZBCRP)
in the price relatives for a given window k. The ZBCRP tries
to buy stocks that are going up and shorts those that are going
The Anti-ZBCRP does the opposite and looks to buy stocks
that have gone down and sell stocks that have gone up
"""

function online_zbcrp(X, k, Anti, abs_exp = False)
    exp_tup = X[:,:]
    M = size(exp_tup)[2]
    onesM = ones(M)
    if Anti == false
        mu = colwise(mean,exp_tup)
    else
        mu = - colwise(mean,exp_tup)
    end

    if abs_exp == true
        if k < M
            sigma = eye(M)
            inv_sigma = sigma
        else
            sigma = cov(exp_tup)
            if cond(sigma) < (1/1e-2) && k > M
                inv_sigma = inv(sigma)
            else
                sigma = diagm(diag(sigma))
                if cond(sigma) > 1/1e-2
                    sigma = eye(M)
                    inv_sigma = inv(sigma)
                elseif any(sigma .== 0) == true
                    sigma = ZeroVar(sigma)
                    inv_sigma = inv(sigma)
                end
            end
        end
        b = OptimSemiLog(mu[:,1], sigma, 1)
    else
        if k < 3
            sigma = eye(M)
            inv_sigma = sigma
        else
            sigma = cov(exp_tup)
            if cond(sigma) < (1/1e-2) && k > M
                inv_sigma = inv(sigma)
            else
                sigma = diagm(diag(sigma))
                if cond(sigma) > 1/1e-2
                    sigma = eye(M)
                    inv_sigma = inv(sigma)
                elseif any(sigma .== 0) == true
                    sigma = zero_var(sigma)
                    inv_sigma = inv(sigma)
                end
            end
        end
        sigma_scale = transpose(onesM)*inv_sigma*onesM
        risk_aver = transpose(onesM) * abs.( inv_sigma * ( mu - onesM * ((transpose(onesM)*inv_sigma*mu)/sigma_scale)))
        if risk_aver==0
            risk_aver = 1
        end
        b = inv_sigma/risk_aver * (mu - onesM*((transpose(onesM)*inv_sigma*mu)/sigma_scale))
    end
    return b
end

function zero_var(S)
    z_ind = (diag(S) .== 0)
    S[z_ind, z_ind] = mean(diag(S))*eye(sum(z_ind))
    return S
end

"""Online ZANTICOR Expert Generating Algorithm
Inputs:
       X: Matrix of price relatives
       K: Window size parameter
       T: Total number of days of price relatives available
       H0: Matrix of previous days experts controls
       absExp: Boolean indicating whether the expert is
       absolute (TRUE) or active(FALSE)
Outputs:
       b: Vector of an individual expert's controls

Workflow:
The online ZANTICOR agent generating algorithm chooses stocks
by searching for correlations between stocks in consecutive
time periods. Essentially the algorithm tries to buy stocks
that will mean revert upwards and short stocks that will mean
revert downwards."""

function online_anticor(X, K, T, H0, abs_exp)
    if T < 2*K
        b = H0
    else
        b = Vector{Float64}(size(X)[2]) .= 0
        exp_tup = X[:,:]
        M = size(exp_tup)[2]
        LX1 = exp_tup[end-2*K+1:end-K,:]-1
        LX2 = exp_tup[end-K+1:end,:]-1
        mu1 = mean(LX1,1)
        mu2 = mean(LX2,1)
        claims = Matrix(M,M) .= 0
        H1 = Vector{Float64}(M)

        if K > 1
            sigma1 = transpose(mapslices(std, LX1, [1]))
            sigma2 = transpose(mapslices(std, LX2, [1]))
            M_cov = (1/(K-1)) * (transpose(LX1-(ones(K)*mu1))) * (LX2-(ones(K)*mu2))
            M_cor = M_cov ./ (sigma1 * transpose(sigma2))
            M_cor[find(isnan.(M_cor))] = 0

            for j in 1:M
                claims[:,j] = claim(j, M_cor[:,j], diag(M_cor), M, claims[:,j], mu2)
            end
        else
            sigma1 = zeros(M)
            sigma2 = zeros(M)
            M_cov = zeros(M,M)
            M_cor = ones(M,M)
            claims = eye(M)
        end
        if abs_exp == true
            transfers = H0 * transpose(ones(M)) .* claims ./ (claims * ones(M) * transpose(ones(M)))
            transfers[find(isnan.(transfers))] = 0
            for i in 1:M
                H1[i] = H0[i] + sum(transfers[:,i] - transfers[i,:])
            end
        else
            for i in 1:M
                H1[i] = H0[i] + 1/3*(sum(claims[:,i] - claims[i,:]))
            end
        end
        b = ReNormAgMix(H1, absExp)
    end
    return b
end

"""Additional function to calculate a claim for column j.
Julia works faster when looping through columns than
through rows."""
function claim(j, Mcorj, McorDiag, M, claim_j, mu)
    for i in 1:M
        if mu[i] >= mu[j] && Mcorj[i] > 0
            claim_j[i] = Mcorj[i] + max(0, -McorDiag[i]) +
            max(0, -McorDiag[j])
        end
    end
    return claim_j
end
