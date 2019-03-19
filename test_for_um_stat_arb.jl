using Distributions, ProgressMeter, Printf, LinearAlgebra

include("ll_um_optim.jl")

function get_inv(x)
    try
        return inv(x)
    catch y
        print(y)
        return false
    end
end

function gen_lambda_theta_pair(val)
    if val < 1
        return [-val-0.5, -1]
    elseif val < (2+sqrt(2))/2
        temp = (val-((2+sqrt(2))/2))/sqrt(2)
        return [temp, temp-0.5]
    elseif val < (3+sqrt(2))/2
        temp = val - (3+sqrt(2))/2
        return [0, temp]
    else
        temp = (val - (3+sqrt(2))/2)/sqrt(2)
        return [temp, temp]
    end
end

function gen_det_seq(leng)
    result = [ones(leng) ones(leng)]'
    vals = collect(range(0,stop=((3+4*sqrt(2))/2), length=leng))
    for i in 1:leng
        result[:,i] = gen_lambda_theta_pair(vals[i])
    end
    return result
end

function t_test_mult(sigma_mat, params)
    t_stats = zeros(Float64, 4)
    t_stats[1] = params[1] / sqrt(sigma_mat[1,1])
    t_stats[2] = (params[4]-params[3] + 0.5)/sqrt(sigma_mat[4,4] + sigma_mat[3,3] - 2*sigma_mat[3,4])
    t_stats[3] = (params[4] + 1)/sqrt(sigma_mat[4,4])
    t_stats[4] = max(-params[3] / sqrt(sigma_mat[3,3]), (params[4]-params[3])/sqrt(sigma_mat[4,4] + sigma_mat[3,3] - 2*sigma_mat[3,4]))
    return t_stats
end

function min_t(sigma_mat, params)
    temp = minimum(t_test_mult(sigma_mat, params))
    return temp
end

function gen_min_t_dist_um(n, nsim, npars, mu, sigma)
    i_n = collect(1:n)

    crit_vals = zeros(npars, 3)
    @showprogress 100 "Computing..." for i in 1:npars
        min_ts = ones(nsim)
        lambda, theta = gen_lambda_theta_pair((i-1)*((3+4*sqrt(2))/2)/npars)
        @printf("\nlambda=%f and theta=%f\n", lambda, theta)
        j = 1
        while j <= nsim

            z_i = rand(Normal(0,1),n)
            delta_v = mu .* i_n.^theta  .+ sigma * i_n.^lambda .* z_i
            res, hess = min_neg_ll_um(delta_v, [mu, sigma, lambda, theta], Newton())
            mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)
            sigma_mat = get_inv(hess)
            if sigma_mat == false || all(diag(sigma_mat) .> 0.0) == false
                res, hess = min_neg_ll_um(delta_v, [mu, sigma, lambda, theta], NewtonTrustRegion())
                mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)
                sigma_mat = get_inv(hess)
                if sigma_mat == false || all(diag(sigma_mat) .> 0.0) == false
                    res, hess = min_neg_ll_um(delta_v, [mu, sigma, lambda, theta], LBFGS())
                    mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)
                    sigma_mat = get_inv(hess)
                    if sigma_mat == false || all(diag(sigma_mat) .> 0.0) == false
                        res, hess = min_neg_ll_um(delta_v, [mu, sigma, lambda, theta], NelderMead())
                        mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)
                        sigma_mat = get_inv(hess)
                        if sigma_mat == false || all(diag(sigma_mat) .> 0.0) == false
                            @printf("\nERROR!!!\n")
                            continue
                        end
                    end
                end
            end
            min_ts[j] = min_t(sigma_mat, [mu_hat, sigma2_hat, lambda_hat, theta_hat])
            j += 1
            # @printf("t-test stat: 4: %.6f, with estimate  %.6f and variance %.6f\n", t_test_stats[4], lambda_hat, sigma_mat[3,3])
            # @printf("min-t test stats: %.6f \n" , min_ts[j])
        end
        crit_vals[i,:] = quantile(min_ts, [0.9,0.95,0.99])
    end
    return crit_vals
end

function get_min_t_um(S_t)
    v = [1.0]
    append!(v, S_t)
    delta_v = v[2:end] - v[1:(end-1)]

    res, hess = min_neg_ll_um(delta_v)
    mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)

    @printf("Maximum LL:  %.6f\n", -res.minimum)
    @printf("True mu: %.6f --- Estimated mu: %.6f\n" , mu, mu_hat)
    @printf("True sigma2: %.6f --- Estimated sigma2: %.6f\n" , sigma^2, sigma2_hat)
    @printf("True lambda: %.6f --- Estimated lambda: %.6f\n" , lambda, lambda_hat)
    @printf("True theta: %.6f --- Estimated theta: %.6f\n" , theta, theta_hat)

    sigma_mat = inv(hess)

    min_t_1 = min_t(sigma_mat, [mu_hat, sigma2_hat, lambda_hat, theta_hat])
    return min_t_1, params, sigma_mat
end

# AIM: 0.4034, 0.6004, and 0.9074 (100,1000)
crit_min_ts = @time gen_min_t_dist_um(500,1000,100,-10e-6,0.01)
maximum(crit_min_ts, dims = 1)

using Gadfly
plot(x=min_ts_um_1, Geom.histogram())
quantile(min_ts_um_1, [0.9,0.95,0.99])

trouble

delta_v1 = trouble[:,2]
res, hess = min_neg_ll_um(delta_v1, [10e-6,0.01,0.0,0.0], LBFGS())

mu_hat, sigma2_hat, lambda_hat, theta_hat = get_minimizers(res)
inv(hess)



     = issuccess(lu(x, check=false))

isinvertible(hess)

a= 1
if a==2 || a<2
    print(a)
end
