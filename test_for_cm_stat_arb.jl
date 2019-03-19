using Distributions, ProgressMeter, Printf
using Gadfly

include("ll_cm_optim.jl")


function t_test_mult(sigma_mat, params)
    t_stats = zeros(Float64, 2)
    t_stats[1] = params[1] / sqrt(sigma_mat[1,1])
    t_stats[2] = -params[3] / sqrt(sigma_mat[3,3])
    return t_stats
end

function min_t(sigma_mat, params)
    temp = minimum(t_test_mult(sigma_mat, params))
    return temp
end




function gen_min_t_dist_cm(n, nsim, mu, sigma, lambda)
    i_n = collect(1:n)

    min_ts = ones(nsim)
    @showprogress 1 "Computing..." for j in 1:nsim
        z_i = rand(Normal(0,1),n)
        delta_v = mu .+ sigma * i_n.^lambda .* z_i
        res, hess = min_neg_ll_cm(delta_v)

        mu_hat, sigma2_hat, lambda_hat = get_minimizers(res)
        if mod(j+1000,1000) == 0
            @printf("Maximum LL:  %.6f\n", -res.minimum)
            @printf("True mu: %.6f --- Estimated mu: %.6f\n" , mu, mu_hat)
            @printf("True sigma2: %.6f --- Estimated sigma2: %.6f\n" , sigma^2, sigma2_hat)
            @printf("True lambda: %.6f --- Estimated lambda: %.6f\n" , lambda, lambda_hat)
        end


        sigma_mat = inv(hess)
        min_ts[j] = min_t(sigma_mat, [mu_hat, sigma2_hat, lambda_hat])
    end
    return min_ts
end

function get_min_t_cm(S_t)
    v = [1.0]
    append!(v, S_t)
    delta_v = v[2:end] - v[1:(end-1)]

    res, hess = min_neg_ll_cm(delta_v)
    mu_hat, sigma2_hat, lambda_hat = get_minimizers(res)

    @printf("Maximum LL:  %.6f\n", -res.minimum)
    @printf("True mu: %.6f --- Estimated mu: %.6f\n" , mu, mu_hat)
    @printf("True sigma2: %.6f --- Estimated sigma2: %.6f\n" , sigma^2, sigma2_hat)
    @printf("True lambda: %.6f --- Estimated lambda: %.6f\n" , lambda, lambda_hat)

    sigma_mat = inv(hess)
    min_t_1 = min_t(sigma_mat, [mu_hat, sigma2_hat, lambda_hat])
    return min_t_1, params, sigma_mat
end


min_ts_1 = gen_min_t_dist_cm(500,5000,0,0.01,0)
plot(x=min_ts_1, Geom.histogram())
quantile(min_ts_1, [0.9,0.95,0.99])

# min_ts_2 = gen_min_t_dist_cm(500,500000,0,0.01,0)
# plot(x=min_ts_2, Geom.histogram())
# quantile(min_ts_2, [0.9,0.95,0.99])

# using DelimitedFiles
# writedlm("data/min_t_dist/min_t_stats.csv", min_ts_1,',')
#
#
# min_ts_2 = gen_min_t_dist_cm(400,5000,0,0.01,0)
# min_ts_3 = gen_min_t_dist_cm(400,5000,0,0.01,0)
# min_ts_4 = gen_min_t_dist_cm(400,5000,0,0.01,0)
# min_ts_5 = gen_min_t_dist_cm(400,5000,0,0.01,0)
# min_ts_6 = gen_min_t_dist_cm(400,5000,0,0.01,0)
#
# cut_off_values_1 = quantile(min_ts_1, [0.9,0.95,0.99])
#
#
#
# cut_off_values_2 = quantile(min_ts_2, [0.9,0.95,0.99])
#
#
#
#
# cut_off_values_3 = quantile(min_ts_3, [0.9,0.95,0.99])
#
#
# cut_off_values_4 = quantile(min_ts_4, [0.9,0.95,0.99])
# cut_off_values_5 = quantile(min_ts_5, [0.9,0.95,0.99])
# cut_off_values_6 = quantile(min_ts_6, [0.9,0.95,0.99])
#
