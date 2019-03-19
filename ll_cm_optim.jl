using Optim


function min_neg_ll_cm(delta_v)
    n = length(delta_v)
    i_n = collect(1:n)
    log_i_n = log.(i_n)

    function g!(grad,x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]
        mu_delta_v_sum = mu .- delta_v
        mu_delta_v_sum_sq = mu_delta_v_sum.^2
        i_2_lambda = i_n.^(-2*lambda)

        grad[1] = -(-(1/exp(gamma)) * sum(i_2_lambda .* mu_delta_v_sum))
        grad[2] = -(-n/2 + (1/(2*exp(gamma))) * sum(i_2_lambda .*  mu_delta_v_sum_sq))
        grad[3] = -(-sum(log_i_n) + (1/exp(gamma)) * sum( i_2_lambda .* log_i_n .* mu_delta_v_sum_sq))
    end

    function h!(H,x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]
        mu_delta_v_sum = mu .- delta_v
        mu_delta_v_sum_sq = mu_delta_v_sum.^2
        i_2_lambda = i_n.^(-2*lambda)


        H[1,1] = -(-1/exp(gamma) * sum(i_2_lambda))
        H[1,2] = -(1/exp(gamma) * sum(i_2_lambda .* mu_delta_v_sum))
        H[1,3] = -(2/exp(gamma) * sum(i_2_lambda .* log_i_n .* mu_delta_v_sum))

        H[2,1] = H[1,2]
        H[2,2] = -(-1/(2 * exp(gamma)) * sum(i_2_lambda .* mu_delta_v_sum_sq))
        H[2,3] = -(-1/exp(gamma) * sum(i_2_lambda .* log_i_n .* mu_delta_v_sum_sq))

        H[3,1] = H[1,3]
        H[3,2] = H[2,3]
        H[3,3] = -(-2/exp(gamma) * sum(i_2_lambda .* (log_i_n .^2) .* mu_delta_v_sum_sq))
    end

    function f(x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]

        return -(-n/2 * gamma - lambda * sum(log(i) for i=1:n)
            -sum(i^(-2*lambda) * (mu - delta_v[i])^2 for i=1:n)/ (2*exp(gamma)))
    end

    x0 = [0.0,log(0.01^2),0.0]
    res = optimize(f, g!, h!, x0)
    hess = zeros(3,3)
    h!(hess, res.minimizer)
    return res, hess

end

function get_minimizers(results)
    mu, gamma, lambda = results.minimizer
    sigma2 = exp(gamma)
    return mu, sigma2, lambda
end
