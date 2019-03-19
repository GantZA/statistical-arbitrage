using Optim

function min_neg_ll_um(delta_v, initial_params, optimizer=Newton())
    n = length(delta_v)
    i_n = collect(1:n)
    log_i_n = log.(i_n)
    log_i_n_sq = log_i_n.^2

    function g!(grad,x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]
        theta = x[4]
        mu_delta_v_sum = mu * i_n .^ theta .- delta_v
        mu_delta_v_sum_sq = mu_delta_v_sum.^2
        i_2_lambda = i_n.^(-2*lambda)
        i_2_lambda_theta = i_n.^(theta-2*lambda)

        grad[1] = -(-(1/exp(gamma)) * sum(i_2_lambda_theta .* mu_delta_v_sum))
        grad[2] = -(-n/2 + (1/(2*exp(gamma))) * sum(i_2_lambda .*  mu_delta_v_sum_sq))
        grad[3] = -(-sum(log_i_n) + (1/exp(gamma)) * sum( i_2_lambda .* log_i_n .* mu_delta_v_sum_sq))
        grad[4]= -(-mu/exp(gamma) * sum(i_2_lambda_theta .*  mu_delta_v_sum .* log_i_n))
    end

    function h!(H,x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]
        theta = x[4]
        mu_delta_v_sum = (mu*(i_n.^theta)) .- delta_v
        mu_delta_v_sum_sq = mu_delta_v_sum.^2
        two_mu_delta_v_sum = 2 * mu * i_n.^theta .- delta_v
        i_2_lambda = i_n.^(-2*lambda)
        i_2_lambda_theta = i_n.^(theta-2*lambda)
        i_2_lambda_2_theta = i_n.^(2*theta-2*lambda)



        H[1,1] = -(-1/exp(gamma) * sum(i_2_lambda_2_theta))
        H[1,2] = -(1/exp(gamma) * sum(i_2_lambda_theta .* mu_delta_v_sum))
        H[1,3] = -(2/exp(gamma) * sum(i_2_lambda_theta .* log_i_n .* mu_delta_v_sum))
        H[1,4] = -(-1/exp(gamma) * sum(i_2_lambda_theta .* log_i_n .* two_mu_delta_v_sum))

        H[2,1] = H[1,2]
        H[2,2] = -(-1/(2 * exp(gamma)) * sum(i_2_lambda .* mu_delta_v_sum_sq))
        H[2,3] = -(-1/exp(gamma) * sum(i_2_lambda .* log_i_n .* mu_delta_v_sum_sq))
        H[2,4] = -( mu/exp(gamma) * sum(i_2_lambda_theta .* log_i_n .* mu_delta_v_sum))

        H[3,1] = H[1,3]
        H[3,2] = H[2,3]
        H[3,3] = -(-2/exp(gamma) * sum(i_2_lambda .* log_i_n_sq .* mu_delta_v_sum_sq))
        H[3,4] = -(2*mu/exp(gamma) * sum(i_2_lambda_theta .* log_i_n_sq .* mu_delta_v_sum))

        H[4,1] = H[1,4]
        H[4,2] = H[2,4]
        H[4,3] = H[3,4]
        H[4,4] = -(-mu/exp(gamma) * sum(i_2_lambda_theta .* log_i_n_sq .* two_mu_delta_v_sum))
    end

    function f(x)
        mu = x[1]
        gamma = x[2]
        lambda = x[3]
        theta = x[4]

        return -(-n/2 * gamma - lambda * sum(log(i) for i=1:n)
            -sum(i^(-2*lambda) * (mu*i^theta - delta_v[i])^2 for i=1:n)/ (2*exp(gamma)))
    end

    x0 = initial_params
    res = optimize(f, g!, h!, x0,optimizer)
    hess = zeros(4,4)
    h!(hess, res.minimizer)
    return res, hess

end

function get_minimizers(results)
    mu, gamma, lambda, theta = results.minimizer
    sigma2 = exp(gamma)
    return mu, sigma2, lambda, theta
end
