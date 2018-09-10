using Optim

function ll_pp(x,mu, log_sigma, lambda, theta)
    i_seq = collect(1:length(x))
    result = -(-1/2 * sum(log.(exp(log_sigma) .* i_seq.^(2*lambda)))
        - 1/(2 * exp(log_sigma)) .* sum(1./i_seq.^(2*lambda) .* (x .- (mu * i_seq .^ theta)).^2))
    return result
end

nvar = 4
func = TwiceDifferentiable(vars -> ll_pp(PLt_JSE2, vars[1], vars[2],
    vars[3], vars[4]), ones(nvar); autodiff=:forward)
opt = optimize(func, ones(nvar))

opt.minimizer
mu_hat,sigma_2_hat, lambda_hat, theta_hat = opt.minimizer[1],exp(opt.minimizer[2]), opt.minimizer[3], opt.minimizer[4]

mu_hat
sigma_2_hat
lambda_hat
theta_hat

mu_hat > 0
lambda_hat < theta_hat
theta_hat > max(lambda_hat-1/2, -1)
