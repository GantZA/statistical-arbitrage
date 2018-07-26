"""Online Learning Algorithm
Inputs:
       q0: Vector of prior expert mixtures
       b0: Vector of prior aggregate portfolio controls
       x: Most recent price relative
       H0: Matrix of prior experts controls
       H1: Matrix of subsequent expert controls
       S_t0: Prior aggregate portfolio wealth
       S_n0: Vector of prior expert portfolio wealth
       Rule_g: Update Rule, either "Uni" or "EWMA"
       lambda: EWMA parameter
       absol: Boolean for Long-only aggregate portfolio (TRUE) and zero-cost (FALSE)
Outputs:
       b1: Vector of latest aggregate portfolio controls
       q1: Vector of latest expert mixtures
       S_t1: Updated aggregate portfolio wealth
       S_n1: Updated expert portfolio wealth

Workflow:
The online learning algorithm updates the aggregate portfolio wealth and the ex-
perts portfolio wealth and determines the expert combinations or mixtures that will
be used to create the aggregate portfolio controls for the next time period.
The OnlineLearn runs through the algorithm, performing the necessary updates,
and calls the relevant functions to renormalize etc."""

function online_learn(q0, b0, x, H0, H1, S_t0, S_n0, Rule_g, eta, absol)

    S_t1 = min(S_t0, 1) * (transpose(b0) * (x-1)+1) + max(0, S_t0-1)
    S_n1 = min.(S_n0, 1) .* (H0 * (x-1)+1) + max.(0, S_n0-1)

    q1 = online_function(Rule_g, q0, S_n1, eta, absol)
    b1 = transpose(H1) * q1

    v = sum(abs.(b1))
    if v > eps()
        b1 = 1/v * b1
        q1 = 1/v *q1
    end
    return b1, q1, S_t1, S_n1
end

function online_function(g, q0, PL_t, eta, absol)
    if g == "Uni"
        q1 = PL_t
        q1 = re_norm_arg_mix(q1, absol)
    elseif g == "EG"
        if all(isapprox.(q0,0)) == true
            q1 = PL_t
            q1 = re_norm_arg_mix(q1, absol)
        else
            temp = q0 .* (exp.((eta*PL_t)/(transpose(q0)*PL_t)))
            q1 = temp / (sum(abs.(temp)))
            a = sum(q1)
            q1 = re_norm_arg_mix(q1, absol)
        end
    elseif g == "EWMA"
        if all(isapprox.(q0,0)) == true
            q1 = PL_t
            q1 = re_norm_arg_mix(q1, absol)
        else
            q1 = eta*q0 + (1-eta)*(q0 .* PL_t)/(transpose(q0)*PL_t)
            q1 = re_norm_arg_mix(q1, absol)
        end
    end
    return q1
end

function re_norm_arg_mix(q0, absol)
    if absol == true
        q1= q0 / sum(q0)
    else
        absSum = sum(abs.(q0 - mean(q0)))
        if absSum > eps()
            q1 = (q0 - mean(q0)) / absSum
        else
            q1 = zeros(size(q0))
        end
        a = sum(abs.(q1))
    end
    return q1
end

function OnlineLearn_Comp(q0, b0, x, H0, H1, S_t0, S_n0, Rule_g, eta, absol)

    S_t1 = S_t0 * ((transpose(b0) * (x-1)) + 1)
    S_n1 = S_n0 .* ((H0 * (x-1)) + 1)

    q1 = OnlineFunction(Rule_g, q0, S_n1, eta, absol)
    b1 = transpose(H1) * q1

    v = sum(abs.(b1))
    if v > eps()
        b1 = 1/v * b1
        q1 = 1/v *q1
    end
    return b1, q1, S_t1, S_n1
end
