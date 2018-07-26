"""Simulation
Inputs:
    dateInd: A boolean indicating whether the user will input the start and end
        dates as date values from the datevec or row indices of the X matrix
    startdate: A string of the start date for the trading period or an index indi-
        cating what time period the trading period begins
    enddate: A string of the last date for the trading period or an index indicating
        what time period the trading period ends
    datevec: If a you have a column vector of dates, ordered from oldest to most
        recent or a vector of integers representing the time periods
    X: The matrix of stock price relatives, with rows representing time and columns
        for each stock. The price relatives rows should match with the datevec and
        each element must be Float64 datatype.
    K: The vector of the expert generating algorithms window parameter, with
        K[1] belonging to ZBCRP, K[2] belonging to the Anti-ZBCRP and K[3] be-
        longing to ZANTICOR.
    W: The vector indicating which algorithms to include. Legacy parameter, set
        to [1,3] for now.
    lambda: The learning algorithm speciﬁc parameter. Used for EWMA only
    Learn_Alg: The update rule, either "Uni" for Universal or "EWMA" for EWMA.
    absPort: Boolean for absolute aggregate portfolio. True indicates absolute
        (Long-only) and False indicates active (zero-cost).
    absExp: Boolean for absolute expert portfolios. True indicates absolute
        (Long-only) and False indicates active (zero-cost).
Outputs:
    b_mat: Matrix of aggregate portfolio controls over time
    q_mat: Matrix of expert mixtures over time
    S_t: Vector of aggregate portfolio wealth over time
    S_nt: Matrix of expert portfolio wealth over time
    H_mat: Collection of Matrices of expert portfolio controls over time
Workflow:
The Simulation function handles all the expert generating algorithms, the
learning algorithms, the data and the parameters. Thus the user is only
required to provide a dataset of price relatives and the parameters for the
expert algorithms and learning algorithm. The Simulation function performs a backtest
on the given dataset and using the given parameters"""

function simulation(date_ind, start_date, end_date,
     X, K, W, lambda, learn_alg, abs_port = false, abs_exp = false)

    if date_ind == true
        start_ind = find(date_vec .== start_date)[1]
        end_ind = find(date_vec .== end_date)[1]
    else
        start_ind = start_date
        end_ind = end_date
    end

    T = end_ind - start_ind
    X = X[start_ind:end_ind, :]
    M = size(X)[2]
    N = sum(K)

    q_mat = Matrix(N,T) .= 0
    b_mat = Matrix(M,T) .= 0

    if abs_exp == true
        H_mat =  [Matrix(N,M) .= 1/M]
    else
        abs_port = false
        H_mat =  [Matrix(N,M) .= 0]
    end

    if abs_port == true
        b_mat[:, 1] = 1/M
        q_mat[:,1] = 1/N
    else
        b_mat[:, 1] = zeros(M)
        q_mat[:,1] = zeros(N)
    end

    S_t = Vector(T+1) .= 0
    S_t[1] = 1

    S_nt = Matrix(N,T+1) .= 1

    for t in 2:T
        H1 = exp_gen(K, W, M, X[1:(t-1),:], H_mat[(t-1)])
        b_mat[: , t], q_mat[:, t], S_t[t], S_nt[:,t] =online_learn(q_mat[:, (t-1)], b_mat[:, (t-1)], X[(t-1),:], H_mat[(t-1)],
        H1, S_t[t-1], S_nt[:, (t-1)] ,
        learn_alg, lambda, abs_port)

        H_mat = push!(H_mat, H1)
    end

    S_t[end] = min(S_t[end-1], 1) *
    (transpose(X.values[end,:]-1) * b_mat[:, end] + 1) +
     max(0, S_t[end-1] - 1)

    S_nt[:,end] = min.(S_nt[:,end-1], 1) .*
    ((H_mat[end] * (X.values[end,:]-1)) + 1) + max.(0,S_nt[:,end-1] - 1)

    return b_mat, q_mat, S_t, S_nt, H_mat
end

"""Returns the expert generating algorithm based on the w value """
function exp_alg(W, X, K, T, H0)
    if W == 1
        return online_zbcrp(X, K, false)
    elseif W == 2
        return online_zbcrp(X, K, true)
    elseif W == 3
        return online_anticor(X, K, T, H0)
    end
end
"""Iterates through all experts in one time period"""
function exp_gen(K, W, M, stocks, H0, abs_exp = false)
    T = size(stocks)[1]
    N = sum(K)
    H_mat = Matrix(N,M) .= 0
    mat_ind = 0
    for w in W[1] : W[2]
        for k in 1:K[w]
            mat_ind += 1
            if w == 3
                H_mat[mat_ind,:] = exp_alg(w,
                stocks[(end - min(2*k, T) +1 : end),:],
                min(k, T), T, H0[mat_ind,:])
            else
                H_mat[mat_ind,:] = exp_alg(w, stocks[(end - min(k, T) +1 : end),:],
                min(k, T), 0, 0)
            end
        end
    end
    return H_mat
end

function Simulation_Comp(dateInd, start_date, end_date,
    date_vec, X, K, W, eta, Learn_Alg, absPort, absExp)

    if dateInd == true
        start_ind = find(date_vec .== start_date)[1]
        end_ind = find(date_vec .== end_date)[1]
    else
        start_ind = start_date
        end_ind = end_date
    end


    T = end_ind - start_ind
    X = X[start_ind:end_ind, :]
    M = size(X)[2]
    N = sum(K)

    q_mat = Matrix(N,T)
    b_mat = Matrix(M,T)

    if absExp == true
        H_mat = [Matrix(N,M)]
        H_mat[1] .= 1/M
    else
        absPort = false
        H_mat = [Matrix(N,M)]
        H_mat[1] .= 0
    end

    if absPort == true
        b_mat[:, 1] = 1/M
        q_mat[:,1] = 1/N
    else
        b_mat[:, 1] = zeros(M)
        q_mat[:,1] = zeros(N)
    end

    S_t = Vector(T+1)
    S_t[1] = 1

    S_nt = Matrix(N,T+1)
    S_nt[:,1] = ones(N)

    for t in 2:T
        H1 = Exp_Gen(K, W, M, X[1:(t-1),:], H_mat[(t-1)],
         absExp)

        b_mat[: , t], q_mat[:, t], S_t[t], S_nt[:,t] =
        OnlineLearn_Comp(q_mat[:, (t-1)],b_mat[:, (t-1)],
         X[(t-1),:], H_mat[(t-1)], H1, S_t[t-1],
         S_nt[:, (t-1)] , Learn_Alg, eta, absPort)

        H_mat = push!(H_mat, H1)
    end

    S_t[end] = S_t[end-1] * (transpose(X[end,:]-1)
    * b_mat[:, end] + 1)

    S_nt[:,end] = S_nt[:,end-1] .* ((H_mat[end]
    * (X[end,:]-1)) + 1)

    return b_mat, q_mat, S_t, S_nt, H_mat
end
