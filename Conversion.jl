using Distributed
using DataFrames
using Combinatorics
using SpecialFunctions
using Dates

function score(parents::Set{Int64}, i::Int64, D::DataFrame, ESS::Float64)
    num_i_states = length(Set(D[:,i]))
    state_counts = combine(DataFrames.groupby(D, names(D)[[j for j in parents]]), names(D)[i] => x -> [sum(x.==s) for s in Set(D[:,i])])
    state_counts = convert(Vector{Int64}, state_counts[:,ncol(state_counts)])
    num_parents_states = Int64(length(state_counts) / num_i_states)
    state_counts = reshape(state_counts, (num_i_states, num_parents_states))
    alpha = ESS / Float64(num_parents_states)
    beta = ESS / (Float64(num_parents_states) * Float64(num_i_states))
    score_ = 0.0
    score_ = score_ - sum([loggamma(state_counts[k,j] + beta) for j in 1:num_parents_states for k in 1:num_i_states])
    score_ = score_ + sum([loggamma(sum(state_counts[:,j]) + alpha) for j in 1:num_parents_states])
    score_ = score_ - Float64(num_parents_states) * loggamma(alpha) 
    score_ = score_ + Float64(num_parents_states) * Float64(num_i_states) * loggamma(beta) 
    return score_
end

function PSI(i::Int64, D::DataFrame, m::Int64, ESS::Float64, greedy::Bool, omega::Int64)
    X_ = Vector{Int64}([j for j in 1:ncol(D) if j != i])
    if greedy == false
        W_list = Vector{Set{Int64}}([Set{Int64}(W) for m_ in 0:m for W in collect(combinations(X_, m_))])
        W_list_ = copy(W_list)
        Score_list = Vector{Float64}([score(W, i, D, ESS) for W in W_list_])
        Score_list_ = copy(Score_list)
        space = length(W_list)
        for (h_, W_) in enumerate(W_list_)
            for (h, W) in Iterators.reverse(enumerate(W_list))
                if (issubset(W_, W) == true) & (issetequal(W_, W) == false) & (Score_list_[h_] <= Score_list[h])
                    deleteat!(W_list, h)
                    deleteat!(Score_list, h)
                end
            end
        end        
    else
        W_list = Vector{Set{Int64}}([Set{Int64}([])])
        Score_list = Vector{Float64}([score(Set{Int64}([]), i, D, ESS)])
        W_old = Vector{Set{Int64}}([Set{Int64}([])])
        space = 1
        for m_ in 1:m
            W_new_ = Set{Set{Int64}}([])
            for W in W_old
                for W1 in setdiff(X_, W)
                    W_new_ = W_new_ ∪ Set{Set{Int64}}([W ∪ Set{Int64}([W1])])
                end
            end
            W_new_ = Vector{Set{Int64}}([W for W in W_new_])
            W_new = Vector{Set{Int64}}([])
            Score_new = Vector{Float64}([])
            space = space + length(W_new_)
            for W in W_new_
                Score = score(W, i, D, ESS)
                for (h, W_) in enumerate(W_list)
                    if (issubset(W_, W) == true) & (Score_list[h] <= Score)
                        break
                    elseif h == length(W_list)
                        push!(W_new, W)
                        push!(Score_new, Score)
                    end
                end
            end
            if 0 < omega < length(Score_new)
                Score = sort(Score_new)[omega + 1]
                W_new = Vector{Set{Int64}}([W for (h, W) in enumerate(W_new) if Score_new[h] < Score])
                Score_new = Vector{Float64}([Score_new[h] for h in 1:length(Score_new) if Score_new[h] < Score])
            end
            append!(W_list, W_new)
            append!(Score_list, Score_new)
            W_old = copy(W_new)
        end
        for (h, W) in Iterators.reverse(enumerate(W_list))
            if length(W) > 1
                for W_ in Vector{Set{Int64}}([Set{Int64}(W_) for m_ in 1:(length(W)-1) for W_ in collect(combinations(Vector{Int64}([W1 for W1 in W]), m_))])
                    if score(W_, i, D, ESS) <= Score_list[h]
                        deleteat!(W_list, h)
                        break
                    end
                end
            end
        end
    end
    return Set{Set{Int64}}(W_list), space
end

function encode(W_list::Set{Set{Int64}}, i::Int64, D::DataFrame, ESS::Float64, mu::Int64)
    Y = Set{Int64}([W_ for W in W_list for W_ in W])
    lambda12_2 = length(W_list) - 1
    U_list = copy(W_list)
    V_list = Set{Set{Int64}}([Set{Int64}([])])
    Z = Set{Int64}([])
    Z_best = Set{Int64}([])
    while length(setdiff(Y, Z)) > 0
        mu_ = min(floor(Int64, length(setdiff(Y, Z))), mu)
        Y_ = Vector{Int64}([W1 for W1 in setdiff(Y, Z)])
        for Z_ in Vector{Set{Int64}}([Set{Int64}(Z_) for m_ in 1:mu_ for Z_ in collect(combinations(Y_, m_))])
            ZZ_ = Z ∪ Z_
            U_list_ = Set{Set{Int64}}([W ∩ ZZ_ for W in W_list])
            V_list_ = Set{Set{Int64}}([W ∩ setdiff(Y, ZZ_) for W in W_list])
            UV_list = Set{Set{Int64}}([U ∪ V for U in U_list_ for V in V_list_])
            if (issubset(W_list, UV_list) == true) & ((length(U_list_) + length(V_list_) - 2) < lambda12_2)
                lambda12_2 = length(U_list_) + length(V_list_) - 2
                U_list = copy(U_list_)
                V_list = copy(V_list_)
                Z_best = copy(ZZ_)
            end
        end
        if Z_best == Z
            break
        else 
            Z = copy(Z_best)
        end
    end
    return U_list, V_list
end

function worker(i::Int64, D::DataFrame, m::Int64, ESS::Float64, greedy::Bool, omega::Int64, mu::Int64)
    time0 = now()
    W_list, space = PSI(i, D, m, ESS, greedy, omega)
    time1 = now() - time0
    time0 = now()
    U_list, V_list = encode(W_list, i, D, ESS, mu)
    time2 = now() - time0
    println("i=", i, " completed (lambda=", length(W_list) - 1, ", lambda_pq=", length(U_list) + length(V_list) - 2, ", space=", space, ")")
    return [W_list, U_list, V_list, space, time1, time2]
end

function EC(D::DataFrame, m::Int64=ncol(D)-1, ESS::Float64=1.0, greedy::Bool=false, omega::Int64=0, mu::Int64=1)
    sort!(D)
    result_list = pmap(i -> worker(i, D, m, ESS, greedy, omega, mu), 1:ncol(D))
    col = vcat(vcat(vcat(["i"], [D_ * "_" * string(ii) for ii in 1:2 for D_ in names(D)]), [D_ for D_ in names(D)]), ["score"])
    output_list = DataFrame([(if (i == (1 + 3 * ncol(D) + 1)) Float64[] else Int64[] end) for i in 1:(1 + 3 * ncol(D) + 1)], col)
    for (i, result) in enumerate(result_list)
        for U in result[2]
            hoge = vcat([i], [(if j in U 1 else 0 end) for j in 1:ncol(D)])
            for V in result[3]
                push!(output_list, vcat(vcat(vcat(hoge, [(if j in V 1 else 0 end) for j in 1:ncol(D)]), [(if j in (U ∪ V) 1 else 0 end) for j in 1:ncol(D)]), [score(U ∪ V, i, D, ESS)]))
            end
        end
    end
    lambda_sum = sum([length(result[1]) - 1 for result in result_list])
    lambda_pq_sum = sum([length(result[2]) - 1 for result in result_list]) + sum([length(result[3]) - 1 for result in result_list])
    space_sum = sum([result[4] for result in result_list])
    time1_sum = sum([result[5] for result in result_list])
    time2_sum = sum([result[6] for result in result_list])
    return output_list, lambda_sum, lambda_pq_sum, space_sum, time1_sum, time2_sum
end