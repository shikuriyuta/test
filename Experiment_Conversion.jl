using Distributed
#addprocs(10)
@everywhere dir_ = "C:/Users/shiku/Work/Paper/ICML2021/" # directory
data = "water" # dataset
m = 5 # maximum size of parent set
N = 10000 # sample size
greedy = false # algorithm1
iter_n = 10 # number of samplings

using Statistics, DataFrames, Feather, CSV, Dates
@everywhere include(dir_ * "Conversion.jl")
lambda_sum_list = Vector{Int64}([])
lambda_pq_sum_list = Vector{Int64}([])
space_sum_list = Vector{Int64}([])
time1_sum_list = Vector{Millisecond}([])
time2_sum_list = Vector{Millisecond}([])
for num in 1:iter_n
    println("num=", string(num), " start")
    D = CSV.read(dir_ * data * "/" * data * "_" * string(N) * "_" * string(num) * ".csv", DataFrame)
    output_list, lambda_sum, lambda_pq_sum, space_sum, time1_sum, time2_sum = EC(D, m, 1.0, greedy)
    CSV.write(dir_ * data * "/" * data * "_" * string(N) * "_" * string(greedy) * "_" * string(num) * ".csv", output_list)
    push!(lambda_sum_list, lambda_sum)
    push!(lambda_pq_sum_list, lambda_pq_sum)
    push!(space_sum_list, space_sum)
    push!(time1_sum_list, time1_sum)
    push!(time2_sum_list, time2_sum)
    println("num=", string(num), " completed (lambda_sum=", lambda_sum, ", lambda_pq_sum=", lambda_pq_sum, ", space_sum=", space_sum, ", time1_sum=", time1_sum, ", time2_sum=", time2_sum, ")")
end
println()
println("Result (", data, ", N=", N, ", greedy=", greedy, ")")
println("lambda_sum: mean ", mean(lambda_sum_list), " std ", std(lambda_sum_list))
println("lambda_pq_sum: mean ", mean(lambda_pq_sum_list), " std ", std(lambda_pq_sum_list))
println("space_sum: mean ", mean(space_sum_list), " std ", std(space_sum_list))
println("time1_sum: mean ", Dates.value(sum(time1_sum_list)) / iter_n)
println("time2_sum: mean ", Dates.value(sum(time2_sum_list)) / iter_n)