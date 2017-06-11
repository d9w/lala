# DA-STDP
# Dennis G Wilson, June 8, 2017
# Adapted from daspnet.m by Eugene M. Izhikevich
# https://www.izhikevich.org/publications/daspnet.m
using Logging
using Distances
@Logging.configure(level=INFO, filename="log")
method = length(ARGS) > 0 ? parse(Int64, ARGS[1]) : 0
method = 2
# method (0) STDP (1) R-STDP (2) DA-STDP (3) A-STDP (4) MS-STDP (5) AMS-STDP
seed = length(ARGS) > 1 ? parse(Int64, ARGS[2]) : 0
srand(seed)

M = 100   # number of post-synaptic connections per neuron
D = 1     # conduction delay (timesteps between pre and post firing)
Ne = 800  # excitatory neurons
Ni = 200  # inhibitory neurons
N = Ne+Ni # total number
a = [0.02*ones(Ne,1); 0.1*ones(Ni,1)] # in u part of model eq
d = [8*ones(Ne,1); 2*ones(Ni,1)] # in u part of model eq
sm = 4    # maximal synaptic strength
absorption = 0.001
const_r_signal = 1.0
stimulus_signal = 15
stimulus_interval = 10

if method == 4 || method == 5
    const_decay = 0.1 # dopamine decay
    max_delay = 10
    poses = randn(N, 3)*0.5
    dist = colwise(Euclidean(), poses', zeros(3, N))
    dist /= maximum(dist)
    decay = exp(-const_decay*dist)
    delay = Int64.(ceil(dist*max_delay))
end

post = [rand(1:N, Ne, M);rand(1:Ne, Ni, M)] # post-synaptic mapping
# s = [0.25*randn(Ne,M)+1.0;0.25*randn(Ni,M)-1.0]      # synaptic weights
s = [ones(Ne,M); -1.0*ones(Ni,M)] # synaptic weights
sd = zeros(N,M)                   # their derivatives
pre = Array{Array{Int64}}(N)
aux = Array{Array{Int64}}(N)

for i = 1:N
    pre[i] = find((post.==i) & (s.>0)) # pre-synaptic mapping
    aux[i] = findn((post.==i) & (s.>0))[1] # neuron of pre-synaptic mapping
end

group_secs = 4000    # time [s] for each group
secs = 2*group_secs   # the duration of experiment [s]
T = 1000          # timesteps per sec
if method == 4 || method == 5
    DA = zeros(N,1)   # level of dopamine above the baseline
else
    DA = 0.0
end
rew = zeros(Int64, 1)          # reward instances

STDP = zeros(N,1)
v = -65*ones(N,1) # membrane potential
u = 0.2.*v        # membrane recovery variable

# TODO: replace this with instrumental coding
groups = shuffle!(collect(1:Ne))
S_group = groups[1:50]
A_group = groups[51:100]
B_group = groups[101:150]

interval = 20     # the coincidence interval for n1 and n2
A_count = 0       # number of A group firings
B_count = 0       # number of B group firings
shist = zeros(Float64, secs, 2) # history of s over time

for sec=0:secs-1                      # simulation of 1 minute
    if sec % stimulus_interval == 0
        stimulus = rand(1:T)
    else
        stimulus = 0
    end
    A_count = 0
    B_count = 0
    for t=1:T                       # simulation of 1 sec
        tinp = 13*(rand(N,1)-0.5)   # random thalamic input
        if method == 3 || method == 5
            tinp += const_r_signal*DA
        end
        if t == stimulus
            tinp[S_group] += stimulus_signal
        end
        fired = find(v.>=30);       # indices of fired neurons
        if length(fired)>0
            v[fired] = -65            # reset v
            u[fired] = u[fired] + d[fired] # reset u
            STDP[fired] = 0.1    # increase STDP
            for k=1:length(fired)
                ind = post[fired[k],:] # post-synaptic indices
                # presynaptic STDP
                sd[pre[fired[k]]] += 1.0*STDP[aux[fired[k]]]
                # postsynaptic STDP
                sd[fired[k],:] -= 1.5*STDP[ind]
                # increase post-synaptic input
                tinp[ind] += s[fired[k], :]
            end
        end
        v = v+0.5*((0.04*v+5).*v+140-u+tinp)    # for numerical stability
        v = v+0.5*((0.04*v+5).*v+140-u+tinp)    # time step is 0.5 ms
        u = u+a.*(0.2*v-u)                      # tau = 20 ms
        STDP *= 0.95

        if method == 1
            DA = 0.0
        else
            DA -= DA*absorption
        end
        if (mod(t,10) == 0)
            if method == 0 || method == 3
                s[1:Ne,:] = max(0,min(sm,s[1:Ne,:]+0.1*sd[1:Ne,:]))
            elseif method == 1 || method == 2
                s[1:Ne,:] = max(0,min(sm,s[1:Ne,:]+(0.002+DA)*sd[1:Ne,:]))
            elseif method == 4 || method == 5
                for i=1:Ne
                    s[i,:] = max(0,min(sm,s[i,:]+(0.002+DA[i])*sd[i,:]))
                end
            end
            sd = 0.99*sd
        end
        if stimulus > 0
            if t >= stimulus && t <= stimulus + interval
                A_count += length(intersect(fired, A_group))
                B_count += length(intersect(fired, B_group))
            end
            if t == stimulus + interval
                if sec < group_secs
                    if A_count > B_count
                        append!(rew, [sec*T+t+ceil(T*rand())])
                    end
                else
                    if B_count > A_count
                        append!(rew, [sec*T+t+ceil(T*rand())])
                    end
                end
            end
        end
        if any(rew .== sec*T+t)
            if method < 4
                DA = DA + 0.5
            end
        end
        if method == 4 || method == 5
            rs = rew[end]+delay .== sec*T+t
            DA[rs] += decay[rs] * 0.5
        end
    end
    if sec % stimulus_interval == 0
        @info("S: ", @sprintf("%d,%d,%d,%d,%0.3f,%0.3f,%0.3f,%d",
                              sec, A_count, B_count, stimulus, mean(s),
                              mean(s[S_group,:][find(indexin(post[S_group, :], A_group))]),
                              mean(s[S_group,:][find(indexin(post[S_group, :], B_group))]),
                              rew[end]))

    end
end
nothing
