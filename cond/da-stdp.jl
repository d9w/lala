# DA-STDP
# Dennis G Wilson, June 8, 2017
# Adapted from daspnet.m by Eugene M. Izhikevich
# https://www.izhikevich.org/publications/daspnet.m

seed = length(ARGS) > 0 ? parse(Int64, ARGS[1]) : 0
srand(seed)

M = 100   # number of post-synaptic connections per neuron
D = 1     # conduction delay (timesteps between pre and post firing)
Ne = 800  # excitatory neurons
Ni = 200  # inhibitory neurons
N = Ne+Ni # total number
a = [0.02*ones(Ne,1); 0.1*ones(Ni,1)] # in u part of model eq
d = [8*ones(Ne,1); 2*ones(Ni,1)] # in u part of model eq
sm = 4    # maximal synaptic strength

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

secs = 120         # the duration of experiment [s]
T = 1000          # timesteps per sec
DA = 0            # level of dopamine above the baseline (use neuron-specific)
rew = Array{Int64}(0)          # reward instances (use array with zero norm)

STDP = zeros(N,1)
v = -65*ones(N,1) # membrane potential
u = 0.2.*v        # membrane recovery variable

# TODO: replace this with instrumental coding
n1 = 1            # presynaptic neuron
syn = 1           # the synapse number to the postsynaptic neuron
n2 = post[n1,syn] # postsynaptic neuron
s[n1,syn] = 1.0   # start with 0 value

interval = 20     # the coincidence interval for n1 and n2
n1f = -100      # the last spike of n1
shist = zeros(Float64, secs, 2) # history of s over time

for sec=0:secs-1                      # simulation of 1 minute
    for t=1:T                       # simulation of 1 sec
        tinp = 13*(rand(N,1)-0.5)   # random thalamic input
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

        DA -= DA*0.001;
        if (mod(t,10) == 0)
            s[1:Ne,:] = max(0,min(sm,s[1:Ne,:]+(0.002+DA)*sd[1:Ne,:]));
            sd = 0.99*sd;
        end;
        if any(fired .== n1)
            n1f = sec*T+t
        end
        if any(fired .== n2)
            if (sec*T+t-n1f<interval) & ((sec*T+t)>n1f)
                append!(rew, [sec*T+t+T+ceil(2*T*rand())])
            end
        end
        if any(rew .== sec*T+t)
            DA = DA+0.5
        end
    end
    shist[sec+1, :] = [s[n1, syn] sd[n1, syn]]
    println("S: ", shist[sec+1, :])
end
println("R: ", rew)
nothing
