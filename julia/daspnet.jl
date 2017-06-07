# daspnet.m: Spiking network with DA-modulated STDP. Based on spnet.m
# Created by Eugene M.tinpzhikevich.                November 3, 2004
# Translated to Julia by Dennis G Wilson, June 6, 2017
#
# This program reproduces the experiment in Fig.1 in
# tinpzhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage
# of STDP and Dopamine Signaling. Cerebral Cortex, 10.1093/cercor/bhl152
#
# n1 - the presynaptic neuron. syn is the synapse to be reinforced.
# Plot: top - spike raster. Bottom left - synaptic strength (blue), the
# eligibility trace (green), and the rewards (red x). Bottom right - the
# distribution of synaptic weights with the chosen synapse marked by red dot.
using Gadfly
Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        key_label_font="Droid Sans",
                        key_label_font_size=16pt,
                        major_label_font_size=18pt,
                        minor_label_font_size=16pt,
                        line_width=0.8mm,
                        point_size=0.5mm,
                        highlight_width=0.0mm,
                        panel_fill=colorant"#ffffff",
                        default_color=colorant"#000000"))
colors = [colorant"#e66101", colorant"#fdb863", colorant"#b2abd2", colorant"#5e3c99",
          colorant"#b35806", colorant"#e08214", colorant"#fdb863", colorant"#fee0b6",
          colorant"#d8daeb", colorant"#b2abd2", colorant"#8073ac", colorant"#542788"]

M=100;                 # number of synapses per neuron
D=1;                   # maximal conduction delay
# excitatory neurons   # inhibitory neurons      # total number
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 # maximal synaptic strength

post=Int64.(ceil([N*rand(Ne,M);Ne*rand(Ni,M)]))
s=[ones(Ne,M);-ones(Ni,M)];         # synaptic weights
sd=zeros(N,M);                      # their derivatives
delays=Array{Any}(N,D)
pre=Array{Any}(N)
aux=Array{Any}(N)
for i=1:N
    if i<=Ne
        for j=1:D
            delays[i,j]=Int64(M/D*(j-1))+(1:Int64(M/D));
        end
    else
        delays[i,1]=1:M;
    end
    pre[i]=find((post.==i) & (s.>0));             # pre excitatory neurons
    aux[i]=Int64.(N*(D-1-ceil(ceil(pre[i]/N)/(M/D)))+1+mod(pre[i]-1,N))
end
STDP = zeros(N,N+1+D);
v = -65*ones(N,1);                      # initial values
u = 0.2.*v;                             # initial values
firings=[-D 0];                         # spike timings

#---------------
# new stuff related to DA-STDP
T=3600;         # the duration of experiment
DA=0;           # level of dopamine above the baseline
rew=[];

n1=1;           # presynaptic neuron
syn=1;          # the synapse number to the postsynaptic neuron
n2=post[n1,syn] # postsynaptic neuron
s[n1,syn]=0;    # start with 0 value

interval = 20;  # the coincidence interval for n1 and n2
n1f=[-100];       # the last spike of n1
n2f=[];         # the last spike of n2
shist=Array{Float64}(0,2)#zeros(1000*T,2);
#--------------

for sec=1:T                             # simulation of 1 day
    for t=1:1000                          # simulation of 1 sec
        tinp=13*(rand(N,1)-0.5);               # random thalamic input
        fired = find(v.>=30);                # indices of fired neurons
        if length(fired)>0
            v[fired]=-65;
            u[fired]=u[fired]+d[fired];
            STDP[fired,t+D]=0.1;
            for k=1:length(fired)
                sd[pre[fired[k]]]+=STDP[N*t+aux[fired[k]]];
            end
            firings=[firings; Int64.([t*ones(length(fired),1) fired])];
            k=size(firings,1);
            while firings[k,1]>t-D
                del = delays[firings[k,2],t-firings[k,1]+1];
                ind = post[firings[k,2],del];
                tinp[ind]=tinp[ind]+s[firings[k,2], del];
                sd[firings[k,2],del]=sd[firings[k,2],del]-1.5*STDP[ind,t+D];
                k=k-1;
            end
        end
        v=v+0.5*((0.04*v+5).*v+140-u+tinp);    # for numerical
        v=v+0.5*((0.04*v+5).*v+140-u+tinp);    # stability time
        u=u+a.*(0.2*v-u);                   # step is 0.5 ms
        STDP[:,t+D+1]=0.95*STDP[:,t+D];     # tau = 20 ms

        DA=DA*0.995;
        if (mod(t,10)==0)
            s[1:Ne,:]=max(0,min(sm,s[1:Ne,:]+(0.002+DA)*sd[1:Ne,:]));
            sd=0.99*sd;
        end;
        if any(fired.==n1)
            append!(n1f, [sec*1000+t])
        end
        if any(fired.==n2)
            append!(n2f, [sec*1000+t])
            if (sec*1000+t-n1f[end]<interval) & (n2f[end]>n1f[end])
                append!(rew, [sec*1000+t+1000+ceil(2000*rand())])
            end
        end
        if any(rew.==sec*1000+t)
            DA=DA+0.5;
        end
        shist = [shist; [s[n1, syn] sd[n1, syn]]]
    end;
    # ---- plot -------
    plt_firings = plot(x=firings[:,1]/1000.0, y=firings[:,2], Geom.point,
                       Coord.cartesian(ymin=0, ymax=1000),
                       Guide.xticks(ticks=collect(0:0.25:1.0)),
                       Guide.yticks(ticks=collect(0:250:1000)),
                       Guide.xlabel("Time [s]"),
                       Guide.ylabel("Neuron"));
    plt_s = plot(x=s[s.>0], Geom.histogram,
                 Guide.xticks(ticks=collect(0.0:0.5:4)),
                 Guide.xlabel("Synaptic weight [mV]"),
                 Guide.ylabel("Count"));
    xrang=1:((sec-1)*1000+t)
    xticks=0.001*collect(0:10000:((sec-1)*1000+t))
    rew_layer = layer(x=0.001*rew, y=0*rew, Geom.point, shape=[Gadfly.cross],
                      style(point_size=1.0mm, default_color=colors[1]))
    plt_shist1 = plot(layer(x=0.001*xrang, y=shist[xrang,1], Geom.line),
                      rew_layer,
                      Guide.xticks(ticks=xticks, label=false),
                      Guide.xlabel(nothing),
                      Guide.ylabel("Synaptic weight [mV]"));
    plt_shist2 = plot(layer(x=0.001*xrang, y=shist[xrang,2], Geom.line),
                      rew_layer,
                      Guide.xticks(ticks=xticks),
                      Guide.xlabel("Time [s]"),
                      Guide.ylabel("Dopamine", orientation=:vertical));
    plt = vstack(hstack(plt_s, plt_firings), vstack(plt_shist1, plt_shist2));
    draw(PNG(string("julia/plot/",sec,".png"), 16inch, 12inch), plt)
    # ---- end plot ------
    STDP[:,1:D+1]=STDP[:,(N+1):(N+1+D)];
    ind = find(firings[:,1] .> N+1-D);
    firings=[-D 0; [firings[ind,1]-N firings[ind,2]]];
end;
