# Created by Eugene M. Izhikevich, February 25, 2003
# Translated to Julia by Dennis G Wilson, June 06, 2017
using Gadfly
Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        key_label_font="Droid Sans",
                        key_label_font_size=16pt,
                        major_label_font_size=18pt,
                        minor_label_font_size=16pt,
                        line_width=0.8mm,
                        point_size=0.3mm,
                        highlight_width=0.0mm,
                        default_color=colorant"#000000"))

# Excitatory neurons    Inhibitory neurons
Ne=800;                 Ni=200;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-65+15*re.^2;        -65*ones(Ni,1)];
d=[8-6*re.^2;           2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne)  -rand(Ne+Ni,Ni)];

v=-65*ones(Ne+Ni,1);    # Initial values of v
u=b.*v;                 # Initial values of u
firings=Array{Int64}(0,2);             # spike timings

for t=1:1000            # simulation of 1000 ms
    tinp=[5*randn(Ne,1);2*randn(Ni,1)]; # thalamic input
    fired=find(v.>=30);    # indices of spikes
    firings=[firings; [t+0*fired fired]];
    if length(fired)>0
        v[fired]=c[fired];
        u[fired]=u[fired]+d[fired];
        tinp=tinp+sum(S[:,fired],2);
    end
    v=v+0.5*(0.04*v.^2+5*v+140-u+tinp); # step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+tinp); # for numerical
    u=u+a.*(b.*v-u);                 # stability
end
plt = plot(x=firings[:,1]/1000.0, y=firings[:,2], Geom.point,
           Coord.cartesian(ymin=0, ymax=1000),
           Guide.xticks(ticks=collect(0:0.25:1.0)),
           Guide.yticks(ticks=collect(0:250:1000)),
           Guide.xlabel("Time [s]"), Guide.ylabel("Neuron"))
draw(PDF("timing.pdf", 8inch, 6inch), plt)
firings
