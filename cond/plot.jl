using Gadfly
using Colors

Gadfly.push_theme(Theme(major_label_font="Droid Sans", minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=16pt))
colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a", colorant"#984ea3",
          colorant"#ff7f00", colorant"#ffff33", colorant"#a65628", colorant"#f781bf"]

function read_res(base::String)

    ress = Array{Array{Float64}}(20)

    for i=0:19
        ress[i+1] = readcsv("$base/s_$i.log")
    end

    allres = zeros(maximum(size.(ress))..., 20)

    for i=1:20
        si = size(ress[i])
        allres[1:si[1], 1:si[2], i] = ress[i]
    end

    allres
end

function plot_single(allres::Array{Float64}, title::String, color_i=0)
    a_gt_b = sum(allres[:,2,:] .> allres[:,3,:], 2)/20.0
    b_gt_a = sum(allres[:,3,:] .> allres[:,2,:], 2)/20.0
    meanres = mean(allres, 3)[:,:,1]
    stdres = std(allres, 3)[:,:,1]
    # x = Int64.(ceil.(linspace(1, size(allres,1), 100)))
    x = 1:size(allres,1)
    layers = [layer(x=meanres[x,1]/10, y=a_gt_b[x], Geom.smooth(method=:loess,smoothing=0.1),
                    style(default_color=colorant"red")),#colors[color_i+1])),
              layer(x=meanres[x,1]/10, y=b_gt_a[x], Geom.smooth(method=:loess,smoothing=0.1),
                    style(default_color=colorant"blue"))]#colors[color_i+2]))]
    plt = plot(layers..., Guide.title(title), Guide.xlabel("Time [s]"),
               Guide.ylabel("Synaptic weight [mV]"))
    draw(PDF(string("plot/",title, ".pdf"), 6inch, 4inch), plt)
    layers
end

function plot_weight(allres::Array{Float64}, title::String, color_i=0)
    meanres = mean(allres, 3)[:,:,1]
    stdres = std(allres, 3)[:,:,1]
    # x = Int64.(ceil.(linspace(1, size(allres,1), 100)))
    x = 1:size(allres,1)
    layers = [layer(x=meanres[x,1]/(10), y=meanres[x,6],
                    ymax=meanres[x,6]+stdres[x,6]*0.5,
                    ymin=meanres[x,6]-stdres[x,6]*0.5,
                    Geom.smooth(method=:loess,smoothing=0.1),
                    Geom.ribbon(),
                    style(default_color=colors[color_i+1])),
              layer(x=meanres[x,1]/(10), y=meanres[x,7],
                    ymax=meanres[x,7]+stdres[x,7]*0.5,
                    ymin=meanres[x,7]-stdres[x,7]*0.5,
                    Geom.smooth(method=:loess,smoothing=0.1),
                    Geom.ribbon(),
                    style(default_color=colors[color_i+2]))]
    plt = plot(layers..., Guide.title(title), Guide.xlabel("Episodes [10s]"),
               Guide.ylabel("Synaptic weight [mV]"))
    draw(PDF(string("plot/",title, ".pdf"), 6inch, 4inch), plt)
    layers
end

function plot_all(pltf::Function)
    stdp_layers = pltf(read_res("results/471926"), "STDP")
    r_stdp_layers = pltf(read_res("results/471926"), "R-STDP")
    da_stdp_layers = pltf(read_res("results/471920"), "DA-STDP")
    a_stdp_layers = pltf(read_res("results/471986"), "IF-STDP")
    dm_stdp_layers = pltf(read_res("results/471936"), "DM-STDP")
    ifdm_stdp_layers = pltf(read_res("results/471940"), "IFDM-STDP")

    ymin = 1.0
    ymax = 1.8

    function plt_layers(layers, xlabel::String)
        plot(layers..., xintercept=[400.0], Guide.xlabel(xlabel), Guide.ylabel(nothing),
             Coord.cartesian(xmin=0, xmax=800, ymin=ymin, ymax=ymax),
             Guide.yticks(ticks=collect(ymin:0.2:ymax), orientation=:horizontal),
             Geom.vline(color=colorant"gray", size=0.01mm),
             Guide.xticks(ticks=[0, 200, 400, 600, 800], orientation=:horizontal))
    end

    plt = hstack(
        plot(r_stdp_layers..., xintercept=[400.0],
             Guide.xlabel("R-STDP"), Guide.ylabel("Synaptic weight [mV]"),
             Coord.cartesian(xmin=0, xmax=800, ymin=ymin, ymax=ymax),
             Guide.yticks(ticks=collect(ymin:0.2:ymax), orientation=:horizontal),
             Geom.vline(color=colorant"gray", size=0.01mm),
             Guide.xticks(ticks=[0, 200, 400, 600, 800], orientation=:horizontal)),
        plt_layers(da_stdp_layers, "DA-STDP"),
        plt_layers(a_stdp_layers, "IF-STDP"),
        plt_layers(dm_stdp_layers, "DM-STDP"),
        plt_layers(ifdm_stdp_layers, "IFDM-STDP"))
    draw(PDF(string("cond_stdp.pdf"), 20inch, 4inch), plt)
end
