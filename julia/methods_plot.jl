using Gadfly
using Colors

Gadfly.push_theme(Theme(major_label_font="Droid Sans", minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=16pt))
colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a", colorant"#984ea3",
          colorant"#ff7f00", colorant"#ffff33", colorant"#a65628", colorant"#f781bf"]
logs = ["rstdp", "dastdp", "ifstdp", "dmstdp", "ifdmstdp"]

function readlog(logname::String)
    res = readlines(open(logname))
    numres = Array{Float64}(length(res),length(split(chomp(split(res[1],"::")[2]),",")))
    for i in eachindex(res)
        numres[i,:] = map(x->parse(Float64,x), split(chomp(split(res[i],"::")[2]),","))
    end
    numres
end

function get_all(base::String)
    res = zeros(length(logs), 1500, 6, 20)
    for i=0:19
        for l in eachindex(logs)
            log = string(base, "/", i, "_", logs[l], ".log")
            numres = readlog(log)
            res[l, :, :, i+1] = numres
        end
    end
    res
end

function plot_all(base::String)
    res = get_all(base)
    meanres = mean(res, 4)[:,:,:,1]
    stdres = std(res, 4)[:,:,:,1]
    layers = []
    x = 1:size(res,2)
    for l in eachindex(logs)
        append!(layers,
                [layer(x=meanres[l,x,3], y=meanres[l,x,2],
                       # ymax=meanres[l,x,2]+stdres[l,x,2]*0.5,
                       # ymin=meanres[l,x,2]-stdres[l,x,2]*0.5,
                       Geom.line,
                       # Geom.ribbon(),
                       style(default_color=colors[l]))])
    end
    plt = plot(layers..., Guide.title("Stingray"),
               Guide.manual_color_key(
                   "",
                   ["R-STDP", "DA-STDP", "IF-STDP", "DM-STDP", "IFDM-STDP"],
                   colors[1:5]),
               Guide.xlabel("Time"),
               Guide.ylabel("Distance"))
    draw(PDF(string("plot/stingray.pdf"), 6inch, 4inch), plt)
    layers
end
