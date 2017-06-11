using Gadfly

Gadfly.push_theme(Theme(major_label_font="Droid Sans", minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        key_label_font_size=16pt))
colors = [colorant"#e66101", colorant"#fdb863", colorant"#b2abd2", colorant"#5e3c99",
          colorant"#b35806", colorant"#e08214", colorant"#fdb863", colorant"#fee0b6",
          colorant"#d8daeb", colorant"#b2abd2", colorant"#8073ac", colorant"#542788"]

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
    meanres = mean(allres, 3)[:,:,1]
    stdres = std(allres, 3)[:,:,1]
    x = Int64.(ceil.(linspace(1, size(allres,1), 100)))
    layers = [layer(x=x, y=meanres[x,1], ymin=meanres[x,1]-stdres[x,1]*0.5,
                    ymax=meanres[x,1]+stdres[x,1]*0.5, Geom.line, Geom.ribbon,
                    style(default_color=colors[color_i+1])),
              layer(x=x, y=meanres[x,3], ymin=meanres[x,3]-stdres[x,3]*0.5,
                    ymax=meanres[x,3]+stdres[x,3]*0.5, Geom.line, Geom.ribbon,
                    style(default_color=colors[color_i+2]))]
    plt = plot(layers..., Guide.title(title), Guide.xlabel("Time [s]"),
               Guide.ylabel("Synaptic weight [mV]"))
    draw(PDF(string(title, ".pdf"), 6inch, 4inch), plt)
    layers
end
