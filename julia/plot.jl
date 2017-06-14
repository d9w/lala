using Gadfly
using Distributions
using DataFrames

Gadfly.push_theme(Theme(major_label_font="Droid Sans", minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=16pt))
colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a", colorant"#984ea3",
          colorant"#ff7f00", colorant"#ffff33", colorant"#a65628", colorant"#f781bf"]

function get_layer(base::String, l::Int64)

    xmax = 50

    allresults = zeros(xmax, 20)
    rmax = zeros(xmax, 20)

    for i=0:19
        res = readcsv("$base/$i.log")
        for j=1:min(size(res,1), xmax)
            allresults[j, i+1] = res[j, 3]
        end
    end

    for i=1:20
        for j=1:xmax
            rmax[j,i] = -minimum(allresults[1:j,i])
        end
    end
    allresults, rmax

    # nrmax = flipdim(sortcols(flipdim(rmax,1), rev=true)[:, 15:20], 1)
    nrmax = rmax

    rmean = mean(nrmax, 2)
    rstd = std(nrmax, 2)

    println(maximum(nrmax))
    layers = Array{Gadfly.Layer,1}()
    append!(layers, layer(x=1:xmax, y=rmean, ymin=rmean-0.5*rstd, ymax=rmean+0.5*rstd,
           Geom.line, Geom.ribbon, style(default_color=colors[l])))
    layers
end

function plot_all(paths::Array{String}, labels::Array{String}, title::String="", color_key::Bool=false)
    layers = Array{Gadfly.Layer,1}()

    for p in eachindex(paths)
        append!(layers, get_layer(paths[p], p))
    end

    plt = plot(layers...,
               Guide.xticks(ticks=[0:10:50;]),
               Guide.xlabel("Generation"), Guide.ylabel(nothing),
               Guide.title(title),
               Guide.manual_color_key("", labels, colors[1:length(paths)]),
               style(key_position=:bottom));
    draw(PDF("cmaes.pdf", 8inch, 6inch), plt)
    plt
end

function plot_cmaes(quad_paths::Array{String}, sting_paths::Array{String})
    labels = ["R-STDP", "DA-STDP", "IFDM-STDP"];

    quad_layers = Array{Gadfly.Layer,1}()
    for p in eachindex(quad_paths)
        append!(quad_layers, get_layer(quad_paths[p], p))
    end

    sting_layers = Array{Gadfly.Layer,1}()
    for p in eachindex(sting_paths)
        append!(sting_layers, get_layer(sting_paths[p], p))
    end

    quad_plt = plot(quad_layers...,
                    Guide.xticks(ticks=[0:10:50;]),
                    Guide.xlabel("Generation"), Guide.ylabel("Fitness"),
                    Guide.title("Quadropus"),
                    # Coord.cartesian(aspect_ratio=2.0),
                    Guide.manual_color_key("", labels, colors[1:length(labels)]),
                    style(key_position=:right));
    sting_plt = plot(sting_layers...,
                     Guide.xticks(ticks=[0:10:50;]),
                     Guide.xlabel("Generation"), Guide.ylabel("Fitness"),
                     Guide.title("Stingray"),
                    # Coord.cartesian(aspect_ratio=2.0),
                     Guide.manual_color_key("", labels, colors[1:length(labels)]),
                    style(key_position=:right));
    plt = hstack(quad_plt, sting_plt)
    draw(PDF("evo.pdf", 14inch, 4inch), plt)
    plt
end

function plot_params(base::String, params::Array{Float64}, title::String)
    labels = ["p<sub>rs</sub>", "p<sub>df</sub>", "p<sub>dd</sub>",
              "p<sub>dat</sub>", "p<sub>dab</sub>"]

    nparams = params[[6, 5, 3, 2, 4],:]

    p = DataFrame()

    n = size(nparams,2)

    p[:values] = nparams[:]
    p[:labels] = repmat(labels, n)
    p[:mean] = repmat(mean(nparams,2), n)[:]
    p[:max] = repmat(nparams[:,1], n)[:]
    p[:rank] = repmat([1:n;]', length(labels))[:]

    plt = plot(
        layer(p, x="labels", y="max", Geom.point,
              style(default_color=colorant"#000000")),
        layer(p, x="labels", y="values", Geom.violin,
              style(default_color=colorant"#000000")),
        Guide.ylabel("Distribution"),
        Guide.yticks(ticks=collect(-0.5:0.5:1.5)),
        Guide.title(title), Guide.xlabel(nothing));
    draw(PDF(string(base, "_params.pdf"), 6inch, 4inch), plt)
    plt
end

function get_params(base::String)

  allfits = zeros(50, 20)
  allparams = zeros(6, 20)

  for i=0:19
    res = readcsv("$base/$i.log")
    for j=1:size(res,1)
      allfits[j, i+1] = res[j, 3]
    end
    maxf = findmin(allfits[:, i+1])
    allparams[1, i+1] = -maxf[1]
    allparams[2:end, i+1] = map(x->mod(x, 1.0), res[maxf[2], 4:end])
  end

  allparams[1, :] /= maximum(allparams[1, :])

  plot_params(base, sortcols(allparams, rev=true)[2:end,1:10])

  allparams
end

function best_params(base::String, title::String="", pbest::Int64=10)

  allres = zeros(40000, 6)

  for i=0:19
    res = readcsv("$base/$i.log")
    for j=1:size(res,1)
      allres[2000*i+j, 1] = -res[j, 3]
      allres[2000*i+j, 2:end] = map(x->mod(x, 1.0), res[j,4:end])
    end
  end

  bestparams = sortrows(allres, rev=true)[1:pbest, :]'

  bestparams[1, :] /= maximum(bestparams[1, :])

  plt = plot_params(base, bestparams, title)

    println(bestparams[:,1])

  allres, bestparams, plt
end

function all_res(base::String)
    ress = Array{Array{Float64}}(20)

    for i=0:19
        ress[i+1] = readcsv("$base/$i.log")
    end

    allres = zeros(50, maximum(size.(ress,2)), 20)
    for i=1:20
        si = size(ress[i])
        allres[1:50, 1:si[2], i] = ress[i][1:50, :]
    end

    allres
end

function best_config(base::String, defaults::String, shape::String, ranges::String, s::Int64=0)
    allres = all_res(base)
    minres = minimum(allres[:,3,:])
    println(minres)
    minind = findn(allres.==minres)
    a = allres[minind[1],4:end,minind[3]]

    defaults = JSON.parsefile(defaults)
    shapedefaults = JSON.parsefile(shape)
    ranges = JSON.parsefile(ranges)

    for i in eachindex(a)
        if (a[i] < 0.0) || (a[i] > 1.0)
            a[i] = mod(a[i], 1.0)
        end
    end
    println(a)
    i = 1
    for (k, v) in ranges
        defaults[k] = v[1]+(v[2]-v[1])*a[i]
        i += 1
    end
    for (k, v) in shapedefaults
        defaults[k] = v
    end
    defaults["seed"] = s
    defaults
end
