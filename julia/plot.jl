using Gadfly
using Distributions
using DataFrames

Gadfly.push_theme(Theme(major_label_font="Droid Sans",minor_label_font="Droid Sans",
                        major_label_font_size=18pt,minor_label_font_size=16pt,line_width=0.8mm,
                        key_label_font="Droid Sans",key_label_font_size=16pt))
# colors = [colorant"#fe9601", colorant"#86269b", colorant"#00d2f1", colorant"#00b796", colorant"#cc0063"]
colors = [colorant"#000000", colorant"#FF1300",colorant"#0086CE",colorant"#B0F700",
          colorant"#000000", colorant"#FF1300",colorant"#0375B2",colorant"#ADF300",
          colorant"#000000", colorant"#E31100",colorant"#066090",colorant"#98D600",
          colorant"#000000", colorant"#B0271C",colorant"#175170",colorant"#7EA61A",
          colorant"#000000", colorant"#75322D",colorant"#1F3B4B",colorant"#5B6F2A"]

function get_layer(base::String, l::Int64)

  xmax = 1000

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
  erange = (50*l):200:xmax

  println(maximum(nrmax))

  layers = Array{Gadfly.Layer,1}()

  append!(layers, layer(x=1:xmax, y=rmean, Geom.line, style(default_color=colors[l])))
  append!(layers, layer(x=erange, y=rmean[erange], ymin=(rmean-rstd./2)[erange], ymax=(rmean+rstd./2)[erange],
                    Geom.point, Geom.errorbar, style(default_color=colors[l])))
  layers
end

function plot_all(paths::Array{String}, labels::Array{String})
  layers = Array{Gadfly.Layer,1}()

  for p in eachindex(paths)
    append!(layers, get_layer(paths[p], p))
  end

  plt = plot(layers...,
             Guide.xlabel("Generation"), Guide.ylabel("Fitness"),
             Guide.manual_color_key("", labels, colors[1:length(paths)]));
  draw(PDF("cmaes.pdf", 8inch, 6inch), plt)
  plt
end

function plot_params(base::String, params::Array{Float64})
  labels = ["p<sub>$i</sub>" for i=0:4]

  nparams = params[[5, 4, 2, 1, 3],:]

  p = DataFrame()

  n = size(nparams,2)

  p[:values] = nparams[:]
  p[:labels] = repmat(labels, n)
  p[:mean] = repmat(mean(nparams,2), n)[:]

  plt = plot(layer(p, x="labels", y="mean", Geom.point, style(default_color=colorant"#000000")),
             layer(p, x="labels", y="values", Geom.violin, style(default_color=colorant"#000000")),
             Coord.cartesian(ymin=0.0, ymax=1.0),
             Guide.yticks(ticks=collect(0.0:0.1:1.0)),
             Guide.title(nothing), Guide.xlabel(nothing), Guide.ylabel(nothing));
  draw(PDF(string(base, "_params.pdf"), 6inch, 4inch), plt)
end

function get_params(base::String)

  allfits = zeros(2000, 20)
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

function best_params(base::String, pbest::Int64=10)

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

  plot_params(base, bestparams[2:end, :])

  allres, bestparams
end

function all_res(base::String)
    ress = Array{Array{Float64}}(20)

    for i=0:19
        ress[i+1] = readcsv("$base/$i.log")
    end

    allres = zeros(maximum(size.(ress))..., 20)
    for i=1:20
        si = size(ress[i])
        allres[1:si[1], 1:si[2], i] = ress[i]
    end

    allres
end

function best_config(base::String, defaults::String, shape::String, ranges::String, s::Int64=0)
    allres = all_res(base)
    minres = minimum(allres[:,3,:])
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
