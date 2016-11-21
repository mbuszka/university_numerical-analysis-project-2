include("spline.jl")
using PlotlyJS

function test_sin()
  n = 10
  xs = linspace(0, 2π, n+1)
  ys = sin(xs)
  ds = cos(xs[2:end])
  a, b, c, g, h, f = build_matrix(xs, ys)
  s, k = spline_cyclic(xs, ys)
  for i in 2 : n-1
    @printf("%e\t%e\t%e\n", c[i-1] * k[i-1] + a[i] * k[i] + b[i] * k[i+1] - g[i], 
            c[i-1] * ds[i-1] + a[i] * ds[i] + b[i] * ds[i+1] - g[i],
            k[i] - ds[i])
  end
  @printf("%e\t%e\t%e\n", c[n-1] * k[n-1] + a[n] * k[n] + b[n] * k[1] - g[n], 
          c[n-1] * ds[n-1] + a[n] * ds[n] + b[n] * ds[1] - g[n],
          k[n] - ds[n])
  a, b = reduce_matrix(a, b, c)
  c = b
  println("-----------------\n")
  for i in 2 : n-1
    @printf("%e\t%e\t%e\n", c[i-1] * k[i-1] + a[i] * k[i] + b[i] * k[i+1] - g[i], 
            c[i-1] * ds[i-1] + a[i] * ds[i] + b[i] * ds[i+1] - g[i],
            k[i] - ds[i])
  end
  @printf("%e\t%e\t%e\n", c[n-1] * k[n-1] + a[n] * k[n] + b[n] * k[1] - g[n], 
          c[n-1] * ds[n-1] + a[n] * ds[n] + b[n] * ds[1] - g[n],
          k[n] - ds[n])
end

function test_k()
  n = 10
  xs = linspace(0, 2π, n+1)
  ys = cos(xs)
  dds = -cos(xs[2:end])
  k = spline_cyclic_k(xs, ys)
  @show k
  plot([scatter(; x = xs[2:end], y = k), scatter(; x = xs[2:end], y = dds)])
end


function test()
  f = cos
  xs = linspace(-π, π, 11)
  ys = f(xs)
  xss = linspace(-π, π, 1000)

  s, k = spline_cyclic_k(xs, ys, debug = true)
  tss = [ s(x) for x = xss ]
  yss = [ x for (x, d, dd) = tss ]
  dss = [ d for (x, d, dd) = tss ]
  ddss = [ dd for (x, d, dd) = tss ]
  plot([ scatter(;x = xss, y = yss), scatter(;x = xss, y = dss)
        ,scatter(;x = xss, y = ddss), scatter(;x = xss, y = f(xss))
        ,scatter(;x = xs[2:end], y = k)])
end
