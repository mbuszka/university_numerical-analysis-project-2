# Maciej Buszka

include("spline.jl")
using Plotly

function test_sin(n = 10)
  test_fun(sin, linspace(0, 2π, n))
end

function test_cos(n = 10)
  test_fun(cos, linspace(0, 2π, n))
end

function randspace(from, to, steps)
  ts = zeros(typeof(to), steps)
  c  = zero(typeof(to))
  for i in 2 : steps
    c += rand(typeof(to))
    ts[i] = c
  end
  for i in 1 : steps
    ts[i] *= (to - from) / c
    ts[i] += from
  end
  ts
end

function test_fun(f, xs)
  ys  = [ f(x) for x in xs ]
  xss = linspace(xs[1], xs[end], size(xs, 1) * 10)
  yss = [ f(x) for x in xss ]
  s   = spline_periodic(xs, ys)
  tss = [ s(x) for x in xss ]
  plot([scatter(; x = xss, y = yss, name = "f")
       ,scatter(; x = xss, y = tss, name = "s")
       ])
end

function sample_ellipse(a, b, points)
  ts = 2π * points
  xs = [ a * cos(t) for t in ts ]
  ys = [ b * sin(t) for t in ts ]
  xs, ys
end

function sample_mysterious()
  xs = [ 3.7, 3.2, 2.7, 2.1, 1.7, 1.1, 0.7, 0.4
       , 0.4, 0.5, 0.3, 0.6, 0.6, 0.7, 0.6, 0.8
       , 0.8, 0.6, 0.8, 0.9, 1.1, 1.4, 1.8, 1.7
       , 1.9, 2.2, 2.1, 2.7, 2.6, 3.3, 3.5, 3.7
       , 3.9, 4.2, 4.3, 4.5, 4.7, 5.0, 5.5, 5.9
       , 6.2, 6.4, 6.3, 6.5, 6.8, 7.2, 7.1, 7.2
       , 6.8, 6.7, 6.8, 6.4, 6.2, 6.9, 6.8, 6.6
       , 6.5, 6.1, 5.5, 5.0, 4.6, 4.1, 3.7, 3.4
       , 3.2, 3.7
       ]

  ys = [ 6.4, 6.7, 6.5, 6.4, 6.0, 5.9, 5.7, 5.7
       , 5.4, 5.0, 4.6, 4.3, 4.0, 3.7, 3.2, 2.9
       , 2.6, 2.4, 2.3, 2.4, 2.2, 2.1, 2.0, 1.8
       , 1.4, 1.5, 1.8, 1.6, 1.4, 1.3, 0.9, 0.6
       , 0.8, 0.7, 0.4, 0.5, 0.7, 0.6, 0.8, 0.6
       , 0.4, 0.3, 0.7, 1.2, 1.7, 2.0, 2.2, 2.4
       , 2.8, 3.2, 3.6, 3.9, 4.2, 4.5, 5.1, 5.6
       , 6.0, 6.2, 6.1, 6.2, 6.2, 6.3, 6.0, 6.1
       , 6.5, 6.4
       ]

  xs, ys
end

function sample_lissajous(a, b, points)
  ts = 2π * points
  xs = [ cos(a * t) for t in ts ]
  ys = [ sin(b * t) for t in ts ]
  xs, ys
end

function sample_hypotrochoid(R, r, d, points)
  ts = (R-r+1) * 2 * π * points
  xs = [ (R - r) * cos(t) + d * cos((R - r) / r * t) for t in ts ]
  ys = [ (R - r) * sin(t) - d * sin((R - r) / r * t) for t in ts ]
  xs, ys
end

function show_curve(curve)
  xs, ys = curve(linspace(0, 1, 1000))
  scatter(; x = xs, y = ys, name = "krzywa")
end

function test_curve(curve, mode, points)
  xs, ys = curve(points)
  n      = size(xs, 1)
  sx, sy, ti = interpolate_curve(xs, ys, mode)
  ts = linspace(0, ti[end], 1000)
  fx = [ sx(t) for t in ts ]
  fy = [ sy(t) for t in ts ]
  scatter(; x = fx, y = fy, name = "interpolant w $n punktach")
end

function interpolate_curve(xs, ys, mode)
  n = size(xs, 1)
  if     mode == "equidistant"
    ts = linspace(0, n-1, n)
  elseif mode == "cumulative"
    d  = similar(xs)
    d[1] = zero(eltype(xs))
    for i in 2 : n
      d[i] = d[i-1] + sqrt((xs[i] - xs[i-1]) ^ 2 + (ys[i] - ys[i-1]) ^ 2)
    end
    ts = [ p / d[n] for p in d ]
  end
  sx, memo = spline_periodic(ts, xs)
  sy, _    = spline_periodic(ts, ys; memoized = memo)
  sx, sy, ts
end
