# Problem Set 3 — Deuda con default (Arellano)
using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()   # instala los paquetes si faltan
cd(@__DIR__)                        # trabaja en la carpeta del archivo

include("simuldef.jl")              # carga simuldef.jl, que a su vez incluye arellano.jl y deuda.jl
mkpath("figures")                   # carpeta de gráficos

# Corrección del spread (mail del profesor, 17/11). simuldef.jl usa (1+r)(1/q - 1), la fórmula del modelo con
# deuda de largo plazo, donde el bono libre de riesgo vale 1. Acá el bono promete pagar 1, así que sin riesgo vale
# 1/(1+r), y el spread trimestral es la tasa implícita 1/q - 1 menos r. Con q = 1/(1+r) da cero, como debe ser.
get_spread(qt, dd::Arellano) = 1 / qt - (1 + dd.pars.r)


# ===================== Constantes de todo el trabajo =====================
const SEMILLA = 1234                                   # semilla de las simulaciones
const TOL     = 1e-6                                   # tolerancia de convergencia (la de mpe!)
const MAXITER = 1000                                   # tope de iteraciones
const BETA0, YH0, CHI0 = 0.953, 0.969, 0.01            # parámetros por defecto del constructor
const BETAS  = [0.90, 0.92, 0.94, 0.953, 0.96, 0.97]   # paciencia (incluye el valor base)
const YHS    = [0.93, 0.95, 0.969, 0.98, 0.99]         # costo de default original (incluye el valor base)
const CHIS   = [0.0001, 0.1]                           # χ que fija la consigna en 1.2
const B_MAX_12 = 0.4                                   # deuda máxima que se muestra en los gráficos de 1.2
const BMAX_ROB = 1.2                                   # grilla de deuda más ancha para el chequeo de 1.3
const T_SIM, T_BURN, N_SIM = 10_000, 1_000, 400        # largo de cada simulación, períodos descartados y repeticiones
const UMBRAL = 0.5                                     # probabilidad que define b*(y)
const ROJO, AZUL, VERDE = "#EF553B", "#636efa", "#00cc96"
const RAMPA  = ["#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]

Random.seed!(SEMILLA)

# Todo lo que se imprime queda también en resultados.txt
LOG = open("resultados.txt", "w")
informar(xs...) = (println(xs...); println(LOG, xs...); flush(LOG))
r4(x) = round(x, digits = 4)


# ===================== Funciones comunes =====================

# Igual que mpe! del curso, pero con la opción de no actualizar el precio de la deuda (para 1.3).
# Devuelve la cantidad de iteraciones y la distancia final.
function resolver!(dd::Default; actualizar_q = true, tol = TOL, maxiter = MAXITER)
    new_v, new_q = similar(dd.v), copy(dd.q)
    dist, iter = 1 + tol, 0
    while dist > tol && iter < maxiter
        iter += 1
        actualizar_q && q_iter!(new_q, dd)                       # precio dada la probabilidad de default
        dist_q = norm(new_q - dd.q) / max(1, norm(dd.q))
        itp_q  = make_itp(dd, new_q)
        vfi_iter!(new_v, itp_q, dd)                              # valor dado el precio
        dist_v = norm(new_v - dd.v) / max(1, norm(dd.v))
        dist   = max(dist_q, dist_v)
        dd.v .= new_v
        dd.q .= new_q
    end
    return iter, dist
end

# Elección de b' robusta para los modelos con default. El opt_value del curso busca un máximo con Optim sobre
# toda la grilla; cuando el objetivo tiene más de un pico (pasa con χ muy chico y con precios fijos) salta de un
# pico a otro entre iteraciones y la función de valor no converge. Con ROBUSTO[] = true primero se evalúa el
# objetivo en una grilla gruesa (uno de cada PASO_ROB puntos de bgrid) y después se refina con Optim alrededor
# del mejor punto. Con ROBUSTO[] = false se usa el opt_value del curso sin cambios.
const ROBUSTO  = Ref(false)
const PASO_ROB = 4
function opt_value(jb, jy, q, itp_v, dd::Arellano)
    ROBUSTO[] || return invoke(opt_value, Tuple{Any, Any, Any, Any, Deuda}, jb, jy, q, itp_v, dd)
    obj_f(bpv) = eval_value(jb, jy, bpv, q, itp_v, dd)[1]
    Nb = length(dd.bgrid)
    jmax, fmax = 1, -Inf
    for j in [1:PASO_ROB:Nb; Nb]
        f = obj_f(dd.bgrid[j])
        f > fmax && ((jmax, fmax) = (j, f))
    end
    lo, hi = dd.bgrid[max(jmax - PASO_ROB, 1)], dd.bgrid[min(jmax + PASO_ROB, Nb)]
    b_star = Optim.maximizer(Optim.maximize(obj_f, lo, hi))
    fmax > obj_f(b_star) && (b_star = dd.bgrid[jmax])       # Optim no evalúa los extremos del intervalo
    vp, c_star = eval_value(jb, jy, b_star, q, itp_v, dd)
    return vp, c_star, b_star
end

# Crea y resuelve un modelo de Arellano con costos originales
function modelo(nombre; actualizar_q = true, robusto = false, kw...)
    dd = Arellano(; kw...)
    t0 = time()
    ROBUSTO[] = robusto
    iter, dist = resolver!(dd; actualizar_q = actualizar_q)
    ROBUSTO[] = false
    informar(rpad(nombre, 34), "| ", iter, " iteraciones | dist = ", round(dist, sigdigits = 3),
             " | ", round(time() - t0, digits = 1), " s")
    return dd
end

# b*(y): nivel de deuda donde la probabilidad de default cruza el umbral. Para cada y se interpola P(b, y) en b y
# se minimiza b -> (P(b, y) - umbral)^2 (opción sugerida en la consigna). La búsqueda se hace entre los dos puntos
# de la grilla donde P cruza el umbral, porque fuera de ahí la función es casi plana (P = 0 o P = 1) y el
# minimizador podría quedarse en esa zona. NaN si P no llega al umbral dentro de la grilla.
function bstar(dd; umbral = UMBRAL)
    bs = fill(NaN, length(dd.ygrid))
    for jy in eachindex(dd.ygrid)
        p  = dd.prob[:, jy]
        jb = findfirst(>=(umbral), p)
        if jb === nothing
            continue
        elseif jb == 1
            bs[jy] = dd.bgrid[1]
        else
            itp = interpolate((dd.bgrid,), p, Gridded(Linear()))
            res = Optim.optimize(b -> (itp(b) - umbral)^2, dd.bgrid[jb-1], dd.bgrid[jb])
            bs[jy] = Optim.minimizer(res)
        end
    end
    return bs
end

# Gráfico de varias series (sin título). La serie número `arriba` se dibuja por encima de las demás.
# Con chico = true usa letra grande y la leyenda debajo, para las figuras que van de a tres por fila.
function grafico(xs, ys, nombres, colores; xt = "", yt = "", arriba = length(xs), chico = false)
    sc = [scatter(x = x, y = y, mode = "lines+markers", name = n, line_color = c, marker_color = c, legendrank = k)
          for (k, (x, y, n, c)) in enumerate(zip(xs, ys, nombres, colores))]
    orden = [setdiff(1:length(sc), arriba); arriba]
    lay = chico ? Layout(xaxis_title = xt, yaxis_title = yt, margin = attr(t = 20), font_size = 22,
                         legend = attr(orientation = "h", x = 0, y = -0.3))  :
                  Layout(xaxis_title = xt, yaxis_title = yt, margin = attr(t = 20))
    return plot(sc[orden], lay)
end
guardar_chico(p, archivo) = savefig(p, archivo, width = 520, height = 520)

# Contour de la probabilidad de default, con la curva b*(y) encima.
# En contour las filas de z van al eje vertical: prob[jb, jy] => x = ygrid, y = bgrid.
function grafico_contour(dd, archivo)
    bs = bstar(dd); m = isfinite.(bs)
    tr = [contour(x = dd.ygrid, y = dd.bgrid, z = dd.prob, colorbar_title = "P(b, y)")]
    any(m) && push!(tr, scatter(x = dd.ygrid[m], y = bs[m], mode = "lines+markers", name = "b*(y)",
                                line_color = "white", marker_color = "white", showlegend = false))
    savefig(plot(tr, Layout(xaxis_title = "y", yaxis_title = "b", margin = attr(t = 20))), archivo)
end

# Gráfico de b*(y) para una familia de modelos que difieren en un parámetro
# Con dif = true grafica la diferencia con el b*(y) del modelo base (sirve cuando las curvas se tapan).
function grafico_bstar(mods, valores, base, simbolo, archivo; dif = false)
    xs, ys, ns, cs = [], [], [], []
    bs0 = bstar(mods[base])
    k = 0
    for val in valores
        dd = mods[val]; bs = bstar(dd)
        if val == base
            dif && continue
            m = isfinite.(bs)
            push!(xs, dd.ygrid[m]); push!(ys, bs[m]); push!(ns, "$simbolo = $val (base)"); push!(cs, ROJO)
        else
            k += 1
            m = dif ? isfinite.(bs .- bs0) : isfinite.(bs)
            push!(xs, dd.ygrid[m]); push!(ys, dif ? (bs .- bs0)[m] : bs[m]); push!(ns, "$simbolo = $val"); push!(cs, RAMPA[k])
        end
    end
    arriba = dif ? length(xs) : findfirst(==(base), valores)
    yt = dif ? "b*(y) − b*(y) con $simbolo = $base" : "b*(y)"
    savefig(grafico(xs, ys, ns, cs, xt = "y", yt = yt, arriba = arriba), archivo)
end

fila(v) = join([lpad(isnan(x) ? "—" : string(round(x, digits = 3)), 7) for x in v])


# ===================== 1.1 Regiones de repago y default =====================
informar("===== 1.1 Regiones de repago y default =====")
dd0 = modelo("Base")
JY = argmin(abs.(dd0.ygrid .- 1))                      # índice del ingreso más cercano a 1
informar("ygrid: ", length(dd0.ygrid), " puntos de ", r4(dd0.ygrid[1]), " a ", r4(dd0.ygrid[end]),
         " | y[", JY, "] = ", r4(dd0.ygrid[JY]))
informar("bgrid: ", length(dd0.bgrid), " puntos de ", dd0.bgrid[1], " a ", dd0.bgrid[end],
         " | paso = ", r4(dd0.bgrid[2] - dd0.bgrid[1]))

grafico_contour(dd0, "figures/prob_contour.pdf")

informar("P creciente en b para todo y: ", all(diff(dd0.prob, dims = 1) .>= -1e-8),
         " | P decreciente en y para todo b: ", all(diff(dd0.prob, dims = 2) .<= 1e-8))
informar("máx de P que sube con y: ", round(maximum(diff(dd0.prob, dims = 2)), sigdigits = 3))
informar("y      = ", fila(dd0.ygrid))
informar("b*(y)  = ", fila(bstar(dd0)))
informar("q(b',y) en y[", JY, "]: ", fila([dd0.q[argmin(abs.(dd0.bgrid .- bv)), JY] for bv in 0:0.1:0.6]))

# Opcional: b*(y) para distintos β y distintos costos de default
informar("--- Opcional 1.1: variar β ---")
mods_beta = Dict(β => (β == BETA0 ? dd0 : modelo("β = $β"; β = β)) for β in BETAS)
informar("--- Opcional 1.1: variar yh ---")
mods_yh = Dict(yh => (yh == YH0 ? dd0 : modelo("yh = $yh"; yh = yh)) for yh in YHS)

informar("b*(y) por β.   y = ", fila(dd0.ygrid))
for β in BETAS
    informar("  β = ", rpad(β, 13), fila(bstar(mods_beta[β])))
end
informar("b*(y) por yh.  y = ", fila(dd0.ygrid))
for yh in YHS
    informar("  yh = ", rpad(yh, 12), fila(bstar(mods_yh[yh])))
end
grafico_bstar(mods_beta, BETAS, BETA0, "β", "figures/bstar_beta.pdf")                     # niveles (las curvas se tapan)
grafico_bstar(mods_beta, BETAS, BETA0, "β", "figures/bstar_beta_dif.pdf", dif = true)     # diferencia con el base
grafico_bstar(mods_yh, YHS, YH0, "yh", "figures/bstar_yh.pdf")


# ===================== 1.2 Mecánica de los shocks de preferencias =====================
informar("===== 1.2 Shocks de preferencias =====")
mods_chi = Dict(χ => modelo("χ = $χ"; χ = χ, robusto = (χ < CHI0)) for χ in CHIS)
mods_chi[CHI0] = dd0

# Chequeo: vR y vD se cruzan dentro del gráfico (b ≤ B_MAX_12) con y = 1
for χ in sort(collect(keys(mods_chi)))
    dd = mods_chi[χ]
    jb = findfirst(dd.vR[:, JY] .< dd.vD[JY])
    informar("χ = ", rpad(χ, 7), "| primer b con vR < vD: ", jb === nothing ? "ninguno" : r4(dd.bgrid[jb]),
             " | máx de V - max(vR, vD) = ", round(maximum(dd.v[:, JY] .- max.(dd.vR[:, JY], dd.vD[JY])), sigdigits = 3),
             " | χ log 2 = ", round(χ * log(2), sigdigits = 3), " | P(0, y=1) = ", r4(dd.prob[1, JY]))
end

# Un solo gráfico para los dos χ: color = función de valor, trazo = χ. vD no depende de b, así que va como una
# línea horizontal (shape de PlotlyJS); las trazas de un solo punto sólo sirven para que aparezca en la leyenda.
function grafico_valores(mods, chis, archivo)
    trazos = ["solid", "dash"]
    tr_vR, tr_V, formas = [], [], []
    for (χ, dash) in zip(chis, trazos)
        dd = mods[χ]; m = dd.bgrid .<= B_MAX_12; b = dd.bgrid[m]
        push!(tr_vR, scatter(x = b, y = dd.vR[m, JY], mode = "lines", name = "v<sup>R</sup>(b, y), χ = $χ",
                             line = attr(color = AZUL, dash = dash, width = 2.5)))
        push!(tr_vR, scatter(x = [b[1]], y = [dd.vD[JY]], mode = "lines", name = "v<sup>D</sup>(y), χ = $χ",
                             line = attr(color = VERDE, dash = dash, width = 2.5)))
        push!(formas, attr(type = "line", layer = "below", xref = "x", x0 = b[1], x1 = b[end], yref = "y", y0 = dd.vD[JY], y1 = dd.vD[JY],
                           line = attr(color = VERDE, dash = dash, width = 2.5)))
        push!(tr_V, scatter(x = b, y = dd.v[m, JY], mode = "lines", name = "V(b, y), χ = $χ",
                            line = attr(color = ROJO, dash = dash, width = 2.5)))
    end
    lay = Layout(xaxis_title = "b", yaxis_title = "valor", margin = attr(t = 20), shapes = formas)
    savefig(plot(map(identity, [tr_vR; tr_V]), lay), archivo)       # V va al final para quedar arriba
end
grafico_valores(mods_chi, CHIS, "figures/valores_chi.pdf")


# ===================== 1.3 Cuántos defaults son por los precios? =====================
informar("===== 1.3 Defaults por precios =====")
# Opción con jerarquía de tipos. SinRiesgo envuelve un costo de default, así que Arellano{SinRiesgo{OG}} sigue siendo
# un Arellano <: Default: usa h, value_default y vfi_iter! del modelo con default. Lo único que cambia, por
# multiple dispatch, es debtprice, que pasa a ser la del modelo sin default (deuda.jl): 1/(1+r) para todo b' e y.
# Así se puede usar mpe! del curso sin tocarlo. mpe! sigue actualizando dd.q, pero el gobierno ya no lo mira.
struct SinRiesgo{C<:DefaultCost} <: DefaultCost
    base::C
end
h(yv, c::SinRiesgo) = h(yv, c.base)
debtprice(dd::Arellano{<:SinRiesgo}, _, _, _) = 1 / (1 + dd.pars.r)

function modelo_sin_riesgo(nombre; kw...)
    base = Arellano(; kw...)
    dd = Arellano{SinRiesgo{typeof(base.cost)}}(base.pars, SinRiesgo(base.cost), base.bgrid, base.ygrid, base.Py,
         base.v, base.vR, base.vD, base.prob, base.gc, base.gb, base.q)
    t0 = time()
    ROBUSTO[] = true
    mpe!(dd; maxiter = MAXITER, verbose = false)
    iter, dist = resolver!(dd)                  # chequeo: si mpe! convergió, una iteración más ya cumple la tolerancia
    ROBUSTO[] = false
    informar(rpad(nombre, 34), "| mpe! del curso | iteración de chequeo: ", iter, " (dist = ", round(dist, sigdigits = 3),
             ") | ", round(time() - t0, digits = 1), " s")
    return dd
end

ddq = modelo_sin_riesgo("Precios de repago seguro")
informar("precio que enfrenta el gobierno: ", r4(debtprice(ddq, 0.3, 1.0, nothing)), " = 1/(1+r)")
informar("y                        = ", fila(dd0.ygrid))
informar("b*(y) equilibrio         = ", fila(bstar(dd0)))
informar("b*(y) precios fijos      = ", fila(bstar(ddq)))
D0, Dq = dd0.prob .> UMBRAL, ddq.prob .> UMBRAL
informar("Puntos (b,y) con P > 0.5: equilibrio = ", sum(D0), " | precios fijos = ", sum(Dq), " | de ", length(D0))
informar("   sólo en equilibrio = ", sum(D0 .& .!Dq), " | sólo con precios fijos = ", sum(Dq .& .!D0), " | en ambos = ", sum(D0 .& Dq))
informar("P_fijos - P_equilibrio: mín = ", round(minimum(ddq.prob - dd0.prob), sigdigits = 3),
         " | máx = ", round(maximum(ddq.prob - dd0.prob), sigdigits = 3))
informar("b'(b, y=1) precios fijos: ", fila([ddq.gb[argmin(abs.(ddq.bgrid .- bv)), JY] for bv in 0:0.1:0.6]))
grafico_contour(ddq, "figures/prob_contour_precios.pdf")
m0, mq = isfinite.(bstar(dd0)), isfinite.(bstar(ddq))
savefig(grafico([ddq.ygrid[mq], dd0.ygrid[m0]], [bstar(ddq)[mq], bstar(dd0)[m0]],
                ["Precios de repago seguro", "Precios de equilibrio"], [AZUL, ROJO], xt = "y", yt = "b*(y)"),
        "figures/bstar_precios.pdf")

# Chequeo: con precios fijos el gobierno se endeuda hasta el tope de la grilla, así que b*(y) puede depender de bmax
ddq2 = modelo_sin_riesgo("Precios fijos, bmax = $BMAX_ROB"; bmax = BMAX_ROB)
informar("b*(y) precios fijos, bmax = ", BMAX_ROB, ": ", fila(bstar(ddq2)))

informar("===== FIN EJERCICIO 1 =====")


# ===================== 2. Calibración del modelo de default =====================
informar("===== 2 Calibración =====")

# Momentos del curso (deuda/PBI, spread promedio y desvío) más la frecuencia anual de default
# (defaults nuevos por cada 100 años con acceso al mercado)
function momentos(pp::SimulPath)
    M = moments(pp)
    d = pp[:d]
    nuevos = count(t -> d[t] == 1 && d[t-1] == 0, 2:length(d))
    M[:deffreq] = 100 * 4 * nuevos / max(1, count(==(0), d[1:end-1]))
    return M
end

# Para un vector de simulaciones: calcula los momentos de cada una y los promedia
function momentos(pps::Vector{SimulPath})
    Ms = [momentos(pp) for pp in pps]
    return Dict(k => mean(M[k] for M in Ms) for k in keys(Ms[1]))
end

# Para un modelo resuelto: N_SIM simulaciones de T_SIM períodos, descartando los primeros T_BURN de cada una
function momentos(dd::Default; semilla = SEMILLA)
    Random.seed!(semilla)
    return momentos([subpath(simul(dd, T = T_SIM + T_BURN), T_BURN + 1, T_SIM + T_BURN) for _ in 1:N_SIM])
end

const CLAVES = [:debtGDP, :avgspr, :stdspr, :deffreq]

# Momentos de una familia de modelos, con dos semillas distintas para chequear que no dependen del azar
function tabla_momentos(mods, valores, simbolo)
    informar(rpad(simbolo, 8), join(lpad.(string.(CLAVES), 11)), "   | dif. máx. entre semillas (%)")
    Ms = []
    for val in valores
        t0 = time()
        M  = momentos(mods[val])
        M2 = momentos(mods[val], semilla = SEMILLA + 1)
        dif = maximum(abs(M2[k] / M[k] - 1) for k in CLAVES[1:3]) * 100
        informar(rpad(val, 8), join([lpad(round(M[k], digits = 2), 11) for k in CLAVES]), "   | ", round(dif, digits = 2),
                 "  (", round(time() - t0, digits = 1), " s)")
        push!(Ms, M)
    end
    return Ms
end

function grafico_momentos(valores, Ms, simbolo, sufijo)
    x = valores
    col(k) = [M[k] for M in Ms]
    guardar_chico(grafico([x], [col(:debtGDP)], ["Deuda / PBI"], [AZUL], xt = simbolo, yt = "Deuda / PBI anual (%)",
                          chico = true), "figures/mom_deuda_$sufijo.pdf")
    guardar_chico(grafico([x, x], [col(:stdspr), col(:avgspr)], ["Desvío", "Promedio"], [AZUL, ROJO],
                          xt = simbolo, yt = "Spread (pb)", chico = true), "figures/mom_spread_$sufijo.pdf")
    guardar_chico(grafico([x], [col(:deffreq)], ["Frecuencia de default"], [VERDE], xt = simbolo,
                          yt = "Defaults por año (%)", chico = true), "figures/mom_def_$sufijo.pdf")
end

informar("--- Momentos en función de β (yh = ", YH0, ") ---")
Ms_beta = tabla_momentos(mods_beta, BETAS, "β")
grafico_momentos(BETAS, Ms_beta, "β", "beta")

informar("--- Momentos en función de yh (β = ", BETA0, ") ---")
Ms_yh = tabla_momentos(mods_yh, YHS, "yh")
grafico_momentos(YHS, Ms_yh, "yh", "yh")

informar("===== FIN =====")
close(LOG)

# Copia las figuras a informe/ para que esa carpeta compile sola (el logo, logo_udesa.jpg, ya está ahí)
mkpath("informe")
for f in readdir("figures")
    cp(joinpath("figures", f), joinpath("informe", f), force = true)
end
