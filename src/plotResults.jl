export plotResults

"Generate plots and pdfs for specified quantities."
function plotResults( )
constants = L"""$\mathrm{Ca}= %$(round(Ca; sigdigits=1))$, $ M= %$(round(M; digits=1))$,
            $ \sigma = %$(round(T; sigdigits=2))$, $\Delta P = %$(round(P; sigdigits=2))$
            """
# Plot vm(Sw)
#     p_vm_Sw = plot(Sw, [vm,  vm_Sw_f_line, vm_linear_fit], label=[L"v_m" L"v_{m,\mathrm{p}}" L"v_{m,\mathrm{lin}}" ], legendfontsize=14, xlabel=L"S_w", ylabel=L"v_m", xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, legend=:bottomright, palette = palette(:tab10) )
#     plot!(Sw, vm_Sw_f_line_c0, ls = :dashdot, label=L"v_{m,p,0}")
# # Optional colorscheme: cgrad(:jet, 3, categorical=:true)
# annotate!(0.5,maximum([vm vm_Sw_f_line])*(1 - 0.20),text(constants, 14))
# display(p_vm_Sw)
# savefig(p_vm_Sw, "../images/article/vmSw-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")


    # Plot vm(v_d)
#         p_vm_vd = scatter(v_d, [vm], label=L"v_m", legendfontsize=14, xlabel=L"v^{\prime}", ylabel=L"v_m", xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, legend=:topleft, palette = palette(:tab10) )
#         plot!(v_d, [vm_Sw_f_line, vm_linear_fit], label=[L"v_{m,\mathrm{p}}" L"v_{m,\mathrm{lin}}"] )
#     plot!(v_d, vm_Sw_f_line_c0, ls = :dashdot, lw=1, label=L"v_{m,p,0}")
# annotate!((0.65, 0.9),text(constants, 14) )
# savefig(p_vm_vd, "../images/article/vmvd-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")
# display(p_vm_vd)

# Recreated velocities
    # p_rec = scatter(Sw, [vw, vn], label=[L"v_w"   L"v_n"  ], ylimits=(-Inf, maximum([vw vn]).*1.1), legendfontsize=10, xlabel=L"S_w", ylabel=L"v_w, v_n",xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, palette = palette(:tab10)  )
    #     plot!(Sw, [vw_rec_hom, vn_rec_hom], label=[L"v_{w,\mathrm{rec}}" L"v_{n,\mathrm{rec}}"],  linecolor = [palette(:tab10)[1] palette(:tab10)[2] ])
    #     annotate!((0.5, 0.9), text(constants, 13))
    # display(p_rec)
    # savefig(p_rec, "../images/article/rec-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")

    # Misc. geometric quantities

    # Saturations
#     p_sat = plot(Sw, [Sw, Sw_tilde, Sn, Sn_tilde], label=[L"S_w" L"\tilde{S}_w" L"S_n" L"\tilde{S}_n"], legendfontsize=12, xlabel=L"S_w", ylabel=L"Saturations", xguidefontsize=17,  yguidefontsize=17, legend=:bottom, palette = palette(:tab10) )
#     p_sat = scatter(Sw, [Sw, Sn], label=[L"S_w"   L"S_n"  ],  legendfontsize=10, xlabel=L"S_w", ylabel=L"\mathrm{Saturations}",xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, legend=:bottom, palette = palette(:tab10)  )
#     plot!(Sw, [Sw_tilde, Sn_tilde], label=[L"\tilde{S}_w" L"\tilde{S}_n"],  linecolor = [palette(:tab10)[1] palette(:tab10)[2] ])
# annotate!((0.5,0.9),text(constants, 13) )
# display(p_sat)
#     savefig(p_sat, "../images/article/sat-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")

    # Cross ratios
#     p_cross = scatter(Sw, -Sw./Sn.*τ, label=L"k" ,  legendfontsize=10, xlabel=L"S_w", ylabel=L"\mathrm{Cross-ratio}",xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, palette = palette(:tab10)  )
#     plot!(Sw, Sw_tilde./Sn_tilde, label=L"k_e",  linecolor = [palette(:tab10)[1] palette(:tab10)[2] ])
# annotate!((0.5,0.9),text(constants, 13) )
# display(p_cross)
#     savefig(p_cross, "../images/article/cross-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")


    # τ
#     p_tau = plot(Sw, τ.^-1, label=:none ,  legendfontsize=10, xlabel=L"S_w", ylabel=L"τ^{-1}",xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, palette = palette(:tab10)  )
# annotate!((0.5,0.9),text(constants, 13) )
# display(p_tau)
# savefig(p_tau, "../images/article/tau-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")


# Plot vm(v_d), linear fits
# p_vm_vd_fits = scatter(v_d, [vm,  vm_Sw_f_line  ], label=[L"v_m" L"v_{m,\mathrm{proj}}" ], legendfontsize=14, xlabel=L"v^{\prime}", ylabel=L"v_m", xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, legend=:topleft, palette = palette(:tab10) )
# plot!(v_d, [vm_linear_fit, vm_hom_linear_fit], label=[ L"v_{m,\mathrm{lin}}" L"v_{m,\mathrm{lin}, \mathrm{p}}"]  )
# annotate!((0.65, 0.9),text(constants, 14) )
# savefig(p_vm_vd_fits, "../images/article/vmvdfits-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")
# display(p_vm_vd_fits)

# Plot vm(v_d), c=0
#     p_vm_vd_fits = plot(v_d, [vm,  vm_Sw_f_line, vm_Sw_f_line_c0 ], label=[L"v_m" L"v_{m,\mathrm{proj}}" L"v_{m,\mathrm{linear}}" ], legendfontsize=14, xlabel=L"v^{\prime}", ylabel=L"v_m", xguidefontsize=17,  yguidefontsize=17, xtickfontsize=14, ytickfontsize=14, legend=:topleft, palette = palette(:tab10) )
#     annotate!(maximum(v_d) - .3maximum(v_d),minimum([vm vm_Sw_f_line]),text(constants, 14) )
# annotate!((0.65, 0.9),text(constants, 14) )
# savefig(p_vm_vd_fits, "../images/article/vmvdfits-M-"*string(M)*"-T-"*string(T)*"-P-"*string(P)*".pdf")
# display(p_vm_vd_fits)
end

    # Plot relperm data
function myplotrelperm(file::String, quants::RelpermQuantities, Asw, Asn)
    filename = basename(file)
    (; σ, μ_w, μ_n, K, ϕ, ΔP) = quants
    darcy = 9.869 * 10^-13
    K = K * darcy
    ΔP = ΔP * 10^6 # To Pa
    v0 = K * abs.(ΔP) / ϕ
    dat, krw, krn = getrelpermdata(file, quants)
    results, data, vmetc, vw_rec_hom_c0, vn_rec_hom_c0 = computeall(dat, quants, Asw, Asn, krw, krn) # Overloaded version
    v_max = maximum(abs.(vmetc.vm[1:end]))
    vd_max = maximum(abs.(vmetc.v_d[1:end]))

    # Plot vm as a function of Sw
    p1 = scatter((dat.Sw), [vmetc.vm[1:end] ./ v_max, results.vm_Sw_f_line ./ v_max, results.vm_linear_fit ./ v_max, (vmetc.v_d .- (1.0./(1 - 0.35)).*dat.v.*(krw ./μ_w .- krn./μ_n)./( krn./μ_n .+ krw./μ_w))./v_max], label=[L"v_m" L"v_{m,\mathrm{p}}" L"v_{m,\mathrm{lin}}"], legendfontsize=14, xlabel=L"S_w", ylabel=L"v_m", xguidefontsize=17, yguidefontsize=17, xtickfontsize=9, ytickfontsize=9, palette=palette(:tab10), shape=[:circle :xcross :square], ma=[nothing nothing 1], ms=[5 8 3.5], msw=[1 2 1])
    # For manual fit of linear vm
    # p1 = scatter((dat.Sw), [vmetc.vm[1:end] ./ v_max, results.vm_Sw_f_line ./ v_max, 0.6.*vmetc.v_d ./ v_max], label=[L"v_m" L"v_{m,\mathrm{p}}" L"v_{m,\mathrm{lin}}"], legendfontsize=14, xlabel=L"S_w", ylabel=L"v_m", xguidefontsize=17, yguidefontsize=17, xtickfontsize=9, ytickfontsize=9, palette=palette(:tab10), shape=[:circle :xcross :square], ma=[nothing nothing 1], ms=[5 8 3.5], msw=[1 2 1])
    # savefig(p1, "../images/article/relperm-vmSw-" * filename * ".pdf")

    # Plot vm as a function of v'
    # p2 = scatter((vmetc.v_d) ./ vd_max, [vmetc.vm[1:end] ./ v_max, results.vm_Sw_f_line ./ v_max, 0.6 .* vmetc.v_d ./ v_max], label=[L"v_m" L"v_{m,\mathrm{p}}" L"v_{m,\mathrm{lin}}"], legendfontsize=14, xlabel=L"v^{\prime}", ylabel=L"v_m", xguidefontsize=17, yguidefontsize=17, xtickfontsize=9, ytickfontsize=9, palette=palette(:tab10), shape=[:circle :xcross :square], ma=[nothing nothing 1], ms=[5 8 3.5], msw=[1 2 1])
    # savefig(p2, "../images/article/relperm-vmvd-" * filename * ".pdf")

    # Plot recreated relperms
    p_rec = scatter(dat.Sw, [krw, krn], label=[L"k_{rw}"   L"k_{rn}"  ], legendfontsize=9, xlabel=L"S_w", ylabel=L"k_{rw}, k_{rn}",xguidefontsize=17,  yguidefontsize=17, xtickfontsize=9, ytickfontsize=9, palette = palette(:tab10)  )
    plot!(dat.Sw, [(results.vw_rec_hom) .* dat.Sw .* μ_w .* 10^-3 ./ v0, (results.vn_rec_hom) .* dat.Sn .* μ_n .* 10^-3 ./ v0], label=[L"k_{rw,\mathrm{rec}}" L"k_{rn,\mathrm{rec}}"], linecolor=[palette(:tab10)[1] palette(:tab10)[2]])
    # Repeat, set viscosities equal.
    file_data = RelpermQuantities(σ, μ_w, μ_w, K./darcy, ϕ, ΔP*10^-6)
    dat, krw, krn = getrelpermdata(file, file_data)
    results, data, vmetc, vm = computeall(dat, file_data, Asw, Asn, krw, krn) # Overloaded version
    plot!(dat.Sw, [(results.vw_rec_hom) .* dat.Sw .* μ_w .* 10^-3 ./ v0, (results.vn_rec_hom) .* dat.Sn .* μ_w .* 10^-3 ./ v0], label=[L"k_{rw,\mathrm{rec}, M=1}" L"k_{rn,\mathrm{rec}, M=1}"], linecolor=[palette(:tab10)[1] palette(:tab10)[2]], ls=:dash)
    ylims!(-0.1, 1.1)
    # savefig(p_rec, "../images/article/relperm-rec-" * filename * ".pdf")
    display(p1)
end


