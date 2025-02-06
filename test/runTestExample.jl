using VmProjective
using Plots

# Parameters for 2d DNM data
# M_vec = ["00.5", "01.0", "02.0", "10.0"]
# T_vec = ["2.0", "3.0", "4.0"]
# P_vec = ["0.20", "0.40", "0.80", "1.00"]
# MTP = ["02.0", "3.0", "0.20"]

# NOTE: check the order of the columns in the data! If necessary, modify the order in 

# Relative permeability dataset
file = "<path_to_data>"
file_data = RelpermQuantities(34.24, 0.65, 0.05, 1810 * 10^-3, 0.222, -30 * 10^3 * 10 - 6) #
dat, krw, krn = getrelpermdata(file, file_data)
Asw = 0.0;
Asn = 0.0;
results, data, vmetc, vm = computeall(dat, file_data, Asw, Asn, krw, krn) # Overloaded version
myplotrelperm(file, file_data, Asw , Asn)
Ca = (10^-3).*(file_data.μ_w  .* dat.vw .+ file_data.μ_n .* dat.vn) ./ (file_data.σ .* 10^-3) 
