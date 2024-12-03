using StructuralIdentifiability

# Original model

ode = @ODEmodel(
    x1'(t) =  -(betae*x2(t) + betan*x3(t) + betas*x4(t))*(x1(t)/(x1(t)+x2(t)+x3(t)+x4(t)+x8(t))) - eta*x1(t),
    x2'(t) = (betae*x2(t) + betan*x3(t) + betas*x4(t))*(x1(t)/(x1(t)+x2(t)+x3(t)+x4(t)+x8(t))) - (1/14)*x2(t),
    x3'(t) = (1-rho)*(1/14)*x2(t) - (gamman)*x3(t),
    x4'(t) = rho*(1/14)*x2(t) - (alpha1+gammas)*x4(t),
    x5'(t) = alpha1*x4(t) - (alpha2+gammah + muh)*x5(t),
    x6'(t) = alpha2*x5(t) - (alpha3+gammaICU+muICU)*x6(t),
    x7'(t) = alpha3*x6(t) - (gammav + muv)*x7(t),
    x8'(t) = gamman*x3(t) + gammas*x4(t) + gammah*x5(t) + gammaICU*x6(t) + gammav*x7(t),
    x9'(t) = eta*x1(t),
    y1(t) = (1/14)*rho*x2(t),
    y2(t) = muh*x5(t) + muv*x7(t) + muICU*x6(t),
    y3(t) = x6(t),
    y4(t) = x7(t),
)

#println(assess_identifiability(ode))
println(assess_identifiability(ode, known_ic=[x1,x2,x3,x4,x5,x6,x7,x8,x9]))


# Transformation #1
# Observe that: x1' + x2' + x9' +  y1 / rho = 0
# Therefore: x1 + x2 + x9 + int_y1 / rho = C1 for some constant C1
# We introduce new state int_y1' = 1/14 * rho * x2 and make it an output
# Then we can use x1 = C1 - x2 - x9 - int_y1 / rho to eliminate x1 from the system

# Transformation #2
# Observe that: x3' + x4' + x5' + x6' + x7' + x8' -  y1 / rho + y2 = 0
# Therefore: x3 + x4 + x5 + x6 + x7 + x8 - int_y1 / rho + int_y2 = C2 for some constant C2
# We introduce new state int_y2' = muh*x5 + muICU*x6 + muv*x7 and make it an output
# Then we can use x4 = C2 - x3 - x5 - x6 - x7 - x8 + int_y1 / rho - int_y2 to eliminate x4 from the system
ode = @ODEmodel(
    x2'(t) = (betae*x2(t) + betan*x3(t) + betas*(C2 - (x3(t)+x5(t)+x6(t)+x7(t)+x8(t)) + int_y1(t) / rho - int_y2 ))*((C1 - x2(t) - x9(t) - int_y1(t) / rho)/((C1 - x2(t) - x9(t) - int_y1(t) / rho)+x2(t)+x3(t)+(C2 - (x3(t)+x5(t)+x6(t)+x7(t)+x8(t)) + int_y1(t) / rho - int_y2 )+x8(t)+x9(t))) - (1/14)*x2(t),
    x3'(t) = (1-rho)*(1/14)*x2(t) - (gamman)*x3(t),
    x5'(t) = alpha1*(C2 - (x3(t)+x5(t)+x6(t)+x7(t)+x8(t)) + int_y1(t) / rho - int_y2 ) - (alpha2+gammah + muh)*x5(t),
    x6'(t) = alpha2*x5(t) - (alpha3+gammaICU+muICU)*x6(t),
    x7'(t) = alpha3*x6(t) - (gammav + muv)*x7(t),
    x8'(t) = gamman*x3(t) + gammas*(C2 - (x3(t)+x5(t)+x6(t)+x7(t)+x8(t)) + int_y1(t) / rho - int_y2 ) + gammah*x5(t) + gammaICU*x6(t) + gammav*x7(t),
    x9'(t) = eta*(C1 - x2(t) - x9(t) - int_y1(t) / rho),
    int_y1'(t) = (1/14) * rho * x2(t),
    int_y2'(t) = muh*x5(t) + muICU*x6(t) + muv*x7(t),
    y1(t) = (1/14)*rho*x2(t),
    y2(t) = muh*x5(t) + muICU*x6(t) + muv*x7(t),
    y3(t) = x6(t),
    y4(t) = x7(t),
    y5(t) = int_y1(t),
    y6(t) = int_y2(t),
)

println(assess_identifiability(ode, known_ic=[x2,x3,x5,x6,x7,x8,x9]))
