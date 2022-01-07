function build_ind_agent!(mod)
# Solver settings
# set_optimizer_attribute(mod, "NumericFocus",3)
set_optimizer_attribute(mod, "OutputFlag",0)

# Extract sets
JY = mod.ext[:sets][:JY]

# Extract parameters
A = mod.ext[:parameters][:A]
e = mod.ext[:parameters][:e]
λ_EUA = mod.ext[:parameters][:λ_EUA]
b_bar = mod.ext[:parameters][:b_bar]
ρ_EUA = mod.ext[:parameters][:ρ_EUA]

# Define variables
b = mod.ext[:variables][:b] = @variable(mod, [jy=JY], lower_bound=0, base_name="EUA") 

# Objective
@objective(mod, Min, sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) +sum(ρ_EUA/2*A[jy]*(b[jy] - b_bar[jy])^2 for jy in JY))
      
# Constraints 
mod.ext[:constraints][:con1]  = @constraint(mod,[jy=JY], 
    sum(b[y2] for y2=1:jy) >= sum(e[y2] for y2=1:jy)
)

optimize!(mod);

return mod
end
