function solve_ind_agent!(mod)
# Extract sets
JY = mod.ext[:sets][:JY]

# Extract parameters
A = mod.ext[:parameters][:A]
e = mod.ext[:parameters][:e]
AC = mod.ext[:parameters][:AC]
λ_EUA = mod.ext[:parameters][:λ_EUA]
b_bar = mod.ext[:parameters][:b_bar]
ρ_EUA = mod.ext[:parameters][:ρ_EUA]

# Define variables
b = mod.ext[:variables][:b] 

# Update expressions
mod.ext[:expressions][:tot_cost] = @expression(mod, 
    sum(A[jy]*AC[jy] for jy in JY)
)

# Update objective
@objective(mod, Min, sum(A[jy]*λ_EUA[jy]*b[jy] for jy in JY) + sum(ρ_EUA/2*(b[jy] - b_bar[jy])^2 for jy in JY))
      
# Update constraints 
for jy in JY
    delete(mod,mod.ext[:constraints][:con1][jy])
end 

mod.ext[:constraints][:con1]  = @constraint(mod,[jy=JY], 
    sum(b[y2] for y2=1:jy) >= sum(e[y2] for y2=1:jy)
)

optimize!(mod);

return mod
end
